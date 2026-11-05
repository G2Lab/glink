// Node matrix construction, in the order build_matrix does it: read rows, callability, inversion
// split, matrix values, samples, filters, complement dropping and pattern collapse.
#include <algorithm>
#include <cmath>
#include <set>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

#include "matrix.hpp"
#include "coordinates.hpp"
#include "conditioning.hpp"
#include "gfa.hpp"

namespace glink {

Values Values::take(const std::vector<int64_t>& keep) const {
    Values out;
    out.rows = keep.size();
    out.cols = cols;
    out.v.resize(out.rows * cols);
    for (size_t i = 0; i < keep.size(); i++) std::memcpy(out.row(i), row(keep[i]), cols);
    return out;
}

double NodeMatrix::density() const {
    size_t cells = values.rows * values.cols, stored = 0;
    for (int8_t x : values.v) stored += x != 0;
    return cells ? double(stored) / double(cells) : 0.0;
}

std::vector<std::string> NodeMatrix::row_labels() const {
    std::vector<std::string> out(node_ids.size());
    for (size_t i = 0; i < node_ids.size(); i++) {
        out[i] = std::to_string(node_ids[i]);
        if (oriented) out[i] += node_orientation[i] == FORWARD ? "+" : node_orientation[i] == REVERSE ? "-" : "";
    }
    return out;
}

namespace {

bool parse_float(const std::string& s, double& out) {
    if (s.empty()) return false;
    char* end;
    out = std::strtod(s.c_str(), &end);
    while (*end && std::isspace((unsigned char)*end)) end++;
    return *end == 0 && !std::isspace((unsigned char)s[0]);
}

}  // namespace

std::optional<Filter> parse_filter(const std::string& spec, std::string& error) {
    Filter f;
    if (spec == "all") { f.kind = "all"; return f; }
    if (spec == "majority") { f.kind = "band"; f.low = 0.5; f.high = 1.0; return f; }
    if (spec == "polymorphic") { f.kind = "band"; f.low = 0.0; f.high = 1.0; return f; }
    if (spec == "common") { f.kind = "common"; f.value = 0.05; f.people = 50; return f; }
    if (starts_with(spec, "@")) {
        f.kind = "list";
        f.path = spec.substr(1);
        if (starts_with(f.path, "~/") && std::getenv("HOME")) f.path = std::getenv("HOME") + f.path.substr(1);
        return f;
    }
    if (starts_with(spec, "minor:")) {
        if (!parse_float(spec.substr(6), f.value)) { error = "'" + spec + "' is not minor:FREQUENCY"; return {}; }
        if (!(0.0 <= f.value && f.value < 0.5)) { error = "minor frequency must be in [0, 0.5), got " + shortest(f.value); return {}; }
        f.kind = "minor";
        return f;
    }
    size_t colon = spec.find(':');
    if (colon != std::string::npos) {
        if (!parse_float(spec.substr(0, colon), f.low) || !parse_float(spec.substr(colon + 1), f.high)) {
            error = "'" + spec + "' is not a LOW:HIGH band";
            return {};
        }
        if (!(0.0 <= f.low && f.low < f.high && f.high <= 1.0)) {
            error = "band must satisfy 0 <= LOW < HIGH <= 1, got " + shortest(f.low) + ":" + shortest(f.high);
            return {};
        }
        f.kind = "band";
        return f;
    }
    error = "'" + spec + "' is not a filter; use one of all, majority, polymorphic, common, LOW:HIGH, minor:F, or @file";
    return {};
}

namespace {

std::vector<uint8_t> listed_mask(const std::vector<int64_t>& ids, const std::string& path) {
    std::unordered_set<int64_t> wanted;
    for (auto& t : split_ws(read_text(path))) {
        size_t end;
        long long v;
        try { v = std::stoll(t, &end); } catch (const std::exception&) { throw Error("invalid node id '" + t + "' in " + path); }
        if (end != t.size()) throw Error("invalid node id '" + t + "' in " + path);
        wanted.insert(v);
    }
    std::vector<uint8_t> mask(ids.size());
    for (size_t i = 0; i < ids.size(); i++) mask[i] = wanted.count(ids[i]) > 0;
    return mask;
}

// Matrix values for a block of packed calls: phased per haplotype, dosage/presence per individual.
struct Columns {
    std::vector<std::pair<size_t, size_t>> pairs;  // dosage/presence
    std::vector<std::string> names;
};

Columns columns_for(const std::vector<std::string>& names, const std::string& matrix_type) {
    Columns c;
    if (matrix_type == "phased") { c.names = names; return c; }
    std::set<std::string> samples;
    for (auto& n : names) {
        size_t hash = n.rfind('#');
        samples.insert(hash == std::string::npos ? n : n.substr(0, hash));
    }
    std::unordered_map<std::string, size_t> lookup;
    for (size_t i = 0; i < names.size(); i++) lookup.emplace(names[i], i);
    for (auto& s : samples) {
        auto a = lookup.find(s + "#1"), b = lookup.find(s + "#2");
        if (a != lookup.end() && b != lookup.end()) {
            c.pairs.emplace_back(a->second, b->second);
            c.names.push_back(s);
        }
    }
    if (c.pairs.empty()) throw Error("no sample has both #1 and #2 haplotypes");
    return c;
}

Values call_values(const Packed& bits, const Packed* called, const Columns& cols, const std::string& matrix_type,
                size_t haplotypes) {
    Values out;
    out.rows = bits.rows;
    out.cols = matrix_type == "phased" ? haplotypes : cols.pairs.size();
    out.v.resize(out.rows * out.cols);
    std::vector<int8_t> calls(haplotypes);
    for (size_t r = 0; r < bits.rows; r++) {
        const uint8_t* b = bits.row(r);
        const uint8_t* c = called ? called->row(r) : nullptr;
        for (size_t h = 0; h < haplotypes; h++) {
            int8_t present = (b[h >> 3] >> (h & 7)) & 1;
            calls[h] = !c || ((c[h >> 3] >> (h & 7)) & 1) ? present : int8_t(-1);
        }
        int8_t* o = out.row(r);
        if (matrix_type == "phased") { std::memcpy(o, calls.data(), haplotypes); continue; }
        for (size_t k = 0; k < cols.pairs.size(); k++) {
            int8_t a = calls[cols.pairs[k].first], bb = calls[cols.pairs[k].second];
            if (a < 0 || bb < 0) o[k] = -1;
            else o[k] = matrix_type == "dosage" ? int8_t(a + bb) : int8_t(a > 0 || bb > 0);
        }
    }
    return out;
}

std::vector<uint8_t> select(const Values& v, const Filter& f, const std::string& matrix_type) {
    std::vector<uint8_t> keep(v.rows);
    const int ploidy = matrix_type == "dosage" ? 2 : 1;
    for (size_t r = 0; r < v.rows; r++) {
        const int8_t* x = v.row(r);
        long long carriers = 0, called = 0, counts = 0;
        for (size_t c = 0; c < v.cols; c++) {
            carriers += x[c] > 0;
            called += x[c] >= 0;
            counts += x[c] > 0 ? x[c] : 0;
        }
        if (f.kind == "band") {
            double s = called > 0 ? double(carriers) / double(called) : NAN;
            keep[r] = s > f.low && s < f.high;
        } else if (f.kind == "minor" || f.kind == "common") {
            long long total = ploidy * called;
            double s = total > 0 ? double(counts) / double(total) : NAN;
            double minor = std::isnan(s) ? NAN : std::min(s, 1.0 - s);
            keep[r] = minor > f.value;
            if (f.kind == "common")  // two phased haplotypes count as one person
                keep[r] = keep[r] && (matrix_type == "phased" ? double(called) / 2 : double(called)) >= f.people;
        } else throw Error("unknown filter '" + f.kind + "'");
    }
    return keep;
}

std::vector<int64_t> indices_of(const std::vector<uint8_t>& mask) {
    std::vector<int64_t> out;
    for (size_t i = 0; i < mask.size(); i++) if (mask[i]) out.push_back(int64_t(i));
    return out;
}

template <class T>
std::vector<T> take(const std::vector<T>& v, const std::vector<int64_t>& keep) {
    std::vector<T> out(keep.size());
    for (size_t i = 0; i < keep.size(); i++) out[i] = v[keep[i]];
    return out;
}

// complement_free: keep a row unless an earlier row is its perfect complement.
// Also sets `partner`: each dropped row's kept complement (-1 for kept rows).
std::vector<uint8_t> complement_free(const Values& v, int ploidy, const Values* compared,
                                     std::vector<int64_t>& partner) {
    const size_t n = v.rows;
    std::vector<std::string> own(n), other(n);
    for (size_t r = 0; r < n; r++) {
        own[r].resize(v.cols);
        other[r].resize(v.cols);
        const int8_t* x = v.row(r);
        const int8_t* m = compared ? compared->row(r) : nullptr;
        for (size_t c = 0; c < v.cols; c++) {
            bool mask = x[c] >= 0 && (!m || m[c]);
            own[r][c] = mask ? char(1 + x[c]) : 0;
            other[r][c] = mask ? char(1 + ploidy - x[c]) : 0;
        }
    }
    std::unordered_map<std::string_view, size_t> first;
    for (size_t r = 0; r < n; r++) first.emplace(own[r], r);  // emplace keeps the earliest
    std::vector<uint8_t> keep(n);
    partner.assign(n, -1);
    for (size_t r = 0; r < n; r++) {
        size_t f_own = first.at(own[r]);
        auto it = first.find(other[r]);
        size_t f_other = it == first.end() ? n : it->second;
        keep[r] = f_other > f_own || (other[r] == own[r] && f_own == r);
        if (!keep[r]) partner[r] = int64_t(f_other);
    }
    return keep;
}

// --collapse-backbone-parents: a parent visited by every haplotype called at it (backbone), or
// skipped by at least broken_min_haplotypes haplotypes that visit its children (broken), is skipped,
// and its children condition on its nearest remaining ancestor. `bits` must hold every graph node,
// in parents' order. `summary` says what moved.
ParentTable collapse_backbone_parents(const ParentTable& parents, const Packed& bits, size_t haplotypes,
                                      int broken_min_haplotypes, Json& summary) {
    const size_t n = parents.node_id.size(), width = bits.width;
    IdIndex index(parents.node_id);
    std::vector<int32_t> parent(n, -1);
    for (size_t i = 0; i < n; i++)
        if (parents.parent_node_id[i] >= 0) {
            const int64_t p = index.find(parents.parent_node_id[i]);
            if (p < 0) throw Error("Conditioning tree names a parent outside the graph");
            parent[i] = int32_t(p);
        }
    // Children of each parent, and every node in tree order, level by level from the roots.
    std::vector<int64_t> first(n + 1, 0);
    for (int32_t p : parent) if (p >= 0) first[p + 1]++;
    for (size_t i = 0; i < n; i++) first[i + 1] += first[i];
    std::vector<int32_t> children(first[n]);
    {
        std::vector<int64_t> at(first.begin(), first.end() - 1);
        for (size_t i = 0; i < n; i++) if (parent[i] >= 0) children[at[parent[i]]++] = int32_t(i);
    }
    std::vector<int32_t> order;
    order.reserve(n);
    std::vector<size_t> level = {0};
    for (size_t i = 0; i < n; i++) if (parent[i] < 0) order.push_back(int32_t(i));
    while (level.back() < order.size()) {
        const size_t begin = level.back(), end = order.size();
        level.push_back(end);
        for (size_t k = begin; k < end; k++) {
            const int32_t p = order[k];
            order.insert(order.end(), children.begin() + first[p], children.begin() + first[p + 1]);
        }
    }
    if (order.size() != n) throw Error("Conditioning tree contains a cycle");

    // Padding bits past the last haplotype are ignored.
    std::vector<uint8_t> pad(width, 0xFF);
    if (haplotypes % 8) pad[width - 1] = uint8_t((1u << (haplotypes % 8)) - 1);
    // Broken: at least N haplotypes reach one of its children while skipping the parent.
    std::vector<uint8_t> broken(n), backbone(n), scratch(width);
    for (size_t p = 0; p < n; p++) {
        if (first[p] == first[p + 1]) continue;
        std::fill(scratch.begin(), scratch.end(), 0);
        for (int64_t k = first[p]; k < first[p + 1]; k++) {
            const uint8_t* row = bits.row(children[k]);
            for (size_t w = 0; w < width; w++) scratch[w] |= row[w];
        }
        const uint8_t* own = bits.row(p);
        long long skipped = 0;
        for (size_t w = 0; w < width; w++) skipped += __builtin_popcount(scratch[w] & ~own[w] & pad[w]);
        broken[p] = skipped >= broken_min_haplotypes;
    }
    // Ancestors resolve first, so a chain of skipped parents collapses to the first informative one.
    std::vector<int32_t> effective = parent;
    for (size_t d = 0; d + 1 < level.size(); d++) {
        for (size_t k = level[d]; k < level[d + 1]; k++) {
            const int32_t p = order[k];
            if (first[p] == first[p + 1]) continue;
            const uint8_t* own = bits.row(p);
            const uint8_t* above = effective[p] >= 0 ? bits.row(effective[p]) : nullptr;
            bool never_absent = true;  // every haplotype called at p visits p
            for (size_t w = 0; w < width && never_absent; w++)
                never_absent = ((above ? above[w] : uint8_t(0xFF)) & ~own[w] & pad[w]) == 0;
            backbone[p] = never_absent;
            for (int64_t c = first[p]; c < first[p + 1]; c++)
                if (backbone[p] || broken[p]) effective[children[c]] = effective[p];
        }
    }
    long long backbone_parents = 0, broken_parents = 0, moved = 0, with_parent = 0;
    ParentTable out;
    out.node_id = parents.node_id;
    out.parent_node_id.resize(n);
    for (size_t i = 0; i < n; i++) {
        backbone_parents += backbone[i] && !broken[i];
        broken_parents += broken[i];
        moved += effective[i] != parent[i];
        with_parent += effective[i] >= 0;
        out.parent_node_id[i] = effective[i] >= 0 ? parents.node_id[effective[i]] : -1;
    }
    Json c = Json::object();
    c["haplotypes"] = (long long)haplotypes;
    c["broken_min_haplotypes"] = broken_min_haplotypes;
    c["backbone_parents"] = backbone_parents;
    c["broken_parents"] = broken_parents;
    c["moved_nodes"] = moved;
    c["with_parent"] = with_parent;
    summary = c;
    return out;
}

// The graph rows a run looks at. A node-id filter needs no genotypes, so it is applied first and
// the other rows are never read.
struct GraphRows {
    std::vector<int64_t> ids;
    Packed bits, forward, reverse;
};

std::optional<std::vector<int64_t>> listed_rows(const std::vector<int64_t>& ids, const std::vector<Filter>& filters) {
    std::optional<std::vector<int64_t>> rows;
    for (const Filter& f : filters) {
        if (f.kind != "list") continue;
        const std::vector<int64_t> current = rows ? take(ids, *rows) : ids;
        const std::vector<int64_t> keep = indices_of(listed_mask(current, f.path));
        if (keep.empty()) throw Error("filter list kept no nodes out of " + thousands(current.size()));
        rows = rows ? take(*rows, keep) : keep;
    }
    return rows;
}

// The rows of each selected node's parent, so that the parents of the other nodes need not be read.
struct ParentRows {
    std::vector<int64_t> ids;
    Packed bits;
};

ParentRows read_parent_rows(const Graph& g, const std::vector<int64_t>& ids, const ParentTable& parents) {
    ParentRows out;
    for (int64_t p : parents.lookup(ids)) if (p >= 0) out.ids.push_back(p);
    std::sort(out.ids.begin(), out.ids.end());
    out.ids.erase(std::unique(out.ids.begin(), out.ids.end()), out.ids.end());
    const IdIndex index(g.ids);
    std::vector<int64_t> rows;
    for (int64_t p : out.ids) rows.push_back(index.find(p));
    out.bits = g.rows(rows);
    return out;
}

// Haplotypes visiting each node's parent; every haplotype when it has none.
Packed parent_visits(const std::vector<int64_t>& ids, const ParentRows& rows, const ParentTable& parents) {
    const std::vector<int64_t> parent = parents.lookup(ids);
    Packed visits(ids.size(), rows.bits.width, 0xFF);
    const IdIndex index(rows.ids);
    for (size_t i = 0; i < ids.size(); i++)
        if (parent[i] >= 0) std::memcpy(visits.row(i), rows.bits.row(size_t(index.find(parent[i]))), rows.bits.width);
    return visits;
}

// A node is inverted when some haplotype reads it forward and another (or the same) in reverse. Each
// walk was oriented to the reference when the graph was read. An inverted node gets one row per
// orientation; other nodes are orientation-free. Returns the orientation of each row.
std::vector<int8_t> split_inversions(GraphRows& rows, std::optional<Packed>& called, long long& inverted_nodes) {
    std::vector<int64_t> keep;
    std::vector<int8_t> orientation;
    for (size_t i = 0; i < rows.ids.size(); i++) {
        const auto any = [&](const Packed& p) {
            return std::any_of(p.row(i), p.row(i) + p.width, [](uint8_t x) { return x != 0; });
        };
        const bool inverted = any(rows.forward) && any(rows.reverse);
        keep.push_back(int64_t(i));
        orientation.push_back(inverted ? FORWARD : UNORIENTED);
        if (inverted) { keep.push_back(int64_t(i)); orientation.push_back(REVERSE); inverted_nodes++; }
    }
    Packed split = rows.bits.take(keep);
    for (size_t k = 0; k < keep.size(); k++) {
        if (orientation[k] == FORWARD) std::memcpy(split.row(k), rows.forward.row(keep[k]), split.width);
        if (orientation[k] == REVERSE) std::memcpy(split.row(k), rows.reverse.row(keep[k]), split.width);
    }
    if (called) called = called->take(keep);
    rows.ids = take(rows.ids, keep);
    rows.bits = std::move(split);
    return orientation;
}

// The columns to keep for --samples. Phased columns are haplotypes, so a sample names two.
std::vector<size_t> select_samples(const std::string& path, const std::string& matrix_type, std::vector<std::string>& columns) {
    const auto selected = split_ws(read_text(path));
    const std::set<std::string> unique(selected.begin(), selected.end());
    if (selected.empty() || unique.size() != selected.size()) throw Error("Sample list must contain unique sample IDs");
    std::vector<std::string> requested;
    if (matrix_type == "phased") {
        for (auto& s : selected) { requested.push_back(s + "#1"); requested.push_back(s + "#2"); }
    } else requested = selected;
    std::unordered_map<std::string, size_t> lookup;
    for (size_t i = 0; i < columns.size(); i++) lookup.emplace(columns[i], i);
    std::vector<std::string> missing;
    std::vector<size_t> pick;
    for (auto& name : requested) {
        const auto it = lookup.find(name);
        if (it == lookup.end()) missing.push_back(name);
        else pick.push_back(it->second);
    }
    if (!missing.empty()) {
        std::string list;
        for (size_t k = 0; k < std::min<size_t>(5, missing.size()); k++) list += (k ? ", '" : "'") + missing[k] + "'";
        throw Error("Requested samples are missing graph calls: [" + list + "]");
    }
    columns = requested;
    return pick;
}

Values take_columns(const Values& v, const std::vector<size_t>& pick) {
    Values out;
    out.rows = v.rows;
    out.cols = pick.size();
    out.v.resize(out.rows * out.cols);
    for (size_t r = 0; r < v.rows; r++)
        for (size_t c = 0; c < pick.size(); c++) out.row(r)[c] = v.row(r)[pick[c]];
    return out;
}

// The matrix rows, and what has to stay in step with them. With split inversions a node id names two rows.
struct MatrixRows {
    Values values;
    std::vector<int64_t> ids;
    std::vector<int8_t> orientation;
    std::optional<Values> compared;  // with drop complements: where complements are compared

    void keep(const std::vector<int64_t>& rows) {
        values = values.take(rows);
        ids = take(ids, rows);
        if (!orientation.empty()) orientation = take(orientation, rows);
        if (compared) compared = compared->take(rows);
    }
};

// Keeps one row per distinct pattern; returns each input row's pattern.
std::vector<int64_t> collapse_patterns(MatrixRows& rows) {
    std::unordered_map<std::string_view, int64_t> seen;
    std::vector<int64_t> first, pattern_of(rows.values.rows);
    for (size_t r = 0; r < rows.values.rows; r++) {
        const std::string_view key(reinterpret_cast<const char*>(rows.values.row(r)), rows.values.cols);
        const auto [it, inserted] = seen.emplace(key, int64_t(first.size()));
        if (inserted) first.push_back(int64_t(r));
        pattern_of[r] = it->second;
    }
    seen.clear();
    rows.keep(first);
    return pattern_of;
}

}  // namespace

NodeMatrix build_matrix(Graph g, const BuildOptions& o) {
    const std::string& matrix_type = o.matrix_type;
    Json model = Json::object();

    std::optional<ParentTable> parents;
    Json parent_metadata;
    if (o.parents) {
        parents = parent_table(*o.parents);
        parent_metadata = o.parents->metadata();
        parent_metadata["path"] = o.parents->path();
    }
    Json backbone_collapse;
    if (o.collapse_backbone) {
        g.load();  // the rule looks at every parent's and child's visits, not only the selected nodes'
        parents = collapse_backbone_parents(*parents, g.bits, g.names.size(), o.broken_min_haplotypes, backbone_collapse);
        stage("  collapse parents");
    }

    // Read the rows.
    GraphRows rows;
    const auto listed = listed_rows(g.ids, o.filters);
    rows.ids = listed ? take(g.ids, *listed) : g.ids;
    ParentRows parent_rows;
    if (o.condition_on_parents) parent_rows = read_parent_rows(g, rows.ids, *parents);
    if (listed) rows.bits = g.rows(*listed);
    else { g.load(); rows.bits = std::move(g.bits); }
    if (g.oriented) {
        rows.forward = listed ? g.forward.take(*listed) : std::move(g.forward);
        rows.reverse = listed ? g.reverse.take(*listed) : std::move(g.reverse);
    }
    g.bits = Packed();
    stage("  read rows");

    // Callability: unknown where a haplotype does not visit the node and is not known to skip it.
    std::optional<Packed> called;
    if (o.condition_on_parents) {
        called = parent_visits(rows.ids, parent_rows, *parents);
        for (size_t k = 0; k < called->data.size(); k++) called->data[k] |= rows.bits.data[k];
        model["absence_model"] = "parent_conditioned";
        if (o.collapse_backbone) model["backbone_collapse"] = backbone_collapse;
    } else {
        if (!o.assume_complete) called = rows.bits;
        model["absence_model"] = o.assume_complete ? "assume_complete" : "observed_only";
    }

    MatrixRows matrix;
    if (g.oriented) {
        long long inverted_nodes = 0;
        matrix.orientation = split_inversions(rows, called, inverted_nodes);
        Json orientation = Json::object();
        orientation["inverted_nodes"] = inverted_nodes;
        orientation["walks_flipped"] = g.flipped;
        orientation["walks"] = g.walks;
        orientation["walks_unanchored"] = g.unanchored;
        orientation["reference"] = o.reference;
        model["orientation"] = orientation;
    }

    // Matrix values.
    const Columns cols = columns_for(g.names, matrix_type);
    stage("  callability");
    matrix.values = call_values(rows.bits, called ? &*called : nullptr, cols, matrix_type, g.names.size());
    matrix.ids = std::move(rows.ids);
    called.reset();
    rows.bits = Packed();
    std::vector<std::string> columns = cols.names;
    if (o.drop_complements && o.condition_on_parents) {
        // Complements are compared only on the haplotypes visiting the node's parent.
        matrix.compared = call_values(parent_visits(matrix.ids, parent_rows, *parents), nullptr, cols, matrix_type, g.names.size());
        const int8_t full = matrix_type == "phased" ? 1 : 2;
        for (auto& x : matrix.compared->v) x = x == full;
    }
    if (!o.samples.empty()) {
        const std::vector<size_t> pick = select_samples(o.samples, matrix_type, columns);
        matrix.values = take_columns(matrix.values, pick);
        if (matrix.compared) matrix.compared = take_columns(*matrix.compared, pick);
    }

    // Filters, complements and patterns.
    for (const Filter& f : o.filters) {
        if (f.kind == "list" || f.kind == "all") continue;
        matrix.keep(indices_of(select(matrix.values, f, matrix_type)));
    }
    NodeMatrix m;
    if (o.drop_complements) {
        std::vector<int64_t> partner;
        const auto keep = complement_free(matrix.values, matrix_type == "phased" ? 1 : 2, matrix.compared ? &*matrix.compared : nullptr, partner);
        // (dropped id, kept complement id) per dropped row
        for (size_t r = 0; r < keep.size(); r++)
            if (!keep[r]) { m.complement_pairs.push_back(matrix.ids[r]); m.complement_pairs.push_back(matrix.ids[partner[r]]); }
        const auto kept = indices_of(keep);
        Json complements = Json::object();
        complements["dropped"] = (long long)(keep.size() - kept.size());
        model["complements"] = complements;
        matrix.keep(kept);
    }
    matrix.compared.reset();
    stage("  filters");

    m.input_node_ids = matrix.ids;
    m.input_node_orientation = matrix.orientation;
    m.oriented = g.oriented;
    if (o.collapse) {
        m.pattern_of = collapse_patterns(matrix);
        m.collapsed = true;
    }
    m.values = std::move(matrix.values);
    m.node_ids = std::move(matrix.ids);
    m.node_orientation = std::move(matrix.orientation);
    m.columns = std::move(columns);
    m.matrix_type = matrix_type;
    m.callability = std::move(model);
    if (o.coordinates) {
        m.coordinates = coordinates_of(*o.coordinates, m.node_ids);
        m.coordinate_chromosomes = o.coordinates->strings("chromosomes");
        m.coordinate_metadata = coordinate_provenance(*o.coordinates);
    }
    if (parents) {
        m.parent_node_ids = parents->lookup(m.input_node_ids);
        m.parent_metadata = std::move(parent_metadata);
    }
    return m;
}

void write_matrix(const NodeMatrix& m, const std::string& path) {
    StoreWriter a(path);
    a.add("values", m.values.rows, m.values.cols, m.values.v.data());
    a.add("node_ids", m.node_ids);
    a.add("columns", m.columns);
    a.add("input_node_ids", m.input_node_ids);
    if (m.collapsed) a.add("pattern_of", m.pattern_of);
    Json metadata = Json::object();
    metadata["matrix_type"] = m.matrix_type;
    metadata["callability"] = m.callability;
    if (!m.coordinates.empty()) {
        for (auto& [k, v] : m.coordinates) a.add("coordinate_" + k, v);
        a.add("coordinate_chromosomes", m.coordinate_chromosomes);
        metadata["coordinates"] = m.coordinate_metadata;
    }
    if (!m.parent_metadata.is_null()) {
        a.add("parent_node_ids", m.parent_node_ids);
        metadata["parents"] = m.parent_metadata;
    }
    if (m.oriented) {
        a.add("node_orientation", m.node_orientation);
        a.add("input_node_orientation", m.input_node_orientation);
    }
    if (m.callability.contains("complements"))
        a.add("complement_pairs", m.complement_pairs.size() / 2, 2, m.complement_pairs.data());
    a.set_metadata(metadata);
    a.publish();
}

}  // namespace glink
