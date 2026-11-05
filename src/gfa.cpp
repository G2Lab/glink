// GFA W/P haplotype paths -> bit-packed node x haplotype table.
//
// The only reader of the GFA. It goes through the file twice: first for everything but the bits,
// then for the haplotype paths once their columns are known.
#include <algorithm>
#include <cstdio>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "gfa.hpp"
#include "store.hpp"

namespace glink {

Packed Packed::take(const std::vector<int64_t>& keep) const {
    Packed out(keep.size(), width);
    for (size_t i = 0; i < keep.size(); i++) std::memcpy(out.row(i), row(keep[i]), width);
    return out;
}

Packed Packed::take(const std::vector<uint8_t>& mask) const {
    std::vector<int64_t> keep;
    for (size_t i = 0; i < mask.size(); i++) if (mask[i]) keep.push_back(int64_t(i));
    return take(keep);
}

IdIndex::IdIndex(const std::vector<int64_t>& ids) : ids_(&ids) {
    if (ids.empty()) return;
    max_ = ids.back();
    if (ids.front() >= 0 && max_ < (int64_t(1) << 31) && max_ <= int64_t(ids.size()) * 16 + 1024) {
        dense_.assign(size_t(max_) + 1, -1);
        for (size_t i = 0; i < ids.size(); i++) dense_[ids[i]] = int32_t(i);
    }
}

int64_t IdIndex::find(int64_t id) const {
    if (!dense_.empty()) return id < 0 || id > max_ ? -1 : dense_[id];
    auto it = std::lower_bound(ids_->begin(), ids_->end(), id);
    return it != ids_->end() && *it == id ? it - ids_->begin() : -1;
}

namespace {

// Tab-separated fields of a line, as views.
struct Fields {
    std::vector<std::pair<size_t, size_t>> spans;
    const std::string* line;

    Fields(const std::string& l, size_t limit = 0) : line(&l) {
        size_t start = 0;
        while (true) {
            if (limit && spans.size() + 1 == limit) { spans.emplace_back(start, l.size()); break; }
            size_t tab = l.find('\t', start);
            if (tab == std::string::npos) { spans.emplace_back(start, l.size()); break; }
            spans.emplace_back(start, tab);
            start = tab + 1;
        }
    }
    size_t size() const { return spans.size(); }
    std::string str(size_t i) const { return line->substr(spans[i].first, spans[i].second - spans[i].first); }
};

int64_t parse_id(const std::string& s, const std::string& path) {
    size_t end = 0;
    long long v;
    try { v = std::stoll(s, &end); } catch (const std::exception&) { throw Error(path + ": invalid node id '" + s + "'"); }
    while (end < s.size() && std::isspace((unsigned char)s[end])) end++;
    if (end != s.size()) throw Error(path + ": invalid node id '" + s + "'");
    return v;
}

// Haplotype name for a W/P line, or "" when it is not a path line. Sets `biological`.
std::string path_name(const std::string& line, bool& biological) {
    biological = false;
    if (starts_with(line, "W\t")) {
        Fields f(line, 5);
        if (f.size() < 3) throw Error("Malformed W line");
        std::string sample = f.str(1);
        biological = sample != "CHM13" && sample != "GRCh38";
        return sample + "#" + f.str(2);
    }
    if (starts_with(line, "P\t")) {
        Fields f(line, 3);
        std::string full = f.str(1);
        size_t a = full.find('#');
        std::string name = full;
        if (a != std::string::npos) {
            size_t b = full.find('#', a + 1);
            name = b == std::string::npos ? full : full.substr(0, b);
        }
        biological = !starts_with(name, "CHM13") && !starts_with(name, "GRCh38");
        return name;
    }
    return "";
}

// Calls visit(node_id, reverse) for each step of a W/P line; returns false when the
// callback rejects a step.
template <class Visit>
bool steps(const std::string& line, bool oriented, const std::string& path, Visit visit) {
    if (line[0] == 'W') {
        Fields f(line);
        if (f.size() < 7) throw Error(path + ": malformed W line");
        size_t i = f.spans[6].first, end = f.spans[6].second;
        while (i < end) {
            char c = line[i];
            if (c >= '0' && c <= '9') {
                bool reverse = false, valid = !oriented;
                if (i > f.spans[6].first) {
                    char before = line[i - 1];
                    reverse = before == '<';
                    if (before == '<' || before == '>') valid = true;
                }
                int64_t v = 0;
                while (i < end && line[i] >= '0' && line[i] <= '9') v = v * 10 + (line[i++] - '0');
                if (valid && !visit(v, reverse)) return false;
            } else i++;
        }
        return true;
    }
    Fields f(line);
    if (f.size() < 3) throw Error(path + ": malformed P line");
    size_t i = f.spans[2].first, end = f.spans[2].second;
    while (i <= end) {
        size_t comma = line.find(',', i);
        if (comma == std::string::npos || comma > end) comma = end;
        std::string step = line.substr(i, comma - i);
        if (step.empty()) throw Error(path + ": empty P line step");
        bool reverse = step.back() == '-';
        if (!visit(parse_id(step.substr(0, step.size() - 1), path), reverse)) return false;
        i = comma + 1;
    }
    return true;
}

void trim_carriage_return(std::string& line) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
}

inline void set_bit(Packed& p, size_t row, size_t col) { p.row(row)[col >> 3] |= uint8_t(1u << (col & 7)); }

// Whether a W/P line is the reference that walks are oriented against.
bool is_reference_line(const std::string& line, const std::string& reference) {
    if (line[0] == 'W') {
        Fields f(line, 5);
        return f.str(1) == reference || f.str(1) + "#" + f.str(2) == reference;
    }
    Fields f(line, 3);
    return f.str(1) == reference || starts_with(f.str(1), reference + "#");
}

// ---- what coordinates need from each line. These are stricter than the graph is: a GFA that
// fails here still gives a graph, and only a run that asks for coordinates reports the problem.

void coordinate_segment(const Fields& f, const std::string& reference, CoordinateInput& input,
                        std::unordered_set<int64_t>& seen) {
    if (f.size() < 3) throw Error("Malformed S record");
    const int64_t id = nonnegative_integer(f.str(1));
    const size_t sequence = f.spans[2].second - f.spans[2].first;
    int64_t length = sequence == 1 && (*f.line)[f.spans[2].first] == '*' ? -1 : int64_t(sequence);
    int64_t offset = -1, rank = -1;
    std::string chrom;
    for (size_t j = 3; j < f.size(); j++) {
        const std::string tag = f.str(j);
        if (starts_with(tag, "LN:i:")) {
            const int64_t ln = nonnegative_integer(tag.substr(5));
            if (length >= 0 && ln != length) throw Error("Sequence/LN length mismatch");
            length = ln;
        }
        if (starts_with(tag, "SO:i:")) offset = nonnegative_integer(tag.substr(5));
        if (starts_with(tag, "SR:i:")) rank = nonnegative_integer(tag.substr(5));
        if (starts_with(tag, "SN:Z:")) chrom = tag.substr(5);
    }
    if (length <= 0) throw Error("Segment needs a positive sequence length or LN tag");
    if (!seen.insert(id).second) throw Error("Duplicate node ID");
    input.ids.push_back(id);
    input.lengths.push_back(length);
    if (offset < 0 || rank != 0 || chrom.empty()) return;
    // PanSN-tagged sequences can name their assembly explicitly.
    const size_t first = chrom.find('#');
    if (first != std::string::npos) {
        if (chrom.substr(0, first) != reference && !starts_with(chrom, reference + "#")) return;
        chrom = chrom.substr(chrom.rfind('#') + 1);
    }
    input.tags.push_back({id, offset, chrom});
}

void coordinate_link(const Fields& f) {
    if (f.size() < 6 || (f.str(2) != "+" && f.str(2) != "-") || (f.str(4) != "+" && f.str(4) != "-"))
        throw Error("Malformed oriented L record");
    if (f.str(5) != "0M") throw Error("Coordinates require zero-overlap L records (0M)");
    nonnegative_integer(f.str(1));
    nonnegative_integer(f.str(3));
}

void coordinate_walk(const std::string& line, const std::string& reference, CoordinateInput& input) {
    if (line[0] == 'W') {
        if (Fields(line, 5).size() < 4) throw Error("Malformed W record");
        if (!is_reference_line(line, reference)) return;
        const Fields f(line);
        if (f.size() < 7) throw Error("Malformed reference W record");
        input.walks.push_back({f.str(3), f.str(6), nonnegative_integer(f.str(4)), nonnegative_integer(f.str(5)), false});
        return;
    }
    if (Fields(line, 3).size() < 3) throw Error("Malformed P record");
    if (!is_reference_line(line, reference)) return;
    const Fields f(line);
    if (f.size() < 4) throw Error("Malformed reference P record");
    const std::string name = f.str(1);
    std::string chrom = name.substr(name.rfind('#') == std::string::npos ? 0 : name.rfind('#') + 1);
    int64_t start = 0, end = -1;
    const size_t bracket = chrom.find('[');
    if (bracket != std::string::npos) {
        const size_t dash = chrom.find('-', bracket), close = chrom.find(']', dash);
        if (dash == std::string::npos || close != chrom.size() - 1) throw Error("Malformed P subrange");
        start = nonnegative_integer(chrom.substr(bracket + 1, dash - bracket - 1));
        end = nonnegative_integer(chrom.substr(dash + 1, close - dash - 1));
        chrom.resize(bracket);
    }
    const std::string overlaps = f.str(3);
    if (overlaps != "*") {
        // P overlaps must not change the sum of segment lengths.
        size_t p = 0;
        while (p < overlaps.size()) {
            const size_t q = overlaps.find(',', p);
            if (overlaps.substr(p, q - p) != "0M") throw Error("Reference P overlaps must be 0M or *");
            if (q == std::string::npos) break;
            p = q + 1;
        }
    }
    input.walks.push_back({chrom, f.str(2), start, end, true});
}

// Whether a W/P line is the backbone that roots the snarl tree.
bool is_backbone(const std::string& line, const std::string& backbone) {
    if (line[0] == 'W') {
        Fields f(line, 5);
        return f.size() > 2 && (backbone == f.str(1) || backbone == f.str(1) + "#" + f.str(2));
    }
    Fields f(line, 3);
    if (f.size() < 2) return false;
    const std::string name = f.str(1);
    const size_t a = name.find('#');
    std::string two = name;
    if (a != std::string::npos) {
        const size_t b = name.find('#', a + 1);
        if (b != std::string::npos) two = name.substr(0, b);
    }
    return backbone == name || backbone == name.substr(0, a) || backbone == two;
}

}  // namespace

Graph read_gfa(const std::string& path, const GraphOptions& options) {
    const bool oriented = options.oriented;
    const std::string& reference = options.reference;
    LineReader reader(path);
    std::string line;
    Graph g;
    g.oriented = oriented;

    // First pass: everything but the bits. A GFA coordinates cannot use is an error only if coordinates are asked for.
    std::vector<int64_t> raw_ids;
    std::set<std::string> nameset;
    std::set<int64_t> tips;
    std::vector<int64_t> reference_nodes;
    std::vector<uint8_t> reference_reverse;
    std::unordered_set<int64_t> segments;
    auto for_coordinates = [&](auto&& take) {
        if (!g.coordinate_error.empty()) return;
        try { take(); } catch (const std::exception& e) { g.coordinate_error = e.what(); }
    };
    while (reader.next(line)) {
        trim_carriage_return(line);
        if (starts_with(line, "S\t")) {
            const Fields f(line);
            raw_ids.push_back(parse_id(f.str(1), path));
            for_coordinates([&] { coordinate_segment(f, reference, g.coordinate_input, segments); });
            continue;
        }
        if (starts_with(line, "L\t")) {
            Fields f(line);
            if (f.size() < 5) throw Error(path + ": malformed L line");
            g.links.push_back({parse_id(f.str(1), path), parse_id(f.str(3), path), f.str(2) == "-", f.str(4) == "-"});
            for_coordinates([&] { coordinate_link(f); });
            continue;
        }
        if (line.size() < 2 || (line[0] != 'W' && line[0] != 'P') || line[1] != '\t') continue;
        bool biological;
        const std::string name = path_name(line, biological);
        for_coordinates([&] { coordinate_walk(line, reference, g.coordinate_input); });
        if (oriented && is_reference_line(line, reference))
            steps(line, true, path, [&](int64_t id, bool back) {
                reference_nodes.push_back(id);
                reference_reverse.push_back(back);
                return true;
            });
        if (is_backbone(line, options.backbone)) {
            int64_t first = -1, last = -1;
            steps(line, false, path, [&](int64_t id, bool) {
                if (first < 0) first = id;
                last = id;
                return true;
            });
            if (first >= 0) { tips.insert(first); tips.insert(last); }
        }
        if (!name.empty() && biological) nameset.insert(name);
    }
    if (raw_ids.empty()) throw Error(path + ": no S lines, so no nodes");
    if (nameset.empty()) throw Error(path + ": no W or P lines for biological haplotypes");
    g.backbone_tips.assign(tips.begin(), tips.end());
    if (!g.coordinate_error.empty()) g.coordinate_input = {};

    g.names.assign(nameset.begin(), nameset.end());
    std::unordered_map<std::string, size_t> lookup;
    for (size_t i = 0; i < g.names.size(); i++) lookup.emplace(g.names[i], i);
    g.ids = std::move(raw_ids);
    std::sort(g.ids.begin(), g.ids.end());
    g.ids.erase(std::unique(g.ids.begin(), g.ids.end()), g.ids.end());
    IdIndex full(g.ids);
    const size_t width = (g.names.size() + 7) / 8;
    g.bits = Packed(g.ids.size(), width);
    std::vector<int8_t> reference_orientation;  // +1 read forward by the reference, -1 reversed, 0 off it or both
    if (oriented) {
        g.forward = Packed(g.ids.size(), width);
        g.reverse = Packed(g.ids.size(), width);
        if (reference_nodes.empty()) throw Error(path + ": no " + reference + " W/P walk to orient haplotypes against");
        std::vector<uint8_t> on_forward(g.ids.size()), on_reverse(g.ids.size());
        for (size_t k = 0; k < reference_nodes.size(); k++) {
            int64_t r = full.find(reference_nodes[k]);
            if (r < 0) throw Error(path + ": path references an undefined node");
            (reference_reverse[k] ? on_reverse : on_forward)[r] = 1;
        }
        reference_orientation.resize(g.ids.size());
        for (size_t r = 0; r < g.ids.size(); r++) reference_orientation[r] = int8_t(on_forward[r]) - int8_t(on_reverse[r]);
    }

    // Second pass: each haplotype path sets its bits.
    std::vector<int64_t> rows;
    std::vector<uint8_t> backwards;
    reader.rewind();
    while (reader.next(line)) {
        trim_carriage_return(line);
        if (!starts_with(line, "W\t") && !starts_with(line, "P\t")) continue;
        bool biological;
        std::string name = path_name(line, biological);
        auto it = lookup.find(name);
        if (it == lookup.end()) continue;
        const size_t column = it->second;
        rows.clear();
        backwards.clear();
        bool ok = steps(line, oriented, path, [&](int64_t id, bool back) {
            int64_t r = full.find(id);
            if (r < 0) return false;
            rows.push_back(r);
            backwards.push_back(back);
            return true;
        });
        if (!ok) throw Error(path + ": path references an undefined node");
        bool flip = false;
        if (oriented) {
            // A walk written on the opposite strand of the reference is read reversed.
            long long agree = 0, disagree = 0;
            for (size_t k = 0; k < rows.size(); k++) {
                const int vote = reference_orientation[rows[k]] * (backwards[k] ? -1 : 1);
                agree += vote > 0;
                disagree += vote < 0;
            }
            g.walks++;
            if (!agree && !disagree) g.unanchored++;
            else if (disagree > agree) { flip = true; g.flipped++; }
        }
        for (size_t k = 0; k < rows.size(); k++) {
            set_bit(g.bits, rows[k], column);
            if (oriented) set_bit(bool(backwards[k]) != flip ? g.reverse : g.forward, rows[k], column);
        }
    }
    return g;
}

void check_graph(const std::string& path, const std::string& gfa, const GraphOptions& options) {
    Json expected = Json::object();
    expected["graph"] = file_sha256(gfa);
    expected["reference"] = options.reference;
    expected["backbone"] = options.backbone;
    if (options.oriented) expected["oriented"] = true;
    if (!Store(path).matches(expected))
        throw Error("Graph store does not match this GFA/reference/backbone: " + path + "; pass --rebuild to replace it");
}

Graph load_graph(const std::string& path, bool oriented) {
    const Store store(path);
    const Json& metadata = store.metadata();
    Graph g;
    g.ids = store.ints("node_id");
    g.names = store.strings("names");
    const auto shape = store.shape("bits");
    const size_t width = (g.names.size() + 7) / 8;
    if (store.item_size("bits") != 1 || shape.size() != 2 || shape[0] != g.ids.size() || shape[1] != width)
        throw Error(path + ": bits do not match " + thousands(g.ids.size()) + " nodes x " +
                    std::to_string(g.names.size()) + " haplotypes");
    g.bits.rows = shape[0];
    g.bits.width = width;
    g.bits_path = path;
    g.on_disk = true;

    const auto a = store.ints("link_a"), b = store.ints("link_b");
    const auto a_reverse = store.ints("link_a_reverse"), b_reverse = store.ints("link_b_reverse");
    if (a.size() != b.size() || a.size() != a_reverse.size() || a.size() != b_reverse.size()) throw Error(path + ": links do not match");
    g.links.reserve(a.size());
    for (size_t i = 0; i < a.size(); i++) g.links.push_back({a[i], b[i], a_reverse[i] != 0, b_reverse[i] != 0});
    g.backbone_tips = store.ints("backbone_tips");
    g.coordinate_error = metadata.at("coordinate_error").get<std::string>();

    CoordinateInput& input = g.coordinate_input;
    input.ids = store.ints("coordinate_node_id");
    input.lengths = store.ints("coordinate_node_length");
    const auto walk_chrom = store.strings("walk_chrom");
    const auto walk_steps_length = store.ints("walk_steps_length");
    const auto walk_start = store.ints("walk_start"), walk_end = store.ints("walk_end"), walk_is_path = store.ints("walk_is_path");
    std::vector<uint8_t> walk_steps(store.shape("walk_steps").at(0));
    if (!walk_steps.empty()) store.read_bytes("walk_steps", walk_steps.data());
    size_t walk_offset = 0;
    for (size_t i = 0; i < walk_chrom.size(); i++) {
        const size_t length = size_t(walk_steps_length.at(i));
        if (walk_offset + length > walk_steps.size()) throw Error(path + ": walks do not match");
        input.walks.push_back({walk_chrom[i], std::string(walk_steps.begin() + walk_offset, walk_steps.begin() + walk_offset + length),
                               walk_start.at(i), walk_end.at(i), walk_is_path.at(i) != 0});
        walk_offset += length;
    }
    const auto tag_chrom = store.strings("tag_chrom");
    const auto tag_id = store.ints("tag_id"), tag_start = store.ints("tag_start");
    for (size_t i = 0; i < tag_chrom.size(); i++) input.tags.push_back({tag_id.at(i), tag_start.at(i), tag_chrom[i]});

    if (oriented) {
        g.oriented = true;
        g.forward = Packed(g.bits.rows, g.bits.width);
        g.reverse = Packed(g.bits.rows, g.bits.width);
        store.read_bytes("forward", g.forward.data.data());
        store.read_bytes("reverse", g.reverse.data.data());
        g.walks = metadata.at("walks").get<long long>();
        g.flipped = metadata.at("flipped").get<long long>();
        g.unanchored = metadata.at("unanchored").get<long long>();
    }
    return g;
}

void save_graph(const Graph& g, const std::string& path, const std::string& identity, const GraphOptions& options) {
    StoreWriter store(path);
    store.add("node_id", g.ids);
    store.add("bits", g.bits.rows, g.bits.width, g.bits.data.data());
    store.add("names", g.names);

    std::vector<int64_t> a, b;
    std::vector<uint8_t> a_reverse, b_reverse;
    for (const Link& l : g.links) {
        a.push_back(l.a); b.push_back(l.b);
        a_reverse.push_back(l.a_reverse); b_reverse.push_back(l.b_reverse);
    }
    store.add("link_a", a);
    store.add("link_b", b);
    store.add("link_a_reverse", a_reverse);
    store.add("link_b_reverse", b_reverse);
    store.add("backbone_tips", g.backbone_tips);

    const CoordinateInput& input = g.coordinate_input;
    store.add("coordinate_node_id", input.ids);
    store.add("coordinate_node_length", input.lengths);
    std::vector<std::string> walk_chrom, tag_chrom;
    std::vector<int64_t> walk_steps_length, walk_start, walk_end, tag_id, tag_start;
    std::vector<uint8_t> walk_steps, walk_is_path;  // steps as bytes: HDF5 cannot hold one string over about 32 MB
    for (const auto& w : input.walks) {
        walk_chrom.push_back(w.chrom);
        walk_steps.insert(walk_steps.end(), w.steps.begin(), w.steps.end());
        walk_steps_length.push_back(int64_t(w.steps.size()));
        walk_start.push_back(w.start); walk_end.push_back(w.end); walk_is_path.push_back(w.is_path);
    }
    for (const auto& t : input.tags) { tag_id.push_back(t.id); tag_start.push_back(t.start); tag_chrom.push_back(t.chrom); }
    store.add("walk_chrom", walk_chrom);
    store.add("walk_steps", walk_steps);
    store.add("walk_steps_length", walk_steps_length);
    store.add("walk_start", walk_start);
    store.add("walk_end", walk_end);
    store.add("walk_is_path", walk_is_path);
    store.add("tag_id", tag_id);
    store.add("tag_start", tag_start);
    store.add("tag_chrom", tag_chrom);

    Json metadata = Json::object();
    metadata["graph"] = identity;
    metadata["reference"] = options.reference;
    metadata["backbone"] = options.backbone;
    metadata["oriented"] = g.oriented;
    metadata["coordinate_error"] = g.coordinate_error;
    metadata["representation"] = "orientation-independent node presence, little-endian packed haplotype bits";
    if (g.oriented) {
        store.add("forward", g.forward.rows, g.forward.width, g.forward.data.data());
        store.add("reverse", g.reverse.rows, g.reverse.width, g.reverse.data.data());
        metadata["walks"] = g.walks;
        metadata["flipped"] = g.flipped;
        metadata["unanchored"] = g.unanchored;
    }
    store.set_metadata(metadata);
    store.publish();
}

std::pair<Graph, bool> prepare_graph(const std::string& gfa, const GraphOptions& options, std::string output, bool rebuild) {
    const std::string& identity = file_sha256(gfa);
    if (output.empty()) {
        Json signature = Json::object();
        signature["graph"] = identity;
        signature["oriented"] = options.oriented;
        signature["reference"] = options.reference;
        signature["backbone"] = options.backbone;
        output = ".glink-cache/graph/" + basename(gfa) + "." + sha256_hex(signature.dump()).substr(0, 16) + ".h5";
    }
    if (resolve(output) == resolve(gfa)) throw Error("Graph store must differ from the GFA");
    make_dirs(dirname(output));
    FileLock lock(output + ".lock");
    const bool reused = exists(output) && !rebuild;
    if (reused) check_graph(output, gfa, options);
    else {
        Graph g = read_gfa(gfa, options);
        save_graph(g, output, identity, options);
    }
    return {load_graph(output, options.oriented), reused};  // a new graph too: its bits need not stay in memory
}

void Graph::release_step_data() {
    links = {};
    backbone_tips = {};
    coordinate_input = {};
}

Packed Graph::rows(const std::vector<int64_t>& keep) const {
    if (!on_disk) return bits.take(keep);
    Packed out(keep.size(), bits.width);
    Store(bits_path).read_rows("bits", keep, out.data.data());
    return out;
}

void Graph::load() {
    if (!on_disk) return;
    Packed all(bits.rows, bits.width);
    Store(bits_path).read_bytes("bits", all.data.data());
    bits = std::move(all);
    on_disk = false;
}

}  // namespace glink
