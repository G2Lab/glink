// Reusable EdgeDepth-style node coordinates: computing them, and caching them in a store.
// Method: directional reference anchors, https://www.biorxiv.org/content/10.64898/2026.07.01.735728v1
//
// Positions are 1-based and closed, with 0 meaning missing; a chromosome index of -1 means missing.
// A node's point is its reference end, else the base after its upstream anchor, else the base
// before its downstream anchor. Unanchored or ambiguous nodes have position 0.
// Upstream is left of the forward node handle and downstream is right; reverse-complement links
// are kept. reference_source=auto prefers matching W/P walks, otherwise rank-0 SN/SO tags are
// assumed to describe the declared reference.
//
// EdgeDepth-style node anchors on an oriented GFA graph.
// Ordered floods compute extrema over ALL reachable first reference nodes,
// not nearest reference nodes by hop count. Each oriented vertex is visited
// once per flood, including graphs containing cycles.
#include <algorithm>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>

#include "coordinates.hpp"
#include "engine.hpp"
#include "gfa.hpp"

namespace glink {

namespace {

constexpr int FROM_UPSTREAM = 2, FROM_DOWNSTREAM = 3;  // position_source
std::string default_path(const std::string& gfa, const std::string& reference, const std::string& source) {
    std::string key = sha256_hex(coordinates_signature(gfa, reference, source).dump()).substr(0, 16);
    std::string safe;
    for (char c : reference) safe += (std::isalnum((unsigned char)c) || c == '-' || c == '_' || c == '.') ? c : '_';
    return ".glink-cache/node_coordinates/" + basename(gfa) + "." + safe + "." + key + ".h5";
}

}  // namespace

Json coordinates_signature(const std::string& gfa, const std::string& reference, const std::string& source) {
    Json j = Json::object();
    j["schema"] = 1;
    j["graph"] = file_sha256(gfa);
    j["reference"] = reference;
    j["requested_reference_source"] = source;
    j["engine_version"] = kCoordinatesVersion;
    return j;
}

Store load_coordinates(const std::string& path) {
    Store store(path);
    const Json* schema = json_find(store.metadata(), "schema");
    if (!schema || *schema != Json(1)) throw Error("Unsupported node coordinate schema in " + path);
    const std::vector<int64_t> node_id = store.ints("node_id");
    if (node_id.empty()) throw Error("Incomplete node coordinate arrays in " + path);
    for (const std::string& name : store.integer_arrays())
        if (store.shape(name).at(0) != node_id.size()) throw Error("Incomplete node coordinate arrays in " + path);
    for (size_t i = 1; i < node_id.size(); i++)
        if (node_id[i] <= node_id[i - 1]) throw Error("Node coordinate IDs must be sorted and unique in " + path);
    return store;
}

ColumnMap coordinates_of(const Store& coordinates, const std::vector<int64_t>& ids) {
    const std::vector<int64_t> node_id = coordinates.ints("node_id");
    const IdIndex index(node_id);
    std::vector<int64_t> rows(ids.size()), missing;
    for (size_t i = 0; i < ids.size(); i++) {
        rows[i] = index.find(ids[i]);
        if (rows[i] < 0 && missing.size() < 10) missing.push_back(ids[i]);
    }
    if (!missing.empty()) {
        std::string list;
        for (size_t k = 0; k < missing.size(); k++) list += (k ? ", " : "") + std::to_string(missing[k]);
        throw Error("Node coordinate map lacks IDs: [" + list + "]");
    }
    ColumnMap out;
    for (const std::string& name : coordinates.integer_arrays()) out[name] = coordinates.ints(name, rows);
    return out;
}

Json coordinate_provenance(const Store& coordinates) {
    Json j = coordinates.metadata();
    j["path"] = coordinates.path();
    return j;
}

namespace {

Coordinates build_coordinates(const Graph& graph, const std::string& reference, const std::string& source,
                              Json& summary) {
    if (!graph.coordinate_error.empty()) throw Error("Node coordinate preprocessing failed: " + graph.coordinate_error);
    Coordinates c;
    try { c = compute_coordinates(graph.coordinate_input, graph.links, reference, source); }
    catch (const std::exception& e) { throw Error(std::string("Node coordinate preprocessing failed: ") + e.what()); }
    summary = Json::object();
    summary["links"] = (long long)c.links;
    summary["reference_nodes"] = (long long)c.reference_nodes;
    summary["reference_source"] = c.reference_source;
    // Off-reference nodes sit on the reference base beside their anchor, not on the anchor.
    for (size_t i = 0; i < c.position.size(); i++) {
        if (c.position[i] > 0 && c.position_source[i] == FROM_UPSTREAM) c.position[i]++;
        if (c.position[i] > 1 && c.position_source[i] == FROM_DOWNSTREAM) c.position[i]--;
    }
    return c;
}

template <class T>
long long count_nonzero(const std::vector<T>& values) {
    long long n = 0;
    for (auto v : values) n += v != 0;
    return n;
}

}  // namespace

std::pair<Store, bool> prepare_coordinates(const Graph& graph, const std::string& gfa, std::string output,
                                           const std::string& reference, const std::string& source, bool rebuild) {
    if ((source != "auto" && source != "walk" && source != "rgfa") || reference.empty())
        throw Error("Choose a reference and reference_source auto/walk/rgfa");
    Json expected = coordinates_signature(gfa, reference, source);
    if (output.empty()) output = default_path(gfa, reference, source);
    if (resolve(output) == resolve(gfa)) throw Error("Node coordinate output must differ from the GFA");
    make_dirs(dirname(output));
    FileLock lock(output + ".lock");
    if (exists(output) && !rebuild) {
        Store coords = load_coordinates(resolve(output));
        if (!coords.matches(expected))
            throw Error("Node coordinate cache does not match graph/reference/method: " + output +
                        "; choose another output or explicitly rebuild it");
        return {std::move(coords), true};
    }
    Json summary;
    Coordinates c = build_coordinates(graph, reference, source, summary);
    Json metadata = expected;
    metadata["nodes"] = (long long)c.node_id.size();
    for (auto& [k, v] : summary.items()) metadata[k] = v;
    long long upstream_only = 0, downstream_only = 0;
    for (size_t i = 0; i < c.upstream.size(); i++) {
        upstream_only += c.upstream[i] > 0 && c.downstream[i] == 0;
        downstream_only += c.downstream[i] > 0 && c.upstream[i] == 0;
    }
    Json s = Json::object();
    s["nodes"] = (long long)c.node_id.size();
    s["placed"] = count_nonzero(c.position);
    s["ambiguous"] = count_nonzero(c.ambiguous);
    s["reversed_anchors"] = count_nonzero(c.reversed_anchors);
    s["upstream_only"] = upstream_only;
    s["downstream_only"] = downstream_only;
    metadata["summary"] = s;
    StoreWriter store(output);
    store.add("node_id", c.node_id);
    store.add("length", c.length);
    store.add("is_reference", c.is_reference);
    store.add("reference_chrom", c.reference_chrom);
    store.add("reference_start", c.reference_start);
    store.add("reference_end", c.reference_end);
    store.add("upstream", c.upstream);
    store.add("downstream", c.downstream);
    store.add("upstream_node", c.upstream_node);
    store.add("downstream_node", c.downstream_node);
    store.add("upstream_chrom", c.upstream_chrom);
    store.add("downstream_chrom", c.downstream_chrom);
    store.add("upstream_ambiguous", c.upstream_ambiguous);
    store.add("downstream_ambiguous", c.downstream_ambiguous);
    store.add("ambiguous", c.ambiguous);
    store.add("reversed_anchors", c.reversed_anchors);
    store.add("position", c.position);
    store.add("chrom", c.chrom);
    store.add("position_source", c.position_source);
    store.add("chromosomes", c.chromosomes);
    store.set_metadata(metadata);
    store.publish();
    return {load_coordinates(resolve(output)), false};
}

// ---- computing them -------------------------------------------------------------

namespace {

using I = int64_t;

[[noreturn]] void fail(const std::string& message) { throw std::invalid_argument(message); }

template <class T>
std::vector<T> reorder(const std::vector<T>& values, const std::vector<size_t>& order) {
    std::vector<T> out;
    out.reserve(order.size());
    for (size_t i : order) out.push_back(values[i]);
    return out;
}

}  // namespace

int64_t nonnegative_integer(const std::string& s) {
    size_t end = 0;
    I n;
    try { n = std::stoll(s, &end); } catch (const std::exception&) { fail("Expected nonnegative integer: " + s); }
    if (end != s.size() || n < 0) fail("Expected nonnegative integer: " + s);
    return n;
}

Coordinates compute_coordinates(const CoordinateInput& input, const std::vector<Link>& links,
                                const std::string& reference, const std::string& source) {
    if (source != "auto" && source != "walk" && source != "rgfa") fail("Invalid reference source");
    const std::vector<I>& ids = input.ids;
    const std::vector<I>& lengths = input.lengths;
    std::unordered_map<I, size_t> index;
    for (size_t i = 0; i < ids.size(); ++i) index.emplace(ids[i], i);
    auto lookup = [&](I id) {
        auto it = index.find(id);
        if (it == index.end()) fail("Walk/link references undefined node " + std::to_string(id));
        return it->second;
    };

    const size_t n = ids.size(), h = 2 * n;
    if (!n || n > size_t(std::numeric_limits<int32_t>::max()) / 2)
        fail("Empty graph or too many nodes for oriented index");
    const bool use_walk = source == "walk" || (source == "auto" && !input.walks.empty());
    if (use_walk && input.walks.empty()) fail("No selected reference W/P walk: " + reference);
    if (!use_walk && input.tags.empty()) fail("No selected reference walk or rank-0 SN/SO tags: " + reference);

    // Reference placements.
    std::vector<int32_t> rchrom(n, -1);
    std::vector<I> rstart(n, 0), rend(n, 0);
    std::vector<uint8_t> is_ref(n, 0), ambiguous(n, 0);
    std::vector<std::string> chromosomes;
    std::unordered_map<std::string, int32_t> chrom_index;
    auto add_ref = [&](size_t i, const std::string& chr, I start) {
        auto pair = chrom_index.emplace(chr, chromosomes.size());
        if (pair.second) chromosomes.push_back(chr);
        const int32_t c = pair.first->second;
        if (start < 0 || lengths[i] > std::numeric_limits<I>::max() - start) fail("Reference coordinate overflow");
        if (is_ref[i] && (rchrom[i] != c || rstart[i] != start + 1 || rend[i] != start + lengths[i]))
            ambiguous[i] = 1;
        is_ref[i] = 1; rchrom[i] = c; rstart[i] = start + 1; rend[i] = start + lengths[i];
    };
    if (use_walk) {
        for (const auto& walk : input.walks) {
            I pos = walk.start;
            size_t p = 0;
            while (p < walk.steps.size()) {
                size_t q;
                I id;
                if (walk.is_path) {
                    q = walk.steps.find(',', p);
                    std::string step = walk.steps.substr(p, q - p);
                    if (step.size() < 2 || (step.back() != '+' && step.back() != '-')) fail("Malformed P step");
                    id = nonnegative_integer(step.substr(0, step.size() - 1));
                } else {
                    if (walk.steps[p] != '>' && walk.steps[p] != '<') fail("Malformed W step");
                    q = walk.steps.find_first_of("><", p + 1);
                    id = nonnegative_integer(walk.steps.substr(p + 1, q - p - 1));
                }
                size_t i = lookup(id);
                add_ref(i, walk.chrom, pos);
                pos += lengths[i];
                if (q == std::string::npos) break;
                p = walk.is_path ? q + 1 : q;
            }
            if (walk.end >= 0 && pos != walk.end)
                fail("Reference walk span differs from sum of segment lengths: " + walk.chrom);
        }
    } else {
        for (const auto& tag : input.tags) add_ref(lookup(tag.id), tag.chrom, tag.start);
    }
    for (size_t i = 0; i < n; ++i) if (ambiguous[i]) { rchrom[i] = -1; rstart[i] = rend[i] = 0; }

    // Incoming CSR of bidirected handles: a->b also implies flip(b)->flip(a).
    std::vector<size_t> offset(h + 1, 0);
    for (const auto& l : links) {
        size_t a = 2 * lookup(l.a) + l.a_reverse, b = 2 * lookup(l.b) + l.b_reverse;
        ++offset[b + 1]; ++offset[(a ^ 1) + 1];
    }
    std::partial_sum(offset.begin(), offset.end(), offset.begin());
    std::vector<int32_t> incoming(offset.back());
    auto cursor = offset;
    for (const auto& l : links) {
        size_t a = 2 * lookup(l.a) + l.a_reverse, b = 2 * lookup(l.b) + l.b_reverse;
        incoming[cursor[b]++] = static_cast<int32_t>(a);
        incoming[cursor[a ^ 1]++] = static_cast<int32_t>(b ^ 1);
    }
    const size_t link_count = links.size();
    std::vector<size_t>().swap(cursor);

    // Flood from reference nodes in coordinate order (ascending for the
    // leftmost anchor, descending for the rightmost); first label wins.
    std::vector<size_t> seeds;
    for (size_t i = 0; i < n; ++i) if (is_ref[i] && !ambiguous[i]) seeds.push_back(i);
    std::sort(seeds.begin(), seeds.end(), [&](size_t a, size_t b) {
        return std::tie(rchrom[a], rstart[a], rend[a], ids[a]) < std::tie(rchrom[b], rstart[b], rend[b], ids[b]);
    });
    auto flood = [&](bool maximum) {
        std::vector<int32_t> labels(h, -1), stack;
        auto spread = [&](size_t seed) {
            for (size_t strand = 0; strand < 2; ++strand) {
                stack.push_back(static_cast<int32_t>(2 * seed + strand));
                while (!stack.empty()) {
                    size_t v = stack.back(); stack.pop_back();
                    for (size_t k = offset[v]; k < offset[v + 1]; ++k) {
                        size_t u = incoming[k];
                        if (is_ref[u / 2] || labels[u] >= 0) continue;
                        labels[u] = static_cast<int32_t>(seed);
                        stack.push_back(static_cast<int32_t>(u));
                    }
                }
            }
        };
        if (maximum) for (auto it = seeds.rbegin(); it != seeds.rend(); ++it) spread(*it);
        else for (size_t seed : seeds) spread(seed);
        return labels;
    };
    auto lo = flood(false), hi = flood(true);

    // A multiply placed reference node remains a stopping boundary; propagate
    // ambiguity rather than passing through it or assigning an arbitrary locus.
    std::vector<uint8_t> bad(h, 0);
    std::vector<int32_t> stack;
    for (size_t i = 0; i < n; ++i) if (is_ref[i] && ambiguous[i]) {
        stack.push_back(static_cast<int32_t>(2 * i));
        stack.push_back(static_cast<int32_t>(2 * i + 1));
    }
    while (!stack.empty()) {
        size_t v = stack.back(); stack.pop_back();
        for (size_t k = offset[v]; k < offset[v + 1]; ++k) {
            size_t u = incoming[k];
            if (is_ref[u / 2] || bad[u]) continue;
            bad[u] = 1;
            stack.push_back(static_cast<int32_t>(u));
        }
    }

    // Per-node anchors and plotting position.
    std::vector<I> upstream(n, 0), downstream(n, 0), upnode(n, -1), downnode(n, -1), position(n, 0);
    std::vector<int32_t> upchrom(n, -1), downchrom(n, -1), chrom(n, -1);
    std::vector<uint8_t> upbad(n, 0), downbad(n, 0), reversed(n, 0), position_source(n, 0);
    for (size_t i = 0; i < n; ++i) {
        if (is_ref[i]) {
            if (!ambiguous[i]) { position[i] = rend[i]; chrom[i] = rchrom[i]; position_source[i] = 1; }
            continue;
        }
        size_t uh = 2 * i + 1, dh = 2 * i;
        auto side_bad = [&](size_t s) {
            return bad[s] || (lo[s] >= 0 && hi[s] >= 0 && rchrom[lo[s]] != rchrom[hi[s]]);
        };
        upbad[i] = side_bad(uh);
        downbad[i] = side_bad(dh);
        if (!upbad[i] && lo[uh] >= 0) {
            size_t j = lo[uh]; upstream[i] = rend[j]; upchrom[i] = rchrom[j]; upnode[i] = ids[j];
        }
        if (!downbad[i] && hi[dh] >= 0) {
            size_t j = hi[dh]; downstream[i] = rstart[j]; downchrom[i] = rchrom[j]; downnode[i] = ids[j];
        }
        ambiguous[i] = upbad[i] || downbad[i] || (upchrom[i] >= 0 && downchrom[i] >= 0 && upchrom[i] != downchrom[i]);
        reversed[i] = upchrom[i] >= 0 && upchrom[i] == downchrom[i] && upstream[i] > downstream[i];
        if (ambiguous[i]) continue;
        if (upstream[i]) { position[i] = upstream[i]; chrom[i] = upchrom[i]; position_source[i] = 2; }
        else if (downstream[i]) { position[i] = downstream[i]; chrom[i] = downchrom[i]; position_source[i] = 3; }
    }

    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return ids[a] < ids[b]; });
    Coordinates out;
    out.node_id = reorder(ids, order);
    out.length = reorder(lengths, order);
    out.is_reference = reorder(is_ref, order);
    out.reference_chrom = reorder(rchrom, order);
    out.reference_start = reorder(rstart, order);
    out.reference_end = reorder(rend, order);
    out.upstream = reorder(upstream, order);
    out.downstream = reorder(downstream, order);
    out.upstream_node = reorder(upnode, order);
    out.downstream_node = reorder(downnode, order);
    out.upstream_chrom = reorder(upchrom, order);
    out.downstream_chrom = reorder(downchrom, order);
    out.upstream_ambiguous = reorder(upbad, order);
    out.downstream_ambiguous = reorder(downbad, order);
    out.ambiguous = reorder(ambiguous, order);
    out.reversed_anchors = reorder(reversed, order);
    out.position = reorder(position, order);
    out.chrom = reorder(chrom, order);
    out.position_source = reorder(position_source, order);
    out.chromosomes = std::move(chromosomes);
    out.links = int64_t(link_count);
    out.reference_nodes = std::count(is_ref.begin(), is_ref.end(), 1);
    out.reference_source = use_walk ? "reference_walk" : "rgfa_SO";
    return out;
}

}  // namespace glink
