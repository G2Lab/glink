// The conditioning tree, in one store. vg computes the snarls (`vg snarls` piped through `vg view -R`),
// and this finds each node's innermost snarl. A node's parent is that snarl's start node.
// Method: vg snarls rooted on backbone tips, with node membership by innermost-snarl flood.
#include <algorithm>
#include <cstdio>
#include <optional>
#include <sys/wait.h>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include "conditioning.hpp"
#include "gfa.hpp"

namespace glink {

namespace {

constexpr int SCHEMA = 1;

std::string shell_quote(const std::string& s) {
    std::string quoted = "'";
    for (char c : s) quoted += c == '\'' ? std::string("'\\''") : std::string(1, c);
    return quoted + "'";
}

// Runs a pipeline of commands in bash, so that a failure anywhere in it counts, and returns the
// last command's stdout. Throws on failure.
std::string run_pipeline(const std::vector<std::vector<std::string>>& commands) {
    std::string pipeline;
    for (size_t k = 0; k < commands.size(); k++) {
        if (k) pipeline += "| ";
        for (const std::string& word : commands[k]) pipeline += shell_quote(word) + " ";
    }
    const std::string errors_path = temporary_file(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp", "glink-vg-", ".err");
    const std::string shell = "bash -o pipefail -c " + shell_quote(pipeline) + " 2> " + shell_quote(errors_path);
    FILE* pipe = ::popen(shell.c_str(), "r");
    if (!pipe) throw Error("Cannot start bash to run vg");
    std::string out;
    char buffer[1 << 16];
    size_t n;
    while ((n = std::fread(buffer, 1, sizeof buffer, pipe)) > 0) out.append(buffer, n);
    const int status = ::pclose(pipe);
    const std::string errors = read_text(errors_path);
    std::remove(errors_path.c_str());
    if (status != 0) {
        const std::string tail = errors.size() > 2000 ? errors.substr(errors.size() - 2000) : errors;
        throw Error("vg failed (" + std::to_string(WIFEXITED(status) ? WEXITSTATUS(status) : 128) + "): " + tail);
    }
    return out;
}

std::string vg_version(const std::string& vg) {
    std::string out = run_pipeline({{vg, "version"}});
    size_t nl = out.find('\n');
    std::string first = out.substr(0, nl);
    while (!first.empty() && std::isspace((unsigned char)first.back())) first.pop_back();
    while (!first.empty() && std::isspace((unsigned char)first.front())) first.erase(first.begin());
    return out.empty() ? "unknown" : first;
}

std::vector<Json> json_lines(const std::string& text, const std::string& source) {
    std::vector<Json> records;
    size_t start = 0, number = 0;
    while (start < text.size()) {
        size_t end = text.find('\n', start);
        if (end == std::string::npos) end = text.size();
        std::string line = text.substr(start, end - start);
        start = end + 1;
        number++;
        if (split_ws(line).empty()) continue;
        Json r;
        try { r = Json::parse(line); }
        catch (const Json::parse_error&) { throw Error(source + ":" + std::to_string(number) + ": expected vg snarl JSON, one record per line"); }
        if (!r.is_object() || !r.contains("start") || !r.contains("end"))
            throw Error(source + ":" + std::to_string(number) + ": expected a snarl with start and end visits");
        records.push_back(std::move(r));
    }
    return records;
}

// ---- snarl tree ------------------------------------------------------------------

// Links as an adjacency list over node sides: 2*i is the left of node index i, 2*i+1 the right.
struct Adjacency {
    std::vector<int64_t> offsets, targets;
};

Adjacency side_adjacency(const Graph& graph, const std::string& gfa) {
    const IdIndex index(graph.ids);
    std::vector<std::pair<int64_t, int64_t>> edges;
    edges.reserve(2 * graph.links.size());
    for (const Link& link : graph.links) {
        const int64_t a = index.find(link.a), b = index.find(link.b);
        if (a < 0 || b < 0) throw Error(gfa + ": link references an undefined node");
        const int64_t side_a = 2 * a + (link.a_reverse ? 0 : 1), side_b = 2 * b + (link.b_reverse ? 1 : 0);
        edges.emplace_back(side_a, side_b);
        edges.emplace_back(side_b, side_a);
    }
    Adjacency adjacency;
    adjacency.offsets.assign(2 * graph.ids.size() + 1, 0);
    for (const auto& edge : edges) adjacency.offsets[edge.first + 1]++;
    for (size_t k = 1; k < adjacency.offsets.size(); k++) adjacency.offsets[k] += adjacency.offsets[k - 1];
    adjacency.targets.resize(edges.size());
    std::vector<int64_t> fill(adjacency.offsets.begin(), adjacency.offsets.end() - 1);
    for (const auto& edge : edges) adjacency.targets[fill[edge.first]++] = edge.second;
    return adjacency;
}

struct SnarlTree {
    std::vector<int64_t> start, end, parent;
    std::vector<uint8_t> start_back, end_back;  // only the flood needs these
};

std::pair<int64_t, bool> visit(const Json& v) {
    const Json* id = json_find(v, "name");
    if (!id) id = json_find(v, "node_id");
    if (!id) throw Error("vg snarl visit has no node");
    int64_t node = id->is_string() ? std::stoll(id->get<std::string>()) : id->get<int64_t>();  // vg writes ids as strings
    const Json* back = json_find(v, "backward");
    return {node, back && back->is_boolean() && back->get<bool>()};
}

// The same key whichever way round a snarl is read.
using Key = std::tuple<int64_t, bool, int64_t, bool>;

Key snarl_key(int64_t s, bool sb, int64_t e, bool eb) {
    return std::min(Key{s, sb, e, eb}, Key{e, !eb, s, !sb});
}

SnarlTree build_snarl_tree(const std::vector<Json>& records, const std::vector<int64_t>& ids) {
    const size_t n = records.size();
    SnarlTree t;
    t.start.resize(n); t.end.resize(n); t.start_back.resize(n); t.end_back.resize(n);
    t.parent.assign(n, -1);
    const IdIndex nodes(ids);
    std::map<Key, int64_t> index;
    std::vector<std::optional<Key>> parent_keys(n);
    for (size_t i = 0; i < n; i++) {
        auto [s, sb] = visit(records[i].at("start"));
        auto [e, eb] = visit(records[i].at("end"));
        if (nodes.find(s) < 0) throw Error("vg snarl start node is not a GFA segment; was the graph chopped?");
        if (nodes.find(e) < 0) throw Error("vg snarl end node is not a GFA segment; was the graph chopped?");
        t.start[i] = s; t.start_back[i] = sb; t.end[i] = e; t.end_back[i] = eb;
        if (!index.emplace(snarl_key(s, sb, e, eb), int64_t(i)).second)
            throw Error("vg reported the snarl from node " + std::to_string(s) + " to " + std::to_string(e) + " twice");
        const Json* p = json_find(records[i], "parent");
        if (p && p->is_object()) {
            auto [ps, psb] = visit(p->at("start"));
            auto [pe, peb] = visit(p->at("end"));
            parent_keys[i] = snarl_key(ps, psb, pe, peb);
        }
    }
    for (size_t i = 0; i < n; i++) {
        if (!parent_keys[i]) continue;
        const auto it = index.find(*parent_keys[i]);
        if (it == index.end()) throw Error("Snarl parent of snarl " + std::to_string(i) + " is not among the reported snarls");
        t.parent[i] = it->second;
    }
    return t;
}

// vg does not list a snarl's nodes, so each node's innermost snarl comes from flooding the graph
// from the snarl's boundaries, jumping over its child snarls.
std::vector<int64_t> assign_nodes(const SnarlTree& t, const std::vector<int64_t>& ids, const Adjacency& adjacency) {
    const IdIndex idx(ids);
    const size_t n = t.start.size();
    std::vector<int64_t> s(n), e(n), s_inner(n), e_inner(n);
    for (size_t i = 0; i < n; i++) {
        s[i] = idx.find(t.start[i]);
        e[i] = idx.find(t.end[i]);
        s_inner[i] = 2 * s[i] + (t.start_back[i] ? 0 : 1);
        e_inner[i] = 2 * e[i] + (t.end_back[i] ? 1 : 0);
    }
    std::vector<std::vector<int64_t>> children(n);
    for (size_t c = 0; c < n; c++) if (t.parent[c] >= 0) children[t.parent[c]].push_back(int64_t(c));
    std::vector<int64_t> node_snarl(ids.size(), -1);
    for (size_t snarl = 0; snarl < n; snarl++) {
        std::unordered_set<int64_t> walls;
        std::unordered_map<int64_t, std::vector<int64_t>> jumps;
        for (int64_t child : children[snarl]) {
            walls.insert(s_inner[child]);
            walls.insert(e_inner[child]);
            jumps[s[child]].push_back(e[child]);
            jumps[e[child]].push_back(s[child]);
        }
        const int64_t b0 = s[snarl], b1 = e[snarl];
        std::vector<int64_t> sides = {s_inner[snarl], e_inner[snarl]};
        while (!sides.empty()) {
            int64_t side = sides.back();
            sides.pop_back();
            for (int64_t k = adjacency.offsets[side]; k < adjacency.offsets[side + 1]; k++) {
                std::vector<int64_t> pending = {adjacency.targets[k] / 2};
                while (!pending.empty()) {
                    int64_t node = pending.back();
                    pending.pop_back();
                    if (node == b0 || node == b1 || node_snarl[node] == int64_t(snarl)) continue;
                    if (node_snarl[node] >= 0)
                        throw Error("Node " + std::to_string(ids[node]) + " lies in two sibling snarls; "
                                    "the vg snarl tree is inconsistent with the GFA");
                    node_snarl[node] = int64_t(snarl);
                    for (int64_t x : {2 * node, 2 * node + 1}) if (!walls.count(x)) sides.push_back(x);
                    auto it = jumps.find(node);
                    if (it != jumps.end()) for (int64_t j : it->second) pending.push_back(j);
                }
            }
        }
    }
    for (size_t i = 0; i < n; i++)
        if (node_snarl[s[i]] != t.parent[i] || node_snarl[e[i]] != t.parent[i])
            throw Error("Snarl boundary nodes must belong to the parent snarl");
    return node_snarl;
}

// ---- running vg --------------------------------------------------------------------

std::vector<std::string> snarl_command(const std::string& vg, bool named, const std::string& backbone, int threads,
                                       const std::vector<int64_t>& tips, const std::string& graph) {
    std::vector<std::string> c = {vg, "snarls"};
    if (named) c.push_back("-n");
    for (auto& x : {std::string("-P"), backbone, std::string("-t"), std::to_string(threads)}) c.push_back(x);
    for (int64_t tip : tips) { c.push_back("-w"); c.push_back(std::to_string(tip)); }
    c.push_back(graph);
    return c;
}

Json command_json(const std::vector<std::vector<std::string>>& commands) {
    Json list = Json::array();
    for (auto& c : commands) {
        Json one = Json::array();
        for (auto& a : c) one.push_back(Json(a));
        list.push_back(one);
    }
    return list;
}

SnarlRun read_saved_snarls(const TreeOptions& o) {
    std::string format = o.snarls_format;
    if (format == "auto") {
        LineReader probe(o.snarls);
        std::string line, prefix;
        while (prefix.size() < 4096 && probe.next(line)) prefix += line + "\n";
        size_t k = 0;
        while (k < prefix.size() && std::isspace((unsigned char)prefix[k])) k++;
        format = k == prefix.size() || prefix[k] == '{' ? "json" : "binary";
    }
    SnarlRun run;
    if (format == "json") {
        LineReader reader(o.snarls);
        std::string line, text;
        while (reader.next(line)) text += line + "\n";
        run.records = json_lines(text, o.snarls);
    } else {
        const std::vector<std::string> view = {o.vg, "view", "-R", o.snarls};
        run.records = json_lines(run_pipeline({view}), o.snarls);
        run.commands = command_json({view});
        run.vg_version = vg_version(o.vg);
    }
    return run;
}

// ---- the store ---------------------------------------------------------------------

// What a tree store built from these inputs records, and what reusing one has to match.
Json tree_signature(const std::string& gfa, const TreeOptions& o) {
    Json j = Json::object();
    j["schema"] = SCHEMA;
    j["graph"] = file_sha256(gfa);
    j["backbone"] = o.backbone;
    j["vg_graph"] = o.vg_graph.empty() ? Json() : Json(file_sha256(o.vg_graph));
    j["snarls"] = o.snarls.empty() ? Json() : Json(file_sha256(o.snarls));
    return j;
}

std::string tree_path(const std::string& gfa, const TreeOptions& o) {
    if (!o.snarl_tree.empty()) return o.snarl_tree;
    if (!o.output.empty()) return o.output;
    const std::string key = sha256_hex(tree_signature(gfa, o).dump()).substr(0, 16);
    return ".glink-cache/tree/" + basename(gfa) + "." + key + ".h5";
}

// Opens a tree store, throwing unless it belongs to these inputs.
Store open_tree(const std::string& path, const std::string& gfa, const TreeOptions& o) {
    Store tree(resolve(path));
    const auto snarls = tree.shape("start_node");
    if (tree.shape("end_node") != snarls || tree.shape("parent") != snarls || tree.shape("node_id") != tree.shape("node_snarl"))
        throw Error("Incompatible conditioning tree: " + path);
    // A tree passed in by name only has to be of this graph; how it was made is its own business.
    const Json expected = o.snarl_tree.empty() ? tree_signature(gfa, o) : Json{{"schema", SCHEMA}, {"graph", file_sha256(gfa)}};
    if (!tree.matches(expected))
        throw Error("Conditioning tree does not match this GFA/backbone/vg graph/snarls: " + path +
                    "; pass --rebuild, or --rebuild-conditioning-tree without an index, to replace it");
    return tree;
}

}  // namespace

// ---- public ------------------------------------------------------------------------

std::vector<int64_t> ParentTable::lookup(const std::vector<int64_t>& ids) const {
    const IdIndex index(node_id);
    std::vector<int64_t> out(ids.size());
    for (size_t i = 0; i < ids.size(); i++) {
        const int64_t row = index.find(ids[i]);
        if (row < 0) throw Error("Node IDs are not in the conditioning tree");
        out[i] = parent_node_id[size_t(row)];
    }
    return out;
}

ParentTable parent_table(const Store& tree) {
    const std::vector<int64_t> starts = tree.ints("start_node"), node_snarl = tree.ints("node_snarl");
    ParentTable table{tree.ints("node_id"), std::vector<int64_t>(node_snarl.size(), -1)};
    for (size_t i = 0; i < node_snarl.size(); i++)
        if (node_snarl[i] >= 0) table.parent_node_id[i] = starts.at(size_t(node_snarl[i]));
    return table;
}

bool needs_snarls(const std::string& gfa, const TreeOptions& o) {
    const std::string path = tree_path(gfa, o);
    if (o.snarl_tree.empty() && (o.rebuild || !exists(path))) return true;
    open_tree(path, gfa, o);
    return false;
}

SnarlRun run_snarls(const Graph& graph, const std::string& gfa, const TreeOptions& o) {
    if (!o.snarls.empty()) return read_saved_snarls(o);
    if (graph.backbone_tips.empty())
        throw Error("No " + o.backbone + " backbone walk in " + gfa + "; the snarl tree must be rooted on it");
    if (!o.vg_graph.empty()) {
        const auto words = split_ws(run_pipeline({{o.vg, "stats", "-N", o.vg_graph}}));
        if (words.empty() || std::stoll(words[0]) != int64_t(graph.ids.size()))
            throw Error(o.vg_graph + " has " + (words.empty() ? "?" : words[0]) + " nodes but " + gfa + " has " +
                        std::to_string(graph.ids.size()) + "; the GFA must be its --no-translation export");
    }
    const auto snarls = snarl_command(o.vg, o.vg_graph.empty(), o.backbone, o.threads, graph.backbone_tips,
                                      o.vg_graph.empty() ? gfa : o.vg_graph);
    const std::vector<std::string> view = {o.vg, "view", "-R", "-"};
    SnarlRun run;
    run.records = json_lines(run_pipeline({snarls, view}), "vg snarls");
    run.commands = command_json({snarls, view});
    run.vg_version = vg_version(o.vg);
    return run;
}

std::pair<Store, bool> prepare_tree(const Graph& graph, const std::string& gfa, const TreeOptions& o, const SnarlRun* run) {
    const std::string path = tree_path(gfa, o);
    if (!o.snarl_tree.empty()) return {open_tree(path, gfa, o), true};
    for (const std::string& input : {gfa, o.snarls, o.vg_graph})
        if (!input.empty() && resolve(path) == resolve(input)) throw Error("Conditioning tree output must differ from the input files");
    make_dirs(dirname(path));
    FileLock lock(path + ".lock");
    if (!needs_snarls(gfa, o)) return {open_tree(path, gfa, o), true};
    if (!run) throw Error("Building the conditioning tree needs the snarls");
    const SnarlTree tree = build_snarl_tree(run->records, graph.ids);
    const std::vector<int64_t> node_snarl = assign_nodes(tree, graph.ids, side_adjacency(graph, gfa));
    Json metadata = tree_signature(gfa, o);
    metadata["with_parent"] = (long long)std::count_if(node_snarl.begin(), node_snarl.end(), [](int64_t s) { return s >= 0; });
    metadata["commands"] = run->commands;
    metadata["vg_version"] = run->vg_version.empty() ? Json() : Json(run->vg_version);
    StoreWriter store(path);
    store.add("start_node", tree.start);
    store.add("end_node", tree.end);
    store.add("parent", tree.parent);
    store.add("node_id", graph.ids);
    store.add("node_snarl", node_snarl);
    store.set_metadata(metadata);
    store.publish();
    return {open_tree(path, gfa, o), false};
}

}  // namespace glink
