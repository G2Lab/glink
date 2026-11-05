// glink index: the graph, vg snarls, node coordinates and the conditioning tree.
#include <cstdio>
#include <optional>

#include "gfa.hpp"
#include "index.hpp"

namespace glink {

namespace {

std::string gfa_stem(const std::string& gfa) {
    std::string name = basename(gfa);
    if (ends_with(name, ".gz")) name = name.substr(0, name.size() - 3);
    if (ends_with(name, ".gfa")) name = name.substr(0, name.size() - 4);
    return name;
}

void report(const char* what, bool reused, const std::string& path) {
    std::printf("index: %s (%s) -> %s\n", what, reused ? "reused" : "built", path.c_str());
    std::fflush(stdout);
}

}  // namespace

IndexPaths index_paths(const std::string& gfa, const std::string& prefix) {
    IndexPaths p;
    p.root = (ends_with(prefix, "/") ? prefix + gfa_stem(gfa) : prefix) + ".glink";
    p.graph = p.root + "/graph.h5";
    p.oriented_graph = p.root + "/graph_oriented.h5";
    p.coordinates = p.root + "/coords.h5";
    p.tree = p.root + "/tree.h5";
    return p;
}

Index prepare_index(const IndexOptions& o) {
    Index index;
    index.paths = index_paths(o.gfa, o.prefix);
    const IndexPaths& p = index.paths;
    // Only building the tree needs vg, and a finished tree does not, so an analysis run can resume
    // where vg is not installed. With a stale tree, stop before reading the GFA.
    TreeOptions tree_options;
    tree_options.output = p.tree; tree_options.backbone = o.backbone; tree_options.vg = o.vg;
    tree_options.vg_graph = o.vg_graph; tree_options.threads = o.threads; tree_options.rebuild = o.rebuild;
    const bool call_vg = needs_snarls(o.gfa, tree_options);

    // An oriented graph also serves a run that does not need orientation.
    GraphOptions graph_options;
    graph_options.oriented = o.oriented;
    graph_options.reference = o.reference;
    graph_options.backbone = o.backbone;
    const bool use_oriented = o.oriented || (!o.rebuild && !exists(p.graph) && exists(p.oriented_graph));
    const std::string graph_path = use_oriented ? p.oriented_graph : p.graph;
    auto [graph, graph_reused] = prepare_graph(o.gfa, graph_options, graph_path, o.rebuild);
    report("graph", graph_reused, graph_path);

    std::optional<SnarlRun> snarls;
    if (call_vg) {
        snarls = run_snarls(graph, o.gfa, tree_options);
        std::printf("index: vg snarls (%zu snarls)\n", snarls->records.size());
        std::fflush(stdout);
    }
    auto [coordinates, coordinates_reused] = prepare_coordinates(graph, o.gfa, p.coordinates, o.reference, o.reference_source, o.rebuild);
    index.coordinates = std::move(coordinates);
    report("node coordinates", coordinates_reused, p.coordinates);
    auto [tree, tree_reused] = prepare_tree(graph, o.gfa, tree_options, snarls ? &*snarls : nullptr);
    index.tree = std::move(tree);
    report("conditioning tree", tree_reused, p.tree);
    index.graph = std::move(graph);
    return index;
}

}  // namespace glink
