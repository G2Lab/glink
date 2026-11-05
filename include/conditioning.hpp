// The conditioning tree: vg's snarls, each node's innermost snarl, and from those each node's parent.
#pragma once

#include "store.hpp"

namespace glink {

struct Graph;

// Each node's parent (-1 for none), sorted by node. Collapsing makes one that is in no store.
struct ParentTable {
    std::vector<int64_t> node_id, parent_node_id;
    std::vector<int64_t> lookup(const std::vector<int64_t>& ids) const;
};

struct TreeOptions {
    std::string output, backbone = "CHM13", vg = "vg", vg_graph, snarls, snarl_tree, snarls_format = "auto";
    int threads = 1;
    bool rebuild = false;
};

// What vg snarls reported, as JSON records, and how it was run.
struct SnarlRun {
    std::vector<Json> records;
    Json commands;
    std::string vg_version;
};

// Whether these options need the snarls: false when a tree store is there to reuse. Throws if
// that store does not belong to these inputs, so this also stops a stale run before it reads the GFA.
bool needs_snarls(const std::string& gfa, const TreeOptions& options);
// The snarls, from vg or from options.snarls.
SnarlRun run_snarls(const Graph& graph, const std::string& gfa, const TreeOptions& options);
// A tree is a store with the snarl arrays start_node, end_node and parent, and the node arrays
// node_id and node_snarl. It reuses or builds one; `run` is needed only to build.
std::pair<Store, bool> prepare_tree(const Graph& graph, const std::string& gfa, const TreeOptions& options, const SnarlRun* run);
// A node's parent is the start node of its innermost snarl.
ParentTable parent_table(const Store& tree);

}  // namespace glink
