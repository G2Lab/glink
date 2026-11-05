// glink index: the graph and everything derived from it, built once under one output prefix.
#pragma once

#include <optional>

#include "coordinates.hpp"
#include "conditioning.hpp"
#include "gfa.hpp"

namespace glink {

struct IndexOptions {
    std::string gfa, prefix, reference = "CHM13", reference_source = "auto", backbone = "CHM13", vg, vg_graph;
    int threads = 1;
    bool rebuild = false, oriented = false;
};

struct IndexPaths {
    std::string root;  // the directory holding the stores below
    std::string graph, oriented_graph, coordinates, tree;
};

// The index is the directory PREFIX.glink. A prefix ending in '/' names it after the GFA.
IndexPaths index_paths(const std::string& gfa, const std::string& prefix);

struct Index {
    IndexPaths paths;
    Graph graph;
    std::optional<Store> coordinates;
    std::optional<Store> tree;
};

// Builds each missing file and reuses each one already there; `rebuild` replaces them all.
// Reports each step on stdout. A file that does not match these inputs is an error, not a rebuild.
Index prepare_index(const IndexOptions& options);

}  // namespace glink
