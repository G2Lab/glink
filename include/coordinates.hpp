// Node coordinates: a store with one integer array per field, each with a row per node, sorted by node ID.
#pragma once

#include "store.hpp"

namespace glink {

struct Graph;

Json coordinates_signature(const std::string& gfa, const std::string& reference, const std::string& source);
// Opens a coordinate store, checking its schema and that its arrays agree.
Store load_coordinates(const std::string& path);
// Every field's values for these nodes. Only their rows are read. *chrom fields index the
// "chromosomes" string list; -1 is none.
ColumnMap coordinates_of(const Store& coordinates, const std::vector<int64_t>& ids);
Json coordinate_provenance(const Store& coordinates);
std::pair<Store, bool> prepare_coordinates(const Graph& graph, const std::string& gfa, std::string output,
                                           const std::string& reference, const std::string& source, bool rebuild);

}  // namespace glink
