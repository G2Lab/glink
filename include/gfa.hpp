// The one reader of the GFA, and the graph store that saves what it read. Every later step works from a Graph.
#pragma once

#include "common.hpp"
#include "engine.hpp"

namespace glink {

// Node x haplotype bits; bit j of row i (little-endian within a byte) = haplotype j visits node i.
struct Packed {
    size_t rows = 0, width = 0;
    std::vector<uint8_t> data;

    Packed() = default;
    Packed(size_t rows, size_t width, uint8_t fill = 0) : rows(rows), width(width), data(rows * width, fill) {}
    uint8_t* row(size_t i) { return data.data() + i * width; }
    const uint8_t* row(size_t i) const { return data.data() + i * width; }
    bool get(size_t i, size_t j) const { return (row(i)[j >> 3] >> (j & 7)) & 1; }
    Packed take(const std::vector<int64_t>& rows) const;
    Packed take(const std::vector<uint8_t>& mask) const;
};

struct Graph {
    std::vector<int64_t> ids;          // sorted unique segment IDs
    std::vector<std::string> names;    // sorted biological haplotype names
    Packed bits;
    bool oriented = false;
    Packed forward, reverse;           // with oriented: reference-oriented reads of each node
    long long walks = 0, flipped = 0, unanchored = 0;  // with oriented
    // For the later steps. coordinate_error is why the GFA cannot give coordinates, if it cannot.
    std::vector<Link> links;
    std::vector<int64_t> backbone_tips;
    CoordinateInput coordinate_input;
    std::string coordinate_error;
    // A loaded graph leaves `bits` on disk until rows() or load() asks for them.
    std::string bits_path;
    bool on_disk = false;

    Packed rows(const std::vector<int64_t>& keep) const;  // selected rows of `bits`
    void load();                                          // all of `bits` into memory
    void release_step_data();                             // the matrix needs none of the later steps' input
};

struct GraphOptions {
    bool oriented = false;           // harmonize walk strands, for --split-inversions
    std::string reference = "CHM13";  // for coordinates, and to orient walks against
    std::string backbone = "CHM13";   // roots the snarl tree
};

Graph read_gfa(const std::string& path, const GraphOptions& options);
// Whether the store at `path` holds this GFA read with these options; throws if not.
void check_graph(const std::string& path, const std::string& gfa, const GraphOptions& options);
Graph load_graph(const std::string& path, bool oriented = false);
void save_graph(const Graph& g, const std::string& path, const std::string& identity, const GraphOptions& options);
// The graph from the store at `output` (default: under .glink-cache/graph/), reading the GFA and
// saving it there only if it is not there yet. The second value says whether it was reused.
std::pair<Graph, bool> prepare_graph(const std::string& gfa, const GraphOptions& options, std::string output = "", bool rebuild = false);

// Sorted-ID lookup: index of each ID, or -1.
class IdIndex {
public:
    explicit IdIndex(const std::vector<int64_t>& sorted_ids);
    int64_t find(int64_t id) const;

private:
    const std::vector<int64_t>* ids_;
    std::vector<int32_t> dense_;
    int64_t max_ = -1;
};

}  // namespace glink
