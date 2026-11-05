// The node matrix: filters, parent conditioning and building and writing the matrix.
#pragma once

#include <optional>

#include "common.hpp"
#include "conditioning.hpp"

namespace glink {

struct Packed;
struct Graph;

constexpr int8_t FORWARD = 1, REVERSE = -1, UNORIENTED = 0;

struct Values {
    size_t rows = 0, cols = 0;
    std::vector<int8_t> v;
    int8_t* row(size_t i) { return v.data() + i * cols; }
    const int8_t* row(size_t i) const { return v.data() + i * cols; }
    Values take(const std::vector<int64_t>& keep) const;
};

struct Filter {
    std::string kind;  // all, band, minor, common, list
    double low = 0, high = 1, value = 0;
    int people = 0;    // common: fewest called people
    std::string path;
};

std::optional<Filter> parse_filter(const std::string& spec, std::string& error);

struct NodeMatrix {
    Values values;
    std::vector<int64_t> node_ids;
    std::vector<std::string> columns;
    std::string matrix_type;
    bool collapsed = false;
    std::vector<int64_t> pattern_of, input_node_ids;
    Json callability;
    ColumnMap coordinates;                 // *chrom fields index coordinate_chromosomes
    std::vector<std::string> coordinate_chromosomes;
    Json coordinate_metadata;
    std::vector<int64_t> parent_node_ids;
    Json parent_metadata;
    bool oriented = false;
    std::vector<int8_t> node_orientation, input_node_orientation;
    std::vector<int64_t> complement_pairs;   // with drop_complements: (dropped id, kept complement id), row-major

    double density() const;
    std::vector<std::string> row_labels() const;
};

struct BuildOptions {
    std::string matrix_type = "dosage", samples;
    std::vector<Filter> filters;
    bool collapse = true, assume_complete = false, condition_on_parents = false, drop_complements = false;
    bool collapse_backbone = false;
    int broken_min_haplotypes = 6;
    std::string reference = "CHM13";  // that walks were oriented against, for the record
    const Store* coordinates = nullptr;
    const Store* parents = nullptr;  // the conditioning tree
};

// Takes the graph, because it frees each part as soon as it can.
NodeMatrix build_matrix(Graph graph, const BuildOptions& options);
void write_matrix(const NodeMatrix& matrix, const std::string& path);

}  // namespace glink
