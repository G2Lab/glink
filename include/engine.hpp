// Public C++ interface of the LD and coordinate engine. Callers validate inputs
// before calling these.
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace glink {

// Two-locus genotype count table, row-major: counts[3*a + b] for states a, b.
using Counts = std::array<uint16_t, 9>;

// A D′ point estimate with its likelihood interval and classification.
struct Interval {
    double dprime, lower, upper;
    int category;
};

// ---- Diploid hard-call D′ (dosage.cpp) ------------------------------------

// Category 1 = strong LD (lower > .70, upper >= .98), 2 = recombination
// (upper < .90), 0 = uninformative. Unsupported tables have NaN estimates.
Interval dosage_confidence(Counts counts);

// Pair categories for k markers. `packed` is k × 3 × width bytes: one bit
// plane per genotype state (0/1/2) per marker. Writes a k × k matrix.
void dosage_classes(const uint8_t* packed, int k, int width, uint8_t* out, int threads);

// Absolute D′ estimates alone, with no confidence interval, under the same
// support rule as dosage_confidence (NaN when unsupported). Writes a k × k matrix.
void dosage_dprime(const uint8_t* packed, int k, int width, float* out, int threads);

// ---- Gabriel blocks (gabriel.cpp) -----------------------------------------

// Category is the PLINK 1.9 CI band 0..6; lower/upper are reporting bounds.
Interval gabriel_confidence(Counts counts, int ploidy);

// Band of a perfectly complementary pair under `neutral_complements`: neither
// strong LD nor recombination, never a block endpoint, but it passes the 3- and
// 4-marker rules' inner-pair checks (as a redundant copy of one marker would).
constexpr uint8_t NEUTRAL = 7;

// CI bands for every pair within `span` bp. `packed` is n × 3 × width bytes
// as for dosage_classes; positions are sorted. Pairs beyond the span and the
// diagonal are band 1 (uninformative). Writes an n × n matrix.
void gabriel_bands(const uint8_t* packed, int n, int width, const int64_t* positions,
                   int64_t span, int ploidy, uint8_t* out, int threads,
                   bool neutral_complements = false);

// Selected blocks as (first, last) marker indices, sorted by last marker.
std::vector<std::pair<int, int>> gabriel_select(const uint8_t* bands, const int64_t* positions,
                                                int n, int64_t span);

// gabriel_bands then gabriel_select without storing the bands: each marker's row of
// pairs within `span` is computed, used and dropped, so memory grows with the widest
// window rather than with the number of pairs.
std::vector<std::pair<int, int>> gabriel_blocks(const uint8_t* packed, int n, int width,
                                                const int64_t* positions, int64_t span,
                                                int ploidy, int threads,
                                                bool neutral_complements = false);

// ---- Node coordinates (coordinates.cpp) -----------------------------------

// Bump when the coordinate algorithm changes, so caches built by an older
// engine are rejected.
constexpr int kCoordinatesVersion = 1;

struct Coordinates {
    // One entry per node, sorted by node ID.
    std::vector<int64_t> node_id, length, reference_start, reference_end, upstream, downstream,
        upstream_node, downstream_node, position;
    std::vector<int32_t> reference_chrom, upstream_chrom, downstream_chrom, chrom;
    std::vector<uint8_t> is_reference, upstream_ambiguous, downstream_ambiguous, ambiguous,
        reversed_anchors, position_source;
    std::vector<std::string> chromosomes;
    int64_t links = 0, reference_nodes = 0;
    std::string reference_source;
};

// An oriented link between two nodes.
struct Link { int64_t a, b; bool a_reverse, b_reverse; };

// What coordinates need from a GFA besides its links, in file order. read_gfa fills this in.
struct CoordinateInput {
    struct Walk { std::string chrom, steps; int64_t start, end; bool is_path; };
    struct Tag { int64_t id, start; std::string chrom; };
    std::vector<int64_t> ids, lengths;
    std::vector<Walk> walks;
    std::vector<Tag> tags;
};

// Throws unless s is a whole nonnegative number. Coordinates are strict about their numbers.
int64_t nonnegative_integer(const std::string& s);

// EdgeDepth-style anchors. source is "auto", "walk" or "rgfa".
Coordinates compute_coordinates(const CoordinateInput& input, const std::vector<Link>& links,
                                const std::string& reference, const std::string& source);

}  // namespace glink
