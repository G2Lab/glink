// Gabriel et al. (2002) LD blocks with PLINK 1.9's rules, with per-block statistics and TSV output.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iterator>
#include <numeric>
#include <random>
#include <unordered_map>

#include "blocks.hpp"
#include "matrix.hpp"
#include "engine.hpp"

namespace glink {

namespace {

constexpr int64_t UNASSIGNED = -1;
constexpr size_t SAMPLE = 4096;  // block statistics use at most this many members

// A square r² matrix held in memory.
class R2 {
public:
    explicit R2(const std::vector<float>& values) : n_(size_t(std::llround(std::sqrt(double(values.size()))))), data_(values.data()) {
        if (n_ * n_ != values.size()) throw Error("r² matrix must be square");
    }
    size_t size() const { return n_; }
    const float* row(size_t i) const { return data_ + i * n_; }

private:
    size_t n_;
    const float* data_;
};


std::vector<size_t> stable_order(const std::vector<int64_t>& positions) {
    std::vector<size_t> order(positions.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return positions[a] < positions[b]; });
    return order;
}

std::vector<int64_t> method_gabriel(const std::vector<int64_t>& positions, const Values& values, int ploidy,
                                    int64_t max_span_bp, int threads, bool neutral_complements) {
    const size_t n = positions.size();
    auto order = stable_order(positions);
    std::vector<int64_t> pos(n);
    for (size_t k = 0; k < n; k++) pos[k] = positions[order[k]];
    if ((ploidy != 1 && ploidy != 2) || values.rows != n) throw Error("Invalid ordered genotype panel");
    for (int8_t x : values.v) if (x < -1 || x > ploidy) throw Error("Invalid ordered genotype panel");
    if (values.cols > 65535 || threads < 1 || max_span_bp < 0) throw Error("Invalid cohort size, thread count or span");
    const size_t width = (values.cols + 7) / 8;
    std::vector<uint8_t> packed(n * 3 * width, 0);
    for (size_t k = 0; k < n; k++) {
        const int8_t* x = values.row(order[k]);
        for (size_t c = 0; c < values.cols; c++)
            if (x[c] >= 0) packed[(k * 3 + size_t(x[c])) * width + c / 8] |= uint8_t(0x80 >> (c % 8));
    }
    auto blocks = gabriel_blocks(packed.data(), int(n), int(width), pos.data(), max_span_bp, ploidy, threads,
                                 neutral_complements);
    std::vector<int64_t> labels(n);
    std::iota(labels.begin(), labels.end(), 0);
    for (auto [first, last] : blocks)
        for (int k = first; k <= last; k++) labels[order[k]] = int64_t(n) + first;
    return labels;
}

std::vector<int64_t> renumber(const std::vector<int64_t>& labels) {
    std::vector<int64_t> out(labels.size(), UNASSIGNED);
    std::unordered_map<int64_t, int64_t> seen;
    for (size_t i = 0; i < labels.size(); i++) {
        if (labels[i] < 0) continue;
        auto [it, inserted] = seen.emplace(labels[i], int64_t(seen.size()));
        out[i] = it->second;
    }
    return out;
}

std::vector<BlockStats> describe(const R2* m, const std::vector<int64_t>& labels, int64_t count,
                                 const std::vector<int64_t>& positions) {
    std::mt19937_64 rng(0);  // fixed seed: the same input gives the same statistics
    std::vector<std::vector<int64_t>> groups(count);
    for (size_t i = 0; i < labels.size(); i++) if (labels[i] >= 0) groups[labels[i]].push_back(int64_t(i));
    std::vector<BlockStats> out;
    for (int64_t b = 0; b < count; b++) {
        auto& members = groups[b];
        BlockStats s;
        s.block = b;
        s.size = int64_t(members.size());
        if (!members.empty()) {
            int64_t lo = INT64_MAX, hi = INT64_MIN;
            for (auto x : members) { lo = std::min(lo, positions[x]); hi = std::max(hi, positions[x]); }
            s.placed = true;
            s.start = lo; s.end = hi; s.span = hi - lo;
        }
        if (members.size() > 1 && m) {
            std::vector<int64_t> chosen;
            if (members.size() > SAMPLE) std::sample(members.begin(), members.end(), std::back_inserter(chosen), SAMPLE, rng);
            else chosen = members;
            std::vector<float> upper;
            upper.reserve(chosen.size() * (chosen.size() - 1) / 2);
            for (size_t a = 0; a < chosen.size(); a++) {
                const float* row = m->row(chosen[a]);
                for (size_t c = a + 1; c < chosen.size(); c++) {
                    float v = row[chosen[c]];
                    if (std::isfinite(v)) upper.push_back(v);
                }
            }
            if (!upper.empty()) {
                float sum = numpy_sum(upper.data(), upper.size());
                s.scored = true;
                s.mean_r2 = double(float(double(sum) / double(upper.size())));
                s.min_r2 = double(*std::min_element(upper.begin(), upper.end()));
            }
        }
        out.push_back(s);
    }
    return out;
}

}  // namespace

int64_t BlockResult::count() const {
    int64_t mx = -1;
    for (auto l : labels) mx = std::max(mx, l);
    return mx + 1;
}

BlockResult build_blocks(const std::vector<float>* r2, const std::vector<int64_t>& positions, const Values& values, int ploidy,
                         int64_t max_span_bp, int threads, bool neutral_complements) {
    std::unique_ptr<R2> matrix;
    if (r2) matrix = std::make_unique<R2>(*r2);
    if (matrix && matrix->size() != values.rows)
        throw Error(std::to_string(values.rows) + " genotype rows for " + std::to_string(matrix->size()) + " features");
    if (positions.size() != values.rows)
        throw Error(std::to_string(positions.size()) + " positions for " + std::to_string(values.rows) + " features");
    BlockResult result;
    result.labels = renumber(method_gabriel(positions, values, ploidy, max_span_bp, threads, neutral_complements));
    result.statistics = describe(matrix.get(), result.labels, result.count(), positions);
    return result;
}

void write_blocks(const BlockResult& result, const std::vector<std::string>& ids, const std::vector<int64_t>& positions,
                  const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw Error("Cannot write " + path);
    std::string out = "feature\tid\tposition\tblock\tblock_size\tblock_start\tblock_end\tblock_span\tblock_mean_r2\tblock_min_r2\n";
    for (size_t i = 0; i < result.labels.size(); i++) {
        const int64_t label = result.labels[i];
        out += std::to_string(i) + "\t" + ids[i] + "\t" + std::to_string(positions[i]) + "\t";
        if (label >= 0) {
            const BlockStats& s = result.statistics[label];
            out += std::to_string(label) + "\t" + std::to_string(s.size) + "\t";
            out += s.placed ? std::to_string(s.start) + "\t" + std::to_string(s.end) + "\t" + std::to_string(s.span) + "\t" : "\t\t\t";
            out += s.scored ? shortest(s.mean_r2) + "\t" + shortest(s.min_r2) : "\t";
        } else out += "\t\t\t\t\t\t";
        out += "\n";
        if (out.size() > (1 << 20)) { std::fwrite(out.data(), 1, out.size(), f); out.clear(); }
    }
    std::fwrite(out.data(), 1, out.size(), f);
    if (std::fclose(f) != 0) throw Error("Cannot write " + path);
}

}  // namespace glink
