// Gabriel et al. (2002) LD blocks.
#pragma once

#include "common.hpp"

namespace glink {

struct Values;

struct BlockStats {
    int64_t block = 0, size = 0;
    bool placed = false;
    int64_t start = 0, end = 0, span = 0;
    bool scored = false;
    double mean_r2 = 0, min_r2 = 0;
};

struct BlockResult {
    std::vector<int64_t> labels;
    std::vector<BlockStats> statistics;
    int64_t count() const;
};

BlockResult build_blocks(const std::vector<float>* r2, const std::vector<int64_t>& positions, const Values& values, int ploidy,
                         int64_t max_span_bp, int threads, bool neutral_complements = false);
void write_blocks(const BlockResult& result, const std::vector<std::string>& ids,
                  const std::vector<int64_t>& positions, const std::string& path);

}  // namespace glink
