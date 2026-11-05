// Block table comparison.
#pragma once

#include "common.hpp"

namespace glink {

// glink-compare-blocks: bases of `region` (CHROM:START-END) that block tables A and B each cover and share, as JSON.
Json compare_block_files(const std::string& a, const std::string& b, const std::string& region);

}  // namespace glink
