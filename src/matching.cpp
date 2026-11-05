// glink-compare-blocks: the shared sequence coverage of Indap et al. (2005).
#include <algorithm>
#include <map>

#include "matching.hpp"

namespace glink {

namespace {

struct Block { int64_t start0, end0, markers; };

std::vector<std::string> split_on(const std::string& line, char separator) {
    std::vector<std::string> out;
    size_t start = 0;
    while (true) {
        const size_t at = line.find(separator, start);
        if (at == std::string::npos) { out.push_back(line.substr(start)); return out; }
        out.push_back(line.substr(start, at - start));
        start = at + 1;
    }
}

// Blocks from PLINK .blocks.det or a glink blocks table.
std::vector<Block> read_blocks(const std::string& path) {
    std::vector<std::string> lines;
    for (const std::string& raw : split_on(read_text(path), '\n'))
        if (raw.find_first_not_of(" \t\r\f\v") != std::string::npos) lines.push_back(raw);
    if (lines.empty()) throw Error(path + ": expected PLINK .blocks.det or a glink blocks table");
    const auto header = split_ws(lines[0]);
    const auto column = [&](const std::vector<std::string>& names, const std::string& name) {
        const auto at = std::find(names.begin(), names.end(), name);
        if (at == names.end()) throw Error(path + ": no column " + name);
        return size_t(at - names.begin());
    };
    std::vector<Block> blocks;
    if (header.size() >= 3 && header[0] == "CHR" && header[1] == "BP1" && header[2] == "BP2") {
        const size_t bp1 = 1, bp2 = 2, snps = column(header, "NSNPS");
        for (size_t i = 1; i < lines.size(); i++) {
            const auto row = split_ws(lines[i]);
            blocks.push_back({std::stoll(row.at(bp1)) - 1, std::stoll(row.at(bp2)), std::stoll(row.at(snps))});
        }
        return blocks;
    }
    if (std::find(header.begin(), header.end(), "block_start") == header.end())
        throw Error(path + ": expected PLINK .blocks.det or a glink blocks table");
    const auto names = split_on(lines[0], '\t');
    const size_t block = column(names, "block"), start = column(names, "block_start"), end = column(names, "block_end"),
                 size = column(names, "block_size");
    // One entry per block label, first-seen order, last row's values.
    std::map<std::string, size_t> seen;
    std::vector<Block> all;
    for (size_t i = 1; i < lines.size(); i++) {
        const auto row = split_on(lines[i], '\t');
        if (row.size() <= std::max({block, start, end, size}) || row[block].empty() || row[start].empty()) continue;
        const Block b{std::stoll(row[start]) - 1, std::stoll(row[end]), std::stoll(row[size])};
        const auto [at, fresh] = seen.emplace(row[block], all.size());
        if (fresh) all.push_back(b);
        else all[at->second] = b;
    }
    // A one-marker group is not a block; PLINK never reports one.
    for (const Block& b : all) if (b.markers > 1) blocks.push_back(b);
    std::stable_sort(blocks.begin(), blocks.end(), [](const Block& x, const Block& y) { return x.start0 < y.start0; });
    return blocks;
}

std::vector<uint8_t> covered(const std::vector<Block>& blocks, int64_t start0, int64_t end0) {
    std::vector<uint8_t> mask(size_t(end0 - start0));
    for (const Block& b : blocks) {
        const int64_t from = std::max(b.start0, start0) - start0, to = std::max<int64_t>(std::min(b.end0, end0) - start0, 0);
        for (int64_t i = from; i < to; i++) mask[size_t(i)] = 1;
    }
    return mask;
}

}  // namespace

Json compare_block_files(const std::string& path_a, const std::string& path_b, const std::string& region) {
    const size_t colon = region.find(':'), dash = region.find('-', colon);
    if (colon == std::string::npos || dash == std::string::npos) throw Error("Region must be CHROM:START-END");
    const int64_t start0 = std::stoll(region.substr(colon + 1, dash - colon - 1)) - 1, end0 = std::stoll(region.substr(dash + 1));
    const std::vector<Block> a = read_blocks(path_a), b = read_blocks(path_b);
    const auto ma = covered(a, start0, end0), mb = covered(b, start0, end0);
    long long a_bp = 0, b_bp = 0, shared = 0;
    for (size_t i = 0; i < ma.size(); i++) { a_bp += ma[i]; b_bp += mb[i]; shared += ma[i] & mb[i]; }
    Json out = Json::object();
    out["a"] = path_a;
    out["b"] = path_b;
    out["a_blocks"] = (long long)a.size();
    out["b_blocks"] = (long long)b.size();
    out["region"] = region;
    out["region_bp"] = (long long)(end0 - start0);
    out["a_bp"] = a_bp;
    out["b_bp"] = b_bp;
    out["shared_bp"] = shared;
    return out;
}

}  // namespace glink
