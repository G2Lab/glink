// glink-compare-blocks: the shared sequence coverage of two block tables, as JSON.
#include <cstdio>
#include <string>
#include <vector>

#include "matching.hpp"

using namespace glink;

namespace {

int run(int argc, char** argv) {
    std::vector<std::string> args;
    std::string region, out;
    for (int i = 1; i < argc; i++) {
        const std::string a = argv[i];
        if (a == "--region" && i + 1 < argc) region = argv[++i];
        else if (a == "--out" && i + 1 < argc) out = argv[++i];
        else args.push_back(a);
    }
    if (args.size() != 2 || region.empty()) {
        std::fprintf(stderr, "usage: glink-compare-blocks A B --region CHROM:START-END [--out scores.json]\n");
        return 2;
    }
    const std::string text = compare_block_files(args[0], args[1], region).dump(2) + "\n";
    if (out.empty()) std::fputs(text.c_str(), stdout);
    else {
        FILE* f = std::fopen(out.c_str(), "w");
        if (!f) throw Error("Cannot write " + out);
        std::fputs(text.c_str(), f);
        std::fclose(f);
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(argc, argv);
    } catch (const std::exception& e) {
        std::fflush(stdout);
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
