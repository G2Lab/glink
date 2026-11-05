// Reading a GFA, against GlinkV2's read_gfa and read_topology, and through a graph store.
#include <zlib.h>

#include "golden.hpp"

using namespace glink;
using golden::TempDir;

namespace {

const Json& expected() {
    static const Json j = golden::json("graph.json");
    return j;
}

void expect_graph_matches(const Graph& g) {
    EXPECT_EQ(g.ids, golden::ints(expected().at("ids")));
    EXPECT_EQ(g.names, golden::strings(expected().at("names")));
    ASSERT_EQ(g.bits.rows, expected().at("rows").get<size_t>());
    ASSERT_EQ(g.bits.width, expected().at("width").get<size_t>());
    EXPECT_EQ(g.bits.data, golden::raw<uint8_t>("graph.bits.u8"));
}

std::string graph_file() { return golden::path("graph.gfa"); }

}  // namespace

TEST(Graph, ReadMatchesGlinkV2) { expect_graph_matches(read_gfa(graph_file(), GraphOptions{})); }

TEST(Graph, OrientedReadMatchesGlinkV2) {
    GraphOptions options;
    options.oriented = true;
    const Graph g = read_gfa(graph_file(), options);
    expect_graph_matches(g);
    EXPECT_EQ(g.forward.data, golden::raw<uint8_t>("graph.forward.u8"));
    EXPECT_EQ(g.reverse.data, golden::raw<uint8_t>("graph.reverse.u8"));
    EXPECT_EQ(g.walks, expected().at("walks").get<long long>());
    EXPECT_EQ(g.flipped, expected().at("flipped").get<long long>());
    EXPECT_EQ(g.unanchored, expected().at("unanchored").get<long long>());
}

TEST(Graph, OrientedReadNeedsItsReference) {
    GraphOptions options;
    options.oriented = true;
    options.reference = "NOBODY";
    EXPECT_THROW(read_gfa(graph_file(), options), Error);
}

TEST(Graph, BackboneTipsMatchGlinkV2) {
    for (auto& [backbone, tips] : expected().at("tips").items()) {
        GraphOptions options;
        options.backbone = backbone;
        EXPECT_EQ(read_gfa(graph_file(), options).backbone_tips, golden::ints(tips)) << backbone;
    }
}

TEST(Graph, GzipReadsTheSame) {
    TempDir tmp;
    const std::string packed = tmp.file("graph.gfa.gz");
    const std::string text = read_text(graph_file());
    gzFile out = gzopen(packed.c_str(), "wb");
    ASSERT_TRUE(out);
    ASSERT_EQ(gzwrite(out, text.data(), unsigned(text.size())), int(text.size()));
    gzclose(out);
    expect_graph_matches(read_gfa(packed, GraphOptions{}));
    EXPECT_EQ(file_sha256(packed), file_sha256(graph_file()));
}

TEST(Graph, StoreHoldsWhatTheGfaHas) {
    TempDir tmp;
    GraphOptions options;
    auto [first, reused] = prepare_graph(graph_file(), options, tmp.file("g.h5"));
    EXPECT_FALSE(reused);
    first.load();
    expect_graph_matches(first);
    auto [loaded, again] = prepare_graph(graph_file(), options, tmp.file("g.h5"));
    EXPECT_TRUE(again);
    loaded.load();
    expect_graph_matches(loaded);
    EXPECT_EQ(loaded.backbone_tips, golden::ints(expected().at("tips").at("CHM13")));
}

TEST(Graph, StoreHoldsOrientation) {
    TempDir tmp;
    GraphOptions options;
    options.oriented = true;
    prepare_graph(graph_file(), options, tmp.file("g.h5"));
    auto [loaded, reused] = prepare_graph(graph_file(), options, tmp.file("g.h5"));
    EXPECT_TRUE(reused);
    loaded.load();
    EXPECT_EQ(loaded.forward.data, golden::raw<uint8_t>("graph.forward.u8"));
    EXPECT_EQ(loaded.reverse.data, golden::raw<uint8_t>("graph.reverse.u8"));
    EXPECT_EQ(loaded.flipped, expected().at("flipped").get<long long>());
}

TEST(Graph, RowsReadsOnlyThose) {
    TempDir tmp;
    GraphOptions options;
    prepare_graph(graph_file(), options, tmp.file("g.h5"));
    auto [loaded, reused] = prepare_graph(graph_file(), options, tmp.file("g.h5"));
    ASSERT_TRUE(loaded.on_disk);
    const std::vector<int64_t> rows = {5, 0, 17};
    const Packed got = loaded.rows(rows);
    const auto bits = golden::raw<uint8_t>("graph.bits.u8");
    const size_t width = expected().at("width").get<size_t>();
    ASSERT_EQ(got.rows, rows.size());
    for (size_t k = 0; k < rows.size(); k++)
        EXPECT_EQ(std::vector<uint8_t>(got.row(k), got.row(k) + width),
                  std::vector<uint8_t>(bits.begin() + rows[k] * width, bits.begin() + (rows[k] + 1) * width));
}

TEST(Graph, PackedAndIdIndexAgreeWithTheBits) {
    const Graph g = read_gfa(graph_file(), GraphOptions{});
    for (size_t j = 0; j < g.names.size(); j++)
        EXPECT_EQ(g.bits.get(3, j), bool((g.bits.row(3)[j >> 3] >> (j & 7)) & 1));
    IdIndex index(g.ids);
    for (size_t i = 0; i < g.ids.size(); i++) EXPECT_EQ(index.find(g.ids[i]), int64_t(i));
    EXPECT_EQ(index.find(-5), -1);
    EXPECT_EQ(index.find(g.ids.back() + 1000), -1);
    const std::vector<int64_t> keep = {2, 0};
    const Packed taken = g.bits.take(keep);
    EXPECT_EQ(std::vector<uint8_t>(taken.row(1), taken.row(1) + g.bits.width),
              std::vector<uint8_t>(g.bits.row(0), g.bits.row(0) + g.bits.width));
}
