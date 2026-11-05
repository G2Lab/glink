// The conditioning tree and node parents, against GlinkV2's build_snarl_tree, assign_nodes and from_tree.
#include "golden.hpp"

using namespace glink;
using golden::TempDir;

namespace {

void expect_tree_matches(const Store& tree, const Json& expected) {
    for (const char* name : {"start_node", "end_node", "parent", "node_id", "node_snarl"})
        EXPECT_EQ(tree.ints(name), golden::ints(expected.at(name))) << name;
}

const Json& expected() {
    static const Json j = golden::json("tree.json");
    return j;
}

std::string gfa() { return golden::path("graph.gfa"); }

Store build(const TreeOptions& o, const std::string& graph_file) {
    const Graph graph = read_gfa(graph_file, GraphOptions{});
    const SnarlRun snarls = run_snarls(graph, graph_file, o);
    return prepare_tree(graph, graph_file, o, &snarls).first;
}

}  // namespace

TEST(Tree, FromSavedJsonMatchesGlinkV2) {
    TempDir tmp;
    TreeOptions o;
    o.output = tmp.file("t.h5");
    o.snarls = golden::path("graph.snarls.json");
    o.snarls_format = "json";
    expect_tree_matches(build(o, gfa()), expected());
}

TEST(Tree, FormatIsDetected) {
    TempDir tmp;
    TreeOptions o;
    o.output = tmp.file("t.h5");
    o.snarls = golden::path("graph.snarls.json");
    expect_tree_matches(build(o, gfa()), expected());
}

TEST(Tree, FromSavedBinaryMatchesGlinkV2) {
    TempDir tmp;
    TreeOptions o;
    o.output = tmp.file("t.h5");
    o.snarls = golden::path("snarls.bin");
    o.vg = GLINK_VG;
    expect_tree_matches(build(o, gfa()), expected());
}

TEST(Tree, BinaryFormatCanBeNamed) {
    TempDir tmp;
    TreeOptions o;
    o.output = tmp.file("t.h5");
    o.snarls = golden::path("snarls.bin");
    o.snarls_format = "binary";
    o.vg = GLINK_VG;
    expect_tree_matches(build(o, gfa()), expected());
}

// vg's own output could change between versions, and GlinkV2's was recorded with one of them.
TEST(Tree, RunningVgMatchesGlinkV2) {
    TempDir tmp;
    TreeOptions o;
    o.output = tmp.file("t.h5");
    o.vg = GLINK_VG;
    expect_tree_matches(build(o, gfa()), expected());
}

TEST(Tree, ThreadedVgMatchesGlinkV2) {
    TempDir tmp;
    TreeOptions o;
    o.output = tmp.file("t.h5");
    o.vg = GLINK_VG;
    o.threads = 4;
    expect_tree_matches(build(o, gfa()), expected());
}

TEST(Tree, TheBackboneRootsItAsInGlinkV2) {
    const std::string forked = golden::path("small/forked.gfa");
    const Json all = golden::json("tree_backbones.json");
    for (auto& [backbone, expected_tree] : all.items()) {
        SCOPED_TRACE(backbone);
        TempDir tmp;
        TreeOptions o;
        o.output = tmp.file("t.h5");
        o.vg = GLINK_VG;
        o.backbone = backbone;
        GraphOptions options;
        options.backbone = backbone;
        const Graph graph = read_gfa(forked, options);
        const SnarlRun snarls = run_snarls(graph, forked, o);
        expect_tree_matches(prepare_tree(graph, forked, o, &snarls).first, expected_tree);
    }
}

TEST(Tree, GbzMatchesGlinkV2) {
    TempDir tmp;
    TreeOptions o;
    o.output = tmp.file("t.h5");
    o.vg = GLINK_VG;
    o.vg_graph = golden::path("graph.gbz");
    expect_tree_matches(build(o, golden::path("graph.gbz.gfa")), golden::json("tree_gbz.json"));
}

TEST(Tree, GbzOfAnotherGraphIsAnError) {
    TempDir tmp;
    TreeOptions o;
    o.output = tmp.file("t.h5");
    o.vg = GLINK_VG;
    o.vg_graph = golden::path("graph.gbz");
    const std::string other = golden::path("small/rgfa.gfa");
    EXPECT_THROW(run_snarls(read_gfa(other, GraphOptions{}), other, o), Error);
}

TEST(Tree, NodeParentsMatchGlinkV2) {
    TempDir tmp;
    TreeOptions o;
    o.output = tmp.file("t.h5");
    o.snarls = golden::path("graph.snarls.json");
    const ParentTable table = parent_table(build(o, gfa()));
    EXPECT_EQ(table.node_id, golden::ints(expected().at("node_id")));
    EXPECT_EQ(table.parent_node_id, golden::ints(expected().at("parent_node_id")));
    EXPECT_THROW(table.lookup({999999}), Error);
}

TEST(TreeStore, SecondRunReusesIt) {
    TempDir tmp;
    TreeOptions o;
    o.output = tmp.file("t.h5");
    o.snarls = golden::path("graph.snarls.json");
    const Graph graph = read_gfa(gfa(), GraphOptions{});
    EXPECT_TRUE(needs_snarls(gfa(), o));
    const SnarlRun snarls = run_snarls(graph, gfa(), o);
    EXPECT_FALSE(prepare_tree(graph, gfa(), o, &snarls).second);
    EXPECT_FALSE(needs_snarls(gfa(), o));
    EXPECT_TRUE(prepare_tree(graph, gfa(), o, nullptr).second);
    o.rebuild = true;
    EXPECT_TRUE(needs_snarls(gfa(), o));
}

TEST(TreeStore, AnExistingTreeCanBeNamed) {
    TempDir tmp;
    TreeOptions o;
    o.output = tmp.file("t.h5");
    o.snarls = golden::path("graph.snarls.json");
    build(o, gfa());
    TreeOptions named;
    named.snarl_tree = tmp.file("t.h5");
    EXPECT_FALSE(needs_snarls(gfa(), named));
    const auto [tree, reused] = prepare_tree(read_gfa(gfa(), GraphOptions{}), gfa(), named, nullptr);
    EXPECT_TRUE(reused);
    expect_tree_matches(tree, expected());
}

TEST(TreeStore, ATreeOfAnotherGraphIsAnError) {
    TempDir tmp;
    TreeOptions o;
    o.output = tmp.file("t.h5");
    o.snarls = golden::path("graph.snarls.json");
    build(o, gfa());
    TreeOptions named;
    named.snarl_tree = tmp.file("t.h5");
    EXPECT_THROW(needs_snarls(golden::path("plain.gfa"), named), Error);
}

TEST(TreeStore, AnotherBackboneIsAnError) {
    TempDir tmp;
    TreeOptions o;
    o.output = tmp.file("t.h5");
    o.snarls = golden::path("graph.snarls.json");
    build(o, gfa());
    o.backbone = "HG000";
    EXPECT_THROW(needs_snarls(gfa(), o), Error);
}

TEST(TreeStore, OutputMustNotBeAnInput) {
    TreeOptions o;
    o.output = gfa();
    o.snarls = golden::path("graph.snarls.json");
    const Graph graph = read_gfa(gfa(), GraphOptions{});
    const SnarlRun snarls = run_snarls(graph, gfa(), o);
    EXPECT_THROW(prepare_tree(graph, gfa(), o, &snarls), Error);
}
