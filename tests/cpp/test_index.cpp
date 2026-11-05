// glink index and --out-prefix. GlinkV2 has no index, so what the results must equal comes from GlinkV2
// and what the index itself does (building, resuming, refusing) from the README.
#include "golden.hpp"

using namespace glink;
using golden::TempDir;

namespace {

std::vector<std::string> steps(const std::string& out) {
    std::vector<std::string> found;
    std::stringstream stream(out);
    std::string line;
    while (std::getline(stream, line)) {
        if (!starts_with(line, "index: ")) continue;
        const size_t open = line.rfind('('), close = line.find(')', open);
        if (open != std::string::npos && line.find(" -> ") != std::string::npos) found.push_back(line.substr(open + 1, close - open - 1));
    }
    return found;
}

using Words = std::vector<std::string>;

class IndexRuns : public ::testing::Test {
protected:
    void SetUp() override {
        golden::write(tmp.file("graph.gfa"), golden::text("graph.gfa"));
        ASSERT_EQ(index({"--vg", GLINK_VG}).exit, 0);
    }
    golden::Run index(const Words& extra, const std::string& prefix = "idx/graph") {
        Words args = {"index", "--gfa", "graph.gfa", "--out-prefix", prefix};
        args.insert(args.end(), extra.begin(), extra.end());
        return golden::run("glink", args, tmp.path());
    }
    std::string store(const std::string& name) { return tmp.file("idx/graph.glink/" + name + ".h5"); }
    void remove(const std::string& name) { ASSERT_EQ(std::remove(store(name).c_str()), 0); }
    std::string stamp(const std::string& name) {
        return golden::run("/usr/bin/stat", {"-c", "%Y %i", store(name)}, tmp.path()).out;
    }
    TempDir tmp;
    const Words built = {"built", "built", "built"}, reused = {"reused", "reused", "reused"};
};

}  // namespace

TEST_F(IndexRuns, HoldsWhatGlinkV2Computes) {
    const Store tree(store("tree"));
    const Json expected_tree = golden::json("tree.json");
    for (const char* name : {"start_node", "end_node", "parent", "node_id", "node_snarl"})
        EXPECT_EQ(tree.ints(name), golden::ints(expected_tree.at(name))) << name;
    const Store coordinates(store("coords"));
    const Json expected = golden::json("coordinates.json").at(0).at("arrays");
    for (auto& [name, values] : expected.items())
        if (name != "chromosomes") {
            EXPECT_EQ(coordinates.ints(name), golden::ints(values)) << name;
        }
    Graph graph = load_graph(store("graph"));
    graph.load();
    EXPECT_EQ(graph.bits.data, golden::raw<uint8_t>("graph.bits.u8"));
}

TEST_F(IndexRuns, SecondRunReusesEverythingWithoutVg) {
    const std::string before = stamp("graph") + stamp("tree") + stamp("coords");
    const golden::Run r = index({"--vg", "/nonexistent/vg"});
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_EQ(steps(r.out), reused);
    EXPECT_EQ(stamp("graph") + stamp("tree") + stamp("coords"), before);
}

TEST_F(IndexRuns, RebuildBuildsEverything) {
    const golden::Run r = index({"--vg", GLINK_VG, "--rebuild"});
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_EQ(steps(r.out), built);
}

TEST_F(IndexRuns, ResumesWhatIsMissing) {
    remove("coords");
    golden::Run r = index({"--vg", "/nonexistent/vg"});
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_EQ(steps(r.out), (Words{"reused", "built", "reused"}));

    remove("graph");
    r = index({"--vg", "/nonexistent/vg"});
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_EQ(steps(r.out), (Words{"built", "reused", "reused"}));

    remove("tree");
    r = index({"--vg", GLINK_VG});
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_EQ(steps(r.out), (Words{"reused", "reused", "built"}));
}

TEST_F(IndexRuns, MissingTreeWithoutVgFails) {
    remove("tree");
    const golden::Run r = index({"--vg", "/nonexistent/vg"});
    EXPECT_NE(r.exit, 0);
    EXPECT_TRUE(golden::clean(r));
    EXPECT_TRUE(steps(r.out).empty() || steps(r.out).back() != "built");
}

TEST_F(IndexRuns, AnotherBackboneNamesRebuild) {
    const golden::Run r = index({"--backbone", "HG000"});
    EXPECT_EQ(r.exit, 1);
    EXPECT_NE(r.err.find("--rebuild"), std::string::npos) << r.err;
}

TEST_F(IndexRuns, ChangedGfaNamesRebuildAndOverwritesNothing) {
    const std::string before = stamp("graph") + stamp("tree") + stamp("coords");
    golden::write(tmp.file("graph.gfa"), golden::text("graph.gfa") + "S\t99999\tA\n");
    golden::Run r = index({"--vg", GLINK_VG});
    EXPECT_EQ(r.exit, 1);
    EXPECT_NE(r.err.find("--rebuild"), std::string::npos) << r.err;
    EXPECT_EQ(stamp("graph") + stamp("tree") + stamp("coords"), before);
    r = golden::run("glink", {"--gfa", "graph.gfa", "--out-prefix", "idx/graph", "--matrix-out", "m", "--assume-complete"}, tmp.path());
    EXPECT_EQ(r.exit, 1) << r.report();
    r = index({"--vg", GLINK_VG, "--rebuild"});
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_EQ(steps(r.out), built);
}

TEST_F(IndexRuns, MovingOrTouchingTheGfaInvalidatesNothing) {
    ASSERT_EQ(std::rename(tmp.file("graph.gfa").c_str(), tmp.file("moved.gfa").c_str()), 0);
    golden::write(tmp.file("graph.gfa"), golden::text("graph.gfa"));
    const golden::Run r = index({"--vg", "/nonexistent/vg"});
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_EQ(steps(r.out), reused);
}

TEST_F(IndexRuns, ACorruptStoreIsAnError) {
    golden::write(store("coords"), "not a store");
    const golden::Run r = index({"--vg", "/nonexistent/vg"});
    EXPECT_EQ(r.exit, 1) << r.report();
    EXPECT_TRUE(golden::clean(r));
}

TEST_F(IndexRuns, OptionsThatConflictAreUsageErrors) {
    const Words base = {"--gfa", "graph.gfa", "--out-prefix", "idx/graph", "--matrix-out", "m"};
    for (const Words& extra : std::vector<Words>{{"--node-coordinates", "c"}, {"--no-node-coordinates"}, {"--conditioning-tree", "t"},
                                                 {"--no-conditioning-tree"}, {"--snarls", "s"}, {"--snarl-tree", "t"},
                                                 {"--rebuild-conditioning-tree"}, {"--presence-bundle", "p"}}) {
        Words args = base;
        args.insert(args.end(), extra.begin(), extra.end());
        const golden::Run r = golden::run("glink", args, tmp.path());
        EXPECT_EQ(r.exit, 2) << extra[0] << "\n" << r.report();
    }
    EXPECT_EQ(golden::run("glink", {"--gfa", "graph.gfa", "--matrix-out", "m", "--rebuild"}, tmp.path()).exit, 2);
    EXPECT_EQ(golden::run("glink", {"--gfa", "graph.gfa", "--out-prefix", "idx/graph", "--matrix-out", "idx/graph.glink/m"}, tmp.path()).exit, 2);
    EXPECT_EQ(golden::run("glink", {"index", "--gfa", "graph.gfa"}, tmp.path()).exit, 2);
    EXPECT_EQ(golden::run("glink", {"index", "--gfa", "absent.gfa", "--out-prefix", "x"}, tmp.path()).exit, 2);
    EXPECT_EQ(index({"--threads", "0"}).exit, 2);
    EXPECT_EQ(index({"--coordinate-reference-source", "nope"}).exit, 2);
}

TEST_F(IndexRuns, PrefixForms) {
    golden::Run r = index({"--vg", GLINK_VG}, "d/");
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_TRUE(exists(tmp.file("d/graph.glink/graph.h5")));
    r = index({"--vg", GLINK_VG}, "bare");
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_TRUE(exists(tmp.file("bare.glink/tree.h5")));
    const IndexPaths p = index_paths("/x/chr1.gfa.gz", "out/");
    EXPECT_EQ(p.root, "out/chr1.glink");
    EXPECT_EQ(p.oriented_graph, "out/chr1.glink/graph_oriented.h5");
    EXPECT_EQ(index_paths("/x/chr1.gfa", "out/pre").coordinates, "out/pre.glink/coords.h5");
}

namespace {

// Each recorded case run from an index, in one step, must still equal GlinkV2.
class ThroughIndex : public ::testing::TestWithParam<std::string> {};

std::vector<std::string> graph_cases() {
    std::vector<std::string> names;
    for (const char* name : {"parents_phased", "backbone_dosage", "gabriel_parents", "gabriel_parents_dosage", "split_parents", "paper_phased", "paper_dosage",
                                    "split_drop", "plain_gabriel", "drop_parents"})
        names.push_back(name);
    return names;
}

}  // namespace

TEST_P(ThroughIndex, ResultsMatchGlinkV2) {
    const Json& spec = golden::spec_of(GetParam());
    const std::string gfa = spec.at("graph").get<std::string>() + ".gfa";
    const std::string name = "cases/" + GetParam() + "/";
    TempDir tmp;
    const Json matrix = golden::json(name + "matrix.json");
    for (const char* run_name : {"first", "again"}) {
        SCOPED_TRACE(run_name);
        Words args = {"--gfa", gfa, "--out-prefix", tmp.file("idx"), "--vg", GLINK_VG, "--matrix-out", tmp.file("m")};
        if (spec.at("r2").get<bool>()) {
            args.push_back("--r2-out");
            args.push_back(tmp.file("r2.ld"));
        }
        if (spec.at("blocks").get<bool>()) {
            args.push_back("--blocks-out");
            args.push_back(tmp.file("b.tsv"));
        }
        for (const std::string& a : golden::analysis_args(spec)) args.push_back(a);
        const golden::Run r = golden::run("glink", args, golden::dir());
        ASSERT_EQ(r.exit, 0) << r.report();
        const Store m(tmp.file("m"));
        EXPECT_TRUE(golden::exactly(golden::flat(tmp.file("m"), "values"), golden::raw<int8_t>(name + "values.i8")));
        EXPECT_EQ(m.ints("node_ids"), golden::ints(matrix.at("node_ids")));
        if (spec.at("blocks").get<bool>()) {
            EXPECT_TRUE(golden::same_blocks(read_text(tmp.file("b.tsv")), golden::text(name + "blocks.tsv")));
        }
    }
}

INSTANTIATE_TEST_SUITE_P(GlinkV2, ThroughIndex, ::testing::ValuesIn(graph_cases()), golden::test_name);
