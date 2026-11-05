// Node coordinates, against GlinkV2's prepare_coordinates, on each graph, reference and source.
#include "golden.hpp"

using namespace glink;
using golden::TempDir;

namespace {

class Coordinate : public ::testing::TestWithParam<size_t> {
protected:
    static const Json& runs() {
        static const Json j = golden::json("coordinates.json");
        return j;
    }
    const Json& run() const { return runs().at(GetParam()); }
    std::string gfa() const { return golden::path(run().at("graph").get<std::string>()); }
    std::string reference() const { return run().at("reference").get<std::string>(); }
    std::string source() const { return run().at("source").get<std::string>(); }
    std::pair<Store, bool> build(const std::string& output) const {
        GraphOptions options;
        options.reference = reference();
        return prepare_coordinates(read_gfa(gfa(), options), gfa(), output, reference(), source(), false);
    }
};

std::string label(const ::testing::TestParamInfo<size_t>& info) { return "run" + std::to_string(info.param); }

std::vector<size_t> run_numbers() {
    std::vector<size_t> numbers(golden::json("coordinates.json").size());
    for (size_t i = 0; i < numbers.size(); i++) numbers[i] = i;
    return numbers;
}

}  // namespace

TEST_P(Coordinate, MatchesGlinkV2) {
    SCOPED_TRACE(run().at("graph").get<std::string>() + " on " + reference() + " from " + source());
    TempDir tmp;
    if (run().contains("error")) {
        try {
            build(tmp.file("c.h5"));
            FAIL() << "GlinkV2 refuses this with: " << run().at("error").get<std::string>();
        } catch (const Error& e) {
            EXPECT_EQ(std::string(e.what()), run().at("error").get<std::string>());
        }
        return;
    }
    auto [store, reused] = build(tmp.file("c.h5"));
    EXPECT_FALSE(reused);
    const Json& arrays = run().at("arrays");
    for (auto& [name, values] : arrays.items()) {
        if (name == "chromosomes") {
            EXPECT_EQ(store.strings(name), golden::strings(values));
        } else {
            EXPECT_EQ(store.ints(name), golden::ints(values)) << name;
        }
    }
    std::vector<std::string> names = store.integer_arrays();
    names.push_back("chromosomes");
    for (const std::string& name : names) EXPECT_TRUE(arrays.contains(name)) << name << " is not in GlinkV2's arrays";
    EXPECT_EQ(store.metadata().at("summary"), run().at("summary"));
    for (const char* key : {"links", "reference_nodes", "reference_source", "nodes"})
        EXPECT_EQ(store.metadata().at(key), run().at(key)) << key;
}

INSTANTIATE_TEST_SUITE_P(GlinkV2, Coordinate, ::testing::ValuesIn(run_numbers()), label);

TEST(CoordinateStore, SecondRunReusesIt) {
    TempDir tmp;
    const std::string gfa = golden::path("graph.gfa");
    const Graph graph = read_gfa(gfa, GraphOptions{});
    EXPECT_FALSE(prepare_coordinates(graph, gfa, tmp.file("c.h5"), "CHM13", "auto", false).second);
    EXPECT_TRUE(prepare_coordinates(graph, gfa, tmp.file("c.h5"), "CHM13", "auto", false).second);
    EXPECT_FALSE(prepare_coordinates(graph, gfa, tmp.file("c.h5"), "CHM13", "auto", true).second);
}

TEST(CoordinateStore, AnotherReferenceIsAnError) {
    TempDir tmp;
    const std::string gfa = golden::path("graph.gfa");
    const Graph graph = read_gfa(gfa, GraphOptions{});
    prepare_coordinates(graph, gfa, tmp.file("c.h5"), "CHM13", "auto", false);
    EXPECT_THROW(prepare_coordinates(graph, gfa, tmp.file("c.h5"), "GRCh38", "auto", false), Error);
}

TEST(CoordinateStore, LookupNamesMissingNodes) {
    TempDir tmp;
    const std::string gfa = golden::path("graph.gfa");
    const Graph graph = read_gfa(gfa, GraphOptions{});
    Store store = prepare_coordinates(graph, gfa, tmp.file("c.h5"), "CHM13", "auto", false).first;
    try {
        coordinates_of(store, {graph.ids[0], 999999});
        FAIL();
    } catch (const Error& e) {
        EXPECT_EQ(std::string(e.what()), "Node coordinate map lacks IDs: [999999]");  // as GlinkV2's lookup
    }
}
