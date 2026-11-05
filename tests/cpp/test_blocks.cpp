// Gabriel block tables and block comparison, against GlinkV2's.
#include "golden.hpp"

using namespace glink;
using golden::TempDir;

namespace {

class Blocks : public ::testing::TestWithParam<std::string> {
protected:
    std::string table(const NodeMatrix& m, const std::vector<float>* r2, int threads) const {
        const Json& spec = golden::spec_of(GetParam());
        TempDir tmp;
        const auto& positions = m.coordinates.at("position");
        const BlockResult result = build_blocks(r2, positions, m.values, m.matrix_type == "phased" ? 1 : 2,
                                                int64_t(spec.at("max_kb").get<double>() * 1000), threads, spec.at("neutral").get<bool>());
        write_blocks(result, m.row_labels(), positions, tmp.file("blocks.tsv"));
        return read_text(tmp.file("blocks.tsv"));
    }
    std::vector<float> r2(const NodeMatrix& m) const {
        const Json& spec = golden::spec_of(GetParam());
        return compute_r2(m, spec.at("min_pairs").get<int>(), spec.at("min_allele_count").get<int>(), 1);
    }
};

}  // namespace

TEST_P(Blocks, TableMatchesGlinkV2) {
    const NodeMatrix m = golden::build_case(golden::spec_of(GetParam()));
    const bool has_r2 = golden::spec_of(GetParam()).at("r2").get<bool>();
    const std::vector<float> values = has_r2 ? r2(m) : std::vector<float>{};
    EXPECT_TRUE(golden::same_blocks(table(m, has_r2 ? &values : nullptr, 1), golden::text("cases/" + GetParam() + "/blocks.tsv")));
}

TEST_P(Blocks, ThreadsChangeNothing) {
    const NodeMatrix m = golden::build_case(golden::spec_of(GetParam()));
    const bool has_r2 = golden::spec_of(GetParam()).at("r2").get<bool>();
    const std::vector<float> values = has_r2 ? r2(m) : std::vector<float>{};
    EXPECT_TRUE(golden::same_blocks(table(m, has_r2 ? &values : nullptr, 8), golden::text("cases/" + GetParam() + "/blocks.tsv")));
}

INSTANTIATE_TEST_SUITE_P(GlinkV2, Blocks, ::testing::ValuesIn(golden::block_case_names()), golden::test_name);

TEST(BuildBlocks, RefusesInconsistentInputs) {
    const NodeMatrix m = golden::build_case(golden::spec_of("gabriel_phased"));
    const auto& positions = m.coordinates.at("position");
    const std::vector<int64_t> too_few(positions.begin(), positions.end() - 1);
    EXPECT_THROW(build_blocks(nullptr, too_few, m.values, 1, 1000, 1), Error);
    const std::vector<float> wrong_size(4, 0.f);
    EXPECT_THROW(build_blocks(&wrong_size, positions, m.values, 1, 1000, 1), Error);
    EXPECT_THROW(build_blocks(nullptr, positions, m.values, 3, 1000, 1), Error);
    EXPECT_THROW(build_blocks(nullptr, positions, m.values, 1, -1, 1), Error);
    EXPECT_THROW(build_blocks(nullptr, positions, m.values, 1, 1000, 0), Error);
}

TEST(CompareBlocks, ScoresMatchGlinkV2) {
    TempDir tmp;
    const std::string a = tmp.file("a"), b = tmp.file("b");
    const Json all = golden::json("compare.json");
    ASSERT_GT(all.at("tables").size(), 10u);
    size_t trial = 0;
    for (const Json& row : all.at("tables")) {
        SCOPED_TRACE("table " + std::to_string(trial++));
        golden::write(a, row.at("a").get<std::string>());
        golden::write(b, row.at("b").get<std::string>());
        const Json got = compare_block_files(a, b, row.at("region").get<std::string>());
        for (const char* key : {"region", "region_bp", "a_blocks", "b_blocks", "a_bp", "b_bp", "shared_bp"})
            EXPECT_EQ(got.at(key), row.at(key)) << key;
    }
}

// The tool itself, including its verdict on what GlinkV2 refuses.
TEST(CompareBlocksTool, WritesScores) {
    TempDir tmp;
    const Json row = golden::json("compare.json").at("tables").at(1);
    golden::write(tmp.file("a"), row.at("a").get<std::string>());
    golden::write(tmp.file("b"), row.at("b").get<std::string>());
    const auto region = row.at("region").get<std::string>();
    const golden::Run printed = golden::run("glink-compare-blocks", {"a", "b", "--region", region}, tmp.path());
    ASSERT_EQ(printed.exit, 0) << printed.report();
    const golden::Run saved = golden::run("glink-compare-blocks", {"a", "b", "--region", region, "--out", "s.json"}, tmp.path());
    ASSERT_EQ(saved.exit, 0) << saved.report();
    EXPECT_EQ(read_text(tmp.file("s.json")), printed.out);
    const Json scores = Json::parse(printed.out);
    EXPECT_EQ(scores.at("a"), "a");
    EXPECT_EQ(scores.at("b"), "b");
    for (const char* key : {"region", "region_bp", "a_blocks", "b_blocks", "a_bp", "b_bp", "shared_bp"})
        EXPECT_EQ(scores.at(key), row.at(key)) << key;
}

TEST(CompareBlocksTool, BadTablesGetGlinkV2sVerdict) {
    // Known bug, pinned so that fixing it fails here: GlinkV3 skips a row with too few columns and
    // carries on, where GlinkV2 raises KeyError. A truncated table is scored as if it were whole.
    const std::map<std::string, int> known_bugs = {{"glink_short_row", 0}};
    const Json all = golden::json("compare.json");
    for (const Json& row : all.at("bad_tables")) {
        const std::string name = row.at("name").get<std::string>();
        SCOPED_TRACE(name);
        TempDir tmp;
        golden::write(tmp.file("bad"), row.at("text").get<std::string>());
        golden::write(tmp.file("good"), row.at("good").get<std::string>());
        const golden::Run r = golden::run("glink-compare-blocks", {"bad", "good", "--region", "chr1:1-200"}, tmp.path());
        EXPECT_TRUE(golden::clean(r));
        const auto known = known_bugs.find(name);
        EXPECT_EQ(r.exit, known != known_bugs.end() ? known->second : row.at("exit").get<int>()) << r.report();
    }
}

TEST(CompareBlocksTool, BadRegionsGetGlinkV2sVerdict) {
    TempDir tmp;
    golden::write(tmp.file("good"), " CHR BP1 BP2 KB NSNPS SNPS\n 1 1 100 0.1 2 a|b\n");
    const Json all = golden::json("compare.json");
    for (const Json& row : all.at("bad_regions")) {
        const std::string region = row.at("region").get<std::string>();
        SCOPED_TRACE(region);
        const golden::Run r = golden::run("glink-compare-blocks", {"good", "good", "--region", region}, tmp.path());
        EXPECT_TRUE(golden::clean(r));
        EXPECT_EQ(r.exit, row.at("exit").get<int>()) << r.report();
    }
}

TEST(CompareBlocksTool, MissingRegionIsAUsageError) {
    TempDir tmp;
    golden::write(tmp.file("good"), " CHR BP1 BP2 KB NSNPS SNPS\n 1 1 100 0.1 2 a|b\n");
    EXPECT_EQ(golden::run("glink-compare-blocks", {}, tmp.path()).exit, 2);
    EXPECT_EQ(golden::run("glink-compare-blocks", {"good", "good"}, tmp.path()).exit, 2);
}
