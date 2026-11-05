// r², D′ and the engines behind them, against GlinkV2, exactly, at any thread count.
#include "dprime.hpp"
#include "golden.hpp"

using namespace glink;

namespace {

std::vector<float> r2_of(const std::string& name, const NodeMatrix& m, int threads) {
    const Json& spec = golden::spec_of(name);
    return compute_r2(m, spec.at("min_pairs").get<int>(), spec.at("min_allele_count").get<int>(), threads);
}

std::vector<uint8_t> unhex(const std::string& hex) {
    std::vector<uint8_t> out(hex.size() / 2);
    for (size_t i = 0; i < out.size(); i++) out[i] = uint8_t(std::stoi(hex.substr(2 * i, 2), nullptr, 16));
    return out;
}

Counts counts_of(const Json& row) {
    Counts c{};
    for (size_t i = 0; i < 9; i++) c[i] = row.at(0).at(i).get<uint16_t>();
    return c;
}

void expect_interval(const Interval& got, const Json& row) {
    EXPECT_TRUE(golden::same(got.dprime, golden::real(row.at(1).get<std::string>()))) << row.at(0) << " D′ " << got.dprime;
    EXPECT_TRUE(golden::same(got.lower, golden::real(row.at(2).get<std::string>()))) << row.at(0) << " lower " << got.lower;
    EXPECT_TRUE(golden::same(got.upper, golden::real(row.at(3).get<std::string>()))) << row.at(0) << " upper " << got.upper;
    EXPECT_EQ(got.category, row.at(4).get<int>()) << row.at(0);
}

const Json& engine() {
    static const Json j = golden::json("engine.json");
    return j;
}

class Ld : public ::testing::TestWithParam<std::string> {};

}  // namespace

TEST_P(Ld, R2MatchesGlinkV2Exactly) {
    if (!golden::spec_of(GetParam()).at("r2").get<bool>()) GTEST_SKIP() << "GlinkV2 wrote no r² for this case";
    const NodeMatrix m = golden::build_case(golden::spec_of(GetParam()));
    const auto expected = golden::symmetric(golden::raw<float>("cases/" + GetParam() + "/r2.f32"), m.values.rows);
    EXPECT_TRUE(golden::exactly(r2_of(GetParam(), m, 1), expected));
}

TEST_P(Ld, DprimeMatchesGlinkV2Exactly) {
    const NodeMatrix m = golden::build_case(golden::spec_of(GetParam()));
    const auto expected = golden::symmetric(golden::raw<double>("cases/" + GetParam() + "/dprime.f64"), m.values.rows);
    EXPECT_TRUE(golden::exactly(compute_dprime(m, 1), expected));
}

TEST_P(Ld, ThreadsChangeNothing) {
    const NodeMatrix m = golden::build_case(golden::spec_of(GetParam()));
    EXPECT_TRUE(golden::exactly(r2_of(GetParam(), m, 8), r2_of(GetParam(), m, 1)));
    EXPECT_TRUE(golden::exactly(compute_dprime(m, 8), compute_dprime(m, 1)));
}

INSTANTIATE_TEST_SUITE_P(GlinkV2, Ld, ::testing::ValuesIn(golden::case_names()), golden::test_name);

// Duplicates and allele flips share a fit inside compute_dprime. Each must still get its own marker's value, wherever it sits.
TEST(Dprime, DuplicatesAndAlleleFlipsKeepGlinkV2sValues) {
    const std::string name = "gabriel_dosage";
    NodeMatrix m = golden::build_case(golden::spec_of(name));
    const size_t n = m.values.rows, cols = m.values.cols;
    const auto v2 = golden::symmetric(golden::raw<double>("cases/" + name + "/dprime.f64"), n);
    size_t defined = 0;
    for (double x : v2) defined += !std::isnan(x);
    ASSERT_GT(defined, n) << "this case has no defined D′ between markers, so it proves nothing";
    // After the n markers: a copy of each, then each flipped, with its missing calls kept.
    std::vector<size_t> source;
    for (size_t k = 0; k < 3 * n; k++) source.push_back(k % n);
    Values grown;
    grown.rows = 3 * n;
    grown.cols = cols;
    grown.v.resize(grown.rows * cols);
    for (size_t r = 0; r < grown.rows; r++)
        for (size_t c = 0; c < cols; c++) {
            const int8_t x = m.values.row(source[r])[c];
            grown.row(r)[c] = r >= 2 * n && x >= 0 ? int8_t(2 - x) : x;
        }
    m.values = std::move(grown);
    std::vector<double> expected(9 * n * n);
    for (size_t i = 0; i < 3 * n; i++)
        for (size_t j = 0; j < 3 * n; j++) expected[i * 3 * n + j] = v2[source[i] * n + source[j]];
    for (int threads : {1, 8}) EXPECT_TRUE(golden::exactly(compute_dprime(m, threads), expected)) << "threads " << threads;
}

TEST(R2, RefusesImpossibleSupport) {
    const NodeMatrix m = golden::build_case(golden::spec_of("minor"));
    EXPECT_THROW(compute_r2(m, 1, 0, 1), Error);
    EXPECT_THROW(compute_r2(m, 2, -1, 1), Error);
}

TEST(Dprime, RefusesAPhasedMatrixOfDosages) {
    NodeMatrix m = golden::build_case(golden::spec_of("minor_dosage"));
    ASSERT_NE(std::find(m.values.v.begin(), m.values.v.end(), int8_t(2)), m.values.v.end());
    m.matrix_type = "phased";
    EXPECT_THROW(compute_dprime(m, 1), Error);
}

TEST(Engine, DosageConfidenceMatchesGlinkV2) {
    const Json& rows = engine().at("tables").at("dosage");
    ASSERT_GT(rows.size(), 1000u);
    for (const Json& row : rows) expect_interval(dosage_confidence(counts_of(row)), row);
}

TEST(Engine, GabrielConfidenceMatchesGlinkV2) {
    for (const Json& row : engine().at("tables").at("gabriel_diploid")) expect_interval(gabriel_confidence(counts_of(row), 2), row);
    for (const Json& row : engine().at("tables").at("gabriel_haploid")) expect_interval(gabriel_confidence(counts_of(row), 1), row);
}

TEST(Engine, DosageClassesMatchGlinkV2) {
    size_t tested = 0;
    for (const Json& p : engine().at("panels")) {
        if (!p.contains("classes")) continue;
        const int k = p.at("markers").get<int>(), width = p.at("width").get<int>();
        const auto packed = unhex(p.at("packed").get<std::string>());
        for (int threads : {1, 8}) {
            std::vector<uint8_t> got(size_t(k) * k);
            dosage_classes(packed.data(), k, width, got.data(), threads);
            EXPECT_EQ(got, p.at("classes").get<std::vector<uint8_t>>()) << p.at("name") << " threads " << threads;
        }
        tested++;
    }
    EXPECT_GT(tested, 0u);
}

TEST(Engine, DosageDprimeMatchesGlinkV2) {
    size_t tested = 0;
    for (const Json& p : engine().at("panels")) {
        if (!p.contains("dprime")) continue;
        const int k = p.at("markers").get<int>(), width = p.at("width").get<int>();
        const auto packed = unhex(p.at("packed").get<std::string>());
        std::vector<double> expected;
        for (const Json& v : p.at("dprime")) expected.push_back(golden::real(v.get<std::string>()));
        for (int threads : {1, 8}) {
            std::vector<float> got(size_t(k) * k);
            dosage_dprime(packed.data(), k, width, got.data(), threads);
            EXPECT_TRUE(golden::exactly(got, expected)) << p.at("name") << " threads " << threads;
        }
        tested++;
    }
    EXPECT_GT(tested, 0u);
}

TEST(Engine, GabrielBandsAndBlocksMatchGlinkV2) {
    for (const Json& p : engine().at("panels")) {
        const int n = p.at("markers").get<int>(), width = p.at("width").get<int>(), ploidy = p.at("ploidy").get<int>();
        const auto packed = unhex(p.at("packed").get<std::string>());
        const auto positions = golden::ints(p.at("positions"));
        for (const Json& s : p.at("spans")) {
            const int64_t span = s.at("span").get<int64_t>();
            const bool neutral = s.at("neutral").get<bool>();
            SCOPED_TRACE(p.at("name").get<std::string>() + " span " + std::to_string(span) + (neutral ? " neutral" : ""));
            std::vector<std::pair<int, int>> expected;
            for (const Json& b : s.at("blocks")) expected.emplace_back(b.at(0).get<int>(), b.at(1).get<int>());
            const auto bands_expected = unhex(s.at("bands").get<std::string>());
            for (int threads : {1, 8}) {
                std::vector<uint8_t> bands(size_t(n) * n);
                gabriel_bands(packed.data(), n, width, positions.data(), span, ploidy, bands.data(), threads, neutral);
                EXPECT_EQ(bands, bands_expected) << "threads " << threads;
                EXPECT_EQ(gabriel_select(bands.data(), positions.data(), n, span), expected) << "threads " << threads;
                EXPECT_EQ(gabriel_blocks(packed.data(), n, width, positions.data(), span, ploidy, threads, neutral), expected)
                    << "threads " << threads;
            }
        }
    }
}
