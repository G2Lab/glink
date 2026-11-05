// The .ld report layouts. GlinkV2 writes .npy, so here the reference is PLINK 1.9 itself: these
// expected texts are what PLINK v1.90b6.21 printed for
//   plink --vcf t.vcf --double-id --allow-extra-chr --r2 square --write-snplist
//   plink --vcf t.vcf --double-id --allow-extra-chr --r2 inter-chr dprime --ld-window-r2 0
// on three markers. GlinkV2's values in these layouts are compared in test_cli.cpp.
#include <limits>
#include <random>

#include "golden.hpp"

using namespace glink;
using golden::TempDir;

namespace {

const std::string PLINK_SQUARE =
    "1\t0.631821\t0.147929\n"
    "0.631821\t1\t0.347798\n"
    "0.147929\t0.347798\t1\n";
const std::string PLINK_SNPLIST = "a\nlongsnpname\nthird\n";
const std::string PLINK_TABLE =
    " CHR_A         BP_A         SNP_A  CHR_B         BP_B         SNP_B           R2           DP \n"
    "    21          100             a     21          205   longsnpname     0.539615     0.734585 \n"
    "    21          100             a     23          300         third     0.176704     0.420362 \n"
    "    21          205   longsnpname     23          300         third     0.213589     0.462157 \n";

NodeMatrix three_markers() {
    NodeMatrix m;
    m.node_ids = {1, 2, 3};
    m.coordinates["chrom"] = {0, 0, 1};
    m.coordinates["position"] = {100, 205, 300};
    m.coordinate_chromosomes = {"chr21", "chrX"};
    return m;
}

// A matrix whose row labels are these, with no coordinates: node ids of the wanted length.
NodeMatrix labelled(std::vector<int64_t> ids) {
    NodeMatrix m;
    m.node_ids = std::move(ids);
    return m;
}

std::string table_of(const NodeMatrix& m, const std::vector<float>& r2, const std::vector<double>& d) {
    TempDir tmp;
    write_dprime_ld(tmp.file("d.ld"), m, r2, d);
    return read_text(tmp.file("d.ld"));
}

}  // namespace

TEST(PlinkFormat, SquareMatchesPlink) {
    TempDir tmp;
    write_r2_ld(tmp.file("r2.ld"), {1.f, 0.631821f, 0.147929f, 0.631821f, 1.f, 0.347798f, 0.147929f, 0.347798f, 1.f}, {"a", "longsnpname", "third"});
    EXPECT_EQ(read_text(tmp.file("r2.ld")), PLINK_SQUARE);
    EXPECT_EQ(read_text(tmp.file("r2.snplist")), PLINK_SNPLIST);
}

// Numbers are written with to_chars, and PLINK's format is defined by %g.
TEST(PlinkFormat, NumbersAreExactlyPrintfG) {
    std::vector<float> values = {0.f, -0.f, 1.f, -1.f, 0.1f, 0.5f, 1e-4f, 1e-5f, 9.99999e-5f, 0.000123456f, 0.999999f, 0.9999995f, 0.99999949f,
                                 1.5e-10f, 999999.f, 999999.5f, 1e6f, 1234567.f, 1e20f, 1e-40f, 1.17549435e-38f, 3.4028235e38f,
                                 std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity()};
    std::mt19937 random(5);
    std::uniform_real_distribution<float> unit(0.f, 1.f), exponent(-44.f, 38.f);
    while (values.size() < 300 * 300) {
        values.push_back(unit(random));
        values.push_back(std::pow(10.f, exponent(random)));
    }
    values.resize(300 * 300);
    TempDir tmp;
    write_r2_ld(tmp.file("r2.ld"), values, std::vector<std::string>(300, "x"));
    std::string expected;
    char text[64];
    for (size_t i = 0; i < 300; i++)
        for (size_t j = 0; j < 300; j++) {
            std::snprintf(text, sizeof text, "%g", double(values[i * 300 + j]));
            expected += text;
            expected += j == 299 ? '\n' : '\t';
        }
    EXPECT_EQ(read_text(tmp.file("r2.ld")), expected);
    // and padded to the table's width
    const std::string table = table_of(labelled({1, 2}), {0.f, 1e-5f, 1e-5f, 0.f}, {0.0, 0.123456789, 0.123456789, 0.0});
    EXPECT_EQ(table.substr(table.find('\n') + 1), "     0            0    1      0            0    2        1e-05     0.123457 \n");
}

TEST(PlinkFormat, SquareWritesNan) {
    TempDir tmp;
    const float nan = std::nanf("");
    write_r2_ld(tmp.file("r2.ld"), {nan, -nan, -nan, 0.25f}, {"1", "2"});
    EXPECT_EQ(read_text(tmp.file("r2.ld")), "nan\tnan\nnan\t0.25\n");
}

TEST(PlinkFormat, TableMatchesPlink) {
    const NodeMatrix m = three_markers();
    const std::vector<float> r2 = {1, 0.539615f, 0.176704f, 0.539615f, 1, 0.213589f, 0.176704f, 0.213589f, 1};
    const std::vector<double> d = {1, 0.734585, 0.420362, 0.734585, 1, 0.462157, 0.420362, 0.462157, 1};
    // Node ids stand for PLINK's SNP names, so the names and the width of their columns differ
    // (IdColumnWidthsMatchPlink has those); every other field must be PLINK's.
    const std::vector<std::string> plink = split_ws(PLINK_TABLE), ours = split_ws(table_of(m, r2, d));
    ASSERT_EQ(ours.size(), plink.size());
    for (size_t i = 0; i < ours.size(); i++) {
        if (i >= 8 && ((i - 8) % 8 == 2 || (i - 8) % 8 == 5)) continue;  // the SNP names
        EXPECT_EQ(ours[i], plink[i]) << "field " << i;
    }
}

TEST(PlinkFormat, IdColumnWidthsMatchPlink) {
    // PLINK: 4 wide for names up to 4 characters, otherwise the longest name plus 2.
    const std::vector<float> r2(4, 0.5f);
    const std::vector<double> d(4, 0.5);
    const std::string narrow = table_of(labelled({1, 2}), r2, d);
    EXPECT_EQ(narrow.substr(0, narrow.find('\n')), " CHR_A         BP_A SNP_A  CHR_B         BP_B SNP_B           R2           DP ");
    EXPECT_EQ(narrow.substr(narrow.find('\n') + 1),
              "     0            0    1      0            0    2          0.5          0.5 \n");
    // The eleven-character name of PLINK's table gives its thirteen-wide column.
    const std::string wide = table_of(labelled({12345678901, 2}), r2, d);
    EXPECT_EQ(wide.substr(0, wide.find('\n')), PLINK_TABLE.substr(0, PLINK_TABLE.find('\n')));
    EXPECT_EQ(wide.substr(wide.find('\n') + 1),
              "     0            0   12345678901      0            0             2          0.5          0.5 \n");
}

TEST(PlinkFormat, ChromosomesAreNamedAsPlinkDoes) {
    const std::vector<std::pair<std::string, std::string>> names = {
        {"chr1", "1"}, {"CHR2", "2"}, {"Chr3", "3"}, {"21", "21"}, {"chrX", "23"}, {"x", "23"}, {"chrY", "24"}, {"XY", "25"},
        {"chrM", "26"}, {"MT", "26"}, {"chr99", "26"}, {"KI270", "KI270"}, {"hs1", "hs1"}, {"chr", "chr"}, {"chrUn_x", "chrUn_x"}};
    for (auto& [name, plink] : names) {
        NodeMatrix m;
        m.node_ids = {1, 2};
        m.coordinates["chrom"] = {0, -1};
        m.coordinates["position"] = {7, 0};
        m.coordinate_chromosomes = {name};
        const std::vector<std::string> fields = split_ws(table_of(m, std::vector<float>(4, 0.5f), std::vector<double>(4, 0.5)));
        ASSERT_EQ(fields.size(), 16u);
        EXPECT_EQ(fields[8], plink) << name;
        EXPECT_EQ(fields[11], "0") << "a node with no chromosome";
    }
}

TEST(PlinkFormat, UndefinedPairsStayInTheTable) {
    const float nan = std::nanf("");
    const double dnan = std::nan("");
    const std::string text = table_of(labelled({1, 2}), {nan, nan, nan, nan}, {dnan, 0.25, 0.25, dnan});
    EXPECT_EQ(text.substr(text.find('\n') + 1), "     0            0    1      0            0    2          nan         0.25 \n");
}

TEST(BinaryFormat, IsRawFloat32RowByRow) {
    TempDir tmp;
    const float nan = std::nanf("");
    write_ld_binary(tmp.file("f.ld.bin"), std::vector<float>{1.f, 0.25f, nan, 0.5f}, {"a", "b"});
    write_ld_binary(tmp.file("d.ld.bin"), std::vector<double>{1.0, -0.25, std::nan(""), 0.5}, {"a", "b"});
    EXPECT_EQ(read_text(tmp.file("f.snplist")), "a\nb\n");
    EXPECT_EQ(read_text(tmp.file("d.snplist")), "a\nb\n");
    const auto floats = [](const std::string& path) {
        const std::string bytes = read_text(path);
        std::vector<float> out(bytes.size() / sizeof(float));
        std::memcpy(out.data(), bytes.data(), bytes.size());
        return out;
    };
    EXPECT_TRUE(golden::exactly(floats(tmp.file("f.ld.bin")), std::vector<float>{1.f, 0.25f, nan, 0.5f}));
    EXPECT_TRUE(golden::exactly(floats(tmp.file("d.ld.bin")), std::vector<float>{1.f, -0.25f, nan, 0.5f}));
    EXPECT_THROW(write_ld_binary(tmp.file("x.ld"), std::vector<float>{1.f}, {"a"}), Error);
    EXPECT_THROW(write_ld_binary(tmp.file("x.ld.bin"), std::vector<float>{1.f, 1.f}, {"a"}), Error);
    EXPECT_THROW(write_ld_binary(tmp.file("x.ld.bin"), std::vector<float>{1.f}, {"a", "b"}), Error);
}

TEST(PlinkFormat, RefusesWhatItCannotWrite) {
    TempDir tmp;
    EXPECT_THROW(write_r2_ld(tmp.file("r2.txt"), {1.f}, {"a"}), Error);
    EXPECT_THROW(write_r2_ld(tmp.file("r2.ld"), {1.f, 1.f, 1.f}, {"a"}), Error);
    EXPECT_THROW(write_r2_ld(tmp.file("r2.ld"), {1.f}, {"a", "b"}), Error);
    const NodeMatrix one = labelled({1});
    EXPECT_THROW(write_dprime_ld(tmp.file("d.txt"), one, {1.f}, {1.0}), Error);
    EXPECT_THROW(write_dprime_ld(tmp.file("d.ld"), one, {1.f}, {1.0, 1.0}), Error);
    EXPECT_THROW(write_dprime_ld(tmp.file("d.ld"), labelled({1, 2}), {1.f}, {1.0}), Error);
    EXPECT_THROW(write_r2_ld(tmp.file("absent/r2.ld"), {1.f}, {"a"}), Error);
}

TEST(PlinkFormat, CommandLineRefusesOtherFileNames) {
    TempDir tmp;
    for (const char* option : {"--r2-out", "--dprime-out"}) {
        const golden::Run r = golden::run("glink", {"--gfa", "graph.gfa", "--no-conditioning-tree", "--no-node-coordinates", "--assume-complete",
                                                    "--matrix-out", tmp.file("m"), option, tmp.file("x.txt")}, golden::dir());
        EXPECT_EQ(r.exit, 2) << option << r.report();
    }
    EXPECT_EQ(golden::run("glink", {"--gfa", "graph.gfa", "--no-conditioning-tree", "--no-node-coordinates", "--assume-complete",
                                    "--r2-out", tmp.file("x.ld"), "--dprime-out", tmp.file("x.ld")}, golden::dir()).exit, 2);
}
