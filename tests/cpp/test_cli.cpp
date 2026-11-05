// The glink command line from end to end, against GlinkV2's command line on the same inputs: what it
// writes, what it reports, and what it refuses.
#include "golden.hpp"

using namespace glink;
using golden::TempDir;

namespace {

// PLINK's numbers: %g, six significant digits.
std::string number(double value, int width = 0) {
    char text[64];
    if (std::isnan(value)) std::snprintf(text, sizeof text, "%*s", width, "nan");
    else std::snprintf(text, sizeof text, "%*g", width, value);
    return text;
}

std::string pad(const std::string& s, size_t width) { return s.size() >= width ? s : std::string(width - s.size(), ' ') + s; }

std::string plink_chromosome(const std::string& name) { return starts_with(name, "chr") ? name.substr(3) : name; }

// r2.ld as GlinkV2's values would read in PLINK's format.
std::string expected_r2_text(const std::vector<float>& r2, size_t n) {
    std::string out;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) out += (j ? "\t" : "") + number(r2[i * n + j]);
        out += "\n";
    }
    return out;
}

std::string expected_dprime_text(const std::vector<float>& r2, const std::vector<double>& d, const Json& matrix) {
    const auto labels = golden::strings(matrix.at("labels"));
    const size_t n = labels.size();
    size_t longest = 0;
    for (auto& l : labels) longest = std::max(longest, l.size());
    const size_t w = longest <= 4 ? 4 : longest + 2;
    std::vector<std::string> chrom(n, "0"), bp(n, "0");
    if (matrix.contains("coordinates")) {
        const auto names = golden::strings(matrix.at("coordinates").at("chrom"));
        const auto positions = golden::ints(matrix.at("coordinates").at("position"));
        for (size_t i = 0; i < n; i++) {
            if (!names[i].empty()) chrom[i] = plink_chromosome(names[i]);
            bp[i] = std::to_string(positions[i]);
        }
    }
    std::string out = pad("CHR_A", 6) + " " + pad("BP_A", 12) + " " + pad("SNP_A", w) + " " + pad("CHR_B", 6) + " " + pad("BP_B", 12) + " " +
                      pad("SNP_B", w) + " " + pad("R2", 12) + " " + pad("DP", 12) + " \n";
    for (size_t i = 0; i < n; i++)
        for (size_t j = i + 1; j < n; j++)
            out += pad(chrom[i], 6) + " " + pad(bp[i], 12) + " " + pad(labels[i], w) + " " + pad(chrom[j], 6) + " " + pad(bp[j], 12) + " " +
                   pad(labels[j], w) + " " + number(r2[i * n + j], 12) + " " + number(d[i * n + j], 12) + " \n";
    return out;
}

// What a run's report says about the results: paths, and how caches were found, are not results.
std::vector<std::string> report_lines(const std::string& out) {
    std::vector<std::string> lines;
    std::stringstream stream(out);
    std::string line;
    while (std::getline(stream, line)) {
        if (line.empty() || starts_with(line, "graph (") || starts_with(line, "index: ") || starts_with(line, "D′ ")) continue;
        line = line.substr(0, line.find(" -> "));
        lines.push_back(line);
    }
    return lines;
}

// GlinkV2 says "node parents (built)" where GlinkV3 says "conditioning tree (built)", and the same
// of coordinates that may have been reused.
std::string neutral_step(std::string line) {
    for (const char* step : {"node parents", "conditioning tree", "node coordinates"}) {
        if (!starts_with(line, step)) continue;
        const size_t colon = line.find(": ");
        return (std::string(step) == "node coordinates" ? "coordinates" : "parents") + line.substr(colon);
    }
    return line;
}

// GlinkV2 prepares parents before coordinates and GlinkV3 after, so their lines come in either order.
std::vector<std::string> setup_steps_sorted(std::vector<std::string> lines) {
    const auto setup = [](const std::string& line) { return neutral_step(line) != line; };
    const auto end = std::stable_partition(lines.begin(), lines.end(), setup);
    std::sort(lines.begin(), end, [](const std::string& a, const std::string& b) { return neutral_step(a) < neutral_step(b); });
    return lines;
}

// The within-block mean as a number: GlinkV2 prints round(x, 3), GlinkV3 three decimals.
void expect_same_report_line(const std::string& got, const std::string& v2) {
    const std::string mean = "mean within-block r² ";
    const size_t a = got.find(mean), b = v2.find(mean);
    if (a == std::string::npos || b == std::string::npos) {
        EXPECT_EQ(neutral_step(got), neutral_step(v2));
        return;
    }
    EXPECT_EQ(got.substr(0, a), v2.substr(0, b));
    const std::string x = got.substr(a + mean.size(), got.size() - a - mean.size() - 1), y = v2.substr(b + mean.size(), v2.size() - b - mean.size() - 1);
    if (y == "None") {
        EXPECT_EQ(x, "None");
    } else {
        EXPECT_EQ(std::strtod(x.c_str(), nullptr), std::strtod(y.c_str(), nullptr)) << got << "\n" << v2;
    }
}

class Cli : public ::testing::TestWithParam<std::string> {
protected:
    const Json& spec() const { return golden::spec_of(GetParam()); }
    std::string graph() const { return spec().at("graph").get<std::string>(); }

    // A whole run of this case: `steps` says how the graph, coordinates and tree are found.
    golden::Run glink_run(const TempDir& tmp, const std::vector<std::string>& steps, const std::vector<std::string>& extra = {}) const {
        std::vector<std::string> args = steps;
        for (const std::string& a : golden::analysis_args(spec())) args.push_back(a);
        for (const std::string& a : std::vector<std::string>{"--matrix-out", tmp.file("matrix")}) args.push_back(a);
        if (spec().at("r2").get<bool>())
            for (const std::string& a : std::vector<std::string>{"--r2-out", tmp.file("r2.ld"), "--dprime-out", tmp.file("dprime.ld")}) args.push_back(a);
        if (spec().at("blocks").get<bool>())
            for (const std::string& a : std::vector<std::string>{"--blocks-out", tmp.file("blocks.tsv")}) args.push_back(a);
        args.insert(args.end(), extra.begin(), extra.end());
        return golden::run("glink", args, golden::dir());
    }
    std::vector<std::string> explicit_steps(const TempDir& tmp) const {
        std::vector<std::string> args = {"--gfa", graph() + ".gfa"};
        if (spec().at("coords").get<bool>()) {
            args.push_back("--node-coordinates");
            args.push_back(tmp.file("coords"));
        } else args.push_back("--no-node-coordinates");
        if (spec().at("tree").get<bool>()) {
            for (const std::string& a : std::vector<std::string>{"--conditioning-tree", tmp.file("tree"), "--snarls", graph() + ".snarls.json"}) args.push_back(a);
        } else args.push_back("--no-conditioning-tree");
        return args;
    }

    void expect_outputs_match(const TempDir& tmp) const {
        const std::string name = "cases/" + GetParam() + "/";
        const Json matrix = golden::json(name + "matrix.json");
        const size_t n = matrix.at("rows").get<size_t>();
        const Store store(tmp.file("matrix"));
        EXPECT_TRUE(golden::exactly(golden::flat(tmp.file("matrix"), "values"), golden::raw<int8_t>(name + "values.i8")));
        EXPECT_EQ(store.shape("values"), (std::vector<size_t>{n, matrix.at("cols").get<size_t>()}));
        EXPECT_EQ(store.ints("node_ids"), golden::ints(matrix.at("node_ids")));
        EXPECT_EQ(store.ints("input_node_ids"), golden::ints(matrix.at("input_node_ids")));
        EXPECT_EQ(store.strings("columns"), golden::strings(matrix.at("columns")));
        if (matrix.contains("pattern_of")) {
            EXPECT_EQ(store.ints("pattern_of"), golden::ints(matrix.at("pattern_of")));
        }
        if (matrix.contains("parent_node_ids")) {
            EXPECT_EQ(store.ints("parent_node_ids"), golden::ints(matrix.at("parent_node_ids")));
        }
        if (matrix.contains("node_orientation")) {
            EXPECT_EQ(store.ints("node_orientation"), golden::ints(matrix.at("node_orientation")));
            EXPECT_EQ(store.ints("input_node_orientation"), golden::ints(matrix.at("input_node_orientation")));
        }
        if (matrix.contains("complement_pairs")) {
            EXPECT_EQ(golden::flat(tmp.file("matrix"), "complement_pairs"), golden::ints(matrix.at("complement_pairs")));
        }
        if (matrix.contains("coordinates")) {
            EXPECT_EQ(store.ints("coordinate_position"), golden::ints(matrix.at("coordinates").at("position")));
        }
        if (spec().at("r2").get<bool>()) {
            const auto r2 = golden::symmetric(golden::raw<float>(name + "r2.f32"), n);
            const auto d = golden::symmetric(golden::raw<double>(name + "dprime.f64"), n);
            EXPECT_EQ(read_text(tmp.file("r2.ld")), expected_r2_text(r2, n));
            std::string ids;
            for (const std::string& label : golden::strings(matrix.at("labels"))) ids += label + "\n";
            EXPECT_EQ(read_text(tmp.file("r2.snplist")), ids);
            EXPECT_EQ(read_text(tmp.file("dprime.ld")), expected_dprime_text(r2, d, matrix));
        }
        if (spec().at("blocks").get<bool>()) {
            EXPECT_TRUE(golden::same_blocks(read_text(tmp.file("blocks.tsv")), golden::text(name + "blocks.tsv")));
        }
    }
};

}  // namespace

TEST_P(Cli, OutputsMatchGlinkV2) {
    TempDir tmp;
    const golden::Run r = glink_run(tmp, explicit_steps(tmp));
    ASSERT_EQ(r.exit, 0) << r.report();
    ASSERT_TRUE(golden::clean(r));
    expect_outputs_match(tmp);
}

TEST_P(Cli, ReportMatchesGlinkV2) {
    TempDir tmp;
    const golden::Run r = glink_run(tmp, explicit_steps(tmp));
    ASSERT_EQ(r.exit, 0) << r.report();
    const auto got = setup_steps_sorted(report_lines(r.out));
    const auto v2 = setup_steps_sorted(golden::strings(golden::json("cases/" + GetParam() + "/matrix.json").at("console")));
    ASSERT_EQ(got.size(), v2.size()) << r.out;
    for (size_t i = 0; i < got.size(); i++) expect_same_report_line(got[i], v2[i]);
}

TEST_P(Cli, ThreadsChangeNothing) {
    TempDir tmp;
    const golden::Run r = glink_run(tmp, explicit_steps(tmp), {"--threads", "4"});
    ASSERT_EQ(r.exit, 0) << r.report();
    expect_outputs_match(tmp);
}

// .ld.bin keeps every digit, so here GlinkV2's values have to match exactly rather than at six digits.
TEST_P(Cli, BinaryOutputsMatchGlinkV2Exactly) {
    if (!spec().at("r2").get<bool>()) GTEST_SKIP() << "GlinkV2 wrote no r² for this case";
    TempDir tmp;
    std::vector<std::string> args = explicit_steps(tmp);
    for (const std::string& a : golden::analysis_args(spec())) args.push_back(a);
    for (const std::string& a : std::vector<std::string>{"--r2-out", tmp.file("r2.ld.bin"), "--dprime-out", tmp.file("d.ld.bin"), "--blocks-out", tmp.file("b.tsv")}) args.push_back(a);
    const golden::Run r = golden::run("glink", args, golden::dir());
    ASSERT_EQ(r.exit, 0) << r.report();
    const auto read_floats = [](const std::string& path) {
        const std::string bytes = read_text(path);
        std::vector<float> out(bytes.size() / sizeof(float));
        if (!out.empty()) std::memcpy(out.data(), bytes.data(), bytes.size());
        return out;
    };
    const std::string name = "cases/" + GetParam() + "/";
    const Json matrix = golden::json(name + "matrix.json");
    const size_t n = matrix.at("rows").get<size_t>();
    EXPECT_TRUE(golden::exactly(read_floats(tmp.file("r2.ld.bin")), golden::symmetric(golden::raw<float>(name + "r2.f32"), n)));
    std::vector<float> d;
    for (double v : golden::symmetric(golden::raw<double>(name + "dprime.f64"), n)) d.push_back(float(v));
    EXPECT_TRUE(golden::exactly(read_floats(tmp.file("d.ld.bin")), d));
    std::string ids;
    for (const std::string& label : golden::strings(matrix.at("labels"))) ids += label + "\n";
    EXPECT_EQ(read_text(tmp.file("r2.snplist")), ids);
    EXPECT_EQ(read_text(tmp.file("d.snplist")), ids);
}

TEST(CliOutputs, BinaryDprimeAloneNeedsNoR2) {
    TempDir tmp;
    const golden::Run r = golden::run("glink", {"--gfa", "graph.gfa", "--no-conditioning-tree", "--no-node-coordinates", "--assume-complete", "--matrix-type", "phased",
                                                "--matrix-out", tmp.file("m"), "--dprime-out", tmp.file("d.ld.bin")}, golden::dir());
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_TRUE(exists(tmp.file("d.ld.bin")));
    EXPECT_FALSE(exists(tmp.file("r2.ld.bin")));
}

// The caches a first run leaves give the same results as building them.
TEST_P(Cli, ReusedCachesGiveTheSameOutputs) {
    TempDir tmp;
    ASSERT_EQ(glink_run(tmp, explicit_steps(tmp)).exit, 0);
    const golden::Run r = glink_run(tmp, explicit_steps(tmp));
    ASSERT_EQ(r.exit, 0) << r.report();
    expect_outputs_match(tmp);
}

INSTANTIATE_TEST_SUITE_P(GlinkV2, Cli, ::testing::ValuesIn(golden::case_names()), golden::test_name);

TEST(CliInputs, AGraphStoreGivesTheSameOutputs) {
    TempDir tmp;
    const Json& spec = golden::spec_of("gabriel_parents");
    ASSERT_EQ(golden::run("glink", {"index", "--gfa", "graph.gfa", "--out-prefix", tmp.file("i"), "--vg", GLINK_VG}, golden::dir()).exit, 0);
    std::vector<std::string> args = {"--gfa", "graph.gfa", "--presence-bundle", tmp.file("i.glink/graph.h5"), "--node-coordinates", tmp.file("c"),
                                     "--conditioning-tree", tmp.file("t"), "--snarls", "graph.snarls.json",
                                     "--r2-out", tmp.file("r2.ld"), "--blocks-out", tmp.file("b.tsv")};
    for (const std::string& a : golden::analysis_args(spec)) args.push_back(a);
    const golden::Run r = golden::run("glink", args, golden::dir());
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_TRUE(golden::same_blocks(read_text(tmp.file("b.tsv")), golden::text("cases/gabriel_parents/blocks.tsv")));
}

TEST(CliInputs, ATreeStoreGivesTheSameOutputs) {
    TempDir tmp;
    const Json& spec = golden::spec_of("gabriel_parents");
    const std::vector<std::string> first = {"--gfa", "graph.gfa", "--no-node-coordinates", "--conditioning-tree", tmp.file("t"), "--snarls", "graph.snarls.json",
                                            "--assume-complete", "--matrix-out", tmp.file("m")};
    ASSERT_EQ(golden::run("glink", first, golden::dir()).exit, 0);
    std::vector<std::string> args = {"--gfa", "graph.gfa", "--node-coordinates", tmp.file("c"), "--snarl-tree", tmp.file("t"),
                                     "--r2-out", tmp.file("r2.ld"), "--blocks-out", tmp.file("b.tsv")};
    for (const std::string& a : golden::analysis_args(spec)) args.push_back(a);
    const golden::Run r = golden::run("glink", args, golden::dir());
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_TRUE(golden::same_blocks(read_text(tmp.file("b.tsv")), golden::text("cases/gabriel_parents/blocks.tsv")));
}

TEST(CliInputs, AGbzGivesTheSameParents) {
    TempDir tmp;
    const golden::Run r = golden::run("glink", {"--gfa", "graph.gbz.gfa", "--vg-graph", "graph.gbz", "--vg", GLINK_VG, "--no-node-coordinates", "--assume-complete",
                                                "--conditioning-tree", tmp.file("t"), "--matrix-out", tmp.file("m")}, golden::dir());
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_EQ(Store(tmp.file("t")).ints("node_snarl"), golden::ints(golden::json("tree_gbz.json").at("node_snarl")));
}

TEST(CliInputs, RebuildReplacesATreeOfOtherInputs) {
    TempDir tmp;
    const std::vector<std::string> common = {"--gfa", "graph.gfa", "--no-node-coordinates", "--assume-complete", "--conditioning-tree", tmp.file("t"),
                                             "--matrix-out", tmp.file("m"), "--vg", GLINK_VG};
    auto with = [&](const std::vector<std::string>& extra) {
        std::vector<std::string> args = common;
        args.insert(args.end(), extra.begin(), extra.end());
        return golden::run("glink", args, golden::dir());
    };
    ASSERT_EQ(with({"--backbone", "HG000"}).exit, 0);
    EXPECT_EQ(with({}).exit, 1) << "a tree rooted elsewhere must not be reused";
    const golden::Run rebuilt = with({"--rebuild-conditioning-tree"});
    ASSERT_EQ(rebuilt.exit, 0) << rebuilt.report();
    EXPECT_EQ(Store(tmp.file("t")).ints("node_snarl"), golden::ints(golden::json("tree.json").at("node_snarl")));
}

TEST(CliInputs, DefaultCachesAreUsed) {
    TempDir tmp;
    golden::write(tmp.file("graph.gfa"), golden::text("graph.gfa"));
    const std::vector<std::string> args = {"--gfa", "graph.gfa", "--vg", GLINK_VG, "--assume-complete", "--matrix-out", "m"};
    golden::Run r = golden::run("glink", args, tmp.path());
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_NE(r.out.find("graph (built)"), std::string::npos) << r.out;
    EXPECT_NE(r.out.find("node coordinates (built)"), std::string::npos) << r.out;
    EXPECT_NE(r.out.find("conditioning tree (built)"), std::string::npos) << r.out;
    r = golden::run("glink", args, tmp.path());
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_NE(r.out.find("graph (reused)"), std::string::npos) << r.out;
    EXPECT_NE(r.out.find("node coordinates (reused)"), std::string::npos) << r.out;
    EXPECT_NE(r.out.find("conditioning tree (reused)"), std::string::npos) << r.out;
}

TEST(CliInputs, DefaultBlocksFileSitsBesideR2) {
    TempDir tmp;
    const Json& spec = golden::spec_of("gabriel_phased");
    std::vector<std::string> args = {"--gfa", "graph.gfa", "--no-conditioning-tree", "--node-coordinates", tmp.file("c"), "--r2-out", tmp.file("r2.ld")};
    for (const std::string& a : golden::analysis_args(spec)) args.push_back(a);
    const golden::Run r = golden::run("glink", args, golden::dir());
    ASSERT_EQ(r.exit, 0) << r.report();
    EXPECT_TRUE(golden::same_blocks(read_text(tmp.file("r2.blocks.tsv")), golden::text("cases/gabriel_phased/blocks.tsv")));
}

TEST(CliInputs, CollapseFlagsTakeTheLastGiven) {
    TempDir tmp;
    const std::vector<std::string> base = {"--gfa", "graph.gfa", "--no-conditioning-tree", "--no-node-coordinates", "--assume-complete", "--matrix-type", "phased", "--matrix-out", tmp.file("m")};
    auto rows = [&](std::vector<std::string> flags) {
        std::vector<std::string> args = base;
        args.insert(args.end(), flags.begin(), flags.end());
        EXPECT_EQ(golden::run("glink", args, golden::dir()).exit, 0);
        return Store(tmp.file("m")).shape("values").at(0);
    };
    const size_t collapsed = golden::json("cases/phased_complete/matrix.json").at("rows").get<size_t>();
    const size_t all = golden::json("cases/phased_complete/matrix.json").at("input_node_ids").size();
    ASSERT_LT(collapsed, all);
    EXPECT_EQ(rows({}), collapsed);
    EXPECT_EQ(rows({"--collapse"}), collapsed);
    EXPECT_EQ(rows({"--no-collapse"}), all);
    EXPECT_EQ(rows({"--collapse", "--no-collapse"}), all);
    EXPECT_EQ(rows({"--no-collapse", "--collapse"}), collapsed);
}

TEST(CliInputs, CoordinateReferenceMatchesGlinkV2) {
    for (const Json& entry : golden::json("coordinates.json")) {
        if (entry.contains("error") || entry.at("graph") != "graph.gfa") continue;
        const std::string reference = entry.at("reference").get<std::string>();
        SCOPED_TRACE(reference);
        TempDir tmp;
        const golden::Run r = golden::run("glink", {"--gfa", "graph.gfa", "--no-conditioning-tree", "--assume-complete",
                                                    "--coordinate-reference", reference, "--coordinate-reference-source", entry.at("source").get<std::string>(),
                                                    "--node-coordinates", tmp.file("c"), "--matrix-out", tmp.file("m")}, golden::dir());
        ASSERT_EQ(r.exit, 0) << r.report();
        EXPECT_EQ(Store(tmp.file("c")).ints("position"), golden::ints(entry.at("arrays").at("position")));
    }
}

namespace {

// `OUT` in a recorded command line is a directory of the test's own. The recording was made in
// V3's spelling, apart from outputs, which only name files.
std::vector<std::string> spelled(const Json& args, const TempDir& tmp) {
    std::vector<std::string> out;
    for (const Json& a : args) {
        std::string word = a.get<std::string>();
        const size_t at = word.find("OUT");
        if (at != std::string::npos) word.replace(at, 3, tmp.path());
        out.push_back(word);
    }
    return out;
}

// What GlinkV2 said, without Python's name for it.
std::string v2_words(std::string message) {
    for (const char* prefix : {"ValueError: ", "glink: error: "})
        if (starts_with(message, prefix)) return message.substr(std::strlen(prefix));
    return starts_with(message, "--gfa ") ? message : "";  // the rest are Python's own, FileNotFoundError and so on
}

// Refusals GlinkV3 words differently. The result is the same; these are not results.
const std::map<std::string, std::string>& different_words() {
    static const std::map<std::string, std::string> words = {
        {"unknown_option", "argparse reports the missing --gfa first"},
        {"bad_filter", "does not list the v2 preset, which GlinkV3 dropped"},
        {"bad_type", "argparse lists the choices"},
        {"bad_reference_source", "argparse lists the choices"},
        {"condition_without_tree", "--coverage is gone, and parents are the conditioning tree"},
        {"condition_and_complete", "--coverage is gone, and parents are the conditioning tree"},
        {"no_tree_and_input", "parents are the conditioning tree"},
        {"blocks_without_coordinates", "only one block method, so it is not named"},
        {"blocks_unplaced", "only one block method, so it is not named"},
        {"blocks_presence", "only one block method, so it is not named"},
        {"two_snarl_sources", "GlinkV3 names all three sources"},
        {"list_not_numbers", "says which ID, not Python's int()"},
        {"snarls_twice", "names the snarl's nodes, not Python's tuple"},
        {"snarls_orphan", "names the snarl's number, not Python's tuple"},
    };
    return words;
}

class Refused : public ::testing::TestWithParam<std::string> {
protected:
    static const Json& lines() {
        static const Json j = golden::json("errors.json").at("command_lines");
        return j;
    }
    const Json& entry() const {
        for (const Json& e : lines())
            if (e.at("name") == GetParam()) return e;
        throw std::runtime_error(GetParam());
    }
};

std::vector<std::string> command_line_names() {
    std::vector<std::string> names;
    const Json all = golden::json("errors.json");
    for (const Json& e : all.at("command_lines")) names.push_back(e.at("name").get<std::string>());
    return names;
}

class BadGfa : public ::testing::TestWithParam<std::string> {};

std::vector<std::string> bad_gfa_names() {
    std::vector<std::string> names;
    const Json all = golden::json("errors.json");
    for (const Json& e : all.at("gfas")) names.push_back(e.at("name").get<std::string>());
    return names;
}

}  // namespace

// The same verdict as GlinkV2: success, a reported failure (1), or a usage error (2), and in the same words.
TEST_P(Refused, MatchesGlinkV2) {
    // Known differences, each pinned so that changing it fails here. GlinkV3 checks this while parsing the command line.
    // argparse also takes any unambiguous prefix of an option's name, and GlinkV3's parser does not.
    const std::map<std::string, int> other_exit = {{"drop_presence", 2}, {"abbreviated_option", 2}};
    TempDir tmp;
    const golden::Run r = golden::run("glink", spelled(entry().at("args"), tmp), golden::dir());
    const auto known = other_exit.find(GetParam());
    const int expected = known != other_exit.end() ? known->second : entry().at("exit").get<int>();
    EXPECT_TRUE(golden::clean(r));
    EXPECT_EQ(r.exit, expected) << "GlinkV2: " << entry().at("message").get<std::string>() << "\nGlinkV3: " << r.report();
    const std::string words = v2_words(entry().at("message").get<std::string>());
    if (words.empty()) return;
    if (different_words().count(GetParam())) {
        EXPECT_EQ(r.err.find(words), std::string::npos) << "the same words now, so it need not be listed as different: " << words;
    } else {
        EXPECT_NE(r.err.find(words), std::string::npos) << "GlinkV2: " << words << "\nGlinkV3: " << r.err;
    }
}

INSTANTIATE_TEST_SUITE_P(GlinkV2, Refused, ::testing::ValuesIn(command_line_names()), golden::test_name);

TEST_P(BadGfa, MatchesGlinkV2) {
    int expected = -1;
    std::string message;
    const Json all = golden::json("errors.json");
    for (const Json& e : all.at("gfas"))
        if (e.at("name") == GetParam()) {
            expected = e.at("exit").get<int>();
            message = e.at("message").get<std::string>();
        }
    TempDir tmp;
    const golden::Run r = golden::run("glink", {"--gfa", golden::path("bad/" + GetParam() + ".gfa"), "--no-conditioning-tree", "--no-node-coordinates", "--assume-complete",
                                                "--matrix-type", "phased", "--matrix-out", tmp.file("m")}, tmp.path());
    EXPECT_TRUE(golden::clean(r));
    EXPECT_EQ(r.exit, expected) << "GlinkV2: " << message << "\nGlinkV3: " << r.report();
}

INSTANTIATE_TEST_SUITE_P(GlinkV2, BadGfa, ::testing::ValuesIn(bad_gfa_names()), golden::test_name);
