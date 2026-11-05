// What the tests share: reading what tests/make_golden.py recorded from GlinkV2, running glink, and
// the fixtures each test builds on.
#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#include <optional>

#include <highfive/highfive.hpp>

#include "core.hpp"

#ifndef GLINK_TEST_DIR
#error "GLINK_TEST_DIR must name tests/"
#endif
#ifndef GLINK_BIN_DIR
#error "GLINK_BIN_DIR must name the directory holding glink"
#endif
#ifndef GLINK_VG
#error "GLINK_VG must name the vg executable"
#endif

namespace golden {

using glink::Json;

// GLINK_GOLDEN_DIR points the tests at a larger recording made with make_golden.py --out.
inline std::string dir() {
    const char* chosen = std::getenv("GLINK_GOLDEN_DIR");
    return chosen ? chosen : std::string(GLINK_TEST_DIR) + "/golden";
}

inline std::string path(const std::string& relative) { return dir() + "/" + relative; }

inline std::string text(const std::string& relative) { return glink::read_text(path(relative)); }

inline Json json(const std::string& relative) { return Json::parse(text(relative)); }

template <class T>
std::vector<T> raw(const std::string& relative) {
    const std::string bytes = text(relative);
    if (bytes.size() % sizeof(T)) throw std::runtime_error(relative + " is not a whole number of items");
    std::vector<T> out(bytes.size() / sizeof(T));
    if (!out.empty()) std::memcpy(out.data(), bytes.data(), bytes.size());
    return out;
}

// Rows i, columns i.. of a symmetric n x n matrix, as the whole matrix.
template <class T>
std::vector<T> symmetric(const std::vector<T>& upper, size_t n) {
    if (upper.size() != n * (n + 1) / 2) throw std::runtime_error("a triangle of the wrong size");
    std::vector<T> out(n * n);
    size_t k = 0;
    for (size_t i = 0; i < n; i++)
        for (size_t j = i; j < n; j++) out[i * n + j] = out[j * n + i] = upper[k++];
    return out;
}

// A float as float.hex wrote it.
inline double real(const std::string& hex) { return std::strtod(hex.c_str(), nullptr); }

// Exact equality, with NaN equal to NaN: these outputs are meant to be bit-identical.
inline bool same(double a, double b) { return (std::isnan(a) && std::isnan(b)) || a == b; }

// Whether every value is exactly GlinkV2's, naming the first that is not and how many differ.
template <class A, class B>
::testing::AssertionResult exactly(const std::vector<A>& got, const std::vector<B>& expected) {
    if (got.size() != expected.size())
        return ::testing::AssertionFailure() << got.size() << " values, GlinkV2 has " << expected.size();
    size_t wrong = 0, first = 0;
    for (size_t i = 0; i < got.size(); i++)
        if (!same(double(got[i]), double(expected[i])) && !wrong++) first = i;
    if (!wrong) return ::testing::AssertionSuccess();
    return ::testing::AssertionFailure() << wrong << " of " << got.size() << " values differ; at " << first << " got "
                                         << glink::shortest(double(got[first])) << ", GlinkV2 has "
                                         << glink::shortest(double(expected[first]));
}

// Every value of a store's array, whatever its shape, row by row.
inline std::vector<int64_t> flat(const std::string& store, const std::string& name) {
    const HighFive::File file(store, HighFive::File::ReadOnly);
    const HighFive::DataSet data = file.getDataSet(name);
    std::vector<int64_t> out(data.getElementCount());
    if (!out.empty()) data.read_raw(out.data());
    return out;
}

// Whether a block table says what GlinkV2's does. Numbers GlinkV2 prints as Python floats ("1.0")
// GlinkV3 prints with the fewest digits that read back the same ("1"), so its two r² columns are
// compared as numbers, exactly, and everything else as text.
inline ::testing::AssertionResult same_blocks(const std::string& got, const std::string& v2) {
    const auto rows = [](const std::string& text) {
        std::vector<std::vector<std::string>> out;
        std::stringstream lines(text);
        std::string line;
        while (std::getline(lines, line)) {
            out.emplace_back();
            std::stringstream cells(line);
            std::string cell;
            while (std::getline(cells, cell, '\t')) out.back().push_back(cell);
            if (!line.empty() && line.back() == '\t') out.back().push_back("");
        }
        return out;
    };
    if (!got.empty() && got.back() != '\n') return ::testing::AssertionFailure() << "no final newline";
    const auto a = rows(got), b = rows(v2);
    if (a.size() != b.size()) return ::testing::AssertionFailure() << a.size() << " lines, GlinkV2 has " << b.size();
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i].size() != b[i].size()) return ::testing::AssertionFailure() << "line " << i << ": " << a[i].size() << " columns, GlinkV2 has " << b[i].size();
        for (size_t c = 0; c < a[i].size(); c++) {
            const bool numeric = i > 0 && c >= 8 && !a[i][c].empty() && !b[i][c].empty();
            const bool equal = numeric ? std::strtod(a[i][c].c_str(), nullptr) == std::strtod(b[i][c].c_str(), nullptr) : a[i][c] == b[i][c];
            if (!equal) return ::testing::AssertionFailure() << "line " << i << " column " << c << ": '" << a[i][c] << "', GlinkV2 has '" << b[i][c] << "'";
        }
    }
    return ::testing::AssertionSuccess();
}

inline std::vector<int64_t> ints(const Json& array) { return array.get<std::vector<int64_t>>(); }
inline std::vector<std::string> strings(const Json& array) { return array.get<std::vector<std::string>>(); }

// A fresh, empty directory removed with the test.
class TempDir {
public:
    TempDir() {
        std::string pattern = ::testing::TempDir() + "glink-XXXXXX";
        std::vector<char> name(pattern.begin(), pattern.end());
        name.push_back(0);
        if (!mkdtemp(name.data())) throw std::runtime_error("mkdtemp failed");
        path_ = name.data();
    }
    ~TempDir() { std::system(("rm -rf '" + path_ + "'").c_str()); }
    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;
    const std::string& path() const { return path_; }
    std::string file(const std::string& name) const { return path_ + "/" + name; }

private:
    std::string path_;
};

inline void write(const std::string& path, const std::string& content) {
    std::ofstream out(path, std::ios::binary);
    if (!(out << content)) throw std::runtime_error("cannot write " + path);
}

inline std::string quote(const std::string& s) {
    std::string out = "'";
    for (char c : s) out += c == '\'' ? std::string("'\\''") : std::string(1, c);
    return out + "'";
}

struct Run {
    int exit = -1;       // -1 when killed by a signal
    int signal = 0;
    std::string out, err;
    std::string report() const { return "exit " + std::to_string(exit) + "\n" + out + err; }
};

// Runs a program from build/bin, or `vg`, from `cwd`. Its own process, so a crash fails one test only.
inline Run run(const std::string& program, const std::vector<std::string>& args, const std::string& cwd) {
    TempDir capture;
    std::string command = "cd " + quote(cwd) + " && exec " + quote(program[0] == '/' ? program : std::string(GLINK_BIN_DIR) + "/" + program);
    for (const std::string& a : args) command += " " + quote(a);
    command += " > " + quote(capture.file("out")) + " 2> " + quote(capture.file("err")) + " < /dev/null";
    const int status = std::system(command.c_str());
    Run r;
    if (WIFEXITED(status)) r.exit = WEXITSTATUS(status);
    if (WIFSIGNALED(status)) r.signal = WTERMSIG(status);
    r.out = glink::read_text(capture.file("out"));
    r.err = glink::read_text(capture.file("err"));
    return r;
}

// A crash, a sanitizer report or an uncaught exception, as opposed to an error glink reports.
inline ::testing::AssertionResult clean(const Run& r) {
    if (r.signal) return ::testing::AssertionFailure() << "killed by signal " << r.signal << "\n" << r.err;
    for (const char* mark : {"AddressSanitizer", "LeakSanitizer", "runtime error:", "terminate called"})
        if (r.err.find(mark) != std::string::npos) return ::testing::AssertionFailure() << mark << "\n" << r.err;
    if (r.exit && r.err.empty()) return ::testing::AssertionFailure() << "exit " << r.exit << " with no message";
    return ::testing::AssertionSuccess();
}

// ---- the recorded cases ----------------------------------------------------------------------

inline const Json& manifest() {
    static const Json m = json("manifest.json");
    return m;
}

inline std::vector<std::string> case_names() {
    std::vector<std::string> names;
    for (auto& [name, spec] : manifest().at("cases").items()) names.push_back(name);
    return names;
}

inline std::vector<std::string> block_case_names() {
    std::vector<std::string> names;
    for (auto& [name, spec] : manifest().at("cases").items())
        if (spec.at("blocks").get<bool>()) names.push_back(name);
    return names;
}

inline const Json& spec_of(const std::string& name) { return manifest().at("cases").at(name); }

// Test names cannot hold anything but these.
inline std::string test_name(const ::testing::TestParamInfo<std::string>& info) { return info.param; }

// ---- what a case runs on ---------------------------------------------------------------------

// One graph with its coordinates and tree, built by glink once per test program.
class Built {
public:
    explicit Built(const std::string& name) {
        gfa_ = path(name + ".gfa");
        glink::GraphOptions options;
        graph_ = glink::read_gfa(gfa_, options);
        options.oriented = true;
        oriented_ = glink::read_gfa(gfa_, options);
        coordinates_.emplace(glink::prepare_coordinates(graph_, gfa_, dir_.file("coords.h5"), "CHM13", "auto", false).first);
        glink::TreeOptions tree;
        tree.output = dir_.file("tree.h5");
        tree.snarls = path(name + ".snarls.json");
        const glink::SnarlRun snarls = glink::run_snarls(graph_, gfa_, tree);
        tree_.emplace(glink::prepare_tree(graph_, gfa_, tree, &snarls).first);
    }
    const glink::Graph& graph(bool oriented = false) const { return oriented ? oriented_ : graph_; }
    const glink::Store& coordinates() const { return *coordinates_; }
    const glink::Store& tree() const { return *tree_; }

private:
    TempDir dir_;
    std::string gfa_;
    glink::Graph graph_, oriented_;
    std::optional<glink::Store> coordinates_, tree_;
};

inline const Built& built(const std::string& graph) {
    static std::map<std::string, std::unique_ptr<Built>> all;
    auto& entry = all[graph];
    if (!entry) entry = std::make_unique<Built>(graph);
    return *entry;
}

// A golden-relative path in a filter: "@nodes.txt" -> "@/golden/nodes.txt".
inline std::string absolute_filter(const std::string& spec) { return spec[0] == '@' ? "@" + path(spec.substr(1)) : spec; }

// The node matrix of a recorded case, through the library.
inline glink::NodeMatrix build_case(const Json& spec) {
    const Built& b = built(spec.at("graph").get<std::string>());
    glink::BuildOptions o;
    o.matrix_type = spec.at("matrix_type").get<std::string>();
    for (const Json& f : spec.at("filters")) {
        std::string error;
        const auto filter = glink::parse_filter(absolute_filter(f.get<std::string>()), error);
        if (!filter) throw std::runtime_error(error);
        o.filters.push_back(*filter);
    }
    if (o.filters.empty()) o.filters.push_back(glink::Filter{"all"});
    o.collapse = spec.at("collapse").get<bool>();
    o.assume_complete = spec.at("assume_complete").get<bool>();
    o.condition_on_parents = spec.at("condition").get<bool>();
    o.collapse_backbone = spec.at("collapse_backbone").get<bool>();
    o.broken_min_haplotypes = spec.at("broken").get<int>();
    o.drop_complements = spec.at("drop").get<bool>();
    if (!spec.at("samples").is_null()) o.samples = path(spec.at("samples").get<std::string>());
    if (spec.at("coords").get<bool>()) o.coordinates = &b.coordinates();
    if (spec.at("tree").get<bool>()) o.parents = &b.tree();
    return glink::build_matrix(b.graph(spec.at("split").get<bool>()), o);
}

// The analysis options of a recorded case on V3's command line, from the golden directory.
inline std::vector<std::string> analysis_args(const Json& spec) {
    std::vector<std::string> a = {"--matrix-type", spec.at("alias").is_null() ? spec.at("matrix_type").get<std::string>()
                                                                              : spec.at("alias").get<std::string>()};
    for (const Json& f : spec.at("filters")) { a.push_back("--filter"); a.push_back(f.get<std::string>()); }
    a.push_back(spec.at("collapse").get<bool>() ? "--collapse" : "--no-collapse");
    for (auto [flag, key] : std::vector<std::pair<const char*, const char*>>{
             {"--assume-complete", "assume_complete"}, {"--condition-on-parents", "condition"},
             {"--collapse-backbone-parents", "collapse_backbone"}, {"--split-inversions", "split"},
             {"--drop-complements", "drop"}, {"--neutral-complements", "neutral"}})
        if (spec.at(key).get<bool>()) a.push_back(flag);
    a.push_back("--broken-parent-min-haplotypes");
    a.push_back(std::to_string(spec.at("broken").get<int>()));
    a.push_back("--min-pairs");
    a.push_back(std::to_string(spec.at("min_pairs").get<int>()));
    a.push_back("--min-allele-count");
    a.push_back(std::to_string(spec.at("min_allele_count").get<int>()));
    if (!spec.at("samples").is_null()) { a.push_back("--samples"); a.push_back(spec.at("samples").get<std::string>()); }
    if (spec.at("blocks").get<bool>()) {
        a.push_back("--blocks");
        a.push_back("--blocks-max-kb");
        a.push_back(glink::shortest(spec.at("max_kb").get<double>()));
    }
    return a;
}

// The chromosome name of each row of a matrix.
inline std::vector<std::string> chromosome_names(const glink::NodeMatrix& m, const std::string& field) {
    std::vector<std::string> out;
    for (int64_t code : m.coordinates.at(field)) out.push_back(code < 0 ? "" : m.coordinate_chromosomes.at(size_t(code)));
    return out;
}

}  // namespace golden
