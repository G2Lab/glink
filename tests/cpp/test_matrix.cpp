// The node matrix of every recorded case, and filter parsing, against GlinkV2's build_matrix and node_filter.
#include "golden.hpp"

using namespace glink;

namespace {

class Matrix : public ::testing::TestWithParam<std::string> {
protected:
    Json expected() const { return golden::json("cases/" + GetParam() + "/matrix.json"); }
    NodeMatrix build() const { return golden::build_case(golden::spec_of(GetParam())); }
};

}  // namespace

TEST_P(Matrix, ValuesMatchGlinkV2) {
    const Json e = expected();
    const NodeMatrix m = build();
    ASSERT_EQ(m.values.rows, e.at("rows").get<size_t>());
    ASSERT_EQ(m.values.cols, e.at("cols").get<size_t>());
    EXPECT_TRUE(golden::exactly(m.values.v, golden::raw<int8_t>("cases/" + GetParam() + "/values.i8")));
    EXPECT_EQ(m.matrix_type, e.at("matrix_type").get<std::string>());
    EXPECT_EQ(m.columns, golden::strings(e.at("columns")));
    EXPECT_EQ(m.density(), golden::real(e.at("density").get<std::string>()));
}

TEST_P(Matrix, RowsMatchGlinkV2) {
    const Json e = expected();
    const NodeMatrix m = build();
    EXPECT_EQ(m.node_ids, golden::ints(e.at("node_ids")));
    EXPECT_EQ(m.input_node_ids, golden::ints(e.at("input_node_ids")));
    EXPECT_EQ(m.row_labels(), golden::strings(e.at("labels")));
    EXPECT_EQ(m.collapsed, e.contains("pattern_of"));
    if (e.contains("pattern_of")) {
        EXPECT_EQ(m.pattern_of, golden::ints(e.at("pattern_of")));
    }
    EXPECT_EQ(m.oriented, e.contains("node_orientation"));
    if (e.contains("node_orientation")) {
        EXPECT_EQ(std::vector<int64_t>(m.node_orientation.begin(), m.node_orientation.end()), golden::ints(e.at("node_orientation")));
        EXPECT_EQ(std::vector<int64_t>(m.input_node_orientation.begin(), m.input_node_orientation.end()),
                  golden::ints(e.at("input_node_orientation")));
    }
    EXPECT_EQ(m.complement_pairs, e.contains("complement_pairs") ? golden::ints(e.at("complement_pairs")) : std::vector<int64_t>{});
}

TEST_P(Matrix, ParentsMatchGlinkV2) {
    const Json e = expected();
    const NodeMatrix m = build();
    EXPECT_EQ(m.parent_node_ids, e.contains("parent_node_ids") ? golden::ints(e.at("parent_node_ids")) : std::vector<int64_t>{});
}

TEST_P(Matrix, CoordinatesMatchGlinkV2) {
    const Json e = expected();
    const NodeMatrix m = build();
    EXPECT_EQ(m.coordinates.empty(), !e.contains("coordinates"));
    if (!e.contains("coordinates")) return;
    for (auto& [name, values] : e.at("coordinates").items()) {
        if (ends_with(name, "chrom")) {
            EXPECT_EQ(golden::chromosome_names(m, name), golden::strings(values)) << name;
        } else {
            ASSERT_TRUE(m.coordinates.count(name)) << name;
            EXPECT_EQ(m.coordinates.at(name), golden::ints(values)) << name;
        }
    }
}

TEST_P(Matrix, CallabilityMatchesGlinkV2) {
    const Json e = expected().at("callability");
    const NodeMatrix m = build();
    EXPECT_EQ(m.callability.at("absence_model"), e.at("absence_model"));
    EXPECT_EQ(m.callability.contains("complements"), e.contains("complements"));
    if (e.contains("complements")) {
        EXPECT_EQ(m.callability.at("complements").at("dropped"), e.at("complements").at("dropped"));
    }
    EXPECT_EQ(m.callability.contains("orientation"), e.contains("orientation"));
    if (e.contains("orientation")) {
        for (auto& [key, value] : e.at("orientation").items()) EXPECT_EQ(m.callability.at("orientation").at(key), value) << key;
    }
    EXPECT_EQ(m.callability.contains("backbone_collapse"), e.contains("backbone_collapse"));
    if (e.contains("backbone_collapse")) {
        for (auto& [key, value] : e.at("backbone_collapse").items())
            if (key != "rule") {
                EXPECT_EQ(m.callability.at("backbone_collapse").at(key), value) << key;
            }
    }
}

INSTANTIATE_TEST_SUITE_P(GlinkV2, Matrix, ::testing::ValuesIn(golden::case_names()), golden::test_name);

TEST(Filters, AcceptAndRejectAsGlinkV2) {
    // Known differences, each pinned so that changing it fails here.
    const std::map<std::string, Json> known = {
        // Python's float() allows spaces around a number; GlinkV3 allows them only after it.
        {"minor: 0.1", {{"error", "'minor: 0.1' is not minor:FREQUENCY"}}},
        {" 0.1:0.9", {{"error", "' 0.1:0.9' is not a LOW:HIGH band"}}},
        // Python prints -1.0 and 2.0; GlinkV3 prints the fewest digits, as in its block tables.
        {"-1:0.5", {{"error", "band must satisfy 0 <= LOW < HIGH <= 1, got -1:0.5"}}},
        {"0:2", {{"error", "band must satisfy 0 <= LOW < HIGH <= 1, got 0:2"}}},
        // Path('') is '.'. Either way the run then fails, reading it.
        {"@", {{"kind", "list"}, {"path", ""}}},
    };
    const Json all = golden::json("filters.json");
    for (const Json& recorded : all) {
        const std::string spec = recorded.at("spec").get<std::string>();
        SCOPED_TRACE("--filter '" + spec + "'");
        const Json& row = known.count(spec) ? known.at(spec) : recorded;
        std::string error;
        const auto filter = parse_filter(spec, error);
        if (row.contains("error")) {
            EXPECT_FALSE(filter.has_value());
            std::string message = row.at("error").get<std::string>();
            // GlinkV3 dropped the v2 preset, on purpose, so does not list it.
            const size_t preset = message.find("all, v2, ");
            if (preset != std::string::npos) message.erase(preset + 5, 4);
            EXPECT_EQ(error, message);
            continue;
        }
        EXPECT_TRUE(filter.has_value()) << error;
        if (!filter) continue;
        const std::string kind = row.at("kind").get<std::string>();
        EXPECT_EQ(filter->kind, kind);
        if (kind == "band") {
            EXPECT_EQ(filter->low, golden::real(row.at("low").get<std::string>()));
            EXPECT_EQ(filter->high, golden::real(row.at("high").get<std::string>()));
        } else if (kind == "minor") {
            EXPECT_EQ(filter->value, golden::real(row.at("value").get<std::string>()));
        } else if (kind == "common") {
            EXPECT_EQ(filter->value, golden::real(row.at("value").get<std::string>()));
            EXPECT_EQ(filter->people, row.at("people").get<int>());
        } else if (kind == "list") {
            EXPECT_EQ(filter->path, row.at("path").get<std::string>());
        }
    }
}

TEST(MatrixStore, HoldsTheMatrix) {
    golden::TempDir tmp;
    const NodeMatrix m = golden::build_case(golden::spec_of("split_drop_plain"));
    write_matrix(m, tmp.file("m"));
    const Store store(tmp.file("m"));
    EXPECT_EQ(store.ints("node_ids"), m.node_ids);
    EXPECT_EQ(store.ints("input_node_ids"), m.input_node_ids);
    EXPECT_EQ(store.ints("pattern_of"), m.pattern_of);
    EXPECT_EQ(store.strings("columns"), m.columns);
    EXPECT_EQ(store.shape("values"), (std::vector<size_t>{m.values.rows, m.values.cols}));
    EXPECT_EQ(golden::flat(tmp.file("m"), "values"), std::vector<int64_t>(m.values.v.begin(), m.values.v.end()));
    EXPECT_EQ(store.metadata().at("matrix_type"), m.matrix_type);
    EXPECT_EQ(store.metadata().at("callability"), m.callability);
}
