// Stores and paths: GlinkV3's own, with no GlinkV2 counterpart, so tested for what their headers promise.
#include "golden.hpp"

using namespace glink;
using golden::TempDir;

TEST(Store, RoundTrip) {
    TempDir tmp;
    const std::string path = tmp.file("sub/s.h5");
    {
        StoreWriter w(path);
        w.add("ints", std::vector<int64_t>{5, -1, 7});
        w.add("small", std::vector<uint8_t>{1, 0, 255});
        w.add("names", std::vector<std::string>{"a", "b#1"});
        const std::vector<uint8_t> table = {1, 2, 3, 4, 5, 6};
        w.add("table", size_t(2), size_t(3), table.data());
        w.set_metadata(Json{{"schema", 1}, {"note", "hi"}});
        w.publish();
    }
    const Store s(path);
    EXPECT_EQ(s.path(), path);
    EXPECT_EQ(s.ints("ints"), (std::vector<int64_t>{5, -1, 7}));
    EXPECT_EQ(s.ints("ints", {2, 0}), (std::vector<int64_t>{7, 5}));
    EXPECT_TRUE(s.ints("ints", {}).empty());
    EXPECT_EQ(s.ints("small"), (std::vector<int64_t>{1, 0, 255}));
    EXPECT_EQ(s.strings("names"), (std::vector<std::string>{"a", "b#1"}));
    EXPECT_EQ(s.shape("table"), (std::vector<size_t>{2, 3}));
    EXPECT_EQ(s.item_size("small"), 1u);
    EXPECT_EQ(s.item_size("ints"), 8u);
    std::vector<uint8_t> all(6), row(3);
    s.read_bytes("table", all.data());
    EXPECT_EQ(all, (std::vector<uint8_t>{1, 2, 3, 4, 5, 6}));
    s.read_rows("table", {1}, row.data());
    EXPECT_EQ(row, (std::vector<uint8_t>{4, 5, 6}));
    EXPECT_EQ(s.metadata().at("note"), "hi");
    EXPECT_TRUE(s.matches(Json{{"schema", 1}}));
    EXPECT_FALSE(s.matches(Json{{"schema", 2}}));
    EXPECT_FALSE(s.matches(Json{{"absent", 1}}));
    EXPECT_THROW(s.ints("absent"), Error);
}

TEST(Store, EmptyArraysRoundTrip) {
    TempDir tmp;
    {
        StoreWriter w(tmp.file("s.h5"));
        w.add("ints", std::vector<int64_t>{});
        w.add("names", std::vector<std::string>{});
        const std::vector<uint8_t> none;
        w.add("table", size_t(0), size_t(3), none.data());
        w.publish();
    }
    const Store s(tmp.file("s.h5"));
    EXPECT_TRUE(s.ints("ints").empty());
    EXPECT_TRUE(s.strings("names").empty());
    EXPECT_EQ(s.shape("table"), (std::vector<size_t>{0, 3}));
}

TEST(Store, NothingAppearsUntilPublished) {
    TempDir tmp;
    {
        StoreWriter w(tmp.file("s.h5"));
        w.add("ints", std::vector<int64_t>{1});
    }
    EXPECT_FALSE(exists(tmp.file("s.h5")));
    EXPECT_EQ(golden::run("/bin/ls", {"-A", tmp.path()}, tmp.path()).out, "");
}

TEST(Store, ReplacesAStoreButNothingElse) {
    TempDir tmp;
    for (int generation : {1, 2}) {
        StoreWriter w(tmp.file("s.h5"));
        w.add("ints", std::vector<int64_t>{generation});
        w.publish();
    }
    EXPECT_EQ(Store(tmp.file("s.h5")).ints("ints"), (std::vector<int64_t>{2}));
    golden::write(tmp.file("data.txt"), "precious");
    {
        StoreWriter w(tmp.file("data.txt"));
        w.add("ints", std::vector<int64_t>{1});
        EXPECT_THROW(w.publish(), Error);
    }
    EXPECT_EQ(read_text(tmp.file("data.txt")), "precious");
}

TEST(Store, AnythingElseIsRefused) {
    TempDir tmp;
    EXPECT_THROW(Store(tmp.file("absent.h5")), Error);
    EXPECT_THROW(Store(tmp.path()), Error);
    golden::write(tmp.file("text.h5"), "not a store");
    EXPECT_THROW(Store(tmp.file("text.h5")), Error);
    golden::write(tmp.file("empty.h5"), "");
    EXPECT_THROW(Store(tmp.file("empty.h5")), Error);
    // An HDF5 file, but not one of ours.
    {
        HighFive::File other(tmp.file("other.h5"), HighFive::File::Truncate);
        other.createDataSet("x", std::vector<int>{1});
    }
    EXPECT_THROW(Store(tmp.file("other.h5")), Error);
    // A store cut short.
    {
        StoreWriter w(tmp.file("whole.h5"));
        w.add("ints", std::vector<int64_t>(100000, 7));
        w.publish();
    }
    const std::string whole = read_text(tmp.file("whole.h5"));
    golden::write(tmp.file("cut.h5"), whole.substr(0, whole.size() / 2));
    EXPECT_THROW(Store(tmp.file("cut.h5")), Error);
}

TEST(CorruptInputs, AreCleanErrors) {
    for (const char* option : {"--node-coordinates", "--conditioning-tree", "--snarl-tree", "--presence-bundle"}) {
        SCOPED_TRACE(option);
        TempDir tmp;
        golden::write(tmp.file("bad"), "not a store");
        std::vector<std::string> args = {"--gfa", "graph.gfa", option, tmp.file("bad"), "--assume-complete", "--matrix-out", tmp.file("m")};
        if (std::string(option) != "--node-coordinates") args.push_back("--no-node-coordinates");
        if (std::string(option) == "--node-coordinates" || std::string(option) == "--presence-bundle") args.push_back("--no-conditioning-tree");
        const golden::Run r = golden::run("glink", args, golden::dir());
        EXPECT_EQ(r.exit, 1) << r.report();
        EXPECT_TRUE(golden::clean(r));
        EXPECT_FALSE(exists(tmp.file("m")));
    }
}

TEST(GraphStore, OfAnotherGraphIsRefused) {
    TempDir tmp;
    const std::string gfa = golden::path("graph.gfa");
    prepare_graph(gfa, GraphOptions{}, tmp.file("g.h5"));
    EXPECT_NO_THROW(check_graph(tmp.file("g.h5"), gfa, GraphOptions{}));
    EXPECT_THROW(check_graph(tmp.file("g.h5"), golden::path("plain.gfa"), GraphOptions{}), Error);
    GraphOptions other;
    other.backbone = "HG000";
    EXPECT_THROW(check_graph(tmp.file("g.h5"), gfa, other), Error);
    EXPECT_THROW(prepare_graph(gfa, GraphOptions{}, gfa), Error);
}

TEST(Paths, BehaveAsDocumented) {
    TempDir tmp;
    EXPECT_EQ(glink::basename("/a/b.txt"), "b.txt");
    EXPECT_EQ(glink::basename("b.txt"), "b.txt");
    EXPECT_EQ(glink::dirname("/a/b.txt"), "/a");
    EXPECT_EQ(glink::dirname("b.txt"), ".");
    EXPECT_EQ(glink::dirname("/b.txt"), "/");
    EXPECT_TRUE(ends_with("a.ld", ".ld"));
    EXPECT_FALSE(ends_with("ld", ".ld"));
    EXPECT_TRUE(starts_with("@x", "@"));
    EXPECT_EQ(split_ws("  a\tb\n c  "), (std::vector<std::string>{"a", "b", "c"}));
    EXPECT_TRUE(split_ws(" \n").empty());
    make_dirs(tmp.file("x/y/z"));
    EXPECT_TRUE(exists(tmp.file("x/y/z")));
    EXPECT_EQ(resolve(tmp.file("x/y/../y/absent")), resolve(tmp.file("x")) + "/y/absent");
    EXPECT_THROW(read_text(tmp.file("absent")), Error);
}

TEST(Numbers, FormatAsPython) {
    EXPECT_EQ(thousands(0), "0");
    EXPECT_EQ(thousands(999), "999");
    EXPECT_EQ(thousands(1000), "1,000");
    EXPECT_EQ(thousands(-1234567), "-1,234,567");
    EXPECT_EQ(fixed3(0.0005), "0.001");
    EXPECT_EQ(fixed3(1.0), "1.000");
    EXPECT_EQ(shortest(0.1), "0.1");
    EXPECT_EQ(shortest(double(0.1f)), "0.10000000149011612");
}

TEST(Sha256, IsOfTheUncompressedText) {
    EXPECT_EQ(sha256_hex("abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    TempDir tmp;
    golden::write(tmp.file("abc"), "abc");
    EXPECT_EQ(file_sha256(tmp.file("abc")), sha256_hex("abc"));
    EXPECT_THROW(file_sha256(tmp.file("absent")), Error);
}

TEST(Lines, ReadPlainTextWithOrWithoutAFinalNewline) {
    TempDir tmp;
    golden::write(tmp.file("a"), "one\ntwo");
    LineReader reader(tmp.file("a"));
    std::string line;
    std::vector<std::string> lines;
    while (reader.next(line)) lines.push_back(line);
    EXPECT_EQ(lines, (std::vector<std::string>{"one", "two"}));
    reader.rewind();
    ASSERT_TRUE(reader.next(line));
    EXPECT_EQ(line, "one");
    EXPECT_THROW(LineReader(tmp.file("absent")), Error);
}
