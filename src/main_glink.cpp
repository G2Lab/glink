// glink: GFA -> node matrix -> r² matrix (and D′, blocks).
#include <cmath>
#include <cstdio>
#include <functional>
#include <set>

#include <tuple>

#include "args.hpp"
#include "core.hpp"

using namespace glink;

namespace {

struct Settings {
    std::string gfa, samples, presence_bundle, node_coordinates, coordinate_reference = "CHM13",
        coordinate_reference_source = "auto", conditioning_tree, snarls, snarl_tree, vg_graph, snarls_format = "auto",
        backbone = "CHM13", vg, matrix_type = "dosage", matrix_out, r2_out,
        dprime_out, blocks_out, out_prefix;
    bool no_node_coordinates = false, no_conditioning_tree = false, rebuild_conditioning_tree = false, assume_complete = false,
         condition_on_parents = false, split_inversions = false, drop_complements = false, collapse = true,
         neutral_complements = false, collapse_backbone_parents = false, rebuild = false;
    int threads = 1, min_pairs = 2, min_allele_count = 0, broken_parent_min_haplotypes = 6;
    bool blocks = false;
    double blocks_max_kb = 200;
    std::vector<Filter> filters;
};

Settings parse(int argc, char** argv) {
    Options p("glink", "Build a node matrix from a GFA and write its r² matrix (C++ GlinkV3).");
    p.option("gfa", "Input GFA; haplotype paths are read from its W/P lines");
    p.option("out-prefix", "Run glink index first and keep its files under this prefix, reusing any already there");
    p.flag("rebuild", "With --out-prefix: replace the index files instead of reusing them");
    p.option("samples", "Analyze these sample IDs, one per line, before frequency filtering");
    p.option("presence-bundle", "Load the graph from this store (from glink index) instead of reading the GFA");
    p.option("node-coordinates", "Build/reuse node coordinates here (default: .glink-cache/node_coordinates/)");
    p.option("coordinate-reference", "Reference W sample, sample#haplotype, or exact P name (default: CHM13)");
    p.option("coordinate-reference-source", "auto, walk or rgfa (default: auto)");
    p.flag("no-node-coordinates", "Skip coordinate preprocessing");
    p.option("conditioning-tree", "Build/reuse the conditioning tree here (default: .glink-cache/tree/)");
    p.flag("no-conditioning-tree", "Skip the conditioning tree");
    p.option("snarls", "Use saved vg snarls output instead of running vg snarls");
    p.option("snarl-tree", "Use this existing conditioning tree store");
    p.option("vg-graph", "Run vg snarls on this GBZ; GFA must be its --no-translation export");
    p.option("snarls-format", "auto, json or binary (default: auto)");
    p.option("backbone", "Backbone W sample or P name used to root newly computed snarls (default: CHM13)");
    p.option("vg", "vg executable used on a parent-cache miss or to decode binary snarls");
    p.option("threads", "Threads for vg snarls, LD and blocks (default: 1)");
    p.flag("rebuild-conditioning-tree", "Explicitly replace the conditioning tree");
    p.flag("assume-complete", "Explicitly treat every untraversed node as absent");
    p.flag("condition-on-parents", "Call each node only on haplotypes that visit its parent");
    p.flag("collapse-backbone-parents", "With --condition-on-parents, skip backbone and broken parents: their "
                                        "children condition on the nearest remaining ancestor");
    p.option("broken-parent-min-haplotypes", "A parent is broken when at least N haplotypes visit one of its "
                                             "children while skipping it (default: 6)");
    p.flag("split-inversions", "Give each inverted node one row per traversal orientation");
    p.flag("drop-complements", "Drop each row that is the perfect complement of an earlier row");
    p.option("min-pairs", "Minimum jointly called observations per LD pair (default: 2)");
    p.option("min-allele-count", "Minimum copies of each allele in the jointly called observations (default: 0)");
    p.option("matrix-type", "dosage, phased (or haplotype) or presence (default: dosage)");
    p.option("filter", "Repeatable: all, majority, polymorphic, common, LOW:HIGH, minor:F, or @file");
    p.flag("collapse", "Correlate distinct genotype patterns rather than every node (default)");
    p.flag("no-collapse", "Correlate every node");
    p.option("matrix-out", "Write the node matrix as a store (a directory) here; omit to not keep it");
    p.option("r2-out", "Write the r² matrix to this .ld file, as plink --r2 square does, with row IDs in the .snplist beside it; "
                       "a name ending in .ld.bin writes raw float32 instead");
    p.option("dprime-out", "Write every pair's r² and D′ to this .ld file, as plink --r2 inter-chr dprime does "
                           "(phased/presence: signed D′; dosage: genotype-likelihood |D′|); "
                           "a name ending in .ld.bin writes the D′ matrix alone, as raw float32");
    p.flag("blocks", "Also call Gabriel et al. (2002) D′ confidence-interval blocks with PLINK 1.9's rules");
    p.option("blocks-max-kb", "Longest block called (default: 200)");
    p.flag("neutral-complements", "--blocks: a pair of perfect complements (no haplotype carrying neither or "
                                  "both) counts as neither strong LD nor recombination and cannot end a block");
    p.option("blocks-out", "Write the block table here (default: alongside --r2-out)");
    p.parse(argc, argv);

    Settings s;
    if (!p.has("gfa")) p.error("the following arguments are required: --gfa");
    s.gfa = p.get("gfa");
    s.out_prefix = p.get("out-prefix");
    s.rebuild = p.has("rebuild");
    if (s.rebuild && s.out_prefix.empty()) p.error("--rebuild needs --out-prefix");
    if (p.has("out-prefix"))
        for (const char* name : {"node-coordinates", "no-node-coordinates", "conditioning-tree", "no-conditioning-tree", "snarls",
                                 "snarl-tree", "rebuild-conditioning-tree", "presence-bundle"})
            if (p.has(name)) p.error(std::string("--out-prefix cannot be combined with --") + name);
    s.samples = p.get("samples");
    s.presence_bundle = p.get("presence-bundle");
    s.node_coordinates = p.get("node-coordinates");
    s.coordinate_reference = p.get("coordinate-reference", "CHM13");
    s.coordinate_reference_source = p.get("coordinate-reference-source", "auto");
    if (s.coordinate_reference_source != "auto" && s.coordinate_reference_source != "walk" &&
        s.coordinate_reference_source != "rgfa")
        p.error("argument --coordinate-reference-source: invalid choice: '" + s.coordinate_reference_source + "'");
    s.no_node_coordinates = p.has("no-node-coordinates");
    s.conditioning_tree = p.get("conditioning-tree");
    s.no_conditioning_tree = p.has("no-conditioning-tree");
    int sources = p.has("snarls") + p.has("snarl-tree") + p.has("vg-graph");
    if (p.has("snarl-tree") && p.has("conditioning-tree")) p.error("--snarl-tree is the tree, so --conditioning-tree has nothing to build");
    if (sources > 1) p.error("arguments --snarls, --snarl-tree and --vg-graph are mutually exclusive");
    s.snarls = p.get("snarls");
    s.snarl_tree = p.get("snarl-tree");
    s.vg_graph = p.get("vg-graph");
    s.snarls_format = p.get("snarls-format", "auto");
    if (s.snarls_format != "auto" && s.snarls_format != "json" && s.snarls_format != "binary")
        p.error("argument --snarls-format: invalid choice: '" + s.snarls_format + "'");
    s.backbone = p.get("backbone", "CHM13");
    s.vg = p.get("vg", "vg");
    s.threads = p.integer("threads", 1);
    s.rebuild_conditioning_tree = p.has("rebuild-conditioning-tree");
    s.assume_complete = p.has("assume-complete");
    s.condition_on_parents = p.has("condition-on-parents");
    s.collapse_backbone_parents = p.has("collapse-backbone-parents");
    s.broken_parent_min_haplotypes = p.integer("broken-parent-min-haplotypes", 6);
    s.split_inversions = p.has("split-inversions");
    s.drop_complements = p.has("drop-complements");
    s.min_pairs = p.integer("min-pairs", 2);
    s.min_allele_count = p.integer("min-allele-count", 0);
    s.matrix_type = p.get("matrix-type", "dosage");
    if (s.matrix_type == "haplotype") s.matrix_type = "phased";
    if (s.matrix_type != "dosage" && s.matrix_type != "phased" && s.matrix_type != "presence")
        p.error("argument --matrix-type: invalid choice: '" + s.matrix_type + "'");
    for (auto& spec : p.all("filter")) {
        std::string error;
        auto f = parse_filter(spec, error);
        if (!f) p.error("argument --filter: " + error);
        s.filters.push_back(*f);
    }
    s.collapse = p.last_index({"collapse", "no-collapse"}, argv, argc) != 1;
    s.matrix_out = p.get("matrix-out");
    s.r2_out = p.get("r2-out");
    s.dprime_out = p.get("dprime-out");
    s.blocks = p.has("blocks");
    s.blocks_max_kb = p.number("blocks-max-kb", 200);
    s.blocks_out = p.get("blocks-out");
    s.neutral_complements = p.has("neutral-complements");

    if (s.no_node_coordinates && !s.node_coordinates.empty())
        p.error("--node-coordinates cannot be combined with --no-node-coordinates");
    if (s.no_conditioning_tree && (!s.conditioning_tree.empty() || sources))
        p.error("--no-conditioning-tree cannot be combined with conditioning tree inputs");
    if (s.no_conditioning_tree && s.rebuild_conditioning_tree)
        p.error("--no-conditioning-tree cannot be combined with --rebuild-conditioning-tree");
    if (s.condition_on_parents && (s.no_conditioning_tree || s.assume_complete))
        p.error("--condition-on-parents needs a conditioning tree and excludes --assume-complete");
    if (s.drop_complements && s.matrix_type == "presence")
        p.error("Presence calls have no complement; drop complements on phased or dosage");
    if (s.split_inversions && !s.presence_bundle.empty())
        p.error("--split-inversions reads orientation from the GFA; drop --presence-bundle");
    if (s.collapse_backbone_parents && !s.condition_on_parents)
        p.error("--collapse-backbone-parents needs --condition-on-parents");
    if (s.broken_parent_min_haplotypes < 1) p.error("--broken-parent-min-haplotypes must be positive");
    if (s.threads < 1) p.error("--threads must be positive");
    if (s.r2_out.empty() && s.matrix_out.empty()) p.error("Give --r2-out, --matrix-out or both");
    for (auto& [name, path] : {std::pair<std::string, std::string>{"--r2-out", s.r2_out}, {"--dprime-out", s.dprime_out}})
        if (!path.empty() && !ends_with(path, ".ld") && !ends_with(path, ".ld.bin"))
            p.error(name + " must end in .ld, as PLINK's reports do, or in .ld.bin for raw float32");
    if (!s.r2_out.empty() && s.r2_out == s.dprime_out) p.error("--r2-out and --dprime-out must differ");
    if (s.r2_out.empty() && s.blocks && s.blocks_out.empty())
        p.error("--blocks without --r2-out needs --blocks-out");
    std::vector<std::string> tree_paths = {s.conditioning_tree, s.snarls, s.snarl_tree, s.vg_graph};
    std::vector<std::string> other_paths = {s.gfa, s.samples, s.node_coordinates, s.matrix_out, s.r2_out,
                                            s.blocks_out};
    for (size_t i = 0; i < tree_paths.size(); i++) {
        if (tree_paths[i].empty()) continue;
        std::vector<std::string> others = other_paths;
        others.insert(others.end(), tree_paths.begin() + i + 1, tree_paths.end());
        for (auto& o : others)
            if (!o.empty() && resolve(o) == resolve(tree_paths[i]))
                p.error("Conditioning tree inputs and outputs must differ from other input/output files");
    }
    if (!s.node_coordinates.empty())
        for (auto& o : {s.gfa, s.matrix_out, s.r2_out, s.blocks_out})
            if (!o.empty() && resolve(o) == resolve(s.node_coordinates))
                p.error("--node-coordinates must differ from other input/output files");
    if (s.min_pairs < 2 || s.min_allele_count < 0)
        p.error("--min-pairs must be >=2 and --min-allele-count must be >=0");
    if (s.filters.empty()) s.filters.push_back(Filter{"all"});
    if (!s.out_prefix.empty()) {
        const std::string root = resolve(index_paths(s.gfa, s.out_prefix).root) + "/";
        for (auto& o : {s.matrix_out, s.r2_out, s.dprime_out, s.blocks_out})
            if (!o.empty() && starts_with(resolve(o) + "/", root)) p.error("Output " + o + " would overwrite the index");
    }
    return s;
}

// pathlib.Path.with_suffix
std::string with_suffix(const std::string& path, const std::string& suffix) {
    std::string name = basename(path);
    size_t dot = name.rfind('.');
    std::string base = dot == std::string::npos || dot == 0 ? path : path.substr(0, path.size() - (name.size() - dot));
    return base + suffix;
}

std::optional<std::vector<int64_t>> block_positions(const Settings& s, const NodeMatrix& m) {
    if (!m.coordinates.empty()) {
        const auto& position = m.coordinates.at("position");
        const auto& chrom = m.coordinates.at("chrom");
        bool placed = true;
        for (auto p : position) placed &= p > 0;
        std::set<int64_t> chroms(chrom.begin(), chrom.end());
        if (placed && chroms.size() == 1) return position;
    }
    if (s.blocks)
        throw Error("--blocks needs an unambiguous coordinate for every row on one chromosome; "
                    "use the coordinate cache to select placed node IDs with --filter @nodes.txt");
    return std::nullopt;
}

void run_blocks(const Settings& s, const NodeMatrix& m, const std::vector<float>* r2, const std::vector<int64_t>& positions) {
    if (m.matrix_type != "phased" && m.matrix_type != "dosage") throw Error("--blocks needs a phased or dosage matrix");
    const int ploidy = m.matrix_type == "phased" ? 1 : 2;
    BlockResult result = build_blocks(r2, positions, m.values, ploidy, int64_t(s.blocks_max_kb * 1000), s.threads, s.neutral_complements);
    std::string out = s.blocks_out.empty() ? with_suffix(s.r2_out, ".blocks.tsv") : s.blocks_out;
    make_dirs(dirname(out));
    write_blocks(result, m.row_labels(), positions, out);
    long long singletons = 0, largest = 0;
    std::vector<double> within;
    for (auto& st : result.statistics) {
        singletons += st.size == 1;
        largest = std::max<long long>(largest, st.size);
        if (st.size > 1 && st.scored) within.push_back(st.mean_r2);
    }
    std::string mean = "None";
    if (!within.empty()) {
        // np.mean over Python floats: float64 NumPy summation.
        mean = fixed3(numpy_sum(within.data(), within.size()) / double(within.size()));
    }
    std::printf("blocks: %s from %s features (%s singletons, largest %s, mean within-block r² %s) -> %s\n",
                thousands(result.count()).c_str(), thousands(result.labels.size()).c_str(), thousands(singletons).c_str(),
                thousands(largest).c_str(), mean.c_str(), out.c_str());
}

int run(const Settings& s) {
    stage("start");
    if (!exists(s.gfa)) {
        std::fprintf(stderr, "--gfa %s does not exist\n", s.gfa.c_str());
        return 2;
    }
    // read the GFA -> vg snarls -> coordinates -> conditioning tree -> matrix -> r², D′, blocks
    Graph graph;
    std::optional<Store> parents;
    std::optional<Store> coordinates;
    if (!s.out_prefix.empty()) {
        IndexOptions io;
        io.gfa = s.gfa; io.prefix = s.out_prefix; io.reference = s.coordinate_reference;
        io.reference_source = s.coordinate_reference_source; io.backbone = s.backbone; io.vg = s.vg;
        io.vg_graph = s.vg_graph; io.threads = s.threads; io.rebuild = s.rebuild; io.oriented = s.split_inversions;
        Index index = prepare_index(io);
        stage("index");
        parents = std::move(index.tree);
        coordinates = std::move(index.coordinates);
        graph = std::move(index.graph);
    } else {
        TreeOptions o;
        o.output = s.conditioning_tree; o.backbone = s.backbone; o.vg = s.vg; o.threads = s.threads; o.vg_graph = s.vg_graph;
        o.snarls = s.snarls; o.snarl_tree = s.snarl_tree; o.snarls_format = s.snarls_format; o.rebuild = s.rebuild_conditioning_tree;
        const bool call_vg = !s.no_conditioning_tree && needs_snarls(s.gfa, o);

        GraphOptions graph_options;
        graph_options.oriented = s.split_inversions;
        graph_options.reference = s.coordinate_reference;
        graph_options.backbone = s.backbone;
        bool graph_reused = true;
        if (s.presence_bundle.empty()) std::tie(graph, graph_reused) = prepare_graph(s.gfa, graph_options);
        else {
            check_graph(s.presence_bundle, s.gfa, graph_options);
            graph = load_graph(s.presence_bundle);
        }
        stage("read graph");
        std::printf("graph (%s) -> %s\n", graph_reused ? "reused" : "built", graph.bits_path.c_str());
        std::fflush(stdout);

        std::optional<SnarlRun> snarls;
        if (call_vg) {
            snarls = run_snarls(graph, s.gfa, o);
            stage("vg snarls");
        }
        if (!s.no_node_coordinates) {
            auto [coords, reused] = prepare_coordinates(graph, s.gfa, s.node_coordinates, s.coordinate_reference,
                                                        s.coordinate_reference_source, false);
            stage("node coordinates");
            coordinates = std::move(coords);
            std::printf("node coordinates (%s): %s/%s placed on %s -> %s\n", reused ? "reused" : "built",
                        thousands(coordinates->metadata().at("summary").at("placed").get<long long>()).c_str(),
                        thousands(coordinates->shape("node_id").at(0)).c_str(), s.coordinate_reference.c_str(), coordinates->path().c_str());
            std::fflush(stdout);
        }
        if (!s.no_conditioning_tree) {
            auto [tree, reused] = prepare_tree(graph, s.gfa, o, snarls ? &*snarls : nullptr);
            stage("conditioning tree");
            parents = std::move(tree);
            std::printf("conditioning tree (%s): %s/%s nodes have a parent -> %s\n", reused ? "reused" : "built",
                        thousands(parents->metadata().at("with_parent").get<long long>()).c_str(),
                        thousands(parents->shape("node_id").at(0)).c_str(), parents->path().c_str());
            std::fflush(stdout);
        }
    }
    graph.release_step_data();
    BuildOptions b;
    b.matrix_type = s.matrix_type; b.filters = s.filters; b.collapse = s.collapse; b.samples = s.samples;
    b.assume_complete = s.assume_complete; b.condition_on_parents = s.condition_on_parents;
    b.drop_complements = s.drop_complements;
    b.collapse_backbone = s.collapse_backbone_parents; b.broken_min_haplotypes = s.broken_parent_min_haplotypes;
    b.reference = s.coordinate_reference;
    b.coordinates = coordinates ? &*coordinates : nullptr;
    b.parents = parents ? &*parents : nullptr;
    NodeMatrix matrix = build_matrix(std::move(graph), b);
    stage("build matrix");
    if (matrix.values.rows == 0 && (!s.r2_out.empty() || s.blocks || !s.dprime_out.empty()))
        throw Error("The filters kept no nodes; there is no LD to compute");
    coordinates.reset();
    parents.reset();
    auto positions = block_positions(s, matrix);
    std::printf("callability: %s; unknown calls are -1\n", matrix.callability.at("absence_model").get<std::string>().c_str());
    if (const Json* c = json_find(matrix.callability, "complements"))
        std::printf("complements: dropped %s rows\n", thousands(c->at("dropped").get<long long>()).c_str());
    if (const Json* o = json_find(matrix.callability, "orientation"))
        std::printf("orientation: %s inverted nodes split into forward/reverse rows\n",
                    thousands(o->at("inverted_nodes").get<long long>()).c_str());
    const char* rows = matrix.collapsed ? "patterns" : "nodes";
    if (matrix.collapsed)
        std::printf("%s nodes -> %s distinct patterns\n", thousands(matrix.pattern_of.size()).c_str(),
                    thousands(matrix.values.rows).c_str());
    if (!s.matrix_out.empty()) write_matrix(matrix, s.matrix_out);
    stage("write matrix");
    std::printf("%s %s x %s %s columns (density %s)", thousands(matrix.values.rows).c_str(), rows,
                thousands(matrix.values.cols).c_str(), matrix.matrix_type.c_str(), fixed3(matrix.density()).c_str());
    std::printf(s.matrix_out.empty() ? "\n" : " -> %s\n", s.matrix_out.c_str());
    std::fflush(stdout);
    // The D′ table reports r² too, so it needs r² even when --r2-out was not asked for. A binary D′ matrix does not.
    std::vector<float> r2;
    const bool binary_r2 = ends_with(s.r2_out, ".ld.bin"), binary_dprime = ends_with(s.dprime_out, ".ld.bin");
    if (!s.r2_out.empty() || (!s.dprime_out.empty() && !binary_dprime)) {
        r2 = compute_r2(matrix, s.min_pairs, s.min_allele_count, s.threads);
        stage("r2");
    }
    if (!s.r2_out.empty()) {
        make_dirs(dirname(s.r2_out));
        if (binary_r2) write_ld_binary(s.r2_out, r2, matrix.row_labels());
        else write_r2_ld(s.r2_out, r2, matrix.row_labels());
        stage("write r2");
        std::printf("r² %s x %s -> %s\n", thousands(matrix.values.rows).c_str(), thousands(matrix.values.rows).c_str(),
                    s.r2_out.c_str());
        std::fflush(stdout);
    }
    if (!s.dprime_out.empty()) {
        make_dirs(dirname(s.dprime_out));
        const std::vector<double> dprime = compute_dprime(matrix, s.threads);
        stage("dprime");
        if (binary_dprime) write_ld_binary(s.dprime_out, dprime, matrix.row_labels());
        else write_dprime_ld(s.dprime_out, matrix, r2, dprime);
        stage("write dprime");
        std::printf("D′ %s x %s -> %s\n", thousands(matrix.values.rows).c_str(), thousands(matrix.values.rows).c_str(),
                    s.dprime_out.c_str());
        std::fflush(stdout);
    }
    if (s.blocks) run_blocks(s, matrix, s.r2_out.empty() ? nullptr : &r2, *positions);
    stage("blocks");
    return 0;
}

int index_main(int argc, char** argv) {
    Options p("glink index", "Build a graph's store, node coordinates and conditioning tree "
                             "under one prefix. Files already there are reused.");
    p.option("gfa", "Input GFA; haplotype paths are read from its W/P lines");
    p.option("out-prefix", "Write the index to the directory PREFIX.glink; a PREFIX ending in / names it after the GFA");
    p.option("coordinate-reference", "Reference W sample, sample#haplotype, or exact P name (default: CHM13)");
    p.option("coordinate-reference-source", "auto, walk or rgfa (default: auto)");
    p.option("backbone", "Backbone W sample or P name that roots the tree (default: CHM13)");
    p.option("vg", "vg executable (default: vg on PATH); not needed once the tree exists");
    p.option("vg-graph", "Run vg snarls on this GBZ; GFA must be its --no-translation export");
    p.option("threads", "Threads for vg snarls (default: 1)");
    p.flag("rebuild", "Replace the index files instead of reusing them");
    p.parse(argc, argv);
    if (!p.has("gfa") || !p.has("out-prefix")) p.error("the following arguments are required: --gfa, --out-prefix");
    IndexOptions o;
    o.gfa = p.get("gfa");
    o.prefix = p.get("out-prefix");
    o.reference = p.get("coordinate-reference", "CHM13");
    o.reference_source = p.get("coordinate-reference-source", "auto");
    if (o.reference_source != "auto" && o.reference_source != "walk" && o.reference_source != "rgfa")
        p.error("argument --coordinate-reference-source: invalid choice: '" + o.reference_source + "'");
    o.backbone = p.get("backbone", "CHM13");
    o.vg = p.get("vg", "vg");
    o.vg_graph = p.get("vg-graph");
    o.threads = p.integer("threads", 1);
    if (o.threads < 1) p.error("--threads must be positive");
    o.rebuild = p.has("rebuild");
    if (!exists(o.gfa)) {
        std::fprintf(stderr, "--gfa %s does not exist\n", o.gfa.c_str());
        return 2;
    }
    prepare_index(o);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    const bool indexing = argc > 1 && std::string(argv[1]) == "index";
    std::optional<Settings> s;
    if (!indexing) s = parse(argc, argv);
    try {
        return indexing ? index_main(argc - 1, argv + 1) : run(*s);
    } catch (const std::exception& e) {
        std::fflush(stdout);
        std::fprintf(stderr, "glink: error: %s\n", e.what());
        return 1;
    }
}
