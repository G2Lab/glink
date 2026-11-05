// Pairwise r² and D′ from bit planes.
//
// The integer sufficient statistics come from popcounts of called / allele bit
// planes, and are combined with float64 operations in a fixed order, so every output
// value is reproducible bit for bit. Both matrices are symmetric, so
// only one triangle is computed and mirrored.
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <string_view>
#include <sys/mman.h>
#include <type_traits>
#include <unistd.h>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "ld.hpp"
#include "matrix.hpp"
#include "engine.hpp"

namespace glink {

namespace {

constexpr size_t TILE = 64;

// Per row: called mask M, allele planes P1 (x >= 1) and P2 (x == 2), 64 columns per word.
struct Planes {
    size_t rows, words;
    bool has_two = false, missing = false;
    std::vector<uint64_t> m, p1, p2;

    explicit Planes(const Values& v) : rows(v.rows), words((v.cols + 63) / 64) {
        m.assign(rows * words, 0);
        p1.assign(rows * words, 0);
        p2.assign(rows * words, 0);
        for (size_t r = 0; r < rows; r++) {
            const int8_t* x = v.row(r);
            for (size_t c = 0; c < v.cols; c++) {
                const uint64_t bit = uint64_t(1) << (c & 63);
                const size_t w = r * words + (c >> 6);
                if (x[c] < 0) { missing = true; continue; }
                m[w] |= bit;
                if (x[c] >= 1) p1[w] |= bit;
                if (x[c] == 2) { p2[w] |= bit; has_two = true; }
            }
        }
    }
    const uint64_t* M(size_t r) const { return &m[r * words]; }
    const uint64_t* P1(size_t r) const { return &p1[r * words]; }
    const uint64_t* P2(size_t r) const { return &p2[r * words]; }
};

inline long long pc(uint64_t x) { return __builtin_popcountll(x); }

// Joint statistics of rows a and b over their jointly called columns.
struct Pair { long long n, sa, sb, aa, bb, xy; };

inline Pair pair_stats(const Planes& p, size_t a, size_t b) {
    const uint64_t *ma = p.M(a), *mb = p.M(b), *xa = p.P1(a), *xb = p.P1(b);
    long long n = 0, s1a = 0, s1b = 0, xy = 0;
    for (size_t w = 0; w < p.words; w++) {
        n += pc(ma[w] & mb[w]);
        s1a += pc(xa[w] & mb[w]);
        s1b += pc(ma[w] & xb[w]);
        xy += pc(xa[w] & xb[w]);
    }
    Pair r{n, s1a, s1b, s1a, s1b, xy};
    if (p.has_two) {
        const uint64_t *ya = p.P2(a), *yb = p.P2(b);
        long long s2a = 0, s2b = 0, x12 = 0, x21 = 0, x22 = 0;
        for (size_t w = 0; w < p.words; w++) {
            s2a += pc(ya[w] & mb[w]);
            s2b += pc(ma[w] & yb[w]);
            x12 += pc(xa[w] & yb[w]);
            x21 += pc(ya[w] & xb[w]);
            x22 += pc(ya[w] & yb[w]);
        }
        r.sa += s2a; r.sb += s2b;
        r.aa = s1a + 3 * s2a; r.bb = s1b + 3 * s2b;
        r.xy += x12 + x21 + x22;
    }
    return r;
}

// Upper-triangle tiles in parallel; f(i, j) for every i <= j.
template <class F>
void triangle(size_t n, int threads, F f) {
    const size_t tiles = (n + TILE - 1) / TILE;
    std::vector<std::pair<size_t, size_t>> work;
    for (size_t a = 0; a < tiles; a++) for (size_t b = a; b < tiles; b++) work.emplace_back(a, b);
    #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (size_t k = 0; k < work.size(); k++) {
        const size_t a0 = work[k].first * TILE, b0 = work[k].second * TILE;
        const size_t a1 = std::min(n, a0 + TILE), b1 = std::min(n, b0 + TILE);
        for (size_t i = a0; i < a1; i++)
            for (size_t j = std::max(b0, i); j < b1; j++) f(i, j);
    }
}

}  // namespace

std::vector<float> compute_r2(const NodeMatrix& matrix, int min_pairs, int min_allele_count, int threads) {
    const Values& v = matrix.values;
    const size_t rows = v.rows;
    const double total = double(rows) * double(rows) * 4;
    if (total > 8.0 * 1024 * 1024 * 1024)
        std::fprintf(stderr, "warning: r² is %s x %s = %.1f GiB in memory; narrow --filter to shrink it\n",
                     thousands(rows).c_str(), thousands(rows).c_str(), total / (1024.0 * 1024 * 1024));
    if (min_pairs < 2 || min_allele_count < 0) throw Error("Invalid pair support requirements");
    const Planes p(v);
    const bool missing = p.missing || min_pairs > 2 || min_allele_count > 0;
    const int ploidy = matrix.matrix_type == "dosage" ? 2 : 1;
    std::vector<float> values(rows * rows);
    float* r2 = values.data();
    {
        if (missing) {
            triangle(rows, threads, [&](size_t i, size_t j) {
                const Pair s = pair_stats(p, i, j);
                const double n = double(s.n), sa = double(s.sa), sb = double(s.sb);
                const double denominator = (n * double(s.aa) - sa * sa) * (n * double(s.bb) - sb * sb);
                bool valid = s.n >= min_pairs && denominator > 0;
                if (min_allele_count) {
                    double low = std::min({sa, ploidy * n - sa, sb, ploidy * n - sb});
                    valid = valid && low >= min_allele_count;
                }
                float value = NAN;
                if (valid) {
                    const double num = n * double(s.xy) - sa * sb;
                    value = float(std::clamp(num * num / denominator, 0.0, 1.0));
                }
                r2[i * rows + j] = r2[j * rows + i] = value;
            });
        } else {
            const double n = double(v.cols);
            std::vector<double> sx(rows), scale(rows);
            std::vector<uint8_t> ok(rows);
            for (size_t r = 0; r < rows; r++) {
                double s = 0, ss = 0;
                for (size_t c = 0; c < v.cols; c++) { s += v.row(r)[c]; ss += double(v.row(r)[c]) * v.row(r)[c]; }
                const double variance = n * ss - s * s;
                sx[r] = s;
                ok[r] = variance > 0;
                scale[r] = std::sqrt(ok[r] ? variance : 1.0);
            }
            triangle(rows, threads, [&](size_t i, size_t j) {
                float value = NAN;
                if (ok[i] && ok[j]) {
                    const Pair s = pair_stats(p, i, j);
                    const double covariance = n * double(s.xy) - sx[i] * sx[j];
                    const double r = std::clamp(covariance / (scale[i] * scale[j]), -1.0, 1.0);
                    value = float(r * r);
                }
                r2[i * rows + j] = r2[j * rows + i] = value;
            });
        }
    }
    return values;
}

std::vector<double> compute_dprime(const NodeMatrix& matrix, int threads) {
    const Values& v = matrix.values;
    const size_t rows = v.rows;
    const Planes p(v);
    std::vector<double> values(rows * rows);
    double* d = values.data();
    if (matrix.matrix_type != "dosage") {
        if (p.has_two)
            throw Error("D′ needs phase or 0/1 calls; this matrix holds values above 1");
        triangle(rows, threads, [&](size_t i, size_t j) {
            const Pair s = pair_stats(p, i, j);
            const double total = double(s.n);
            const double pa = double(s.sa) / total, pb = double(s.sb) / total;
            const double dd = double(s.xy) / total - pa * pb;
            auto nanmin = [](double x, double y) { return std::isnan(x) || std::isnan(y) ? NAN : std::min(x, y); };
            const double limit = dd >= 0 ? nanmin(pa * (1 - pb), (1 - pa) * pb) : nanmin(pa * pb, (1 - pa) * (1 - pb));
            double value = limit > 0 ? dd / limit : NAN;
            if (!std::isnan(value)) value = std::clamp(value, -1.0, 1.0);
            d[i * rows + j] = d[j * rows + i] = value;
        });
    } else {
        // One genotype-state plane per marker and state, and each marker's distinct pattern. A marker and
        // its allele flip have the same |D′|, so they share a pattern, as do exact duplicates. This is only
        // to save work: every marker gets its value back below.
        const size_t width = (v.cols + 7) / 8, stride = 3 * width;
        std::vector<uint8_t> planes(rows * stride, 0), flipped(rows * stride, 0);
        for (size_t r = 0; r < rows; r++) {
            const int8_t* x = v.row(r);
            for (size_t c = 0; c < v.cols; c++) {
                if (x[c] < 0) continue;
                const uint8_t bit = uint8_t(0x80 >> (c % 8));
                planes[r * stride + size_t(x[c]) * width + c / 8] |= bit;
                flipped[r * stride + size_t(2 - x[c]) * width + c / 8] |= bit;
            }
        }
        std::unordered_map<std::string_view, int> seen;
        std::vector<int> pattern(rows);
        std::vector<uint8_t> distinct;
        for (size_t r = 0; r < rows; r++) {
            const uint8_t* a = &planes[r * stride];
            const uint8_t* b = &flipped[r * stride];
            const uint8_t* chosen = std::memcmp(b, a, stride) < 0 ? b : a;
            const std::string_view key(reinterpret_cast<const char*>(chosen), stride);
            const auto [it, fresh] = seen.emplace(key, int(seen.size()));
            if (fresh) distinct.insert(distinct.end(), chosen, chosen + stride);
            pattern[r] = it->second;
        }
        const size_t k = seen.size();
        std::vector<float> fitted(k * k);
        dosage_dprime(distinct.data(), int(k), int(width), fitted.data(), threads);
        #pragma omp parallel for schedule(static) num_threads(threads)
        for (size_t i = 0; i < rows; i++)
            for (size_t j = 0; j < rows; j++) d[i * rows + j] = fitted[size_t(pattern[i]) * k + size_t(pattern[j])];
    }
    return values;
}

namespace {

// A text file written whole to PATH.partial, then renamed, so that PATH is always complete.
class TextOut {
public:
    explicit TextOut(const std::string& path) : path_(path), partial_(path + ".partial"), file_(std::fopen(partial_.c_str(), "w")) {
        if (!file_) throw Error("Cannot write " + path);
    }
    ~TextOut() {
        if (file_) { std::fclose(file_); std::remove(partial_.c_str()); }
    }
    TextOut(const TextOut&) = delete;
    TextOut& operator=(const TextOut&) = delete;

    std::string buffer;
    void flush() {
        if (std::fwrite(buffer.data(), 1, buffer.size(), file_) != buffer.size()) throw Error("Cannot write " + path_);
        buffer.clear();
    }
    void maybe_flush() { if (buffer.size() > (size_t(1) << 20)) flush(); }
    void raw(const void* data, size_t bytes) {
        if (std::fwrite(data, 1, bytes, file_) != bytes) throw Error("Cannot write " + path_);
    }
    void close() {
        flush();
        FILE* f = file_;
        file_ = nullptr;
        if (std::fclose(f) != 0 || std::rename(partial_.c_str(), path_.c_str()) != 0) {
            std::remove(partial_.c_str());
            throw Error("Cannot write " + path_);
        }
    }

private:
    std::string path_, partial_;
    FILE* file_;
};

// PLINK's numbers: as %g, six significant digits, and always "nan", whatever a NaN's sign. to_chars
// writes the same text as printf, about three times as fast, and these files hold billions of numbers.
void put_number(std::string& out, double value, int width = 0) {
    char text[32];
    const char* end;
    if (std::isnan(value)) {
        std::memcpy(text, "nan", 3);
        end = text + 3;
    } else {
        end = std::to_chars(text, text + sizeof text, value, std::chars_format::general, 6).ptr;
    }
    const int n = int(end - text);
    if (n < width) out.append(size_t(width - n), ' ');
    out.append(text, size_t(n));
}

// A chromosome as PLINK's reports name it: any case of chr dropped, human chromosomes as numbers
// (X 23, Y 24, XY 25, MT 26, and any larger number 26). Other contigs keep their names.
std::string plink_chromosome(const std::string& name) {
    std::string rest = name;
    if (rest.size() > 3 && std::tolower((unsigned char)rest[0]) == 'c' && std::tolower((unsigned char)rest[1]) == 'h' &&
        std::tolower((unsigned char)rest[2]) == 'r')
        rest = rest.substr(3);
    std::string upper;
    for (unsigned char c : rest) upper += char(std::toupper(c));
    if (upper == "X") return "23";
    if (upper == "Y") return "24";
    if (upper == "XY") return "25";
    if (upper == "M" || upper == "MT") return "26";
    bool digits = !rest.empty() && rest.size() < 10;
    for (unsigned char c : rest) digits = digits && std::isdigit(c);
    return digits ? std::to_string(std::min(std::stoi(rest), 26)) : name;
}

size_t square_size(size_t cells) {
    const size_t n = size_t(std::llround(std::sqrt(double(cells))));
    if (n * n != cells) throw Error("LD matrix is not square");
    return n;
}

void write_snplist(const std::string& stem, const std::vector<std::string>& ids) {
    TextOut snplist(stem + ".snplist");
    for (auto& id : ids) { snplist.buffer += id; snplist.buffer += '\n'; }
    snplist.close();
}

// The row IDs, then the matrix as float32. Only D′ on phased calls has to narrow, a row at a time.
template <class T>
void write_binary(const std::string& path, const std::vector<T>& values, const std::vector<std::string>& ids) {
    const size_t n = square_size(values.size());
    if (ids.size() != n) throw Error("One ID per row is needed for " + path);
    if (!ends_with(path, ".ld.bin")) throw Error(path + " must end in .ld.bin");
    write_snplist(path.substr(0, path.size() - 7), ids);
    TextOut out(path);
    if constexpr (std::is_same_v<T, float>) {
        if (!values.empty()) out.raw(values.data(), values.size() * sizeof(float));
    } else {
        std::vector<float> row(n);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) row[j] = float(values[i * n + j]);
            out.raw(row.data(), n * sizeof(float));
        }
    }
    out.close();
}

}  // namespace

void write_ld_binary(const std::string& path, const std::vector<float>& values, const std::vector<std::string>& ids) {
    write_binary(path, values, ids);
}

void write_ld_binary(const std::string& path, const std::vector<double>& values, const std::vector<std::string>& ids) {
    write_binary(path, values, ids);
}

void write_r2_ld(const std::string& path, const std::vector<float>& r2, const std::vector<std::string>& ids) {
    const size_t n = square_size(r2.size());
    if (ids.size() != n) throw Error("One ID per row is needed for " + path);
    if (!ends_with(path, ".ld")) throw Error(path + " must end in .ld");
    write_snplist(path.substr(0, path.size() - 3), ids);
    TextOut out(path);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (j) out.buffer += '\t';
            put_number(out.buffer, r2[i * n + j]);
        }
        out.buffer += '\n';
        out.maybe_flush();
    }
    out.close();
}

void write_dprime_ld(const std::string& path, const NodeMatrix& m, const std::vector<float>& r2, const std::vector<double>& dprime) {
    const size_t n = square_size(r2.size());
    if (dprime.size() != r2.size()) throw Error("r² and D′ differ in size");
    if (!ends_with(path, ".ld")) throw Error(path + " must end in .ld");
    const std::vector<std::string> ids = m.row_labels();
    if (ids.size() != n) throw Error("One ID per row is needed for " + path);
    std::vector<std::string> chrom(n, "0");
    std::vector<int64_t> position(n, 0);
    if (!m.coordinates.empty()) {
        const auto codes = m.coordinates.at("chrom");
        position = m.coordinates.at("position");
        for (size_t i = 0; i < n; i++)
            if (codes[i] >= 0) chrom[i] = plink_chromosome(m.coordinate_chromosomes[size_t(codes[i])]);
    }
    // PLINK's ID column: 4 wide for IDs up to 4 characters, otherwise the longest ID plus 2.
    size_t longest = 0;
    for (auto& id : ids) longest = std::max(longest, id.size());
    const int w = longest <= 4 ? 4 : int(longest) + 2;
    auto field = [](std::string& out, const char* format, int width, const char* value) {
        char text[512];
        const int k = std::snprintf(text, sizeof text, format, width, value);
        if (k < 0 || size_t(k) >= sizeof text) throw Error("An ID is too long for an LD table");
        out.append(text, size_t(k));
    };
    TextOut out(path);
    field(out.buffer, "%*s", 6, "CHR_A");
    field(out.buffer, " %*s", 12, "BP_A");
    field(out.buffer, " %*s", w, "SNP_A");
    field(out.buffer, " %*s", 6, "CHR_B");
    field(out.buffer, " %*s", 12, "BP_B");
    field(out.buffer, " %*s", w, "SNP_B");
    field(out.buffer, " %*s", 12, "R2");
    field(out.buffer, " %*s", 12, "DP");
    out.buffer += " \n";
    for (size_t i = 0; i < n; i++) {
        const std::string bp_a = std::to_string(position[i]);
        for (size_t j = i + 1; j < n; j++) {
            field(out.buffer, "%*s", 6, chrom[i].c_str());
            field(out.buffer, " %*s", 12, bp_a.c_str());
            field(out.buffer, " %*s", w, ids[i].c_str());
            field(out.buffer, " %*s", 6, chrom[j].c_str());
            field(out.buffer, " %*s", 12, std::to_string(position[j]).c_str());
            field(out.buffer, " %*s", w, ids[j].c_str());
            out.buffer += ' ';
            put_number(out.buffer, r2[i * n + j], 12);
            out.buffer += ' ';
            put_number(out.buffer, dprime[i * n + j], 12);
            out.buffer += " \n";
            out.maybe_flush();
        }
    }
    out.close();
}

}  // namespace glink
