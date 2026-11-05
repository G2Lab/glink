// Independent implementation of the PLINK/Haploview Gabriel block rules.
// The D′ point fit remains Glink's own likelihood calculation.
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "dprime.hpp"

namespace glink {
namespace {

// Likelihood surface for the CI, with the association oriented positive.
struct Surface {
    std::array<double, 4> base;
    double ambiguous, p, q, dprime;
    bool valid;
};

Surface phased_surface(const Counts& c) {
    Surface out = {{}, 0, 0, 0, NAN, false};
    double n = c[0] + c[1] + c[3] + c[4];
    double a = c[3] + c[4], b = c[1] + c[4];
    if (n <= 0 || std::min({a, n - a, b, n - b}) <= 0) return out;
    out.base = {double(c[0]), double(c[1]), double(c[3]), double(c[4])};
    out.p = a / n; out.q = b / n;
    double d = c[4] / n - out.p * out.q;
    double dm = d >= 0 ? std::min(out.p * (1 - out.q), (1 - out.p) * out.q)
                       : std::min(out.p * out.q, (1 - out.p) * (1 - out.q));
    out.dprime = std::clamp(std::abs(d) / dm, 0., 1.);
    if (d < 0) {
        std::swap(out.base[0], out.base[1]); std::swap(out.base[2], out.base[3]);
        out.q = 1 - out.q;
    }
    out.valid = true;
    return out;
}

Surface diploid_surface(Counts c) {
    Surface out = {{}, 0, 0, 0, NAN, false};
    c = canonical(c);
    const Totals t = totals(c);
    if (t.n <= 0 || std::min({t.a, 2 * t.n - t.a, t.b, 2 * t.n - t.b}) <= 0) return out;
    const DiploidFit fit = fit_diploid(c, t);
    out.base = fit.base; out.ambiguous = fit.ambiguous;
    out.p = fit.p; out.q = fit.q; out.dprime = fit.dprime;
    out.valid = true;
    return out;
}

Surface surface(const Counts& c, int ploidy) {
    return ploidy == 1 ? phased_surface(c) : diploid_surface(c);
}

// PLINK 1.9 CI band 0..6; optionally the reporting bounds.
uint8_t band(const Surface& s, double* bounds = nullptr) {
    if (!s.valid) {
        if (bounds) { bounds[0] = NAN; bounds[1] = NAN; }
        return 1;
    }
    double dm = std::min(s.p * (1 - s.q), (1 - s.p) * s.q), logs[101], mx = -INFINITY;
    for (int k = 0; k <= 100; k++) {
        double d = dm * k * .01;
        std::array<double, 4> f = {(1 - s.p) * (1 - s.q) + d, (1 - s.p) * s.q - d,
                                   s.p * (1 - s.q) - d, s.p * s.q + d};
        // Match the Haploview/PLINK endpoint convention at D′ = 1.
        if (k == 100) for (double& x : f) x = std::max(x, 1e-10);
        logs[k] = likelihood(s.base, s.ambiguous, f);
        mx = std::max(mx, logs[k]);
    }
    if (!std::isfinite(mx)) return 1;
    double tail[102] = {};
    for (int k = 100; k >= 0; k--) tail[k] = tail[k + 1] + std::exp(logs[k] - mx);
    double total = tail[0], cut = total * .95;
    if (bounds) {
        double cdf = 0;
        int lo = 0, hi = 100;
        bool found = false;
        for (int k = 0; k <= 100; k++) {
            cdf += std::exp(logs[k] - mx);
            if (!found && cdf >= total * .05) { lo = k; found = true; }
            if (cdf >= total * .95) { hi = k; break; }
        }
        bounds[0] = std::max(0, lo - 1) * .01;
        bounds[1] = std::min(100, hi + 1) * .01;
    }
    // Seven CI categories used by PLINK 1.9's default block rules. In
    // particular, grid-tail boundary indices differ from reporting quantiles.
    if (tail[97] * 20 <= total) return tail[89] * 20 <= total ? 0 : 1;
    if (cut < tail[72]) return cut < tail[82] ? 6 : 5;
    if (cut >= tail[71]) return cut >= tail[52] ? 2 : 3;
    return 4;
}

struct Candidate { int64_t span; int left, right; };

}  // namespace

Interval gabriel_confidence(Counts counts, int ploidy) {
    const Surface s = surface(counts, ploidy);
    double bounds[2];
    const int category = band(s, bounds);
    return {s.dprime, bounds[0], bounds[1], category};
}

namespace {

// One pair's CI band, memoised per thread on its count table.
using BandCache = std::unordered_map<Counts, uint8_t, CountsHash>;

// A pair of perfect complements: no observation carries neither or both
// (phased: no 0/0 or 1/1; diploid: only 0/2, 1/1 and 2/0). Checked before
// canonicalisation, whose allele flips turn complements into identical rows.
// Both markers must still be polymorphic, else the pair stays uninformative.
bool complementary(const Counts& c, int ploidy) {
    if (ploidy == 1) return c[0] == 0 && c[4] == 0 && c[1] > 0 && c[3] > 0;
    if (c[0] || c[1] || c[3] || c[5] || c[7] || c[8]) return false;
    const int n = c[2] + c[4] + c[6], a = c[4] + 2 * c[6];  // marker i's allele count
    return a > 0 && a < 2 * n;
}

uint8_t pair_band(const PackedStates& states, int i, int j, int ploidy, BandCache& cache,
                  bool neutral_complements) {
    Counts c = states.pair(i, j, ploidy);
    if (neutral_complements && complementary(c, ploidy)) return NEUTRAL;
    if (ploidy == 2) c = canonical(c);
    else c = std::min(c, Counts{c[0], c[3], 0, c[1], c[4], 0, 0, 0, 0});
    auto it = cache.find(c);
    if (it != cache.end()) return it->second;
    uint8_t value = band(surface(c, ploidy));
    if (cache.size() < 1500000) cache.emplace(c, value);
    return value;
}

// First marker within `span` bp before each marker.
std::vector<int> window_starts(const int64_t* pos, int n, int64_t span) {
    std::vector<int> first(n);
    for (int j = 0; j < n; j++) first[j] = std::lower_bound(pos, pos + j, pos[j] - span) - pos;
    return first;
}

// Bands stored for every pair, n × n.
struct DenseBands {
    const uint8_t* bands;
    int n;
    uint8_t operator()(int i, int j) const { return bands[size_t(i) * n + j]; }
};

// PLINK 1.9's block rules, fed one marker j at a time with the bands of every pair
// (i, j) within the span. `bands(a, b)` (a < b) is only asked for b == j or for
// pairs at most three markers apart.
struct Selector {
    const int64_t* pos;
    std::vector<int64_t> strong, recomb;
    std::vector<Candidate> candidates;

    Selector(const int64_t* pos, int n) : pos(pos), strong(n, 0), recomb(n, 0) {}

    template <class Bands>
    void add(int j, int first, const Bands& bands) {
        const double fraction = .95 + 0.00000000002910383045673370361328125;
        int64_t added_s = 0, added_r = 0;
        for (int i = j - 1; i >= first; i--) {
            uint8_t endpoint = bands(i, j);
            // A neutral pair is neither strong nor recombinant and cannot end a
            // block, but passes the 3- and 4-marker rules' inner-pair checks.
            added_s += endpoint >= 5 && endpoint != NEUTRAL; added_r += endpoint == 0;
            strong[i] += added_s; recomb[i] += added_r;
            if (endpoint < 4 || endpoint == NEUTRAL) continue;
            int markers = j - i + 1;
            int64_t length = pos[j] - pos[i];
            bool pass = false;
            if (markers == 2) pass = endpoint == 6 && length <= 20000;
            else if (markers == 3)
                pass = length <= 30000 && bands(i, i + 1) >= 3 && bands(i + 1, j) >= 3;
            else if (markers == 4) {
                pass = true;
                for (int a = i; a <= j && pass; a++) for (int b = a + 1; b <= j; b++)
                    if ((a != i || b != j) && bands(a, b) < 3) { pass = false; break; }
            } else pass = strong[i] + recomb[i] >= 6 && strong[i] > (strong[i] + recomb[i]) * fraction;
            if (pass) candidates.push_back({length, i, j});
        }
    }

    // Longest candidates first, never sharing a marker.
    std::vector<std::pair<int, int>> choose(int n) {
        std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
            if (a.span != b.span) return a.span > b.span;
            if (a.left != b.left) return a.left > b.left;
            return a.right > b.right;
        });
        std::vector<uint8_t> assigned(n, 0);
        std::vector<std::pair<int, int>> selected;
        for (auto c : candidates) {
            if (assigned[c.left] || assigned[c.right]) continue;
            selected.emplace_back(c.left, c.right);
            std::fill(assigned.begin() + c.left, assigned.begin() + c.right + 1, 1);
        }
        std::sort(selected.begin(), selected.end(), [](auto a, auto b) {
            return a.second != b.second ? a.second < b.second : a.first < b.first;
        });
        return selected;
    }
};

}  // namespace

void gabriel_bands(const uint8_t* packed, int n, int width, const int64_t* pos,
                   int64_t span, int ploidy, uint8_t* out, int threads, bool neutral_complements) {
    const PackedStates states(packed, n, width);
    const std::vector<int> first = window_starts(pos, n, span);
    std::fill(out, out + size_t(n) * n, uint8_t(1));
    #pragma omp parallel num_threads(threads)
    {
        BandCache cache;
        #pragma omp for schedule(dynamic, 2)
        for (int i = 0; i < n; i++)
            for (int j = first[i]; j < i; j++)
                out[size_t(i) * n + j] = out[size_t(j) * n + i] =
                    pair_band(states, i, j, ploidy, cache, neutral_complements);
    }
}

std::vector<std::pair<int, int>> gabriel_select(const uint8_t* bands, const int64_t* pos, int n, int64_t span) {
    const std::vector<int> first = window_starts(pos, n, span);
    const DenseBands dense{bands, n};
    Selector selector(pos, n);
    for (int j = 1; j < n; j++) selector.add(j, first[j], dense);
    return selector.choose(n);
}

std::vector<std::pair<int, int>> gabriel_blocks(const uint8_t* packed, int n, int width, const int64_t* pos,
                                                int64_t span, int ploidy, int threads,
                                                bool neutral_complements) {
    // Rows of bands are computed a chunk at a time in parallel, handed to the
    // selector in order and dropped; only bands of markers at most three apart
    // are kept, for the 3- and 4-marker rules.
    constexpr int CHUNK = 256;
    const PackedStates states(packed, n, width);
    const std::vector<int> first = window_starts(pos, n, span);
    std::vector<std::array<uint8_t, 3>> near(n, {1, 1, 1});   // near[b][d - 1] = band(b - d, b)
    std::vector<std::vector<uint8_t>> rows(CHUNK);
    std::vector<BandCache> caches(std::max(threads, 1));
    Selector selector(pos, n);
    for (int start = 1; start < n; start += CHUNK) {
        const int stop = std::min(n, start + CHUNK);
        for (int j = start; j < stop; j++) rows[j - start].assign(j - first[j], uint8_t(1));
        #pragma omp parallel num_threads(threads)
        {
#ifdef _OPENMP
            BandCache& cache = caches[omp_get_thread_num()];
#else
            BandCache& cache = caches[0];
#endif
            #pragma omp for schedule(dynamic, 1)
            for (int j = start; j < stop; j++)
                for (int i = first[j]; i < j; i++)
                    rows[j - start][i - first[j]] = pair_band(states, j, i, ploidy, cache, neutral_complements);
        }
        for (int j = start; j < stop; j++) {
            const std::vector<uint8_t>& row = rows[j - start];
            const int f = first[j];
            for (int d = 1; d <= 3 && j - d >= f; d++) near[j][d - 1] = row[j - d - f];
            selector.add(j, f, [&](int a, int b) -> uint8_t {
                if (b == j) return a >= f ? row[a - f] : 1;
                const int d = b - a;
                return d >= 1 && d <= 3 && a >= first[b] ? near[b][d - 1] : 1;
            });
        }
    }
    return selector.choose(n);
}

}  // namespace glink
