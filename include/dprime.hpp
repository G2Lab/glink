// Two-locus likelihood pieces shared by the dosage and Gabriel D′ engines.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

#include "engine.hpp"

namespace glink {

struct CountsHash {
    size_t operator()(const Counts& c) const {
        uint64_t h = 1469598103934665603ULL;
        for (auto x : c) { h ^= x; h *= 1099511628211ULL; }
        return h;
    }
};

// Smallest equivalent table under locus exchange and allele flips at either
// locus, so symmetric tables share cache entries and fit orientation.
inline Counts canonical(const Counts& c) {
    Counts best = c;
    for (int tr = 0; tr < 2; tr++) for (int fx = 0; fx < 2; fx++) for (int fy = 0; fy < 2; fy++) {
        Counts candidate;
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
            int a = tr ? j : i, b = tr ? i : j;
            if (fx) a = 2 - a;
            if (fy) b = 2 - b;
            candidate[3 * i + j] = c[3 * a + b];
        }
        best = std::min(best, candidate);
    }
    return best;
}

// Log likelihood of haplotype frequencies f given known haplotype counts and
// h ambiguous double heterozygotes.
inline double likelihood(const std::array<double, 4>& base, double h, const std::array<double, 4>& f) {
    double value = 0;
    for (int k = 0; k < 4; k++) if (base[k]) {
        if (f[k] <= 0) return -INFINITY;
        value += base[k] * std::log(f[k]);
    }
    if (h) {
        double ambiguous = f[0] * f[3] + f[1] * f[2];
        if (ambiguous <= 0) return -INFINITY;
        value += h * std::log(ambiguous);
    }
    return value;
}

// Haplotype counts resolved without phase (all but the double heterozygotes).
inline std::array<double, 4> known(const Counts& c) {
    return {double(2 * c[0] + c[1] + c[3]), double(2 * c[2] + c[1] + c[5]),
            double(2 * c[6] + c[3] + c[7]), double(2 * c[8] + c[5] + c[7])};
}

inline double poly(double z, const std::array<double, 4>& v) {
    return ((v[0] * z + v[1]) * z + v[2]) * z + v[3];
}

struct Totals { int n, a, b; };  // individuals and alt-allele counts per locus

inline Totals totals(const Counts& c) {
    Totals t = {0, 0, 0};
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
        t.n += c[3 * i + j]; t.a += i * c[3 * i + j]; t.b += j * c[3 * i + j];
    }
    return t;
}

// Maximum-likelihood diploid fit, oriented so the association is positive:
// for negative D the second locus is flipped and q, base follow.
struct DiploidFit {
    std::array<double, 4> base;
    double ambiguous, p, q, dprime;
};

// c must be canonical, with n > 0 and both alleles observed at both loci.
inline DiploidFit fit_diploid(Counts c, const Totals& t) {
    const int n = t.n;
    auto base = known(c);
    double h = c[4], best = -INFINITY, best_z = 0;
    // All stationary solutions of the ambiguous-double-heterozygote cubic,
    // plus its boundaries. Bracket roots between derivative turning points.
    std::vector<double> cuts = {0., 1.}, candidates = {0., 1.};
    std::array<double, 4> v = {2 * h * h,
        h * (base[0] + base[3] - base[1] - base[2] - 2 * h) - h * h,
        base[0] * base[3] + (base[1] + h) * (base[2] + h) - h * (base[0] + base[3]),
        -base[0] * base[3]};
    if (h) {
        double disc = v[1] * v[1] - 3 * v[0] * v[2];
        if (disc >= 0) for (double s : {-1., 1.}) {
            double z = (-v[1] + s * std::sqrt(disc)) / (3 * v[0]);
            if (z > 0 && z < 1) cuts.push_back(z);
        }
        std::sort(cuts.begin(), cuts.end());
        for (double z : cuts) candidates.push_back(z);
        for (size_t i = 1; i < cuts.size(); i++) {
            double left = cuts[i - 1], right = cuts[i], fl = poly(left, v), fr = poly(right, v);
            if ((fl < 0 && fr > 0) || (fl > 0 && fr < 0)) {
                for (int step = 0; step < 55; step++) {
                    double mid = (left + right) / 2, fm = poly(mid, v);
                    if ((fl < 0 && fm < 0) || (fl > 0 && fm > 0)) { left = mid; fl = fm; }
                    else right = mid;
                }
                candidates.push_back((left + right) / 2);
            }
        }
    }
    for (double z : candidates) {
        std::array<double, 4> f = {(base[0] + h * z) / (2 * n), (base[1] + h * (1 - z)) / (2 * n),
                                   (base[2] + h * (1 - z)) / (2 * n), (base[3] + h * z) / (2 * n)};
        double ll = likelihood(base, h, f);
        if (ll > best) { best = ll; best_z = z; }
    }
    double p = double(t.a) / (2 * n), q = double(t.b) / (2 * n);
    double d = (base[3] + h * best_z) / (2 * n) - p * q;
    if (std::abs(d) < 1e-12) d = 0;
    double limit = d >= 0 ? std::min(p * (1 - q), (1 - p) * q) : std::min(p * q, (1 - p) * (1 - q));
    DiploidFit fit;
    fit.dprime = std::clamp(std::abs(d) / limit, 0., 1.);
    if (d < 0) {
        for (int row = 0; row < 3; row++) std::swap(c[3 * row], c[3 * row + 2]);
        q = 1 - q; base = known(c);
    }
    fit.base = base; fit.ambiguous = h; fit.p = p; fit.q = q;
    return fit;
}

// Genotype bit planes for many markers, for popcount pair tables.
class PackedStates {
public:
    // packed: markers × 3 × width bytes (np.packbits of each state, per marker).
    PackedStates(const uint8_t* packed, int markers, int width)
        : words_((width + 7) / 8), bits_(size_t(markers) * 3 * words_, 0) {
        for (int i = 0; i < markers; i++) for (int g = 0; g < 3; g++)
            std::memcpy(&bits_[(size_t(i) * 3 + g) * words_], packed + (size_t(i) * 3 + g) * width, width);
    }

    // Joint state counts of markers i and j over states 0..max_state.
    Counts pair(int i, int j, int max_state) const {
        Counts c = {};
        for (int a = 0; a <= max_state; a++) for (int b = 0; b <= max_state; b++) {
            const uint64_t* x = &bits_[(size_t(i) * 3 + a) * words_];
            const uint64_t* y = &bits_[(size_t(j) * 3 + b) * words_];
            for (int w = 0; w < words_; w++) c[3 * a + b] += __builtin_popcountll(x[w] & y[w]);
        }
        return c;
    }

private:
    int words_;
    std::vector<uint64_t> bits_;
};

}  // namespace glink
