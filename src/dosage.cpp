// Diploid two-locus likelihood D′ with a normalized grid confidence interval.
// Genotypes are hard calls 0/1/2; missing samples never enter a count table.
#include <unordered_map>

#include "dprime.hpp"

namespace glink {

Interval dosage_confidence(Counts c) {
    // At exactly zero association the positive orientation is unidentified.
    // Choose the same table under allele flips/locus exchange before fitting,
    // so rounding at the independence solution cannot change CI direction.
    c = canonical(c);
    const Totals t = totals(c);
    Interval result = {NAN, NAN, NAN, 0};
    if (t.n < 50 || std::min({t.a, 2 * t.n - t.a, t.b, 2 * t.n - t.b}) < 5) return result;
    const DiploidFit fit = fit_diploid(c, t);
    result.dprime = fit.dprime;
    // Normalize likelihoods on D′ = 0, .01, ..., 1 at the observed allele
    // frequencies, with the MLE association oriented positive.
    const double p = fit.p, q = fit.q;
    double dm = std::min(p * (1 - q), (1 - p) * q), ll[101], mx = -INFINITY;
    for (int k = 0; k <= 100; k++) {
        double z = dm * k / 100.;
        std::array<double, 4> f = {(1 - p) * (1 - q) + z, (1 - p) * q - z, p * (1 - q) - z, p * q + z};
        // Suppress floating point residue at mathematically zero boundaries.
        if (k == 100) {
            if (p * (1 - q) <= (1 - p) * q) f[2] = 0;
            if ((1 - p) * q <= p * (1 - q)) f[1] = 0;
        }
        ll[k] = likelihood(fit.base, fit.ambiguous, f);
        mx = std::max(mx, ll[k]);
    }
    if (!std::isfinite(mx)) return result;
    double total = 0;
    for (auto& x : ll) { x = std::exp(x - mx); total += x; }
    double cumulative = 0;
    int lo = 0, hi = 100;
    bool found = false;
    for (int k = 0; k <= 100; k++) {
        cumulative += ll[k] / total;
        if (!found && cumulative >= .05) { lo = k; found = true; }
        if (cumulative >= .95) { hi = k; break; }
    }
    result.lower = lo / 100.;
    result.upper = hi / 100.;
    result.category = (lo > 70 && hi >= 98) ? 1 : (hi < 90 ? 2 : 0);
    return result;
}

void dosage_classes(const uint8_t* packed, int k, int width, uint8_t* out, int threads) {
    const PackedStates states(packed, k, width);
    #pragma omp parallel num_threads(threads)
    {
        std::unordered_map<Counts, uint8_t, CountsHash> cache;
        #pragma omp for schedule(dynamic, 4)
        for (int i = 0; i < k; i++) for (int j = 0; j <= i; j++) {
            const Counts c = canonical(states.pair(i, j, 2));
            auto it = cache.find(c);
            uint8_t category;
            if (it != cache.end()) category = it->second;
            else {
                category = uint8_t(dosage_confidence(c).category);
                if (cache.size() < 1000000) cache.emplace(c, category);
            }
            out[size_t(i) * k + j] = out[size_t(j) * k + i] = category;
        }
    }
}

void dosage_dprime(const uint8_t* packed, int k, int width, float* out, int threads) {
    const PackedStates states(packed, k, width);
    #pragma omp parallel num_threads(threads)
    {
        std::unordered_map<Counts, float, CountsHash> cache;
        #pragma omp for schedule(dynamic, 4)
        for (int i = 0; i < k; i++) for (int j = 0; j <= i; j++) {
            const Counts c = canonical(states.pair(i, j, 2));
            auto it = cache.find(c);
            float value;
            if (it != cache.end()) value = it->second;
            else {
                const Totals t = totals(c);
                value = (t.n < 50 || std::min({t.a, 2 * t.n - t.a, t.b, 2 * t.n - t.b}) < 5)
                    ? NAN : float(fit_diploid(c, t).dprime);
                if (cache.size() < 1000000) cache.emplace(c, value);
            }
            out[size_t(i) * k + j] = out[size_t(j) * k + i] = value;
        }
    }
}

}  // namespace glink
