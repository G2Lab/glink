// r² and D′ between every pair of rows of a node matrix, and PLINK's reports of them.
#pragma once

#include "common.hpp"

namespace glink {

struct NodeMatrix;

// An n x n matrix, row by row; NaN where a pair is undefined.
std::vector<float> compute_r2(const NodeMatrix& matrix, int min_pairs, int min_allele_count, int threads);
// Phased/presence: signed D′. Dosage: genotype-likelihood |D′|, computed in single precision.
std::vector<double> compute_dprime(const NodeMatrix& matrix, int threads);

// What plink --r2 square --write-snplist writes: PATH (.ld), the matrix as tab-separated text, and
// the row IDs in order beside it, in the same name with .snplist. Numbers have six significant digits.
void write_r2_ld(const std::string& path, const std::vector<float>& r2, const std::vector<std::string>& ids);

// The same n x n matrix as raw little-endian float32, row by row, in PATH (.ld.bin), with the row IDs
// beside it in the same name with .snplist. Much faster than text, and keeps every digit of r².
void write_ld_binary(const std::string& path, const std::vector<float>& values, const std::vector<std::string>& ids);
void write_ld_binary(const std::string& path, const std::vector<double>& values, const std::vector<std::string>& ids);

// What plink --r2 inter-chr dprime writes: a table of every pair, CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2 DP.
// Unlike PLINK, a pair with undefined r² stays in the table as nan, because its D′ may be defined.
// A row without a coordinate has CHR and BP 0.
void write_dprime_ld(const std::string& path, const NodeMatrix& matrix, const std::vector<float>& r2,
                     const std::vector<double>& dprime);

}  // namespace glink
