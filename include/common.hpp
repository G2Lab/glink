// Shared support for the C++ glink: errors, JSON and number formatting.
#pragma once

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

namespace glink {

// A user-facing failure; the command line prints it and exits nonzero.
struct Error : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// ---- formatting (util.cpp) ---------------------------------------------------

std::string shortest(double value);            // the fewest digits that read back as exactly this value
std::string thousands(long long value);        // 1234567 -> 1,234,567
std::string fixed3(double value);              // three decimals
// With GLINK_TIMING set, print seconds since start and since the last stage to stderr.
void stage(const char* name);

// np.add.reduce of a contiguous float array: pairwise summation within each
// 8192-element buffer, the buffers accumulated in order.
template <class T>
T numpy_pairwise(const T* a, size_t n) {
    if (n < 8) {
        T res = 0;
        for (size_t i = 0; i < n; i++) res += a[i];
        return res;
    }
    if (n <= 128) {
        T r[8];
        for (int j = 0; j < 8; j++) r[j] = a[j];
        size_t i = 8;
        for (; i < n - (n % 8); i += 8)
            for (int j = 0; j < 8; j++) r[j] += a[i + j];
        T res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; i++) res += a[i];
        return res;
    }
    size_t n2 = n / 2;
    n2 -= n2 % 8;
    return numpy_pairwise(a, n2) + numpy_pairwise(a + n2, n - n2);
}

template <class T>
T numpy_sum(const T* a, size_t n) {
    T res = 0;
    for (size_t i = 0; i < n; i += 8192) res += numpy_pairwise(a + i, std::min<size_t>(8192, n - i));
    return res;
}

// ---- paths and files (util.cpp) ---------------------------------------------

bool exists(const std::string& path);
std::string resolve(const std::string& path);  // pathlib.Path.resolve(), non-strict
std::string basename(const std::string& path);
std::string dirname(const std::string& path);
void make_dirs(const std::string& path);
std::string read_text(const std::string& path);
std::vector<std::string> split_ws(const std::string& text);
bool ends_with(const std::string& s, const std::string& suffix);
bool starts_with(const std::string& s, const std::string& prefix);

// An exclusive flock on `path`, released on destruction.
class FileLock {
public:
    explicit FileLock(const std::string& path);
    ~FileLock();
    FileLock(const FileLock&) = delete;
    FileLock& operator=(const FileLock&) = delete;

private:
    int fd_;
};

// A temporary file name in `directory` (created empty), as tempfile.NamedTemporaryFile.
std::string temporary_file(const std::string& directory, const std::string& prefix, const std::string& suffix);

// Lines of a plain or gzip text file, read in large chunks.
class LineReader {
public:
    explicit LineReader(const std::string& path);
    ~LineReader();
    bool next(std::string& line);  // line without its '\n'
    void rewind();

private:
    bool fill();
    void* gz_;
    std::vector<char> buffer_;
    size_t start_ = 0, end_ = 0;
    bool eof_ = false;
};

// ---- JSON --------------------------------------------------------------------

using Json = nlohmann::json;

// The value at key, or nullptr when it is absent or j is not an object.
inline const Json* json_find(const Json& j, const std::string& key) {
    const auto it = j.find(key);
    return it == j.end() ? nullptr : &*it;
}

// Named columns of a store, widened to int64_t.
using ColumnMap = std::map<std::string, std::vector<int64_t>>;

// ---- graph identity (util.cpp) ----------------------------------------------

// SHA-256 of a file's uncompressed contents: what caches compare, so a copied or moved file still
// matches. Each file is hashed once per run; do not ask about a file this run will write again.
const std::string& file_sha256(const std::string& path);
std::string sha256_hex(const std::string& text);

}  // namespace glink
