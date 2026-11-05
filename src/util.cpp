// Formatting, hashing, path and file helpers.
#include <charconv>
#include <chrono>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fcntl.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <zlib.h>

#include <openssl/evp.h>


#include "common.hpp"

namespace glink {

// ---- number formatting -----------------------------------------------------

std::string shortest(double v) {
    char buf[64];
    return std::string(buf, std::to_chars(buf, buf + sizeof buf, v).ptr);
}

std::string thousands(long long v) {
    std::string digits = std::to_string(v < 0 ? -v : v), out;
    for (size_t i = 0; i < digits.size(); i++) {
        if (i && (digits.size() - i) % 3 == 0) out += ',';
        out += digits[i];
    }
    return (v < 0 ? "-" : "") + out;
}

std::string fixed3(double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.3f", v);
    return buf;
}

void stage(const char* name) {
    static const bool on = std::getenv("GLINK_TIMING") != nullptr;
    if (!on) return;
    static const auto start = std::chrono::steady_clock::now();
    static auto last = start;
    auto now = std::chrono::steady_clock::now();
    std::fprintf(stderr, "[timing] %-28s %8.2f s (+%.2f)\n", name,
                 std::chrono::duration<double>(now - start).count(), std::chrono::duration<double>(now - last).count());
    last = now;
}

// ---- paths -----------------------------------------------------------------

bool exists(const std::string& path) {
    struct stat st;
    return ::stat(path.c_str(), &st) == 0;
}

std::string basename(const std::string& path) {
    size_t slash = path.find_last_of('/');
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

std::string dirname(const std::string& path) {
    size_t slash = path.find_last_of('/');
    if (slash == std::string::npos) return ".";
    return slash == 0 ? "/" : path.substr(0, slash);
}

std::string resolve(const std::string& path) {
    char buf[PATH_MAX];
    if (::realpath(path.c_str(), buf)) return buf;
    // Non-strict: resolve the longest existing prefix, append the rest.
    std::string absolute = path;
    if (absolute.empty() || absolute[0] != '/') {
        if (!::getcwd(buf, sizeof buf)) throw Error("getcwd failed");
        absolute = std::string(buf) + "/" + path;
    }
    std::vector<std::string> parts, kept;
    size_t start = 0;
    while (start <= absolute.size()) {
        size_t end = absolute.find('/', start);
        if (end == std::string::npos) end = absolute.size();
        std::string part = absolute.substr(start, end - start);
        if (!part.empty() && part != ".") {
            if (part == "..") { if (!kept.empty()) kept.pop_back(); }
            else kept.push_back(part);
        }
        start = end + 1;
    }
    std::string prefix;
    size_t i = 0;
    for (; i < kept.size(); i++) {
        std::string next = prefix + "/" + kept[i];
        if (!::realpath(next.c_str(), buf)) break;
        prefix = buf;
    }
    for (; i < kept.size(); i++) prefix += "/" + kept[i];
    return prefix.empty() ? "/" : prefix;
}

void make_dirs(const std::string& path) {
    if (path.empty() || path == "." || path == "/" || exists(path)) return;
    make_dirs(dirname(path));
    if (::mkdir(path.c_str(), 0777) != 0 && errno != EEXIST) throw Error("Cannot create directory " + path);
}

std::string read_text(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) throw Error("Cannot read " + path);
    std::string out;
    char buf[1 << 16];
    size_t n;
    while ((n = std::fread(buf, 1, sizeof buf, f)) > 0) out.append(buf, n);
    std::fclose(f);
    return out;
}

std::vector<std::string> split_ws(const std::string& text) {
    std::vector<std::string> out;
    size_t i = 0;
    while (i < text.size()) {
        while (i < text.size() && std::isspace((unsigned char)text[i])) i++;
        size_t j = i;
        while (j < text.size() && !std::isspace((unsigned char)text[j])) j++;
        if (j > i) out.push_back(text.substr(i, j - i));
        i = j;
    }
    return out;
}

bool ends_with(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool starts_with(const std::string& s, const std::string& prefix) {
    return s.compare(0, prefix.size(), prefix) == 0;
}

FileLock::FileLock(const std::string& path) {
    fd_ = ::open(path.c_str(), O_WRONLY | O_CREAT | O_APPEND, 0644);
    if (fd_ < 0) throw Error("Cannot open lock " + path);
    if (::flock(fd_, LOCK_EX) != 0) throw Error("Cannot lock " + path);
}

FileLock::~FileLock() {
    ::flock(fd_, LOCK_UN);
    ::close(fd_);
}

std::string temporary_file(const std::string& directory, const std::string& prefix, const std::string& suffix) {
    std::string pattern = directory + "/" + prefix + "XXXXXX" + suffix;
    std::vector<char> buf(pattern.begin(), pattern.end());
    buf.push_back(0);
    int fd = ::mkstemps(buf.data(), int(suffix.size()));
    if (fd < 0) throw Error("Cannot create a temporary file in " + directory);
    ::close(fd);
    return buf.data();
}

namespace {

std::string hex(const unsigned char* digest, unsigned int size) {
    static const char digits[] = "0123456789abcdef";
    std::string out;
    for (unsigned int i = 0; i < size; i++) { out += digits[digest[i] >> 4]; out += digits[digest[i] & 15]; }
    return out;
}

}  // namespace

const std::string& file_sha256(const std::string& path) {
    static std::unordered_map<std::string, std::string> hashes;
    const std::string resolved = resolve(path);
    const auto known = hashes.find(resolved);
    if (known != hashes.end()) return known->second;
    if (!exists(resolved)) throw Error("No such file: " + path);
    gzFile f = gzopen(resolved.c_str(), "rb");
    if (!f) throw Error("Cannot open " + path);
    EVP_MD_CTX* hash = EVP_MD_CTX_new();
    if (!hash || !EVP_DigestInit_ex(hash, EVP_sha256(), nullptr)) throw Error("Cannot start SHA-256");
    std::vector<unsigned char> buffer(1 << 22);
    int n;
    while ((n = gzread(f, buffer.data(), unsigned(buffer.size()))) > 0) EVP_DigestUpdate(hash, buffer.data(), size_t(n));
    gzclose(f);
    unsigned char digest[EVP_MAX_MD_SIZE];
    unsigned int size = 0;
    EVP_DigestFinal_ex(hash, digest, &size);
    EVP_MD_CTX_free(hash);
    if (n < 0) throw Error("Read error: " + path);
    return hashes[resolved] = hex(digest, size);
}

std::string sha256_hex(const std::string& text) {
    unsigned char digest[EVP_MAX_MD_SIZE];
    unsigned int size = 0;
    if (!EVP_Digest(text.data(), text.size(), digest, &size, EVP_sha256(), nullptr)) throw Error("Cannot compute SHA-256");
    return hex(digest, size);
}

// ---- line reader -------------------------------------------------------------

LineReader::LineReader(const std::string& path) : buffer_(1 << 22) {
    gzFile f = gzopen(path.c_str(), "rb");
    if (!f) throw Error("Cannot open " + path);
    gzbuffer(f, 1 << 20);
    gz_ = f;
}

LineReader::~LineReader() { gzclose(static_cast<gzFile>(gz_)); }

void LineReader::rewind() {
    gzrewind(static_cast<gzFile>(gz_));
    start_ = end_ = 0;
    eof_ = false;
}

bool LineReader::fill() {
    if (eof_) return false;
    if (start_ > 0) {
        std::memmove(buffer_.data(), buffer_.data() + start_, end_ - start_);
        end_ -= start_;
        start_ = 0;
    }
    if (end_ == buffer_.size()) buffer_.resize(buffer_.size() * 2);
    int n = gzread(static_cast<gzFile>(gz_), buffer_.data() + end_, unsigned(buffer_.size() - end_));
    if (n < 0) throw Error("Read error");
    if (n == 0) { eof_ = true; return false; }
    end_ += n;
    return true;
}

bool LineReader::next(std::string& line) {
    size_t scanned = start_;
    while (true) {
        const void* nl = std::memchr(buffer_.data() + scanned, '\n', end_ - scanned);
        if (nl) {
            size_t pos = static_cast<const char*>(nl) - buffer_.data();
            line.assign(buffer_.data() + start_, pos - start_);
            start_ = pos + 1;
            return true;
        }
        size_t offset = end_ - start_;
        if (!fill()) {
            if (start_ == end_) return false;
            line.assign(buffer_.data() + start_, end_ - start_);
            start_ = end_;
            return true;
        }
        scanned = start_ + offset;
    }
}

}  // namespace glink
