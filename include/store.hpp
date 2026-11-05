// A store: one HDF5 file holding glink's arrays, string lists and JSON metadata. It is written
// beside its path and renamed into place, so a store that exists is complete.
#pragma once

#include <optional>

#include <highfive/highfive.hpp>

#include "common.hpp"

namespace glink {

class Store {
public:
    explicit Store(const std::string& path);  // throws Error if this is not a glink store

    const std::string& path() const { return path_; }
    const Json& metadata() const { return metadata_; }
    bool matches(const Json& expected) const;  // whether the metadata holds each of these keys with these values
    std::vector<std::string> integer_arrays() const;
    std::vector<size_t> shape(const std::string& name) const;
    size_t item_size(const std::string& name) const;  // bytes per element
    std::vector<int64_t> ints(const std::string& name) const;  // any integer type, widened
    std::vector<int64_t> ints(const std::string& name, const std::vector<int64_t>& rows) const;  // these rows, in order
    std::vector<std::string> strings(const std::string& name) const;
    void read_bytes(const std::string& name, uint8_t* out) const;  // a whole byte array
    void read_rows(const std::string& name, const std::vector<int64_t>& rows, uint8_t* out) const;  // of a 2-D byte array

private:
    HighFive::DataSet dataset(const std::string& name) const;
    std::string path_;
    std::optional<HighFive::File> file_;
    Json metadata_;
};

// Builds a store in a temporary file beside `path`; publish() moves it into place, and without
// publish() nothing is left behind. An existing store at `path` is replaced; anything else there is an error.
class StoreWriter {
public:
    explicit StoreWriter(const std::string& path);
    ~StoreWriter();
    StoreWriter(const StoreWriter&) = delete;
    StoreWriter& operator=(const StoreWriter&) = delete;

    template <class T>
    void add(const std::string& name, const std::vector<T>& values) { file_->createDataSet(name, values); }
    template <class T>
    void add(const std::string& name, size_t rows, size_t cols, const T* data) {
        auto dataset = file_->createDataSet<T>(name, HighFive::DataSpace({rows, cols}));
        if (rows && cols) dataset.write_raw(data);
    }
    void set_metadata(const Json& metadata) { metadata_ = metadata; }
    void publish();

private:
    std::string path_, partial_;
    std::optional<HighFive::File> file_;
    Json metadata_;
    bool published_ = false;
};

}  // namespace glink
