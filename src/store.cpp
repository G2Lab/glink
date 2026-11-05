#include <filesystem>
#include <unistd.h>

#include "store.hpp"

namespace fs = std::filesystem;

namespace glink {

namespace {

constexpr const char* FORMAT = "glink-store";
constexpr int VERSION = 1;
constexpr const char* METADATA = "_metadata";
const HighFive::SilenceHDF5 silence;  // failures surface as exceptions; keep HDF5 from also printing its stack

}  // namespace

Store::Store(const std::string& path) : path_(path) {
    const Error invalid(path_ + " is not a version " + std::to_string(VERSION) + " glink store; rebuild it");
    if (!fs::is_regular_file(path_)) throw invalid;
    try {
        file_.emplace(path_, HighFive::File::ReadOnly);
        if (!file_->hasAttribute("format") || file_->getAttribute("format").read<std::string>() != FORMAT ||
            !file_->hasAttribute("version") || file_->getAttribute("version").read<int>() != VERSION)
            throw invalid;
        metadata_ = Json::parse(file_->getDataSet(METADATA).read<std::string>());
    } catch (const HighFive::Exception&) {
        throw invalid;
    }
}

HighFive::DataSet Store::dataset(const std::string& name) const {
    if (!file_->exist(name)) throw Error(path_ + " has no array '" + name + "'");
    return file_->getDataSet(name);
}

bool Store::matches(const Json& expected) const {
    for (auto& [key, value] : expected.items()) {
        const Json* have = json_find(metadata_, key);
        if (!have || *have != value) return false;
    }
    return true;
}

std::vector<std::string> Store::integer_arrays() const {
    std::vector<std::string> names;
    for (const std::string& name : file_->listObjectNames())
        if (name != METADATA && file_->getDataSet(name).getDataType().getClass() == HighFive::DataTypeClass::Integer)
            names.push_back(name);
    return names;
}

std::vector<size_t> Store::shape(const std::string& name) const { return dataset(name).getDimensions(); }

size_t Store::item_size(const std::string& name) const { return dataset(name).getDataType().getSize(); }

std::vector<int64_t> Store::ints(const std::string& name) const { return dataset(name).read<std::vector<int64_t>>(); }

std::vector<int64_t> Store::ints(const std::string& name, const std::vector<int64_t>& rows) const {
    if (rows.empty()) return {};
    return dataset(name).select(HighFive::ElementSet(std::vector<size_t>(rows.begin(), rows.end()))).read<std::vector<int64_t>>();
}

std::vector<std::string> Store::strings(const std::string& name) const {
    const auto data = dataset(name);
    if (!data.getElementCount()) return {};  // HighFive cannot read an empty string dataset
    return data.read<std::vector<std::string>>();
}

void Store::read_bytes(const std::string& name, uint8_t* out) const { dataset(name).read_raw(out); }

void Store::read_rows(const std::string& name, const std::vector<int64_t>& rows, uint8_t* out) const {
    const auto data = dataset(name);
    const size_t width = data.getDimensions().at(1);
    for (size_t i = 0; i < rows.size(); i++)
        data.select({size_t(rows[i]), 0}, {1, width}).read_raw(out + i * width);
}

StoreWriter::StoreWriter(const std::string& path) : path_(path), metadata_(Json::object()) {
    make_dirs(dirname(path_));
    partial_ = path_ + ".partial-" + std::to_string(::getpid());
    file_.emplace(partial_, HighFive::File::Truncate);
}

StoreWriter::~StoreWriter() {
    file_.reset();
    std::error_code ignored;
    if (!published_) fs::remove(partial_, ignored);
}

void StoreWriter::publish() {
    file_->createAttribute("format", std::string(FORMAT));
    file_->createAttribute("version", VERSION);
    file_->createDataSet(METADATA, metadata_.dump());
    file_->flush();
    file_.reset();
    const Store written(partial_);  // what was written reads back
    if (fs::exists(path_)) Store existing(path_);  // refuses to replace anything but a store
    fs::rename(partial_, path_);
    published_ = true;
}

}  // namespace glink
