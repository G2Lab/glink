// A small argparse-style option parser shared by the glink command lines.
#pragma once

#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

namespace glink {

class Options {
public:
    struct Spec {
        std::string name;  // without leading dashes
        bool value;        // takes a value
        std::string help;
    };

    Options(std::string prog, std::string description) : prog_(std::move(prog)), description_(std::move(description)) {}

    void flag(const std::string& name, const std::string& help) { specs_.push_back({name, false, help}); }
    void option(const std::string& name, const std::string& help) { specs_.push_back({name, true, help}); }

    [[noreturn]] void error(const std::string& message) const {
        std::fprintf(stderr, "usage: %s [options]\n%s: error: %s\n", prog_.c_str(), prog_.c_str(), message.c_str());
        std::exit(2);
    }

    void parse(int argc, char** argv) {
        for (int i = 1; i < argc; i++) {
            std::string a = argv[i];
            if (a == "-h" || a == "--help") { help(); std::exit(0); }
            if (a.rfind("--", 0) != 0) error("unrecognized arguments: " + a);
            std::string name = a.substr(2), value;
            bool inline_value = false;
            size_t eq = name.find('=');
            if (eq != std::string::npos) { value = name.substr(eq + 1); name = name.substr(0, eq); inline_value = true; }
            const Spec* spec = nullptr;
            for (auto& s : specs_) if (s.name == name) spec = &s;
            if (!spec) error("unrecognized arguments: " + a);
            if (spec->value) {
                if (!inline_value) {
                    if (i + 1 >= argc) error("argument --" + name + ": expected one argument");
                    value = argv[++i];
                }
                values_[name].push_back(value);
            } else {
                if (inline_value) error("argument --" + name + ": ignored explicit argument '" + value + "'");
                values_[name].push_back("");
            }
        }
    }

    bool has(const std::string& name) const { return values_.count(name) > 0; }
    std::string get(const std::string& name, const std::string& fallback = "") const {
        auto it = values_.find(name);
        return it == values_.end() ? fallback : it->second.back();
    }
    std::vector<std::string> all(const std::string& name) const {
        auto it = values_.find(name);
        return it == values_.end() ? std::vector<std::string>{} : it->second;
    }
    // The last of several mutually exclusive flags given, e.g. --collapse/--no-collapse.
    int last_index(const std::vector<std::string>& names, char** argv, int argc) const {
        int found = -1;
        for (int i = 1; i < argc; i++)
            for (size_t k = 0; k < names.size(); k++)
                if (std::string(argv[i]) == "--" + names[k]) found = int(k);
        return found;
    }
    int integer(const std::string& name, int fallback) const {
        if (!has(name)) return fallback;
        char* end;
        const std::string text = get(name);  // `end` points into this, so it must outlive the check
        long v = std::strtol(text.c_str(), &end, 10);
        if (*end || get(name).empty()) error("argument --" + name + ": invalid int value: '" + get(name) + "'");
        return int(v);
    }
    double number(const std::string& name, double fallback) const {
        if (!has(name)) return fallback;
        char* end;
        const std::string text = get(name);  // `end` points into this, so it must outlive the check
        double v = std::strtod(text.c_str(), &end);
        if (*end || get(name).empty()) error("argument --" + name + ": invalid float value: '" + get(name) + "'");
        return v;
    }

    void help() const {
        std::printf("usage: %s [options]\n\n%s\n\noptions:\n", prog_.c_str(), description_.c_str());
        for (auto& s : specs_)
            std::printf("  --%s%s\n      %s\n", s.name.c_str(), s.value ? " VALUE" : "", s.help.c_str());
    }

private:
    std::string prog_, description_;
    std::vector<Spec> specs_;
    std::map<std::string, std::vector<std::string>> values_;
};

}  // namespace glink
