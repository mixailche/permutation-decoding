#pragma once

#include <stdexcept>
#include <unordered_map>
#include <string>

namespace running {
    class ArgsReader {
    public:
        using key_t = std::string;
        using value_t = std::string;

        ArgsReader(int argc, char** argv);
        ArgsReader(int argc, char** argv, std::unordered_map<key_t, value_t> defaultValues);

        bool HasArg(const key_t& key) const;

        template <typename Parse>
        decltype(auto) GetArg(const key_t& key, Parse parse) const;

        std::string GetStringArg(const key_t& key) const;
        size_t GetNumberArg(const key_t& key) const;
        double GetDoubleArg(const key_t& key) const;
        std::vector<size_t> GetNumberListArg(const key_t& key) const;

    private:
        std::unordered_map<key_t, value_t> mArgs;

        constexpr static char KEY_PREFIX = '-';
        constexpr static char LIST_SEPARATOR = ',';
    };
}

template <typename Parse>
decltype(auto) running::ArgsReader::GetArg(const key_t& key, Parse parse) const
{
    auto it = mArgs.find(key);
    if (it == mArgs.end()) {
        throw std::invalid_argument("Argument not present: `" + key + "`");
    }
    return parse(it->second.c_str());
}
