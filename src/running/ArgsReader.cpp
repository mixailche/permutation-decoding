#include "ArgsReader.hpp"
#include "ArgsReader.hpp"
#include "running/ArgsReader.hpp"

using running::ArgsReader;

ArgsReader::ArgsReader(int argc, char** argv)
    : ArgsReader(argc, argv, std::unordered_map<key_t, value_t>())
{}

ArgsReader::ArgsReader(int argc, char** argv, std::unordered_map<key_t, value_t> defaultValues)
    : mArgs(std::move(defaultValues))
{
    for (size_t index = 1; index + 1 < argc; index++) {
        std::string key = argv[index];
        if (key.empty() || key[0] != KEY_PREFIX) {
            continue;
        }
        mArgs.insert_or_assign(key.substr(1), std::string(argv[++index]));
    }
}

bool ArgsReader::HasArg(const key_t& key) const
{
    return mArgs.contains(key);
}

std::string ArgsReader::GetStringArg(const key_t& key) const
{
    return GetArg(key, [](const char* str) { return std::string(str); });
}

size_t ArgsReader::GetNumberArg(const key_t& key) const
{
    return GetArg(key, [](const char* str) { return std::stoull(str); });
}

double ArgsReader::GetDoubleArg(const key_t& key) const
{
    return GetArg(key, [](const char* str) { return std::stod(str); });
}

std::vector<size_t> ArgsReader::GetNumberListArg(const key_t& key) const
{
    auto raw = GetStringArg(key);
    auto total = raw.size();

    std::vector<size_t> result;
    size_t offset = 0;
    while (offset < total) {
        size_t i;
        for (i = 0; offset + i < total && std::isdigit(raw[offset + i]); i++);
        if (raw[i] == LIST_SEPARATOR) {
            result.push_back(std::stoull(raw.substr(offset, i + 1)));
            offset += i + 1;
        }
        else if (offset + i == total) {
            break;
        }
        else {
            throw std::invalid_argument("Wrong number format");
        }
    }

    return result;
}
