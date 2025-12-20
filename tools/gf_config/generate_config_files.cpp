#include "fourdst/config/config.h"
#include "gridfire/config/config.h"
#include <source_location>
#include <filesystem>

#include  "CLI/CLI.hpp"

consteval std::string_view strip_namespaces(const std::string_view fullName) {
    const size_t pos = fullName.rfind("::");
    if (pos == std::string_view::npos) {
        return fullName;
    }
    return fullName.substr(pos + 2);
}

template <typename T>
consteval std::string_view get_type_name() {
    constexpr std::string_view name = std::source_location::current().function_name();
    const auto pos = name.find("T = ");
    if (pos == std::string_view::npos) return name;
    const auto start = pos + 4;
    const auto end = name.rfind(']');
    return name.substr(start, end - start);
}

int main(int argc, char** argv) {
    CLI::App app{"GridFire Sandbox Application."};

    std::string outputPath = ".";

    app.add_option("-p,--path", outputPath, "path to save generated config files (default: current directory)");

    CLI11_PARSE(app, argc, argv);

    const std::filesystem::path outPath(outputPath);
    if (!std::filesystem::exists(outPath)) {
        std::cerr << "Error: The specified path does not exist: " << outputPath << std::endl;
        return 1;
    }

    fourdst::config::Config<gridfire::config::GridFireConfig> configConfig;
    const std::string_view name = strip_namespaces(get_type_name<gridfire::config::GridFireConfig>());

    const std::string defaultConfigFilePath = (outPath / (std::string(name) + ".toml")).string();
    const std::string schemaFilePath = (outPath / (std::string(name) + ".schema.json")).string();
    configConfig.save(defaultConfigFilePath);
    configConfig.save_schema(schemaFilePath);
}