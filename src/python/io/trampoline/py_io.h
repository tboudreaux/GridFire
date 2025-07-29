#pragma once

#include "gridfire/io/io.h"

class PyNetworkFileParser final : public gridfire::io::NetworkFileParser {
    gridfire::io::ParsedNetworkData parse(const std::string &filename) const override;
};