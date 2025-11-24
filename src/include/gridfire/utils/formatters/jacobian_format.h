#pragma once

#include <format>

#include "gridfire/engine/types/jacobian.h"



template <>
struct std::formatter<gridfire::engine::NetworkJacobian> {
    static constexpr auto parse(const std::format_parse_context& ctx) {
        return ctx.begin();
    }

    static auto format(const gridfire::engine::NetworkJacobian& obj, std::format_context& ctx) {
        return std::format_to(ctx.out(), "NetworkJacobian(shape=({}, {}), nnz={}, rank={}, singular={})",
            std::get<0>(obj.shape()),
            std::get<1>(obj.shape()),
            obj.nnz(),
            obj.rank(),
            obj.singular() ? "True" : "False"
        );
    }
};
