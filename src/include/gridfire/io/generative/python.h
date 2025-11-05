#pragma once

#include <string>
#include <vector>

#include "gridfire/reaction/reaction.h"
#include "gridfire/engine/engine_abstract.h"

namespace gridfire::io::gen {
    struct PyFunctionDef {
        std::string func_name;
        std::string func_code;
        std::vector<std::string> module_req;
    };

    PyFunctionDef exportReactionToPy(const reaction::Reaction& reaction);

    std::string exportEngineToPy(const DynamicEngine& engine);

    void exportEngineToPy(const DynamicEngine& engine, const std::string& fileName);
}
