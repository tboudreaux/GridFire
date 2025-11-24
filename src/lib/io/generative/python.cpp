#include "gridfire/io/generative/python.h"
#include "fourdst/atomic/atomicSpecies.h"

#include <string>
#include <vector>
#include <format>
#include <fstream>

#include <optional>

#include "gridfire/engine/engine_abstract.h"

namespace {
    template <typename T>
    std::string join(std::vector<T> arr, const std::string& delim) {
        if (arr.empty()) return {};

        size_t total = delim.size() * (arr.size() - 1);
        for (const auto& e : arr) total += std::string(e).size();

        std::string out;
        out.reserve(total);

        for (size_t i = 0; i < arr.size(); i++) {
            out += std::string(arr[i]);
            if (i < arr.size() - 1) {
                out += delim;
            }
        }
        return out;
    }

    std::string format_reaction_function_name(const gridfire::reaction::Reaction& reaction) {
        const std::vector<fourdst::atomic::Species>& reactants = reaction.reactants();
        const std::vector<fourdst::atomic::Species>& products = reaction.products();

        auto format_species_names = [](const std::vector<fourdst::atomic::Species>& arr) -> std::vector<std::string> {
            std::vector<std::string> out;
            for (const auto& r : arr) {
                auto name = std::string(r.name());
                std::erase_if(name, [](const char c) {
                    return c == '-';
                });
                out.push_back(name);
            }
            return out;
        };

        const std::vector<std::string> reactant_names = format_species_names(reactants);
        const std::vector<std::string> product_names = format_species_names(products);

        std::string funcName = std::format(
            "{}_to_{}",
            join<std::string>(reactant_names, "_"),
            join<std::string>(product_names, "_")
        );

        return funcName;
    }

    std::string tfFunc = R"(def tf(T9: float) -> Dict[str, float]:
    """
    Calculate T9^(1/3) & T9^(5/3) & ln(T9)
    """
    T913 = T9 ** (1.0 / 3.0)
    T953 = T9 ** (5.0 / 3.0)
    logT9 = np.log(T9)
    return {
        "T913": T913,
        "T953": T953,
        "logT9": logT9
    }
)";
}

namespace gridfire::io::gen {
    PyFunctionDef exportReactionToPy(const reaction::Reaction& reaction) {
        std::string funcName = format_reaction_function_name(reaction);
        std::optional<std::vector<reaction::RateCoefficientSet>> rateSets = reaction.getRateCoefficients();

        std::string signature = std::format(
        "def {}(T9: float) -> float:",
            funcName
        );
        if (!rateSets || rateSets -> empty()) {
            const std::string funcCode = std::format(
                "{}\n#No Reaclib Rates Defined\n    return 0.0",
                signature,
                funcName
            );
            return {funcName, funcCode, {}};
        }

        std::string commentHeader = std::format(
            "    # Calculates the reaction rate for the reaction: {}",
            reaction.id()
        );

        std::string callTf = R"(    tfv = tf(T9)
    T913 = tfv["T913"]
    T953 = tfv["T953"]
    logT9 = tfv["logT9"]
)";
        std::string initRate = "    rate = 0.0\n";

        std::vector<std::string> rateLines;
        for (const auto&[a0, a1, a2, a3, a4, a5, a6] : *rateSets) {
            std::string rateLine = std::format(
                "    rate += np.exp(({}) + ({} / T9) + ({} / T913) + ({} * T913) + ({} * T9) + ({} * T953) + ({} * logT9))",
                a0,
                a1,
                a2,
                a3,
                a4,
                a5,
                a6
            );
            rateLines.push_back(rateLine);
        }

        const std::string rate_body = join<std::string>(rateLines, "\n");

        const std::string returnCall = "    return rate";

        const std::vector<std::string> functionLines = {
            signature,
            commentHeader,
            callTf,
            initRate,
            rate_body,
            returnCall
        };

        const std::string funcCode = join<std::string>(functionLines, "\n");

        return {funcName, funcCode, {"np"}};

    }

    std::string exportEngineToPy(const engine::DynamicEngine& engine) {
        auto reactions = engine.getNetworkReactions();
        std::vector<std::string> functions;
        functions.emplace_back(R"(import numpy as np
from typing import Dict, List, Tuple, Callable)");
        functions.push_back(tfFunc);
        for (const auto& reaction : reactions) {
            PyFunctionDef funcDef = exportReactionToPy(*reaction);
            functions.push_back(funcDef.func_code);
        }
        return join<std::string>(functions, "\n\n");
    }

    void exportEngineToPy(const engine::DynamicEngine &engine, const std::string &fileName) {
        const std::string funcCode = exportEngineToPy(engine);
        std::ofstream outFile(fileName);
        outFile << funcCode;
        outFile.close();
    }

}