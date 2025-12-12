#pragma once

#include <string>
#include <vector>

#include "gridfire/reaction/reaction.h"
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/scratchpads/blob.h"

/**
 * @brief Namespace for generative input/output functionalities.
 *
 * This namespace contains functions and structures related to exporting
 * reactions and engines to Python code. It provides tools to convert
 * internal representations of reactions and engines into Python functions
 * and scripts that can be used for further analysis or integration with
 * Python-based workflows.
 */
namespace gridfire::io::gen {
    /**
     * @brief Structure representing a Python function definition.
     *
     * This structure holds the name, code, and module requirements for
     * a Python function that represents a reaction.
     */
    struct PyFunctionDef {
        std::string func_name;
        std::string func_code;
        std::vector<std::string> module_req;
    };

    /**
     * @brief Exports a reaction to a Python function definition.
     *
     * @param reaction The reaction to export.
     * @return A PyFunctionDef structure containing the Python function
     *         representation of the reaction.
     *
     * This function converts the given reaction into a Python function
     * definition, including the necessary code and module requirements.
     */
    PyFunctionDef exportReactionToPy(const reaction::Reaction& reaction);

    /**
     * @brief Exports a dynamic engine to a Python script.
     *
     * @param engine The dynamic engine to export.
     * @return A string containing the Python script representation of the engine.
     *
     * This function converts the given dynamic engine into a Python script
     * that can be used to recreate the engine's functionality in Python.
     */
    std::string exportEngineToPy(engine::scratch::StateBlob& ctx, engine::DynamicEngine& engine);

    /**
     * @brief Exports a dynamic engine to a Python file.
     *
     * @param ctx
     * @param engine The dynamic engine to export.
     * @param fileName The name of the file to write the Python script to.
     *
     * This function writes the Python script representation of the given
     * dynamic engine to the specified file.
     */
    void exportEngineToPy(engine::scratch::StateBlob &ctx, const engine::DynamicEngine& engine, const std::string& fileName);
}
