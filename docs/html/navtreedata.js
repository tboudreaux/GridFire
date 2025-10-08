/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "GridFire", "index.html", [
    [ "Introduction", "index.html#autotoc_md0", [
      [ "Design Philosophy and Workflow", "index.html#autotoc_md1", null ],
      [ "Funding", "index.html#autotoc_md2", null ]
    ] ],
    [ "Usage", "index.html#autotoc_md3", [
      [ "Python installation", "index.html#autotoc_md4", [
        [ "pypi", "index.html#autotoc_md5", null ],
        [ "source", "index.html#autotoc_md6", null ],
        [ "source for developers", "index.html#autotoc_md7", null ]
      ] ],
      [ "Automatic Build and Installation", "index.html#autotoc_md8", [
        [ "Script Build and Installation Instructions", "index.html#autotoc_md9", null ],
        [ "Currently, known good platforms", "index.html#autotoc_md10", null ]
      ] ],
      [ "Manual Build Instructions", "index.html#autotoc_md11", [
        [ "Prerequisites", "index.html#autotoc_md12", [
          [ "Required", "index.html#autotoc_md13", null ],
          [ "Optional", "index.html#autotoc_md14", null ]
        ] ],
        [ "Install Scripts", "index.html#autotoc_md15", [
          [ "Ease of Installation", "index.html#autotoc_md16", null ],
          [ "Reproducibility", "index.html#autotoc_md17", null ],
          [ "Examples", "index.html#autotoc_md18", [
            [ "TUI config and saving", "index.html#autotoc_md19", null ],
            [ "TUI config loading and meson setup", "index.html#autotoc_md20", null ],
            [ "CLI config loading, setup, and build", "index.html#autotoc_md21", null ]
          ] ]
        ] ],
        [ "Dependency Installation on Common Platforms", "index.html#autotoc_md22", null ],
        [ "Building the C++ Library", "index.html#autotoc_md23", [
          [ "Clang vs. GCC", "index.html#autotoc_md24", null ]
        ] ],
        [ "Installing the Library", "index.html#autotoc_md25", null ],
        [ "Minimum compiler versions", "index.html#autotoc_md26", null ]
      ] ],
      [ "Code Architecture and Logical Flow", "index.html#autotoc_md27", null ],
      [ "Engines", "index.html#autotoc_md28", [
        [ "GraphEngine", "index.html#autotoc_md29", null ],
        [ "GraphEngine Configuration Options", "index.html#autotoc_md30", null ],
        [ "Available Partition Functions", "index.html#autotoc_md31", null ],
        [ "AutoDiff", "index.html#autotoc_md32", null ]
      ] ],
      [ "Reaclib in GridFire", "index.html#autotoc_md33", null ],
      [ "Engine Views", "index.html#autotoc_md34", [
        [ "A Note about composability", "index.html#autotoc_md35", null ]
      ] ],
      [ "Numerical Solver Strategies", "index.html#autotoc_md36", [
        [ "NetworkSolverStrategy<EngineT>:", "index.html#autotoc_md37", null ],
        [ "NetIn and NetOut", "index.html#autotoc_md38", [
          [ "NetIn", "index.html#autotoc_md39", null ],
          [ "NetOut", "index.html#autotoc_md40", null ]
        ] ],
        [ "DirectNetworkSolver (Implicit Rosenbrock Method)", "index.html#autotoc_md41", null ],
        [ "Algorithmic Workflow in DirectNetworkSolver", "index.html#autotoc_md42", null ],
        [ "Future Solver Implementations", "index.html#autotoc_md43", null ]
      ] ],
      [ "Python Extensibility", "index.html#autotoc_md44", null ]
    ] ],
    [ "Usage Examples", "index.html#autotoc_md45", [
      [ "C++", "index.html#autotoc_md46", [
        [ "GraphEngine Initialization", "index.html#autotoc_md47", null ],
        [ "Adaptive Network View", "index.html#autotoc_md48", null ],
        [ "Composition Initialization", "index.html#autotoc_md49", null ],
        [ "Common Workflow Example", "index.html#autotoc_md50", [
          [ "Workflow Components and Effects", "index.html#autotoc_md51", null ]
        ] ],
        [ "Callback Example", "index.html#autotoc_md52", [
          [ "Callback Context", "index.html#autotoc_md53", null ]
        ] ]
      ] ],
      [ "Python", "index.html#autotoc_md54", [
        [ "Common Workflow Example", "index.html#autotoc_md55", null ],
        [ "Python callbacks", "index.html#autotoc_md56", null ]
      ] ]
    ] ],
    [ "Related Projects", "index.html#autotoc_md57", null ],
    [ "GridFire Python Usage Guide", "md_docs_2static_2usage.html", [
      [ "Installation", "md_docs_2static_2usage.html#autotoc_md60", [
        [ "1.1 PyPI Release", "md_docs_2static_2usage.html#autotoc_md61", null ],
        [ "1.2 Development from Source", "md_docs_2static_2usage.html#autotoc_md62", null ]
      ] ],
      [ "Why These Engines and Views?", "md_docs_2static_2usage.html#autotoc_md64", null ],
      [ "Step-by-Step Example", "md_docs_2static_2usage.html#autotoc_md66", null ],
      [ "Visualizing Reaction Networks", "md_docs_2static_2usage.html#autotoc_md68", null ],
      [ "Beyond the Basics", "md_docs_2static_2usage.html#autotoc_md70", null ]
    ] ],
    [ "Namespaces", "namespaces.html", [
      [ "Namespace List", "namespaces.html", "namespaces_dup" ],
      [ "Namespace Members", "namespacemembers.html", [
        [ "All", "namespacemembers.html", null ],
        [ "Functions", "namespacemembers_func.html", null ],
        [ "Variables", "namespacemembers_vars.html", null ],
        [ "Typedefs", "namespacemembers_type.html", null ],
        [ "Enumerations", "namespacemembers_enum.html", null ],
        [ "Enumerator", "namespacemembers_eval.html", null ]
      ] ]
    ] ],
    [ "Concepts", "concepts.html", "concepts" ],
    [ "Classes", "annotated.html", [
      [ "Class List", "annotated.html", "annotated_dup" ],
      [ "Class Index", "classes.html", null ],
      [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
      [ "Class Members", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Functions", "functions_func.html", "functions_func" ],
        [ "Variables", "functions_vars.html", "functions_vars" ],
        [ "Typedefs", "functions_type.html", null ],
        [ "Enumerations", "functions_enum.html", null ],
        [ "Enumerator", "functions_eval.html", null ],
        [ "Related Symbols", "functions_rela.html", null ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "File Members", "globals.html", [
        [ "All", "globals.html", null ],
        [ "Functions", "globals_func.html", null ],
        [ "Variables", "globals_vars.html", null ],
        [ "Typedefs", "globals_type.html", null ],
        [ "Macros", "globals_defs.html", null ]
      ] ]
    ] ],
    [ "Examples", "examples.html", "examples" ]
  ] ]
];

var NAVTREEINDEX =
[
"_2_users_2tboudreaux_2_programming_24_d_s_t_a_r_2_grid_fire_2src_2include_2gridfire_2engine_2engine_approx8_8h-example.html",
"classgridfire_1_1_dynamic_engine.html#a8f66fb76c9f983ee93bf8b95cdc85596",
"classgridfire_1_1_multiscale_partitioning_engine_view.html#a40ce2c18063cd56cffd58419b9ae96bf",
"classgridfire_1_1exceptions_1_1_stale_engine_error.html",
"classgridfire_1_1rates_1_1weak_1_1_weak_reaction.html#aefc75b3149bee2da878177279ad375da",
"classgridfire_1_1reaction_1_1_reaction_set.html#abae9fbc4b4306479a6cd3173ac628987",
"classgridfire_1_1trigger_1_1_or_trigger.html#a6da04b5b15371f7087be03bb337dbdd6",
"engine_2bindings_8h.html#a73a2ec4150b3954550f16b97102ec993",
"namespacegridfire_1_1approx8.html#aa04f907d4ef6a1b6b2a9a28d4bb53882",
"structgridfire_1_1_multiscale_partitioning_engine_view_1_1_cache_stats.html#afc5299ebf09f9b208f65619012902b77",
"structgridfire_1_1partition_1_1record_1_1_rauscher_thielemann_partition_data_record.html#acf5c38826050594aee98d746d0ba40b2"
];

var SYNCONMSG = 'click to disable panel synchronization';
var SYNCOFFMSG = 'click to enable panel synchronization';