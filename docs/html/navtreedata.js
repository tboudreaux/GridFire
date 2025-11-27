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
    [ "Introduction", "index.html#autotoc_md2", [
      [ "Design Philosophy and Workflow", "index.html#autotoc_md3", null ],
      [ "Funding", "index.html#autotoc_md4", null ]
    ] ],
    [ "Usage", "index.html#autotoc_md5", [
      [ "Python installation", "index.html#autotoc_md6", [
        [ "pypi", "index.html#autotoc_md7", null ],
        [ "source", "index.html#autotoc_md8", null ],
        [ "source for developers", "index.html#autotoc_md9", [
          [ "Patching Shared Object Files", "index.html#autotoc_md10", null ]
        ] ]
      ] ],
      [ "Automatic Build and Installation", "index.html#autotoc_md11", [
        [ "Script Build and Installation Instructions", "index.html#autotoc_md12", null ],
        [ "Currently, known good platforms", "index.html#autotoc_md13", null ]
      ] ],
      [ "Manual Build Instructions", "index.html#autotoc_md14", [
        [ "Prerequisites", "index.html#autotoc_md15", [
          [ "Required", "index.html#autotoc_md16", null ],
          [ "Optional", "index.html#autotoc_md17", null ]
        ] ],
        [ "Install Scripts", "index.html#autotoc_md18", [
          [ "Ease of Installation", "index.html#autotoc_md19", null ],
          [ "Reproducibility", "index.html#autotoc_md20", null ],
          [ "Examples", "index.html#autotoc_md21", [
            [ "TUI config and saving", "index.html#autotoc_md22", null ],
            [ "TUI config loading and meson setup", "index.html#autotoc_md23", null ],
            [ "CLI config loading, setup, and build", "index.html#autotoc_md24", null ]
          ] ]
        ] ],
        [ "Dependency Installation on Common Platforms", "index.html#autotoc_md25", null ],
        [ "Building the C++ Library", "index.html#autotoc_md26", [
          [ "Clang vs. GCC", "index.html#autotoc_md27", null ]
        ] ],
        [ "Installing the Library", "index.html#autotoc_md28", null ],
        [ "Minimum compiler versions", "index.html#autotoc_md29", null ]
      ] ],
      [ "Code Architecture and Logical Flow", "index.html#autotoc_md30", null ],
      [ "Engines", "index.html#autotoc_md31", [
        [ "GraphEngine", "index.html#autotoc_md32", null ],
        [ "GraphEngine Configuration Options", "index.html#autotoc_md33", null ],
        [ "Available Partition Functions", "index.html#autotoc_md34", null ],
        [ "AutoDiff", "index.html#autotoc_md35", null ]
      ] ],
      [ "Reaclib in GridFire", "index.html#autotoc_md36", null ],
      [ "Engine Views", "index.html#autotoc_md37", [
        [ "A Note about composability", "index.html#autotoc_md38", null ]
      ] ],
      [ "Numerical Solver Strategies", "index.html#autotoc_md39", [
        [ "NetworkSolverStrategy<EngineT>:", "index.html#autotoc_md40", null ],
        [ "NetIn and NetOut", "index.html#autotoc_md41", [
          [ "NetIn", "index.html#autotoc_md42", null ],
          [ "NetOut", "index.html#autotoc_md43", null ]
        ] ],
        [ "CVODESolverStrategy", "index.html#autotoc_md44", null ]
      ] ],
      [ "Python Extensibility", "index.html#autotoc_md45", null ]
    ] ],
    [ "Usage Examples", "index.html#autotoc_md46", [
      [ "C++", "index.html#autotoc_md47", [
        [ "GraphEngine Initialization", "index.html#autotoc_md48", null ],
        [ "Adaptive Network View", "index.html#autotoc_md49", null ],
        [ "Composition Initialization", "index.html#autotoc_md50", null ],
        [ "Common Workflow Example", "index.html#autotoc_md51", null ],
        [ "Callback and Policy Example", "index.html#autotoc_md52", [
          [ "Callback Context", "index.html#autotoc_md53", null ]
        ] ]
      ] ],
      [ "Python", "index.html#autotoc_md54", [
        [ "Python Example for End Users", "index.html#autotoc_md55", null ]
      ] ]
    ] ],
    [ "External Usage", "index.html#autotoc_md56", [
      [ "C API Overview", "index.html#autotoc_md57", [
        [ "C Example", "index.html#autotoc_md58", null ]
      ] ],
      [ "Fortran API Overview", "index.html#autotoc_md59", [
        [ "Fortran Example", "index.html#autotoc_md60", null ]
      ] ]
    ] ],
    [ "Related Projects", "index.html#autotoc_md61", null ],
    [ "GridFire Python Usage Guide", "md_docs_2static_2usage.html", [
      [ "Installation", "md_docs_2static_2usage.html#autotoc_md64", [
        [ "1.1 PyPI Release", "md_docs_2static_2usage.html#autotoc_md65", null ],
        [ "1.2 Development from Source", "md_docs_2static_2usage.html#autotoc_md66", null ]
      ] ],
      [ "Why These Engines and Views?", "md_docs_2static_2usage.html#autotoc_md68", null ],
      [ "Step-by-Step Example", "md_docs_2static_2usage.html#autotoc_md70", null ],
      [ "Visualizing Reaction Networks", "md_docs_2static_2usage.html#autotoc_md72", null ],
      [ "Beyond the Basics", "md_docs_2static_2usage.html#autotoc_md74", null ]
    ] ],
    [ "Deprecated List", "deprecated.html", null ],
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
        [ "Enumerations", "globals_enum.html", null ],
        [ "Enumerator", "globals_eval.html", null ],
        [ "Macros", "globals_defs.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"____init_____8py.html",
"classgridfire_1_1engine_1_1_dynamic_engine.html#a24103027a38e02d4342161164554d332",
"classgridfire_1_1engine_1_1_graph_engine.html#af7dab9063bb92aac64b6dfa755074739",
"classgridfire_1_1engine_1_1_reaction.html#a0c533797f18c2fe4945f1afb3fb06d5d",
"classgridfire_1_1partition_1_1_composite_partition_function.html#a85aaac230e9de2fd50d4d453f6d5def8",
"classgridfire_1_1policy_1_1_proton_proton_chain_policy.html",
"classgridfire_1_1reaction_1_1_logical_reaclib_reaction.html#a69508797ca0793f6d3f58de4e0b97667",
"classgridfire_1_1reaction_1_1_weak_reaclib_reaction.html#ac3fc463a1cea6854ea989732e1e5f65a",
"classgridfire_1_1trigger_1_1_not_trigger.html#ae453a72df20fdd416443ef7867f86b92",
"classgridfire_1_1utils_1_1_column_base.html#ab6e4cb14c8a0ad7b6b9ad0fbde766cba",
"functions_vars_i.html",
"namespacegridfire_1_1engine.html#ace3a790fb57920fd12328fa9c5a3a1de",
"structgridfire_1_1_net_in.html#a0a8d820cfeaa92ee31f253795c57e0d1",
"structgridfire_1_1rates_1_1weak_1_1_rate_data_row.html#ac08d42071950006ccf3c438253ef9ee9",
"weak__types_8h.html"
];

var SYNCONMSG = 'click to disable panel synchronization';
var SYNCOFFMSG = 'click to enable panel synchronization';