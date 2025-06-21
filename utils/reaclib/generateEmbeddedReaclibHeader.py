import re
import sys
from collections import defaultdict
from typing import List, Tuple
import numpy as np
import hashlib

#import dataclasses
from dataclasses import dataclass

@dataclass
class Reaction:
    reactants: List[str]
    products: List[str]
    label: str
    chapter: int
    qValue: float
    coeffs: List[float]
    projectile: str        # added
    ejectile: str         # added
    reactionType: str     # added (e.g. "(p,γ)")
    reverse: bool
    def format_rp_name(self) -> str:
        # radiative or particle captures: 2 reactants -> 1 product
        if len(self.reactants) == 2 and len(self.products) == 1:
            target = self.reactants[0]
            prod = self.products[0]
            return f"{target}({self.projectile},{self.ejectile}){prod}"
        # fallback: join lists
        react_str = '+'.join(self.reactants)
        prod_str = '+'.join(self.products)
        return f"{react_str}->{prod_str}"


    def __repr__(self):
        return f"Reaction({self.format_rp_name()})"


def evaluate_rate(coeffs: List[float], T9: float) -> float:
    rateExponent: float = coeffs[0] + \
                        coeffs[1] / T9 + \
                        coeffs[2] / (T9 ** (1/3)) + \
                        coeffs[3] * (T9 ** (1/3)) + \
                        coeffs[4] * T9 + \
                        coeffs[5] * (T9 ** (5/3)) + \
                        coeffs[6] * (np.log(T9))
    return np.exp(rateExponent)

class ReaclibParseError(Exception):
    """Custom exception for parsing errors."""
    def __init__(self, message, line_num=None, line_content=None):
        self.line_num = line_num
        self.line_content = line_content
        full_message = f"Error"
        if line_num is not None:
            full_message += f" on line {line_num}"
        full_message += f": {message}"
        if line_content is not None:
            full_message += f"\n  -> Line content: '{line_content}'"
        super().__init__(full_message)


def format_cpp_identifier(name: str) -> str:
    name_map = {'p': 'H_1', 'n': 'n_1', 'd': 'd', 't': 't', 'a': 'a'}
    if name.lower() in name_map:
        return name_map[name.lower()]
    match = re.match(r"([a-zA-Z]+)(\d+)", name)
    if match:
        element, mass = match.groups()
        return f"{element.capitalize()}_{mass}"
    return f"{name.capitalize()}_1"

def parse_reaclib_entry(entry: str) -> Tuple[List[str], str, float, List[float], bool]:
    pattern = re.compile(r"""^([1-9]|1[0-1])\r?\n
[ \t]*
((?:[A-Za-z0-9-*]+[ \t]+)*
[A-Za-z0-9-*]+)
[ \t]+
([A-Za-z0-9+]+)
[ \t]+
([+-]?(?:\d+\.\d*|\.\d+)(?:[eE][+-]?\d+)) 
[ \t\r\n]+
[ \t\r\n]*([+-]?\d+\.\d*e[+-]?\d+)\s*
([+-]?\d+\.\d*e[+-]?\d+)\s*
([+-]?\d+\.\d*e[+-]?\d+)\s*
([+-]?\d+\.\d*e[+-]?\d+)\s*
([+-]?\d+\.\d*e[+-]?\d+)\s*
([+-]?\d+\.\d*e[+-]?\d+)\s*
([+-]?\d+\.\d*e[+-]?\d+)
""", re.MULTILINE | re.VERBOSE)
    match = pattern.match(entry)
    reverse = True if entry.split('\n')[1][48] == 'v' else False
    return match, reverse


def get_rp(group: str, chapter: int) -> Tuple[List[str], List[str]]:
    rpGroupings = {
        1: (1, 1), 2: (1, 2), 3: (1, 3), 4: (2, 1), 5: (2, 2),
        6: (2, 3), 7: (2, 4), 8: (3, 1), 9: (3, 2), 10: (4, 2), 11: (1, 4)
    }
    species = group.split()
    nReact, nProd = rpGroupings[chapter]
    reactants = species[:nReact]
    products = species[nReact:nReact + nProd]
    return reactants, products


def determine_reaction_type(reactants: List[str], products: List[str]) -> Tuple[str, str, str]:
    # assumes no reverse flag applied
    projectile = ''
    ejectile = ''
    # projectile is the lighter reactant (p, n, he4)
    for sp in reactants:
        if sp in ('p', 'n', 'he4'):
            projectile = sp
            break
    # ejectile logic
    if len(products) == 1:
        ejectile = 'g'
    elif 'he4' in products:
        ejectile = 'a'
    elif 'p' in products:
        ejectile = 'p'
    elif 'n' in products:
        ejectile = 'n'
    reactionType = f"({projectile},{ejectile})"
    return projectile, ejectile, reactionType
def extract_groups(match: re.Match, reverse: bool) -> Reaction:
    groups = match.groups()
    chapter = int(groups[0].strip())
    rawGroup = groups[1].strip()
    rList, pList = get_rp(rawGroup, chapter)
    if reverse:
        rList, pList = pList, rList
    proj, ejec, rType = determine_reaction_type(rList, pList)
    reaction = Reaction(
        reactants=rList,
        products=pList,
        label=groups[2].strip(),
        chapter=chapter,
        qValue=float(groups[3].strip()),
        coeffs=[float(groups[i].strip()) for i in range(4, 11)],
        projectile=proj,
        ejectile=ejec,
        reactionType=rType,
        reverse=reverse
    )
    return reaction
def format_emplacment(reaction: Reaction) -> str:
    reactantNames = [f'{format_cpp_identifier(r)}' for r in reaction.reactants]
    productNames = [f'{format_cpp_identifier(p)}' for p in reaction.products]

    reactants_cpp = [f'fourdst::atomic::{format_cpp_identifier(r)}' for r in reaction.reactants]
    products_cpp = [f'fourdst::atomic::{format_cpp_identifier(p)}' for p in reaction.products]

    label = f"{'_'.join(reactantNames)}_to_{'_'.join(productNames)}_{reaction.label.upper()}"


    reactants_str = ', '.join(reactants_cpp)
    products_str = ', '.join(products_cpp)

    q_value_str = f"{reaction.qValue:.6e}"
    chapter_str = reaction.chapter

    rate_sets_str = ', '.join([str(x) for x in reaction.coeffs])
    emplacment: str = f"s_all_reaclib_reactions.emplace(\"{label}\", REACLIBReaction(\"{label}\", {chapter_str}, {{{reactants_str}}}, {{{products_str}}}, {q_value_str}, \"{reaction.label}\", {{{rate_sets_str}}}, {"true" if reaction.reverse else "false"}));"

    return emplacment




def generate_reaclib_header(reaclib_filepath: str, culling: float, T9: float, verbose: bool) -> str:
    """
    Parses a JINA REACLIB file using regular expressions and generates a C++ header file string.

    Args:
        reaclib_filepath: The path to the REACLIB data file.

    Returns:
        A string containing the complete C++ header content.
    """
    with open(reaclib_filepath, 'r') as file:
        content = file.read()
        fileHash = hashlib.sha256(content.encode('utf-8')).hexdigest()
    # split the file into blocks of 4 lines each
    lines = content.split('\n')
    entries = ['\n'.join(lines[i:i+4]) for i in range(0, len(lines), 4) if len(lines[i:i+4]) == 4]
    reactions = list()
    for entry in entries:
        m, r = parse_reaclib_entry(entry)
        if m is not None:
            reac = extract_groups(m, r)
            reactions.append(reac)

    # --- Generate the C++ Header String ---
    cpp_lines = [
        "// This file is automatically generated. Do not edit!",
        "// Generated on: " + str(np.datetime64('now')),
        "// REACLIB file hash (sha256): " + fileHash,
        "// Generated from REACLIB data file: " + reaclib_filepath,
        "// Culling threshold: rate >" + str(culling) + " at T9 = " + str(T9),
        "// Note that if the culling threshold is set to 0.0, no reactions are culled.",
        "// Includes %%TOTAL%% reactions.",
        "// Note: Only reactions with species defined in the atomicSpecies.h header will be included at compile time.",
        "#pragma once",
        "#include \"atomicSpecies.h\"",
        "#include \"species.h\"",
        "#include \"reaclib.h\"",
        "\nnamespace gridfire::reaclib {\n",
        """    
    inline void initializeAllReaclibReactions() {
        if (s_initialized) return; // already initialized
        s_initialized = true;
        s_all_reaclib_reactions.clear();
        s_all_reaclib_reactions.reserve(%%TOTAL%%); // reserve space for total reactions
    """
    ]
    totalSkipped = 0
    totalIncluded = 0
    for reaction in reactions:
        reactantNames = [f'{format_cpp_identifier(r)}' for r in reaction.reactants]
        productNames = [f'{format_cpp_identifier(p)}' for p in reaction.products]
        reactionName = f"{'_'.join(reactantNames)}_to_{'_'.join(productNames)}_{reaction.label.upper()}"
        if culling > 0.0:
            rate = evaluate_rate(reaction.coeffs, T9)
            if rate < culling:
                if verbose:
                    print(f"Skipping reaction {reactionName} with rate {rate:.6e} at T9={T9} (culling threshold: {culling} at T9={T9})")
                totalSkipped += 1
                continue
            else:
                totalIncluded += 1

        defines = ' && '.join(set([f"defined(SERIF_SPECIES_{name.upper().replace('-', '_min_').replace('+', '_add_').replace('*', '_mult_')})" for name in reactantNames + productNames]))
        cpp_lines.append(f"         #if {defines}")
        emplacment = format_emplacment(reaction)
        cpp_lines.append(f"             {emplacment}")
        cpp_lines.append(f"         #endif // {defines}")
    cpp_lines.append("\n    }\n} // namespace gridfire::reaclib\n")
    return "\n".join(cpp_lines), totalSkipped, totalIncluded

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Generate a C++ header from a REACLIB file.")
    parser.add_argument("reaclib_file", type=str, help="Path to the REACLIB data file.")
    parser.add_argument("-o", "--output", type=str, default=None, help="Output file path (default: stdout).")
    parser.add_argument('-c', "--culling", type=float, default=0.0, help="Culling threshold for reaction rates at T9 (when 0.0, no culling is applied).")
    parser.add_argument('-T', '--T9', type=float, default=0.01, help="Temperature in billions of Kelvin (default: 0.01) to evaluate the reaction rates for culling.")
    parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose output.")
    args = parser.parse_args()

    try:
        cpp_header_string, skipped, included = generate_reaclib_header(args.reaclib_file, args.culling, args.T9, args.verbose)
        cpp_header_string = cpp_header_string.replace("%%TOTAL%%", str(included))
        print("--- Generated C++ Header (Success!) ---")
        if args.output:
            with open(args.output, 'w') as f:
                f.write(cpp_header_string)
            print(f"Header written to {args.output}")
            print(f"Total reactions included: {included}, Total reactions skipped: {skipped}")
        else:
            print(cpp_header_string)
    except ReaclibParseError as e:
        print(f"\n--- PARSING FAILED ---")
        print(e, file=sys.stderr)
        sys.exit(1)
