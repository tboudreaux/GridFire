import re
import sys
from collections import defaultdict
from re import Match
from typing import List, Tuple, Any, LiteralString
import numpy as np
from serif.atomic import species
from serif.atomic import Species
from serif.constants import Constants
import hashlib
from collections import Counter
import math

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
    rpName: str           # added
    reactionType: str     # added (e.g. "(p,γ)")
    reverse: bool
    def format_rp_name(self) -> str:
        return self.rpName


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
    name_map = {'p': 'H_1', 'd': 'H_2', 't': 'H_3', 'n': 'n_1', 'a': 'He_4'}
    if name.lower() in name_map:
        return name_map[name.lower()]
    match = re.match(r"([a-zA-Z]+)(\d+)", name)
    if match:
        element, mass = match.groups()
        return f"{element.capitalize()}_{mass}"
    return f"{name.capitalize()}_1"

def parse_reaclib_entry(entry: str) -> tuple[Match[str] | None, bool]:
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

def translate_names_to_species(names: List[str]) -> List[Species]:
    sp = list()
    split_alpha_digits = lambda inputString: re.match(r'([A-Za-z]+)[-+*]?(\d+)$', inputString).groups()
    for name in names:
        if name in ('t', 'a', 'd', 'n', 'p'):
            name = {'t': 'H-3', 'a': 'He-4', 'd': 'H-2', 'n': 'n-1', 'p': 'H-1'}[name]
        else:
            name = '-'.join(split_alpha_digits(name)).capitalize()
        try:
            sp.append(species[name])
        except Exception as e:
            print("Error: Species not found in database:", name, e)
            raise ReaclibParseError(f"Species '{name}' not found in species database.", line_content=name)
    return sp



def determine_reaction_type(reactants: List[str],
                            products: List[str],
                            qValue: float
                            ) -> Tuple[str, List[str], List[str], str]:
    """
    Return (targetToken, projectiles, ejectiles, residualToken)

    • targetToken  – the nucleus that appears before the parenthesis (A)
    • projectiles  – every explicit projectile that must be written inside ‘( … )’
    • ejectiles    – every explicit ejectile that must be written after the comma
    • residualToken– the main heavy product that appears after the parenthesis (D)

    Photons and neutrinos are added / omitted exactly the way JINA REACLIB expects:
      – γ is explicit only when it is a **projectile** (photodisintegration)
      – ν/ν̄ are never explicit
    """

    if abs(qValue - 4.621) < 1e-6:
        print("Looking at he3(he3, 2p)he4")
    # --- helper look-ups ----------------------------------------------------
    reactantSpecies = translate_names_to_species(reactants)
    productSpecies  = translate_names_to_species(products)

    # Heaviest reactant → target (A); heaviest product → residual (D)
    targetSpecies   = max(reactantSpecies,  key=lambda s: s.mass())
    residualSpecies = max(productSpecies,   key=lambda s: s.mass())

    # Any other nuclear reactant is the normal projectile candidate
    nuclearProjectiles = [x for x in reactantSpecies]
    nuclearProjectiles.remove(targetSpecies)

    nuclearEjectiles = [x for x in productSpecies]
    nuclearEjectiles.remove(residualSpecies)

    # --- bulk bookkeeping (nuclei only) -------------------------------------
    aReact = sum(s.a() for s in reactantSpecies)
    zReact = sum(s.z() for s in reactantSpecies)
    nReact = len(reactantSpecies)

    aProd  = sum(s.a() for s in productSpecies)
    zProd  = sum(s.z() for s in productSpecies)
    nProd  = len(productSpecies)

    dA = aReact - aProd          # must be 0 – abort if not
    dZ = zReact - zProd          # ≠0 ⇒ leptons needed
    dN = nReact - nProd          # ±1 ⇒ photon candidate

    assert dA == 0, (
        f"Baryon number mismatch: A₍react₎={aReact}, A₍prod₎={aProd}"
    )

    projectiles: List[str] = []
    ejectiles:   List[str] = []

    # -----------------------------------------------------------------------
    # 1.  Charged-lepton bookkeeping (|ΔZ| = 1) ------------------------------
    # -----------------------------------------------------------------------
    if abs(dZ) == 1:
        # Proton → neutron  (β⁻ / e- capture)
        if dZ == -1:
            # Electron capture when (i) exo-thermic and (ii) nucleus count unchanged
            if qValue > 0 and dN == 0:
                projectiles.append("e-")      # write e- as projectile
            else:
                ejectiles.append("e-")        # β⁻ decay: e- is an ejectile

        # Neutron → proton  (β⁺ / positron capture – capture is vanishingly rare)
        elif dZ == 1:
            ejectiles.append("e+")            # β⁺ / weak-proton capture

        # Neutrino companion is implicit – never written
        # (dL is automatically fixed by hiding ν or ν̄)

    # -----------------------------------------------------------------------
    # 2.  Photon bookkeeping (ΔZ = 0) ----------------------------------------
    # -----------------------------------------------------------------------
    if dZ == 0:
        # Two → one nucleus  and exothermic  ⇒ radiative capture (γ ejectile, implicit)
        if dN == 1 and qValue > 0:
            ejectiles.append("g")
            pass  # γ is implicit; nothing to write

        # One → two nuclei  and endothermic ⇒ photodisintegration (γ projectile, explicit)
        elif dN == -1 and qValue < 0:
            projectiles.append("g")

    # -----------------------------------------------------------------------
    # 3.  Add the ordinary nuclear projectile (if any) -----------------------
    # -----------------------------------------------------------------------
    if nuclearProjectiles:
        for nucP in nuclearProjectiles:
            name = nucP.name().replace("-", "").lower()
            if name in ('h1', 'h2', 'h3', 'he4', 'n1', 'p'):
                name = name.replace('h1', 'p').replace('h2', 'd').replace('h3', 't').replace('he4', 'a').replace('n1', 'n')
            projectiles.append(name)  # REACLIB allows exactly one

    if nuclearEjectiles:
        for nucE in nuclearEjectiles:
            name = nucE.name().replace("-", "").lower()
            if name in ('h1', 'h2', 'h3', 'he4', 'n1', 'p'):
                name = name.replace('h1', 'p').replace('h2', 'd').replace('h3', 't').replace('he4', 'a').replace('n1', 'n')
            ejectiles.append(name)

    # -----------------------------------------------------------------------
    # 4.  Build return values -------------------------------------------------
    # -----------------------------------------------------------------------
    targetToken   = targetSpecies.name().replace("-", "").lower()
    residualToken = residualSpecies.name().replace("-", "").lower()

    if targetToken in ('h1', 'h2', 'h3', 'n1', 'p'):
        targetToken = targetToken.replace('h1', 'p').replace('h2', 'd').replace('h3', 't').replace('n1', 'n')
    if residualToken in ('h1', 'h2', 'h3', 'n1', 'p'):
        residualToken = residualToken.replace('h1', 'p').replace('h2', 'd').replace('h3', 't').replace('n1', 'n')

    uniqueProjectiles = set(projectiles)
    uniqueEjectiles = set(ejectiles)
    numPerUniqueProjectiles = {x: projectiles.count(x) for x in uniqueProjectiles}
    numPerUniqueEjectiles = {x: ejectiles.count(x) for x in uniqueEjectiles}
    formatedProjectileNames = [f"{numPerUniqueProjectiles[x]}{x}" if numPerUniqueProjectiles[x] > 1 else x for x in uniqueProjectiles]
    formatedEjectileNames = [f"{numPerUniqueEjectiles[x]}{x}" if numPerUniqueEjectiles[x] > 1 else x for x in uniqueEjectiles]
    rType = f"({" ".join(formatedProjectileNames)},{' '.join(formatedEjectileNames)})"
    reactionKey = f"{targetToken}{rType}{residualToken}"


    return targetToken, projectiles, ejectiles, residualToken, reactionKey, rType


def extract_groups(match: re.Match, reverse: bool) -> Reaction:
    groups = match.groups()
    chapter = int(groups[0].strip())
    rawGroup = groups[1].strip()
    rList, pList = get_rp(rawGroup, chapter)
    if 'c12' in rList and 'mg24' in pList:
        print("Found it!")
    if reverse:
        rList, pList = pList, rList
    qValue = float(groups[3].strip())
    target, proj, ejec, residual, key, rType = determine_reaction_type(rList, pList, qValue)
    reaction = Reaction(
        reactants=rList,
        products=pList,
        label=groups[2].strip(),
        chapter=chapter,
        qValue=float(groups[3].strip()),
        coeffs=[float(groups[i].strip()) for i in range(4, 11)],
        projectile=proj,
        ejectile=ejec,
        rpName=key,
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
    emplacment: str = f"s_all_reaclib_reactions.emplace(\"{label}\", REACLIBReaction(\"{label}\", \"{reaction.format_rp_name()}\", {chapter_str}, {{{reactants_str}}}, {{{products_str}}}, {q_value_str}, \"{reaction.label}\", {{{rate_sets_str}}}, {"true" if reaction.reverse else "false"}));"

    return emplacment




def generate_reaclib_header(reaclib_filepath: str, culling: float, T9: float, verbose: bool) -> tuple[
    LiteralString, int | Any, int | Any]:
    """
    Parses a JINA REACLIB file using regular expressions and generates a C++ header file string.

    Args:
        reaclib_filepath: The path to the REACLIB data file.
        culling: The threshold for culling reactions based on their rates at T9.
        T9: The temperature in billions of Kelvin to evaluate the reaction rates for culling.
        verbose: If True, prints additional information about skipped reactions.

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
            try:
                reac = extract_groups(m, r)
            except ReaclibParseError as e:
                continue
            if verbose:
                print(f"Parsed reaction: {reac.format_rp_name()} ({reac.coeffs}) with label {reac.label} (reverse: {reac.reverse})")
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
        "#include \"fourdst/composition/atomicSpecies.h\"",
        "#include \"fourdst/composition/species.h\"",
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
    energy = list()
    energyFile = open('energy.txt', 'w')
    energyFile.write("name;maxEnergy;QValue,reactants;products;a0;a1;a2;a3;a4;a5;a6\n")
    for reaction in reactions:
        maxEnergy = calculate_peak_importance(reaction)
        energyFile.write(f"{reaction.format_rp_name()};{maxEnergy};{reaction.qValue};{' '.join(reaction.reactants)};{' '.join(reaction.products)};{';'.join([str(x) for x in reaction.coeffs])}\n")
        energy.append(maxEnergy)
        reactantNames = [f'{format_cpp_identifier(r)}' for r in reaction.reactants]
        productNames = [f'{format_cpp_identifier(p)}' for p in reaction.products]
        reactionName = f"{'_'.join(reactantNames)}_to_{'_'.join(productNames)}_{reaction.label.upper()}"
        if culling > 0.0:
            rate = evaluate_rate(reaction.coeffs, T9)
            if rate < culling:
                if verbose:
                    print(f"Skipping reaction {reaction.format_rp_name()} ({reactionName}) with rate {rate:.6e} at T9={T9} (culling threshold: {culling} at T9={T9})")
                totalSkipped += 1
                continue
            else:
                totalIncluded += 1
        else:
            totalIncluded += 1

        defines = ' && '.join(set([f"defined(SERIF_SPECIES_{name.upper().replace('-', '_min_').replace('+', '_add_').replace('*', '_mult_')})" for name in reactantNames + productNames]))
        cpp_lines.append(f"         #if {defines}")
        emplacment = format_emplacment(reaction)
        cpp_lines.append(f"             {emplacment}")
        cpp_lines.append(f"         #endif // {defines}")
    cpp_lines.append("\n    }\n} // namespace gridfire::reaclib\n")
    energyFile.close()
    #save energy data to a file
    return "\n".join(cpp_lines), totalSkipped, totalIncluded

def calculate_peak_importance(reaction: Reaction) -> float:
    TGrid = np.logspace(-3, 2, 100)  # Temperature grid from 0.001 to 100 T9
    RhoGrid = np.logspace(0.0, 6.0, 100)  # Density grid from 1e0 to 1e3 g/cm^3
    N_A: float = Constants['N_a'].value
    u: float = Constants['u'].value
    max_energy_proxy: float = 0.0

    if not reaction.reactants:
        return 0.0

    numReactants: int = len(reaction.reactants)
    maxRate: float = 0.0

    reactantCount: Counter = Counter(reaction.reactants)
    factorial_correction: float = 1.0
    for count in reactantCount.values():
        if count > 1:
            factorial_correction *= math.factorial(count)

    molar_correction_factor = 1.0
    if numReactants > 1:
        molar_correction_factor = N_A ** (numReactants - 1)

    Y_ideal = 1.0 / numReactants

    mass_term = 1.0
    split_alpha_digits = lambda inputString: re.match(r'([A-Za-z]+)(\d+)$', inputString).groups()
    for reactant in reaction.reactants:
        try:
            if reactant in ('t', 'a', 'he4', 'd', 'n', 'p'):
                reactant = {'t': 'H-3', 'a': 'He-4', 'he4': 'He-4', 'd': 'H-2', 'n': 'n-1', 'p': 'H-1'}[reactant]
            else:
                reactant = '-'.join(split_alpha_digits(reactant)).capitalize()
                # print(f"Parsing reactant {reactant} using split_alpha_digits")
            reactantMassAMU = species[reactant].mass()
            reactantMassG = reactantMassAMU * u
            mass_term *= (Y_ideal/ reactantMassG)
        except Exception as e:
            # print(f"Error: Reactant {reactant} not found in species database. (what: {e})")
            return 0.0

    for T9 in TGrid:
        k_reaclib = evaluate_rate(reaction.coeffs, T9)
        for rho in RhoGrid:
            n_product_no_rho = mass_term / factorial_correction
            full_rate = (n_product_no_rho *( rho ** numReactants) * k_reaclib) / molar_correction_factor
            energy_proxy = full_rate * abs(reaction.qValue)
            if energy_proxy > max_energy_proxy:
                max_energy_proxy = energy_proxy
    print(f"For reaction {reaction.format_rp_name()}, max energy proxy: {max_energy_proxy:.6e} MeV")
    return max_energy_proxy

# def smart_cull(reactions: List[Reaction], verbose: bool = False):

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
