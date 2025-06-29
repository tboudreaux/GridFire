"""
REACLIB Reaction Formatting and Parsing Utilities
================================================

This module provides utilities for parsing, formatting, and analyzing nuclear reaction
data in the REACLIB format. It includes tools for extracting reaction information,
computing reaction rates, and exporting data for use in C++ or binary formats.

Algorithms and Quantum Number Accounting
----------------------------------------
The core of this module is the quantum number bookkeeping performed in
`determine_reaction_type`, which ensures conservation of baryon number (A),
charge (Z), and handles leptonic and photonic processes:

- **Baryon Number Conservation**: The difference in total mass number (A) between
  reactants and products must be zero. If not, an assertion error is raised.
- **Charge Conservation**: The difference in total charge (Z) is checked:
    - |ΔZ| = 1: Indicates a weak process (beta decay or electron/positron capture).
      The code distinguishes between electron/positron as projectile or ejectile
      based on the sign of ΔZ and the Q-value.
    - ΔZ = 0: If the number of nuclei changes by ±1, photons are involved
      (radiative capture or photodisintegration).
- **Projectiles and Ejectiles**: The heaviest reactant is considered the target,
  and the heaviest product is the residual. Other nuclei are projectiles/ejectiles.
  Special handling is provided for light nuclei (p, d, t, n, a).

Usage Examples
--------------
Parse a REACLIB file and export to CSV:

    python format.py path/to/reaclib/file -o output.csv -f csv

Programmatic usage:

    from utils.reaclib.format import parse_reaclib_entry, extract_groups

    entry = '''1
    h1 he4 c12 mg24    wkb    1.234e+00
    1.234e+01 2.345e+02 3.456e+03 4.567e+04 5.678e+05 6.789e+06 7.890e+07
    '''
    match, reverse = parse_reaclib_entry(entry)
    if match:
        reaction = extract_groups(match, reverse)
        print(reaction)

Functions
---------
- `parse_reaclib_entry(entry)`: Parse a REACLIB entry string.
- `extract_groups(match, reverse)`: Extracts a Reaction object from a regex match.
- `determine_reaction_type(reactants, products, qValue)`: Determines projectiles,
  ejectiles, and reaction type using quantum number accounting.
- `evaluate_rate(coeffs, T9)`: Computes the reaction rate at temperature T9.
- `format_cpp_identifier(name)`: Formats a species name for C++ code.
- `create_reaction_dataframe(reactions)`: Converts a list of Reaction objects to a DataFrame.
- `write_reactions_binary(reactions, output_path)`: Writes reactions to a binary file.

Classes
-------
- `Reaction`: Dataclass representing a nuclear reaction.
- `ReaclibParseError`: Exception for parsing errors.

See function/class docstrings for further details.

"""

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
import argparse
import pandas as pd
import struct

from dataclasses import dataclass

@dataclass
class Reaction:
    reactants: List[str]
    products: List[str]
    label: str
    chapter: int
    qValue: float
    coeffs: List[float]
    projectile: str
    ejectile: str
    rpName: str
    reactionType: str
    reverse: bool
    def format_rp_name(self) -> str:
        return self.rpName


    def __repr__(self):
        return f"Reaction({self.format_rp_name()})"


def evaluate_rate(coeffs: List[float], T9: float) -> float:
    """
    Evaluate the REACLIB reaction rate at a given temperature.

    Parameters
    ----------
    coeffs : list of float
        The 7 REACLIB coefficients (a0..a6) for the reaction.
    T9 : float
        Temperature in units of 10^9 K.

    Returns
    -------
    float
        The reaction rate at the specified temperature.

    Notes
    -----
    The rate is computed as:
        rate = exp(a0 + a1/T9 + a2/T9^{1/3} + a3*T9^{1/3} + a4*T9 + a5*T9^{5/3} + a6*ln(T9))
    """
    rateExponent: float = coeffs[0] + \
                        coeffs[1] / T9 + \
                        coeffs[2] / (T9 ** (1/3)) + \
                        coeffs[3] * (T9 ** (1/3)) + \
                        coeffs[4] * T9 + \
                        coeffs[5] * (T9 ** (5/3)) + \
                        coeffs[6] * (np.log(T9))
    return np.exp(rateExponent)

class ReaclibParseError(Exception):
    """
    Exception raised for errors encountered while parsing REACLIB entries.

    Parameters
    ----------
    message : str
        Description of the error.
    line_num : int, optional
        Line number where the error occurred.
    line_content : str, optional
        Content of the problematic line.
    """
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
    """
    Convert a REACLIB species name to a C++-friendly identifier.

    Parameters
    ----------
    name : str
        The REACLIB species name (e.g., 'h1', 'he4', 'c12', 'p', 'a').

    Returns
    -------
    str
        The formatted C++ identifier (e.g., 'H-1', 'He-4', 'C-12', etc.).
    """
    name_map = {'p': 'H-1', 'd': 'H-2', 't': 'H-3', 'n': 'n-1', 'a': 'He-4'}
    if name.lower() in name_map:
        return name_map[name.lower()]
    match = re.match(r"([a-zA-Z]+)(\d+)", name)
    if match:
        element, mass = match.groups()
        return f"{element.capitalize()}-{mass}"
    return f"{name.capitalize()}-1"

def parse_reaclib_entry(entry: str) -> tuple[Match[str] | None, bool]:
    """
    Parse a single REACLIB entry string using a regular expression.

    Parameters
    ----------
    entry : str
        The REACLIB entry as a string (typically 4 lines).

    Returns
    -------
    match : re.Match or None
        The regex match object if parsing was successful, else None.
    reverse : bool
        True if the entry is marked as a reverse reaction, else False.

    Notes
    -----
    The function uses a regular expression to extract chapter, reactants/products,
    label, Q-value, and coefficients. The 'reverse' flag is determined by the
    character at a fixed position in the entry.
    """
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
    """
    Split a REACLIB group string into reactants and products based on chapter.

    Parameters
    ----------
    group : str
        The whitespace-separated list of species.
    chapter : int
        The REACLIB chapter number (determines reactant/product count).

    Returns
    -------
    reactants : list of str
        List of reactant species names.
    products : list of str
        List of product species names.
    """
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
    """
    Convert a list of REACLIB species names to Species objects.

    Parameters
    ----------
    names : list of str
        List of REACLIB species names (e.g., 'h1', 'he4', 'c12', etc.).

    Returns
    -------
    list of Species
        List of Species objects corresponding to the names.

    Raises
    ------
    ReaclibParseError
        If a species name cannot be found in the species database.
    """
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
                            ) -> Tuple[str, List[str], List[str], str, str, str]:
    """
    Analyze a reaction for quantum number conservation and classify projectiles/ejectiles.

    Parameters
    ----------
    reactants : list of str
        List of reactant species names.
    products : list of str
        List of product species names.
    qValue : float
        Q-value of the reaction (MeV).

    Returns
    -------
    targetToken : str
        Name of the heaviest reactant (target nucleus).
    projectiles : list of str
        List of projectile names (including leptons/photons if present).
    ejectiles : list of str
        List of ejectile names (including leptons/photons if present).
    residualToken : str
        Name of the heaviest product (residual nucleus).
    reactionKey : str
        Unique string key for the reaction.
    rType : str
        String representation of the reaction type.

    Notes
    -----
    This function performs quantum number bookkeeping:
    - Checks baryon and charge conservation.
    - Identifies weak (leptonic) and photonic processes.
    - Determines projectiles/ejectiles based on mass and reaction type.
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
        # Two → one nucleus  and exothermic  ⇒ radiative capture (γ ejectile, (seems to normally be implicit, but I am writing it explicitly))
        if dN == 1 and qValue > 0:
            ejectiles.append("g")

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
    """
    Extract a Reaction object from a regex match of a REACLIB entry.

    Parameters
    ----------
    match : re.Match
        The regex match object from `parse_reaclib_entry`.
    reverse : bool
        Whether the reaction is a reverse reaction.

    Returns
    -------
    Reaction
        The parsed Reaction dataclass instance.
    """
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
    """
    Format a Reaction object as a C++ emplacement statement.

    Parameters
    ----------
    reaction : Reaction
        The Reaction object to format.

    Returns
    -------
    str
        The C++ code string for emplacing the reaction.
    """
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

def calculate_peak_importance(reaction: Reaction) -> float:
    """
    Estimate the peak energy importance of a reaction over a grid of T and rho.

    Parameters
    ----------
    reaction : Reaction
        The Reaction object to analyze.

    Returns
    -------
    float
        The maximum energy proxy (rate * |Q|) found over the grid.

    Notes
    -----
    The function evaluates the reaction rate over a grid of temperature (T9)
    and density (rho), multiplies by |Q|, and returns the maximum value.
    """
    TGrid = np.logspace(-3, 2, 100)
    RhoGrid = np.logspace(0.0, 6.0, 100)
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
            reactantMassAMU = species[reactant].mass()
            reactantMassG = reactantMassAMU * u
            mass_term *= (Y_ideal/ reactantMassG)
        except Exception as e:
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

def create_reaction_dataframe(reactions: List[Reaction]) -> pd.DataFrame:
    """
    Convert a list of Reaction objects into a pandas DataFrame.

    Parameters
    ----------
    reactions : list of Reaction
        List of Reaction objects.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns for all reaction properties.
    """
    reaction_data = []
    for reaction in reactions:
        record = {
            'id': f"{'_'.join(reaction.reactants)}_to_{'_'.join(reaction.products)}_{reaction.label.upper()}",
            'rpName': reaction.rpName,
            'chapter': reaction.chapter,
            'reactants': ' '.join(reaction.reactants),
            'products': ' '.join(reaction.products),
            'qValue': reaction.qValue,
            'is_reverse': reaction.reverse,
            'label': reaction.label,
            'a0': reaction.coeffs[0],
            'a1': reaction.coeffs[1],
            'a2': reaction.coeffs[2],
            'a3': reaction.coeffs[3],
            'a4': reaction.coeffs[4],
            'a5': reaction.coeffs[5],
            'a6': reaction.coeffs[6]
        }
        reaction_data.append(record)
        
    return pd.DataFrame(reaction_data)

def write_reactions_binary(reactions: List[Reaction], output_path: str):
    """
    Write a list of Reaction objects to a binary file.

    Parameters
    ----------
    reactions : list of Reaction
        List of Reaction objects to write.
    output_path : str
        Path to the output binary file.

    Notes
    -----
    Each reaction is packed using struct with a fixed format for chapter, Q-value,
    coefficients, reverse flag, label, rpName, reactants, and products.
    """
    record_format = '<i d 7d ? 8s 64s 128s 128s'
    
    with open(output_path, 'wb') as f:
        for reaction in reactions:
            label_bytes = reaction.label.encode('utf-8')[:7].ljust(8, b'\0')
            rpName_bytes = reaction.rpName.encode('utf-8')[:63].ljust(64, b'\0')
            reactants_str = ' '.join([format_cpp_identifier(x) for x in reaction.reactants])
            products_str = ' '.join([format_cpp_identifier(x) for x in reaction.products])
            reactants_bytes = reactants_str.encode('utf-8')[:127].ljust(128, b'\0')
            products_bytes = products_str.encode('utf-8')[:127].ljust(128, b'\0')

            packed_data = struct.pack(
                record_format,
                reaction.chapter,
                reaction.qValue,
                *reaction.coeffs,
                reaction.reverse,
                label_bytes,
                rpName_bytes,
                reactants_bytes,
                products_bytes
            )
            f.write(packed_data)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate a CSV file from a REACLIB file.")
    parser.add_argument("reaclib_file", type=str, help="Path to the REACLIB data file.")
    parser.add_argument("-o", "--output", type=str, default="reactions.csv", help="Output CSV file path.")
    parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose output.")
    parser.add_argument('-f', '--format', type=str, choices=['csv', 'bin'], default='bin', help="Output format")
    args = parser.parse_args()

    try:
        with open(args.reaclib_file, 'r') as file:
            content = file.read()
        
        lines = content.split('\n')
        entries = ['\n'.join(lines[i:i+4]) for i in range(0, len(lines), 4) if len(lines[i:i+4]) == 4 and lines[i].strip()]

        parsed_reactions = []
        for i, entry in enumerate(entries):
            m, r = parse_reaclib_entry(entry)
            if m is not None:
                try:
                    reac = extract_groups(m, r)
                    parsed_reactions.append(reac)
                except ReaclibParseError as e:
                    if args.verbose:
                        print(f"Skipping entry starting at line {i*4 + 1} due to parsing error: {e}", file=sys.stderr)
                    continue
        
        print(f"Successfully parsed {len(parsed_reactions)} reactions from {args.reaclib_file}")

        reaction_df = create_reaction_dataframe(parsed_reactions)


        if args.format == 'csv':
            reaction_df.to_csv(args.output, index=False)
            print("--- CSV Generation (Success!) ---")
            print(f"Reaction data written to {args.output}")
        else:
            write_reactions_binary(parsed_reactions, args.output)

            print("--- Binary File Generation (Success!) ---")
            print(f"Reaction data written to {args.output}")

    except FileNotFoundError:
        print(f"Error: Input file not found at {args.reaclib_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\n--- AN UNEXPECTED ERROR OCCURRED ---")
        print(e, file=sys.stderr)
        sys.exit(1)
