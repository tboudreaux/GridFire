#!/usr/bin/env python3
"""Generate MESA reference data for the GridFire validation suite by running bbq on GridFire's own networks.

For every case in cases.json this
  1. exports the case's GridFire network to a MESA .net file (gf_mesa_setup),
  2. asks bbq for its isotope ordering and writes the initial composition in that order,
  3. runs a bbq hydrostatic burn with the case's step grid and BBQ_DISABLE_EOS=1, and
  4. converts bbq's output to reference/<case>/reference.csv (GridFire species names) with provenance.

The comparison test (mesa_comparison_test) reads only the files written here, so references can be generated on any
machine with MESA + bbq and committed (or uploaded as CI artifacts) for machines without them.

Requirements: MESA_DIR pointing at MESA r23.05.1 with mesa-r23.05.1-net_burn.patch (and, for unscreened runs,
mesa-r23.05.1-no-screening.patch) applied, a bbq executable built
against it with the paper's bbq patches (validation/ManuscriptFigures/BBQCompare/diffs/lib_*.patch), and the
gf_mesa_setup tool from a GridFire build directory. Standard library only.

    generate_reference.py --bbq /path/to/bbq --setup build/tests/gf_mesa_setup [--case NAME ...]
"""
import argparse
import datetime
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PATCH_DIR = HERE.parents[2] / "validation" / "ManuscriptFigures" / "BBQCompare" / "diffs"

INLIST_TEMPLATE = """&bbq
    net_name = 'network.net'
    max_steps = 1000000
    eps = 1d-8
    odescal = 1d-10
    stptry = 0
    {mode}
    write_iso_list = .true.
    iso_list_filename = 'iso.list'
/

&sampling
/

&profile
/

&hydrostatic
    min_time = {min_log_time}
    max_time = {max_log_time}
    log_time = .true.
    num_times = {num_times}
    logT = {logT}
    logRho = {logRho}
    input_composition_filename = 'comp.txt'
    output_filename = 'bbq_output.txt'
/

&eos
/

&nuclear
    screening_mode = '{screening_mode}'
/
"""

# GridFire's default weak screening is Salpeter (1954), which MESA calls 'salpeter'. It must be set in &nuclear: bbq
# reads no other group for it and silently keeps its default ('chugunov') otherwise. MESA r23.05.1's no-screening
# option ('') zeroes every screened rate unless mesa-r23.05.1-no-screening.patch is applied.
SCREENING_MODE = "salpeter"


def run(cmd, cwd, env=None):
    result = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True)
    (Path(cwd) / f"{Path(cmd[0]).name}.log").write_text(result.stdout + result.stderr)
    if result.returncode != 0:
        sys.exit(f"{' '.join(map(str, cmd))} failed in {cwd} (exit {result.returncode}); see the .log file there")
    return result.stdout


def git_describe(path):
    try:
        return subprocess.run(["git", "-C", str(path), "describe", "--always", "--dirty", "--tags"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def mesa_version(mesa_dir):
    version_file = Path(mesa_dir) / "data" / "version_number"
    return version_file.read_text().strip() if version_file.exists() else "unknown"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_pairs(path):
    return [line.split() for line in Path(path).read_text().splitlines() if line.strip()]


# MESA reports reactions it cannot evaluate (e.g. reverse handles missing from its reaclib index, weak decays missing
# from its tables) as non-fatal messages and carries on with bad rates, or later aborts on a NaN rate table.
MESA_FAILURE = re.compile(r"^\s*(r_\S+) failed in (\S+)|^\s*(r_\S+) rates\s+\d+\s+NaN", re.MULTILINE)
PARTICLE_CODES = {"p": "h1", "n": "neut", "a": "he4", "d": "h2", "t": "h3", "g": None}


def mesa_participants(name):
    """Sorted iso names of everything in a MESA rate name, e.g. r_f17_na_n14 -> f17, he4, n14, neut."""
    participants = []
    for token in name.removeprefix("r_").split("_"):
        if token == "neut" or re.search(r"\d", token):
            participants.append(token)
        elif token in ("wk", "to", "ec"):
            continue
        elif all(c in PARTICLE_CODES for c in token):
            participants += [PARTICLE_CODES[c] for c in token if PARTICLE_CODES[c]]
        else:
            return None
    return ",".join(sorted(participants))


def mesa_failures(log_text):
    failures = {}
    for m in MESA_FAILURE.finditer(log_text):
        name = m.group(1) or m.group(3)
        failures[name] = f"MESA: {name} " + (f"failed in {m.group(2)}" if m.group(1) else "rate table is NaN")
    return failures


def run_bbq(args, out_dir, env, allow_failure=False):
    """Run bbq; return the MESA rate failures it reported (bbq may exit 0 despite them)."""
    result = subprocess.run([str(args.bbq), "inlist"], cwd=out_dir, env=env, capture_output=True, text=True)
    log = result.stdout + result.stderr
    (out_dir / "bbq.log").write_text(log)
    failures = mesa_failures(log)
    if result.returncode != 0 and not failures and not allow_failure:
        sys.exit(f"bbq failed in {out_dir} (exit {result.returncode}); see bbq.log there")
    return failures


def generate_case(case, args, provenance_base):
    out_dir = HERE / "reference" / case["name"]
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)

    env = dict(os.environ, BBQ_DISABLE_EOS="1")
    params = {k: case[k] for k in ("min_log_time", "max_log_time", "num_times", "logT", "logRho")}
    excluded_file = out_dir / "excluded_reactions.txt"
    excluded = {}  # GridFire reaction id -> reason

    # 1. Export GridFire's network and let MESA load it. Reactions MESA cannot evaluate are removed from the case (in
    #    both codes: the comparison test rebuilds the network with the same exclusions) until MESA loads it cleanly.
    for attempt in range(10):
        excluded_file.write_text("".join(f"{gid}  # {why}\n" for gid, why in sorted(excluded.items())))
        run([str(args.setup), str(HERE / "cases.json"), case["name"], str(out_dir), str(excluded_file)], cwd=out_dir)
        (out_dir / "inlist").write_text(INLIST_TEMPLATE.format(mode="just_write_isos = .true.", screening_mode=SCREENING_MODE, **params))
        failures = run_bbq(args, out_dir, env, allow_failure=True)
        if not failures:
            break
        by_participants = {}
        for line in (out_dir / "reaction_map.txt").read_text().splitlines():
            gid, _mesa_name, participants = line.split("\t")
            by_participants.setdefault(participants, []).append(gid)
        for mesa_name, why in failures.items():
            key = mesa_participants(mesa_name)
            matches = by_participants.get(key, []) if key else []
            if not matches:
                sys.exit(f"{case['name']}: cannot map MESA failure '{why}' to a GridFire reaction (participants {key})")
            for gid in matches:
                excluded[gid] = why
        print(f"{case['name']}: excluding {len(excluded)} reaction(s) MESA cannot evaluate, retrying")
    else:
        sys.exit(f"{case['name']}: MESA still reports failures after excluding {len(excluded)} reactions")

    gf_to_mesa = {gf: mesa for gf, mesa in read_pairs(out_dir / "species_map.txt")}
    mesa_to_gf = {mesa: gf for gf, mesa in gf_to_mesa.items()}
    initial = {gf: float(x) for gf, x in read_pairs(out_dir / "initial_composition.txt")}
    isos = [line.split()[-1] for line in (out_dir / "iso.list").read_text().splitlines() if line.strip()]
    unknown = [iso for iso in isos if iso not in mesa_to_gf]
    if unknown:
        sys.exit(f"{case['name']}: bbq isotopes not in GridFire's network: {unknown}")

    # 2. Initial composition in bbq's order
    (out_dir / "comp.txt").write_text("".join(f"{initial[mesa_to_gf[iso]]:.17e}\n" for iso in isos))

    # 3. The burn
    (out_dir / "inlist").write_text(INLIST_TEMPLATE.format(mode="use_hydrostatic = .true.", screening_mode=SCREENING_MODE, **params))
    failures = run_bbq(args, out_dir, env)
    if failures:
        sys.exit(f"{case['name']}: MESA reported rate failures during the burn: {sorted(failures)}")

    # 4. Convert to GridFire naming
    raw = (out_dir / "bbq_output.txt").read_text().split("\n")
    header = raw[0].split()
    if header[:4] != ["age", "dt", "eps_nuc", "eps_neu"]:
        sys.exit(f"{case['name']}: unexpected bbq output header {header[:6]}")
    columns = ["age", "dt", "energy_step", "eps_neu"] + [mesa_to_gf[iso] for iso in header[4:]]

    provenance = dict(provenance_base, case=case, net_file_sha256=sha256(out_dir / "network.net"),
                      excluded_reactions=excluded,
                      columns={"age": "cumulative burn time [s]",
                               "dt": "duration of this step [s]",
                               "energy_step": "nuclear energy released during the step, net of neutrinos [erg/g] (bbq eps_nuc = avg_eps_nuc * dt)",
                               "eps_neu": "bbq eps_neu as reported by MESA net_1_zone_burn",
                               "<species>": "mass fraction at the end of the step"})
    with open(out_dir / "reference.csv", "w") as f:
        for key in ("mesa_version", "bbq_version", "bbq_disable_eos", "generated"):
            f.write(f"# {key}: {provenance[key]}\n")
        f.write(",".join(columns) + "\n")
        for line in raw[1:]:
            values = line.split()
            if values:
                f.write(",".join(values) + "\n")
    (out_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    # Keep the inputs needed to reproduce the run and the converted reference; drop logs and bbq's raw output
    for scratch in ("bbq.log", "gf_mesa_setup.log", "fourdst.log", "bbq_output.txt"):
        (out_dir / scratch).unlink(missing_ok=True)
    print(f"{case['name']}: {len(isos)} isotopes, {len(excluded)} reactions excluded, "
          f"{sum(1 for l in raw[1:] if l.split())} rows -> {out_dir / 'reference.csv'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bbq", type=Path, required=True, help="bbq executable (built with the paper patches)")
    parser.add_argument("--setup", type=Path, required=True, help="gf_mesa_setup from a GridFire build")
    parser.add_argument("--case", action="append", help="only generate these cases (default: all)")
    args = parser.parse_args()
    args.bbq, args.setup = args.bbq.resolve(), args.setup.resolve()

    mesa_dir = os.environ.get("MESA_DIR")
    if not mesa_dir:
        sys.exit("MESA_DIR must be set to the MESA installation bbq was built against")

    provenance_base = {
        "mesa_dir": mesa_dir,
        "mesa_version": mesa_version(mesa_dir),
        "bbq_version": git_describe(args.bbq.parent),
        "patches_sha256": {
            "mesa-r23.05.1-net_burn.patch": sha256(HERE / "mesa-r23.05.1-net_burn.patch"),
            "mesa-r23.05.1-no-screening.patch": sha256(HERE / "mesa-r23.05.1-no-screening.patch"),
            **{p.name: sha256(p) for p in sorted(PATCH_DIR.glob("lib_*.patch"))},
        },
        "bbq_disable_eos": "1",
        "screening_mode": SCREENING_MODE,
        "gridfire_version": git_describe(HERE),
        "generated": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
    }

    cases = json.loads((HERE / "cases.json").read_text())["cases"]
    selected = [c for c in cases if not args.case or c["name"] in args.case]
    if not selected:
        sys.exit(f"no matching cases in cases.json: {args.case}")
    for case in selected:
        generate_case(case, args, provenance_base)


if __name__ == "__main__":
    main()
