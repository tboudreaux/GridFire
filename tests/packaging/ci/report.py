#!/usr/bin/env python3
"""
Combine the packaging CI results (results-*.json, written by lib.sh's ci_finish) into one table of platforms by
checks, and optionally into shields.io endpoint badges, one per platform and check.

    report.py <results dir>                         print the Markdown table
    report.py <results dir> --present               only the platforms with results (one job's own summary)
    report.py <results dir> --badges <dir>          also write <dir>/<platform>-<check>.json and <dir>/status.md
    report.py --readme-table <owner/repo> <branch>  print the README's badge table (badges read from <branch>)

Exits non-zero if any check failed or a platform reported no results.
"""
import argparse
import json
import pathlib
import sys

# Display order and descriptions; checks not listed here are appended after them
CHECKS = {
    "wheel": "wheel builds and repairs as a release does",
    "python-install": "wheel installs next to the fourdst wheel and imports",
    "gf-config": "C++, C and Fortran programs built with gf-config flags",
    "pkg-config": "C++, C and Fortran programs built with the wheel's pkg-config files",
    "threads-libgomp": "concurrent contexts from a GCC/libgomp OpenMP host program",
    "cmake-find": "CMake find_package against the wheel (C++, C, Fortran)",
    "cmake-fetch": "CMake FetchContent building GridFire from source (C++, C, Fortran)",
    "meson-subproject": "Meson subproject from a clean source tree (C++, C, Fortran)",
}
PLATFORMS = ["linux-x86_64", "linux-aarch64", "macos-arm64"]
SYMBOL = {"pass": "✅ pass", "fail": "❌ fail", "skipped": "⚪ skipped", None: "❔ no result"}
COLOUR = {"pass": "brightgreen", "fail": "red", "skipped": "lightgrey", None: "lightgrey"}


def load(results_dir):
    results = {}
    for path in sorted(pathlib.Path(results_dir).rglob("results-*.json")):
        data = json.loads(path.read_text())
        results.setdefault(data["platform"], {}).update(data["checks"])
    return results


def table(results, present_only=False):
    platforms = PLATFORMS + sorted(p for p in results if p not in PLATFORMS)
    if present_only:
        platforms = [p for p in platforms if p in results]
    reported = {c for p in platforms for c in results.get(p, {})}
    checks = [c for c in CHECKS if c in reported or not present_only] + sorted(reported - set(CHECKS))
    lines = ["| Check | " + " | ".join(platforms) + " |", "|---" * (len(platforms) + 1) + "|"]
    for check in checks:
        cells = [SYMBOL[results.get(p, {}).get(check)] for p in platforms]
        lines.append(f"| `{check}` | " + " | ".join(cells) + " |")
    return "\n".join(lines), platforms, checks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir", nargs="?")
    parser.add_argument("--badges", help="write shields.io endpoint JSON files and status.md here")
    parser.add_argument("--present", action="store_true", help="only platforms that reported results")
    parser.add_argument("--readme-table", nargs=2, metavar=("REPO", "BRANCH"))
    args = parser.parse_args()

    if args.readme_table:
        repo, branch = args.readme_table
        print("| Check | " + " | ".join(PLATFORMS) + " |")
        print("|---" * (len(PLATFORMS) + 1) + "|")
        for check, description in CHECKS.items():
            cells = []
            for platform in PLATFORMS:
                url = f"https://raw.githubusercontent.com/{repo}/{branch}/badges/{platform}-{check}.json"
                cells.append(f"![{platform} {check}](https://img.shields.io/endpoint?url={url})")
            print(f"| {description} | " + " | ".join(cells) + " |")
        return 0

    if not args.results_dir:
        parser.error("a results directory is required")
    results = load(args.results_dir)
    text, platforms, checks = table(results, args.present)
    print(text)

    if args.badges:
        out = pathlib.Path(args.badges)
        out.mkdir(parents=True, exist_ok=True)
        for platform in platforms:
            for check in checks:
                status = results.get(platform, {}).get(check)
                badge = {"schemaVersion": 1, "label": check, "message": status or "no result",
                         "color": COLOUR[status]}
                (out / f"{platform}-{check}.json").write_text(json.dumps(badge) + "\n")
        (out / "status.md").write_text("# GridFire packaging checks\n\n" + text + "\n")

    missing = [] if args.present else [p for p in PLATFORMS if p not in results]
    failed = [(p, c) for p, r in results.items() for c, s in r.items() if s == "fail"]
    for p in missing:
        print(f"no results from {p}", file=sys.stderr)
    for p, c in failed:
        print(f"{c} failed on {p}", file=sys.stderr)
    return 1 if missing or failed else 0


if __name__ == "__main__":
    sys.exit(main())
