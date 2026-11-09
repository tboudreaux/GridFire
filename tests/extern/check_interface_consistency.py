#!/usr/bin/env python3
"""
Check that every bind(C) interface in the Fortran module matches the C API it binds to.

For each Fortran `function`/`subroutine ... bind(C, name="gf_x")` the C header must declare gf_x with the same number
of arguments; each argument must agree on pointer vs value passing and, where the C type is a plain scalar, on the
C interop kind (double <-> c_double, int/enum <-> c_int, size_t <-> c_size_t, char <-> c_char); the result types
must agree the same way. A void* parameter accepts any pointer-like Fortran argument.

Usage: check_interface_consistency.py <gridfire_extern.h> <gridfire_mod.f90>. Exits 1 and lists every mismatch.
"""
import re
import sys

C_KINDS = {"double": "c_double", "int": "c_int", "size_t": "c_size_t", "char": "c_char"}


def parse_c_type(text, enums):
    """(base kind or None, is_pointer) for a C parameter or return type."""
    text = re.sub(r"\b(const|struct)\b", " ", text).strip()
    pointer = "*" in text
    base = text.replace("*", " ").split()
    base = base[0] if base else "void"
    if base == "enum":  # 'enum NAME' is an int
        return "c_int", pointer
    if base in enums:
        return "c_int", pointer
    return C_KINDS.get(base), pointer  # None: void or a struct (opaque)


def parse_header(path):
    text = open(path).read()
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    text = re.sub(r"//[^\n]*", " ", text)
    enums = set(re.findall(r"enum\s+(\w+)\s*\{", text))
    functions = {}
    for ret, name, params in re.findall(r"([\w\s\*]+?)\b(gf_\w+)\s*\(([^)]*)\)\s*;", text):
        ret = ret.split(";")[-1].split("}")[-1]
        params = params.strip()
        args = [] if params in ("", "void") else [p.strip() for p in params.split(",")]
        arg_types = []
        for a in args:
            # drop the parameter name (last identifier, unless the declaration ends in '*')
            a = re.sub(r"\[\s*\]$", "*", a)
            a = re.sub(r"\s*\b\w+$", "", a) if not a.rstrip().endswith("*") else a
            arg_types.append(parse_c_type(a, enums))
        functions[name] = (parse_c_type(ret, enums), arg_types)
    return functions


def parse_fortran(path):
    lines = open(path).read().splitlines()
    # join continuation lines and drop comments
    joined, buf = [], ""
    for line in lines:
        line = line.split("!")[0].rstrip()
        if line.endswith("&"):
            buf += line[:-1] + " "
            continue
        joined.append((buf + line.lstrip("& ")).strip())
        buf = ""
    interfaces = []
    i = 0
    header = re.compile(
        r"^(?:(?P<rtype>[\w\(\), ]+?)\s+)?(?P<kind>function|subroutine)\s+(?P<fname>\w+)\s*\((?P<args>[^)]*)\)"
        r"(?:\s*result\s*\((?P<result>\w+)\))?\s*bind\s*\(\s*C\s*,\s*name\s*=\s*\"(?P<cname>\w+)\"\s*\)", re.I)
    while i < len(joined):
        m = header.match(joined[i])
        if not m:
            i += 1
            continue
        args = [a.strip().lower() for a in m.group("args").split(",") if a.strip()]
        result = (m.group("result") or m.group("fname")).lower() if m.group("kind").lower() == "function" else None
        decls = {}
        i += 1
        while i < len(joined) and not re.match(r"^end\s+(function|subroutine)", joined[i], re.I):
            d = re.match(r"^(?P<type>[^:]+?)\s*::\s*(?P<names>.+)$", joined[i])
            if d:
                attrs = d.group("type").lower()
                for name in d.group("names").split(","):
                    decls[name.strip().split("(")[0].strip().lower()] = attrs
            i += 1
        interfaces.append((m.group("fname"), m.group("cname"), args, result, decls))
    return interfaces


def fortran_kind(attrs):
    """(kind, is_pointer_like) for a Fortran dummy argument's declaration."""
    kind = None
    m = re.search(r"\b(c_double|c_int|c_size_t|c_char|c_ptr|c_funptr)\b", attrs)
    if m:
        kind = m.group(1)
    if kind is None and re.match(r"^\s*type\s*\(\s*\w+\s*\)", attrs):
        kind = "derived"  # a bind(C) derived type, i.e. a struct
    by_value = re.search(r"\bvalue\b", attrs) is not None
    pointer_like = kind == "c_ptr" or not by_value
    return kind, pointer_like


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    c_functions = parse_header(sys.argv[1])
    interfaces = parse_fortran(sys.argv[2])
    problems = []
    interfaces = [i for i in interfaces if i[1].startswith("gf_")]  # bindings to libc (e.g. strlen) are not ours
    for fname, cname, args, result, decls in interfaces:
        where = f"{fname} (bind C name {cname})"
        if cname not in c_functions:
            problems.append(f"{where}: no such function in the C header")
            continue
        (c_ret_kind, c_ret_ptr), c_args = c_functions[cname]
        if len(args) != len(c_args):
            problems.append(f"{where}: {len(args)} Fortran arguments, {len(c_args)} in C")
            continue
        for n, (arg, (c_kind, c_ptr)) in enumerate(zip(args, c_args), start=1):
            if arg not in decls:
                problems.append(f"{where}: argument {arg} has no declaration")
                continue
            f_kind, f_ptr = fortran_kind(decls[arg])
            if c_ptr != f_ptr:
                problems.append(f"{where}: argument {n} ({arg}) is {'a pointer' if c_ptr else 'by value'} in C "
                                f"but {'a pointer' if f_ptr else 'by value'} in Fortran")
            elif c_kind is not None and f_kind not in (c_kind, "c_ptr"):
                problems.append(f"{where}: argument {n} ({arg}) is {c_kind} in C but {f_kind} in Fortran")
        if result is not None:
            f_kind, _ = fortran_kind(decls.get(result, ""))
            expected = "c_ptr" if c_ret_ptr else c_ret_kind
            if expected is not None and f_kind != expected:
                problems.append(f"{where}: returns {expected} in C but {f_kind} in Fortran")
        elif c_ret_kind is not None or c_ret_ptr:
            problems.append(f"{where}: a Fortran subroutine, but the C function returns a value")

    print(f"Checked {len(interfaces)} Fortran bind(C) interfaces against {len(c_functions)} C functions")
    for p in problems:
        print("MISMATCH:", p)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
