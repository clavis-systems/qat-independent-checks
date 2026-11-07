"""Minimal OpenQASM 2.0 parser written for this check (no external parser).

Supports the subset used by the SU(2) LSH circuits: one qreg, gates with
optional numeric parameters (expressions with pi allowed), barriers and
measurements ignored. Raises on anything unknown.
"""
import math
import re

_GATE_RE = re.compile(r"^([a-z_][a-z0-9_]*)\s*(\(([^)]*)\))?\s+(.+)$")
_QARG_RE = re.compile(r"^([a-z_][a-z0-9_]*)\[(\d+)\]$")

KNOWN = {"x": (1, 0), "h": (1, 0), "rz": (1, 1), "cx": (2, 0), "swap": (2, 0),
         "y": (1, 0), "z": (1, 0), "s": (1, 0), "sdg": (1, 0), "sx": (1, 0),
         "rx": (1, 1), "ry": (1, 1), "p": (1, 1), "u1": (1, 1), "cz": (2, 0),
         "rzz": (2, 1), "id": (1, 0)}


def _eval_param(expr):
    expr = expr.strip()
    if not re.fullmatch(r"[0-9eE+\-*/(). pi]+", expr):
        raise ValueError("unsupported parameter expression: %r" % expr)
    return float(eval(expr, {"__builtins__": {}}, {"pi": math.pi}))


def parse_qasm(path):
    """Return (nqubits, gates) with gates = list of (name, params tuple, qubits tuple)."""
    with open(path, "r", encoding="utf-8") as f:
        text = f.read()
    text = re.sub(r"//[^\n]*", "", text)
    stmts = [s.strip() for s in text.split(";")]
    nq = None
    regname = None
    gates = []
    for s in stmts:
        if not s:
            continue
        if s.startswith("OPENQASM") or s.startswith("include"):
            continue
        if s.startswith("qreg"):
            m = re.fullmatch(r"qreg\s+([a-z_][a-z0-9_]*)\[(\d+)\]", s)
            if nq is not None or m is None:
                raise ValueError("unexpected qreg: %r" % s)
            regname, nq = m.group(1), int(m.group(2))
            continue
        if s.startswith("creg") or s.startswith("barrier") or s.startswith("measure"):
            continue
        m = _GATE_RE.match(s)
        if m is None:
            raise ValueError("cannot parse statement: %r" % s)
        name = m.group(1)
        params = tuple(_eval_param(p) for p in m.group(3).split(",")) if m.group(3) else ()
        qargs = []
        for qa in m.group(4).split(","):
            mm = _QARG_RE.match(qa.strip())
            if mm is None or mm.group(1) != regname:
                raise ValueError("bad qarg in %r" % s)
            q = int(mm.group(2))
            if not 0 <= q < nq:
                raise ValueError("qubit out of range in %r" % s)
            qargs.append(q)
        if name not in KNOWN:
            raise ValueError("unknown gate %r" % name)
        nqa, npar = KNOWN[name]
        if len(qargs) != nqa or len(params) != npar:
            raise ValueError("arity mismatch in %r" % s)
        if len(set(qargs)) != len(qargs):
            raise ValueError("repeated qubit in %r" % s)
        gates.append((name, params, tuple(qargs)))
    return nq, gates
