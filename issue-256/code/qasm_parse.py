"""Minimal OpenQASM 3 parser for the OLE circuits (own code, no external parser).

Supports: header, include, qubit[N] q;, user `gate name args { body }` definitions
(body made of known gates without parameters), barrier, and the gates
rx, rz, s, sdg, sx, h, cz, x, z (+ user-defined ones). Parameters are evaluated with a
tiny AST evaluator (numbers, pi, + - * / and unary minus only).

Matrices follow the OpenQASM 3 stdgates.inc conventions:
  rx(t) = exp(-i t X/2), rz(l) = exp(-i l Z/2) (stdgates: gphase(-l/2) U(0,0,l)),
  s = diag(1, i), sdg = diag(1, -i), sx = (1/2)[[1+i,1-i],[1-i,1+i]], h, cz = diag(1,1,1,-1).
Global phases never matter for f = 2^-n Tr(O C^dag O C).
"""
import ast
import math
import re
import numpy as np

_ALLOWED_BIN = {ast.Add: lambda a, b: a + b, ast.Sub: lambda a, b: a - b,
                ast.Mult: lambda a, b: a * b, ast.Div: lambda a, b: a / b}


def eval_param(expr):
    tree = ast.parse(expr.strip(), mode="eval")

    def ev(node):
        if isinstance(node, ast.Expression):
            return ev(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.Name) and node.id in ("pi", "π"):
            return math.pi
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -ev(node.operand)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.UAdd):
            return ev(node.operand)
        if isinstance(node, ast.BinOp) and type(node.op) in _ALLOWED_BIN:
            return _ALLOWED_BIN[type(node.op)](ev(node.left), ev(node.right))
        raise ValueError(f"unsupported expression: {expr!r}")

    return ev(tree)


def mat1(name, params=()):
    c = complex
    if name == "rx":
        t = params[0]
        co, si = math.cos(t / 2), math.sin(t / 2)
        return np.array([[co, -1j * si], [-1j * si, co]], dtype=np.complex128)
    if name == "rz":
        l = params[0]
        return np.array([[np.exp(-0.5j * l), 0], [0, np.exp(0.5j * l)]], dtype=np.complex128)
    if name == "s":
        return np.array([[1, 0], [0, 1j]], dtype=np.complex128)
    if name == "sdg":
        return np.array([[1, 0], [0, -1j]], dtype=np.complex128)
    if name == "sx":
        return 0.5 * np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]], dtype=np.complex128)
    if name == "h":
        return np.array([[1, 1], [1, -1]], dtype=np.complex128) / math.sqrt(2)
    if name == "x":
        return np.array([[0, 1], [1, 0]], dtype=np.complex128)
    if name == "z":
        return np.array([[1, 0], [0, -1]], dtype=np.complex128)
    raise KeyError(name)


BUILTIN_1Q = {"rx", "rz", "s", "sdg", "sx", "h", "x", "z"}
BUILTIN_2Q = {"cz"}


def _strip_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


_qref = re.compile(r"^\s*([A-Za-z_]\w*)\s*\[\s*(\d+)\s*\]\s*$")
_stmt = re.compile(r"^([A-Za-z_]\w*)\s*(\((.*)\))?\s*(.*)$", re.S)


def parse_qasm(path):
    """Return dict with 'nreg', 'gates' (list of dicts in time order), 'barriers' (positions)."""
    text = _strip_comments(open(path, encoding="utf-8").read())
    userdefs = {}
    # extract gate definitions
    def repl(m):
        name, args, body = m.group(1), m.group(2), m.group(3)
        argnames = [a.strip() for a in args.split(",") if a.strip()]
        ops = []
        for st in body.split(";"):
            st = st.strip()
            if not st:
                continue
            mm = _stmt.match(st)
            gname, targ = mm.group(1), mm.group(4).strip()
            if mm.group(2):
                raise ValueError("parametrised gates in user definitions unsupported")
            ops.append((gname, [argnames.index(t.strip()) for t in targ.split(",")]))
        userdefs[name] = (len(argnames), ops)
        return ""
    text = re.sub(r"gate\s+([A-Za-z_]\w*)\s+([^{]*)\{([^}]*)\}", repl, text)

    gates, barriers, nreg, regname = [], [], None, None
    for st in text.split(";"):
        st = st.strip()
        if not st:
            continue
        if st.startswith("OPENQASM") or st.startswith("include"):
            continue
        m = re.match(r"^qubit\s*\[\s*(\d+)\s*\]\s*([A-Za-z_]\w*)$", st)
        if m:
            nreg, regname = int(m.group(1)), m.group(2)
            continue
        mm = _stmt.match(st)
        if not mm:
            raise ValueError(f"cannot parse statement {st!r}")
        name, pstr, targ = mm.group(1), mm.group(3), mm.group(4)
        qs = []
        for t in targ.split(","):
            r = _qref.match(t)
            if not r or r.group(1) != regname:
                raise ValueError(f"bad qubit ref {t!r} in {st!r}")
            qs.append(int(r.group(2)))
        if name == "barrier":
            barriers.append(len(gates))
            continue
        params = tuple(eval_param(p) for p in pstr.split(",")) if pstr else ()
        if name in BUILTIN_1Q:
            assert len(qs) == 1, st
            gates.append(dict(name=name, qubits=(qs[0],), params=params,
                              mat=mat1(name, params), src=st))
        elif name in BUILTIN_2Q:
            assert len(qs) == 2 and qs[0] != qs[1], st
            gates.append(dict(name=name, qubits=(qs[0], qs[1]), params=params, mat=None, src=st))
        elif name in userdefs:
            nargs, ops = userdefs[name]
            assert len(qs) == nargs and nargs == 1 and not params, st
            m = np.eye(2, dtype=np.complex128)
            for gname, idx in ops:
                assert gname in BUILTIN_1Q
                m = mat1(gname) @ m          # later op multiplies from the left
            gates.append(dict(name=name, qubits=(qs[0],), params=params, mat=m, src=st))
        else:
            raise ValueError(f"unknown gate {name!r}")
    return dict(nreg=nreg, gates=gates, barriers=barriers, userdefs=userdefs)


if __name__ == "__main__":
    import sys, collections, hashlib
    p = sys.argv[1]
    print("sha256", hashlib.sha256(open(p, "rb").read()).hexdigest())
    c = parse_qasm(p)
    cnt = collections.Counter(g["name"] for g in c["gates"])
    used = sorted({q for g in c["gates"] for q in g["qubits"]})
    print("register", c["nreg"], "used qubits", len(used), "gates", len(c["gates"]), dict(cnt))
    print("barriers", len(c["barriers"]))
