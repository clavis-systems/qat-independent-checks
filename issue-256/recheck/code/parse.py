"""Own QASM parser for the OLE circuit. Returns list of (name, angle_or_None, qubits)."""
import re, math

def ang(s):
    s = s.strip().replace("pi", "(%r)" % math.pi)
    return float(eval(s, {"__builtins__": {}}))

def parse(path):
    gates = []
    in_gate_def = False
    for line in open(path):
        line = line.strip()
        if line.startswith("gate "):
            in_gate_def = True
        if in_gate_def:
            if line.startswith("}"):
                in_gate_def = False
            continue
        if not line or line.startswith(("OPENQASM", "include", "qubit")):
            continue
        if line.startswith("barrier"):
            gates.append(("barrier", None, ()))
            continue
        m = re.match(r"(\w+)(?:\((.*)\))?\s+(.*);", line)
        name, a, qs = m.groups()
        qs = tuple(int(x) for x in re.findall(r"q\[(\d+)\]", qs))
        gates.append((name, ang(a) if a else None, qs))
    return gates

if __name__ == "__main__":
    import sys
    g = parse(sys.argv[1])
    print(len(g), "items")
