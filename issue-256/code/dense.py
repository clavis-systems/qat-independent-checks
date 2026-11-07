"""Independent dense references for small sub-circuits.

own_dense_f : own numpy unitary simulator acting gate by gate on raw parsed gates
              (no merging, no peephole, no light cone), f = 2^-k sum s_a s_b |C_ba|^2.
qiskit_f    : circuit parsed by qiskit.qasm3 (independent parser), restricted to the same qubits,
              unitary from qiskit.quantum_info.Operator.
"""
import numpy as np


def own_dense_f(gates, qubits, obs):
    qubits = list(qubits)
    k = len(qubits)
    pos = {q: i for i, q in enumerate(qubits)}
    U = np.eye(2 ** k, dtype=np.complex128).reshape((2,) * k + (2 ** k,))
    for g in gates:
        if g["name"] == "cz":
            a, b = pos[g["qubits"][0]], pos[g["qubits"][1]]
            idx = [slice(None)] * (k + 1)
            idx[a] = 1; idx[b] = 1
            U[tuple(idx)] *= -1
        else:
            a = pos[g["qubits"][0]]
            U = np.tensordot(g["mat"], U, axes=([1], [a]))
            U = np.moveaxis(U, 0, a)
    C = U.reshape(2 ** k, 2 ** k)          # row = output b (qubit 0 most significant), col = input a
    s = np.ones(2 ** k)
    for q in obs:
        if q in pos:
            bit = (np.arange(2 ** k) >> (k - 1 - pos[q])) & 1
            s = s * (1 - 2 * bit)
    P = np.abs(C) ** 2
    return float(s @ P @ s) / 2 ** k


def qiskit_f(qasm_path, qubits, obs, upto=None, gate_filter=None):
    import qiskit.qasm3
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Operator
    qc = qiskit.qasm3.loads(open(qasm_path).read())
    qubits = list(qubits)
    pos = {q: i for i, q in enumerate(qubits)}
    sub = QuantumCircuit(len(qubits))
    for inst in qc.data:
        if inst.operation.name == "barrier":
            continue
        idx = [qc.find_bit(b).index for b in inst.qubits]
        if all(i in pos for i in idx):
            sub.append(inst.operation, [pos[i] for i in idx])
    C = Operator(sub).data                   # qiskit: little-endian, bit j of index = qubit j
    k = len(qubits)
    s = np.ones(2 ** k)
    for q in obs:
        if q in pos:
            bit = (np.arange(2 ** k) >> pos[q]) & 1
            s = s * (1 - 2 * bit)
    P = np.abs(C) ** 2
    return float(s @ P @ s) / 2 ** k, sub.size()


def qiskit_gate_list(qasm_path):
    import qiskit.qasm3
    qc = qiskit.qasm3.loads(open(qasm_path).read())
    out = []
    for inst in qc.data:
        if inst.operation.name == "barrier":
            continue
        idx = tuple(qc.find_bit(b).index for b in inst.qubits)
        out.append((inst.operation.name, idx, tuple(float(p) for p in inst.operation.params)))
    return out
