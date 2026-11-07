import numpy as np
from qmps import *

def split_file(gates):
    """returns init bit flips (list of x on qubits before first non-x) and per-step gate lists"""
    k = 0
    while gates[k][0] == 'x': k += 1
    init = gates[:k]; body = gates[k:]
    starts = [i for i in range(len(body)-1) if body[i][0]=='cx' and body[i][2]==(1,0) and body[i+1][0]=='h' and body[i+1][2]==(1,)]
    starts.append(len(body))
    steps = [body[starts[i]:starts[i+1]] for i in range(len(starts)-1)]
    assert starts[0]==0
    return init, steps

def bits_from_init(nq, init):
    b = np.zeros(nq, int)
    for n,a,q in init: b[q[0]] ^= 1
    return b
