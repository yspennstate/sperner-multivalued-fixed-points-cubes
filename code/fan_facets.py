"""Labelings of the (d-1)-cube by corners of the d-cube that send adjacent vertices to equal or adjacent corners.

For d = 3 and d = 4 the script lists the distinct label sets and checks that each lies in a proper face of the d-cube
or is an affine cut set, that is, is strictly separated from the remaining corners by a hyperplane.
Requires numpy and scipy.
"""
import numpy as np
from scipy.optimize import linprog


def separable(A, m):
    corners = [np.array([(c >> j) & 1 for j in range(m)], float) * 2 - 1 for c in range(1 << m)]
    rows = []
    for c in range(1 << m):
        s = -1 if c in A else 1          # w.x - t >= 1 on A, w.x - t <= -1 off A
        rows.append(list(s * corners[c]) + [-s])
    res = linprog(np.zeros(m + 1), A_ub=np.array(rows), b_ub=-np.ones(1 << m),
                  bounds=[(None, None)] * (m + 1), method="highs")
    return res.status == 0


def label_sets(n):
    m, N = n + 1, 1 << n
    order = sorted(range(N), key=lambda v: (bin(v).count("1"), v))
    nbrs = {v: [v ^ (1 << i) for i in range(n)] for v in range(N)}
    found, lab = set(), {}

    def extend(k):
        if k == N:
            found.add(frozenset(lab.values()))
            return
        v = order[k]
        for c in range(1 << m):
            if all(bin(c ^ lab[u]).count("1") <= 1 for u in nbrs[v] if u in lab):
                lab[v] = c
                extend(k + 1)
                del lab[v]

    for c in range(1 << m):
        lab[order[0]] = c
        extend(1)
        del lab[order[0]]
    return found


for d in (3, 4):
    sets = label_sets(d - 1)
    bad = [A for A in sets
           if not any(len({(x >> j) & 1 for x in A}) == 1 for j in range(d)) and not separable(A, d)]
    print(f"d = {d}: {len(sets)} label sets, {len(bad)} neither in a proper face nor an affine cut set")
