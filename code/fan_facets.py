"""Labelings of the (d-1)-cube by corners of the d-cube that send adjacent vertices to equal or adjacent corners.

For d = 3 and d = 4 the script lists the distinct label sets and checks that each lies in a proper face of the d-cube
or is an affine cut set, that is, is strictly separated from the remaining corners by a hyperplane. A candidate
hyperplane is found by linear programming in floating point; it is then rounded to rationals and the strict
inequalities w.c > t on the set and w.c < t off it are checked in exact arithmetic, so a set is counted as an affine
cut set only with an exact certificate.
Requires numpy and scipy.
"""
from fractions import Fraction

import numpy as np
from scipy.optimize import linprog


def corner(c, m):
    return [2 * ((c >> j) & 1) - 1 for j in range(m)]


def exact_cut(A, m, w, t):
    for c in range(1 << m):
        s = sum(wi * xi for wi, xi in zip(w, corner(c, m))) - t
        if (c in A and s <= 0) or (c not in A and s >= 0):
            return False
    return True


def separable(A, m):
    rows = []
    for c in range(1 << m):
        s = -1 if c in A else 1          # w.x - t >= 1 on A, w.x - t <= -1 off A
        rows.append([s * x for x in corner(c, m)] + [-s])
    res = linprog(np.zeros(m + 1), A_ub=np.array(rows, float), b_ub=-np.ones(1 << m),
                  bounds=[(None, None)] * (m + 1), method="highs")
    if res.status != 0:
        return False
    for den in (1, 2, 4, 8, 16, 64, 256, 1024, 10 ** 6):
        q = [Fraction(v).limit_denominator(den) for v in res.x]
        if exact_cut(A, m, q[:m], q[m]):
            return True
    return False


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
    in_face = [A for A in sets if any(len({(x >> j) & 1 for x in A}) == 1 for j in range(d))]
    rest = [A for A in sets if A not in in_face]
    cut = [A for A in rest if separable(A, d)]
    print(f"d = {d}: {len(sets)} label sets, {len(in_face)} in a proper face, {len(cut)} affine cut sets "
          f"with an exact certificate, {len(rest) - len(cut)} neither")
