"""Labelings of the vertices of the 3-cube by its corners (as vectors in {-1,1}^3).

Claim checked: if the labels miss a corner and every facet's label set is an affine cut set (strictly separated
from the remaining corners by a hyperplane), then the trilinear interpolant of the labels has no zero in the closed
cube, so the cube is not centrally labeled in Musin's sense.

All arithmetic is exact. Affine cut sets are the sets cut off by a hyperplane w.x > b with w in {-3,...,3}^3 and b a
half-integer; there are 104 of them (empty and full set included), the number of threshold functions of three
variables, so none is missed. Whether 0 lies in the convex hull of a set of points is decided with integers or
fractions through Caratheodory's theorem. A labeling with 0 outside the hull of its labels has no zero, so only the
others are examined: up to the symmetries of the cube acting on the vertices and on the labels, each is checked by
cutting the cube into subcubes, on each of which the interpolant is again trilinear and so takes values in the convex
hull of its values at the corners of the subcube. Requires numpy.
"""
import itertools
from fractions import Fraction as Fr
import numpy as np

d = 3
V = [tuple(((c >> j) & 1) * 2 - 1 for j in range(d)) for c in range(8)]


def det3(a, b, c):
    return a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0])


def zero_in_hull(P):
    """0 in conv(P) for points P of R^3, exactly: 0 lies in a segment, triangle or tetrahedron of points of P."""
    P = list(set(P))
    if (0, 0, 0) in P:
        return True
    for a, b in itertools.combinations(P, 2):
        cr = (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
        if cr == (0, 0, 0) and any(a[j]*b[j] < 0 for j in range(d)):
            return True
    for a, b, c in itertools.combinations(P, 3):
        if det3(a, b, c) != 0:
            continue
        n = ((b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]), (b[2]-a[2])*(c[0]-a[0]) - (b[0]-a[0])*(c[2]-a[2]),
             (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]))
        k = max(range(3), key=lambda j: abs(n[j]))
        if n[k] == 0:
            continue
        i0, i1 = [j for j in range(3) if j != k]

        def ar(p, q, r):
            return (q[i0]-p[i0])*(r[i1]-p[i1]) - (q[i1]-p[i1])*(r[i0]-p[i0])
        o = (0, 0, 0)
        s = (ar(o, b, c), ar(a, o, c), ar(a, b, o))
        if min(s) >= 0 or max(s) <= 0:
            return True
    for a, b, c, e in itertools.combinations(P, 4):
        def sd(p, q, r, t):
            return det3(tuple(q[j]-p[j] for j in range(3)), tuple(r[j]-p[j] for j in range(3)),
                        tuple(t[j]-p[j] for j in range(3)))
        if sd(a, b, c, e) == 0:
            continue
        o = (0, 0, 0)
        t = [sd(o, b, c, e), sd(a, o, c, e), sd(a, b, o, e), sd(a, b, c, o)]
        if min(t) >= 0 or max(t) <= 0:
            return True
    return False


# exact tables over the 256 sets of corners (bit c of m set = corner c in the set)
cutsets = set()
for w in itertools.product(range(-3, 4), repeat=3):
    vals = [sum(a * b for a, b in zip(w, v)) for v in V]
    for b2 in range(-21, 22, 2):
        cutsets.add(sum(1 << c for c in range(8) if 2 * vals[c] > b2))
assert len(cutsets) == 104
cut = np.array([m in cutsets for m in range(256)])
hz = np.array([m > 0 and zero_in_hull([V[c] for c in range(8) if m >> c & 1]) for m in range(256)])
facets = [[v for v in range(8) if (v >> j & 1) == s] for j in range(d) for s in (0, 1)]


def interp(x, L):
    out = [Fr(0)] * d
    for v in range(8):
        w = Fr(1)
        for j in range(d):
            w *= x[j] if v >> j & 1 else 1 - x[j]
        for j in range(d):
            out[j] += w * V[L[v]][j]
    return tuple(out)


def boxes_with_zero(L, maxdepth):
    """Subcubes of [0,1]^3 (side 2^-maxdepth) not excluded by the hull test; [] proves there is no zero."""
    stack, left = [((Fr(0),) * 3, Fr(1), 0)], []
    while stack:
        lo, s, dep = stack.pop()
        vals = [interp(tuple(lo[j] + s * ((c >> j) & 1) for j in range(d)), L) for c in range(8)]
        if not zero_in_hull(vals):
            continue
        if dep >= maxdepth:
            left.append((lo, s))
            continue
        h = s / 2
        for c in range(8):
            stack.append((tuple(lo[j] + h * ((c >> j) & 1) for j in range(d)), h, dep + 1))
    return left


# control: Musin's example (identity except that the vertex (1,1,1) gets the label of (1,1,0)) has a zero at
# (1/2,1/2,2/3) and its facet x_3 = 1 is not an affine cut set
musin = [0, 1, 2, 3, 4, 5, 6, 3]
assert interp((Fr(1, 2), Fr(1, 2), Fr(2, 3)), musin) == (0, 0, 0)
assert [F for F in facets if not cut[sum(1 << musin[v] for v in F)]] == [[4, 5, 6, 7]]
assert boxes_with_zero(musin, 6)

C = []
N = 8 ** 8
for start in range(0, N, 1 << 21):
    idx = np.arange(start, min(N, start + (1 << 21)), dtype=np.int64)
    L = np.empty((len(idx), 8), np.uint8)
    for v in range(8):
        L[:, v] = (idx >> (3 * v)) & 7
    bit = (np.uint16(1) << L.astype(np.uint16))
    allm = np.bitwise_or.reduce(bit, axis=1)
    ok = (allm != 255) & hz[allm]
    for F in facets:
        ok &= cut[np.bitwise_or.reduce(bit[:, F], axis=1)]
    C.append(L[ok])
C = np.concatenate(C)
print("labelings missing a corner, with no problematic facet and 0 in the hull of the labels:", len(C))
syms = []
for perm in itertools.permutations(range(d)):
    for flip in range(8):
        syms.append([sum((((v >> perm[j]) & 1) ^ ((flip >> j) & 1)) << j for j in range(d)) for v in range(8)])
syms = np.array(syms, np.int64)
w8 = (8 ** np.arange(8)).astype(np.int64)
canon = np.full(len(C), np.iinfo(np.int64).max, np.int64)
Ci = C.astype(np.int64)
for g in syms:
    Lg = Ci[:, g]
    for h in syms:
        canon = np.minimum(canon, h[Lg] @ w8)
reps = [[int((c >> (3 * v)) & 7) for v in range(8)] for c in np.unique(canon)]
print("classes under the symmetries of the cube (on vertices and on labels):", len(reps))
bad = 0
for L in reps:
    left = boxes_with_zero(L, 10)
    bad += bool(left)
    print(L, "no zero" if not left else f"{len(left)} subcubes not excluded")
print("classes whose interpolant may vanish:", bad)
