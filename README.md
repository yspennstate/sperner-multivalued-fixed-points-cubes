# Sperner's Lemma and Fixed Points of Multivalued Maps on Cubes

Yitzchak Shmalo

For an upper semicontinuous map $f$ of $[-1,1]^d$ with nonempty compact convex values, label each vertex $v$ of a
cubical subdivision by the coordinate signs of $y_v-v$, where $y_v\in f(v)$, inward at a zero coordinate on the
boundary. Some cube then carries all $2^d$ labels or has a facet whose labels are not the corners on one side of a
hyperplane, and along a subsequence of subdivisions of vanishing mesh some sequence of such cubes, or of such facets,
shrinks to a fixed point of $f$. For Lipschitz maps, Gidea and the author found cubes carrying $d+1$ sign labels with
both signs in every coordinate; here no Lipschitz bound is assumed. In the plane this alternative is a case of Musin's
lemma for quadrangulations. In general it follows from a Sperner lemma for cells labeled by the vertices of a polytope,
which holds when the piecewise affine map determined by the labels has nonzero degree on the boundary and some interior
point of the target is joined to each label missing from a facet by a segment that avoids the convex hull of the labels
on that facet. Kakutani's interpolation argument, with Sperner's lemma in place of Brouwer's theorem, gives on a cube
Kakutani's theorem and a fixed point theorem of Bich for maps that need not be upper semicontinuous.

`main.tex` and `main.pdf` are the paper; `Sperner.png` and `crosssectioncube.png` are its two figures.
`code/fan_facets.py` checks the facets of Fan labelings for d = 3, 4 (numpy and scipy), and
`code/central_cube3.py` checks in exact arithmetic that a 3-cube whose labels miss a corner and which has no
problematic facet is not centrally labeled (numpy). Their outputs are in `code/*_output.txt`.
