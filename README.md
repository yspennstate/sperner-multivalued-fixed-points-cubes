# Sperner's Lemma and Fixed Points of Multivalued Maps on Cubes

Yitzchak Shmalo

For an upper semicontinuous map $f$ of $[-1,1]^d$ with nonempty compact convex values, label each vertex $v$ of a
cubical subdivision by the coordinate signs of $y_v-v$, where $y_v\in f(v)$, inward at a zero coordinate on the
boundary. Some cube then carries all $2^d$ labels or has a facet whose labels are not the corners on one side of a
hyperplane, and as the mesh tends to zero such cubes or facets shrink, along a subsequence, to a fixed point. No
Lipschitz bound is assumed. In the plane this alternative is a case of Musin's lemma for quadrangulations. In general it
follows from a Sperner lemma for polytopal cells with nonzero boundary degree and a point joined to each label missing
from a facet by a segment avoiding the hull of that facet's labels. Kakutani's interpolation argument gives Kakutani's
theorem on the cube and a fixed point theorem of Bich.

`main.tex` and `main.pdf` are the paper; `Sperner.png` and `crosssectioncube.png` are its two figures.
`code/fan_facets.py` is the search for d <= 4 cited in the introduction (numpy and scipy).
