# degen_vert2lcon
 
Find a linearly constrained convex hull of a set of points.
In the case that the convex hull lies in a degenerate affine subspace, identify and return linear equality constraints (hyperplanes) characterizing said subspace.

## Usage

[ $A$, $b$, $A_{eq}$, $b_{eq}$ ] = degen_vert2lcon(Points)

## Output

- [ $A_{eq}$, $b_{eq}$ ]
    - $A_eq \in \mathbb{R}^{n_0 \times d}$
    - $b_eq \in \mathbb{R}^{n_0 \times 1}$

Linear equality constraints that characterize the subspace on which the convex hull lies. x is in the releant subspace iff $A_{eq}x - b_{eq} = 0$.

- [ $A$, $b$ ]
    - $A \in \mathbb{R}^{n_1 \times d}$
    - $b \in \mathbb{R}^{n_1 \times 1}$

Linear inequality constraints that characterize the convex hull. $x$ is in the convex hull iff $Ax < b$ and $A_{eq} x = b_{eq}$.

## Input

$Points$ - $d \times N$ matrix where each column is a point in $\mathbb{R}^d$.