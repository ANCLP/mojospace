# Barnes-Hut Algorithm and Hierarchical Trees (Quadtrees & Octrees)

## 1. Overview and Problem Statement

In gravitational and astrophysical n-body simulations, every particle interacts with every other particle via Newton's law of universal gravitation:

$$\mathbf{F}_i = \sum_{j \neq i} G \frac{m_i m_j}{\|\mathbf{r}_j - \mathbf{r}_i\|^3} (\mathbf{r}_j - \mathbf{r}_i)$$

Unlike short-range contact or fluid forces, gravity is **long-range** ($\propto 1/r^2$), meaning interactions do not diminish to zero past a cutoff radius $r_{\text{cut}}$. Evaluating all interactions directly requires evaluating all $N(N-1)/2$ pairs:
$$\mathcal{O}(N^2) \text{ operations}$$

For $N = 10^5$, $\mathcal{O}(N^2) \approx 5 \times 10^9$ interactions per time step.

The **Barnes-Hut algorithm** (1986) reduces this complexity to:
$$\mathcal{O}(N \log N)$$
by clustering distant groups of particles into single pseudoparticles located at their collective center of mass.

---

## 2. Tree Construction (Quadtree in 2D, Octree in 3D)

Space is partitioned recursively into hierarchical axis-aligned bounding boxes:
- **Quadtree (2D)**: Each internal node subdivides into $2^2 = 4$ quadrants (NW, NE, SW, SE).
- **Octree (3D)**: Each internal node subdivides into $2^3 = 8$ octants:
  - Octant 0: $(x < x_{\text{mid}}, y < y_{\text{mid}}, z < z_{\text{mid}})$
  - Octant 1: $(x \ge x_{\text{mid}}, y < y_{\text{mid}}, z < z_{\text{mid}})$
  - Octant 2: $(x < x_{\text{mid}}, y \ge y_{\text{mid}}, z < z_{\text{mid}})$
  - Octant 3: $(x \ge x_{\text{mid}}, y \ge y_{\text{mid}}, z < z_{\text{mid}})$
  - Octant 4: $(x < x_{\text{mid}}, y < y_{\text{mid}}, z \ge z_{\text{mid}})$
  - Octant 5: $(x \ge x_{\text{mid}}, y < y_{\text{mid}}, z \ge z_{\text{mid}})$
  - Octant 6: $(x < x_{\text{mid}}, y \ge y_{\text{mid}}, z \ge z_{\text{mid}})$
  - Octant 7: $(x \ge x_{\text{mid}}, y \ge y_{\text{mid}}, z \ge z_{\text{mid}})$

Subdivision continues until each leaf contains at most $K$ particles (typically $K = 1$).

---

## 3. Mass and Center-of-Mass Aggregation

During an upward post-order traversal of the tree, each node $C$ aggregates the total mass and center of mass of all particles in its subtree:

### Total Mass:
$$M_C = \sum_{i \in C} m_i = \sum_{k=0}^{7} M_{\text{child}_k}$$

### Center of Mass:
$$\mathbf{R}_C = \frac{1}{M_C} \sum_{i \in C} m_i \mathbf{r}_i = \frac{1}{M_C} \sum_{k=0}^{7} M_{\text{child}_k} \mathbf{R}_{\text{child}_k}$$

This zeroth-order multipole expansion (monopole approximation) treats the entire volume as a single point mass $M_C$ positioned at $\mathbf{R}_C$.

---

## 4. Multipole Acceptance Criterion (MAC)

When computing the force exerted on a target particle $i$ located at $\mathbf{r}_i$, the algorithm traverses the tree starting from the root node:

Let:
- $s$ be the spatial side length of the node's bounding box:
  $$s = \max(x_{\max} - x_{\min}, y_{\max} - y_{\min}, z_{\max} - z_{\min})$$
- $d$ be the Euclidean distance from target $\mathbf{r}_i$ to the node's center of mass $\mathbf{R}_C$:
  $$d = \|\mathbf{R}_C - \mathbf{r}_i\|$$

The opening angle ratio is:
$$\theta = \frac{s}{d}$$

### Decision Rule:
Given an accuracy threshold parameter $\theta_{\text{crit}}$ (typically $\theta_{\text{crit}} \in [0.5, 0.8]$):
1. **If $\theta < \theta_{\text{crit}}$**: The node is sufficiently distant. Approximate the force contribution of all particles in the subtree by treating the entire node as a single particle of mass $M_C$ at $\mathbf{R}_C$. Do not recurse into children.
2. **If $\theta \ge \theta_{\text{crit}}$**:
   - If the node is an **internal node**, recurse into all active children.
   - If the node is a **leaf**, compute direct pairwise forces for all particles in the leaf ($j \neq i$).

```
             Target Particle (r_i)
                     *
                    / \
                   /   \  d = ||R_C - r_i||
                  /     \
                 v       v
            +-----------------+
            |  Node Center    |
            |     of Mass     |
            |      (R_C)      | s (node width)
            |                 |
            +-----------------+
            Condition: s / d < theta_crit
```

---

## 5. Plummer Softening

To prevent numerical divergence when two particles pass arbitrarily close to each other ($\|\mathbf{r}_j - \mathbf{r}_i\| \to 0$), the force law incorporates a **Plummer softening parameter** $\epsilon > 0$:

$$\mathbf{F}_{ij} = G \frac{m_i m_j}{(d^2 + \epsilon^2)^{3/2}} (\mathbf{r}_j - \mathbf{r}_i)$$

where $d^2 = \|\mathbf{r}_j - \mathbf{r}_i\|^2$.

The corresponding potential energy is softened as:
$$\Phi(r) = - \frac{G M}{\sqrt{r^2 + \epsilon^2}}$$

---

## 6. Geometric Radius Neighbor Search with Trees

Beyond force evaluation, Quadtrees and Octrees can also be utilized for geometric radius queries (e.g., finding all particles within distance $R$ of point $\mathbf{p}$):

### Pruning Rule (AABB vs Sphere):
For a node with bounding box $B = [\mathbf{x}_{\min}, \mathbf{x}_{\max}]$:
1. Compute the minimum distance squared from query point $\mathbf{p}$ to the box $B$:
   $$d_{\min}^2 = \sum_{k \in \{x, y, z\}} \begin{cases} (x_{\min, k} - p_k)^2, & p_k < x_{\min, k} \\ (p_k - x_{\max, k})^2, & p_k > x_{\max, k} \\ 0, & \text{otherwise} \end{cases}$$
2. **If $d_{\min}^2 > R^2$**, the sphere does not intersect box $B$ $\implies$ **Prune subtree completely**.
3. **Otherwise**, traverse children or check leaf particles.

---

## 7. Published Sources & Citations

1. **Barnes, J., & Hut, P.** (1986). *A hierarchical $O(N \log N)$ force-calculation algorithm*. **Nature**, 324(6096), 446–449.
   - DOI: [10.1038/324446a0](https://doi.org/10.1038/324446a0)
   - *The seminal paper that introduced tree-code hierarchical n-body approximation.*
2. **Plummer, H. C.** (1911). *On the problem of distribution in globular star clusters: of the particular cases of spherical symmetry, and the general case of ellipsoidal symmetry*. **Monthly Notices of the Royal Astronomical Society**, 71(5), 460–470.
   - DOI: [10.1093/mnras/71.5.460](https://doi.org/10.1093/mnras/71.5.460)
3. **Dehnen, W.** (2000). *A Very Fast and Accurate Barnes-Hut Algorithm*. **The Astrophysical Journal**, 536(1), L39–L42.
   - DOI: [10.1086/312724](https://doi.org/10.1086/312724)
4. **Salmon, J. K., & Warren, M. S.** (1994). *Skeletons from the Treecode Closet: A Performance Analysis of the Barnes-Hut Algorithm*. **Journal of Computational Physics**, 111(1), 136–155.
   - DOI: [10.1006/jcph.1994.1050](https://doi.org/10.1006/jcph.1994.1050)
