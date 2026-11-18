# Spatial Hashing in Particle Simulations

## 1. Overview and Problem Statement

In discrete particle simulations (such as Smoothed Particle Hydrodynamics (SPH), the Discrete Element Method (DEM), and molecular dynamics), each particle interacts only with neighbors lying within a finite cutoff radius $r_{\text{cut}}$.

A naive all-pairs distance check requires evaluating:
$$N_{\text{pairs}} = \frac{N(N - 1)}{2} \implies \mathcal{O}(N^2)$$
evaluations, which rapidly becomes prohibitive for $N > 10^3$.

**Spatial hashing** transforms this $\mathcal{O}(N^2)$ bottleneck into an expected $\mathcal{O}(N)$ operation by discretizing continuous physical space into regular hyper-cubic cells of size $l \ge r_{\text{cut}}$ and indexing them into a finite hash table.

---

## 2. Coordinate Discretization

Given a particle position $\mathbf{x} = (x, y, z) \in \mathbb{R}^3$, the continuous coordinate is mapped to discrete integer cell indices $\mathbf{c} = (c_x, c_y, c_z) \in \mathbb{Z}^3$:

$$c_x = \left\lfloor \frac{x - x_{\min}}{l} \right\rfloor, \quad c_y = \left\lfloor \frac{y - y_{\min}}{l} \right\rfloor, \quad c_z = \left\lfloor \frac{z - z_{\min}}{l} \right\rfloor$$

where:
- $l$ is the cell edge length (typically chosen as $l = r_{\text{cut}}$ or $l = 2 r_{\text{particle}}$).
- $\mathbf{x}_{\min} = (x_{\min}, y_{\min}, z_{\min})$ is an optional domain origin offset.

---

## 3. The Teschner et al. (2003) Prime Hash Function

To support unbounded or large sparse domains without allocating a massive 3D grid array, 3D integer coordinates are mapped to a 1D hash table index using the XOR-prime hash formulation introduced by Teschner et al. (2003):

$$\mathcal{H}(c_x, c_y, c_z) = \Big( (c_x \cdot p_1) \oplus (c_y \cdot p_2) \oplus (c_z \cdot p_3) \Big) \pmod N$$

For 2D domains:
$$\mathcal{H}(c_x, c_y) = \Big( (c_x \cdot p_1) \oplus (c_y \cdot p_2) \Big) \pmod N$$

### Standard Prime Constants
The large prime numbers recommended by Teschner et al. for minimizing hash collisions across Cartesian integer lattices are:
$$p_1 = 73856093$$
$$p_2 = 19349663$$
$$p_3 = 83492791$$

Here:
- $\oplus$ denotes the bitwise XOR operation.
- $N$ is the hash table capacity, ideally chosen as a prime number (e.g., $N \approx 2 \cdot N_{\text{particles}}$) to avoid harmonic patterns with grid coordinates.
- Non-negative modulo: for negative coordinates $c_i < 0$, the mathematical modulo $(a \pmod N + N) \pmod N$ is used to guarantee non-negative indices in $[0, N - 1]$.

---

## 4. Neighbor Search Algorithm (Query Stencils)

Because $l \ge r_{\text{cut}}$, any particle that is within Euclidean distance $\|\mathbf{x}_j - \mathbf{x}_i\| \le r_{\text{cut}}$ from $\mathbf{x}_i$ must reside either in the exact same cell as particle $i$ or in one of the immediately adjacent neighboring cells.

### Query Stencils:
- **In 2D (9 cells)**:
  $$\Delta \mathbf{c} \in \{-1, 0, 1\} \times \{-1, 0, 1\}$$
- **In 3D (27 cells)**:
  $$\Delta \mathbf{c} \in \{-1, 0, 1\} \times \{-1, 0, 1\} \times \{-1, 0, 1\}$$

```
2D Stencil (9 cells)         3D Stencil (27 cells)
+------+------+------+        3x3x3 grid centered
| -1,1 | 0,1  | 1,1  |        around target cell (0,0,0)
+------+------+------+
| -1,0 | target 1,0  |
+------+------+------+
| -1,-1| 0,-1 | 1,-1 |
+------+------+------+
```

### Algorithm Steps:
1. **Clear Grid**: Reset bucket headers without deallocating underlying arrays.
2. **Build / Insert**: For each particle $i = 0 \dots N - 1$:
   - Calculate $\mathbf{c}_i = \lfloor (\mathbf{x}_i - \mathbf{x}_{\min}) / l \rfloor$.
   - Compute hash $h = \mathcal{H}(\mathbf{c}_i)$.
   - Store particle index $i$ in bucket $h$.
3. **Neighbor Query**: For a query point $\mathbf{x}_{\text{query}}$ with radius $r_{\text{cut}}$:
   - Calculate target cell $\mathbf{c}_{\text{query}}$.
   - For each offset in the stencil:
     - Compute cell coordinate $\mathbf{c}' = \mathbf{c}_{\text{query}} + \Delta \mathbf{c}$.
     - Look up candidates in bucket $\mathcal{H}(\mathbf{c}')$.
     - Filter candidates by verifying:
       $$\|\mathbf{x}_j - \mathbf{x}_{\text{query}}\|^2 \le r_{\text{cut}}^2$$
       (Using squared Euclidean distance to eliminate expensive square root calculations).

---

## 5. Performance and Complexity

| Operation | Naive All-Pairs | Spatial Hashing (Expected) | Spatial Hashing (Worst-Case) |
|---|---|---|---|
| Grid Construction | N/A | $\mathcal{O}(N)$ | $\mathcal{O}(N)$ |
| Single Radius Query | $\mathcal{O}(N)$ | $\mathcal{O}(k)$ ($k$: local neighbors) | $\mathcal{O}(N)$ (all hash to same cell) |
| All-Pairs Neighbor Search | $\mathcal{O}(N^2)$ | $\mathcal{O}(N)$ | $\mathcal{O}(N^2)$ |
| Memory Overhead | $\mathcal{O}(1)$ | $\mathcal{O}(N + \text{TableSize})$ | $\mathcal{O}(N + \text{TableSize})$ |

---

## 6. Published Sources & Citations

1. **Teschner, M., Heidelberger, B., Müller, M., Pomeranets, D., & Gross, M.** (2003). *Optimized Spatial Hashing for Collision Detection of Deformable Objects*. In **Proceedings of Vision, Modeling, Visualization (VMV '03)**, Munich, Germany, pp. 47–54.
   - [Paper Link (ETH Zürich)](https://graphics.ethz.ch/publications/papers/paperTes03.php)
2. **Green, S.** (2008). *Particle Simulation using CUDA*. NVIDIA Whitepaper.
   - Describes compact cell-linked lists and sorted key hash techniques for multi-core GPU/CPU architectures.
3. **Ihmsen, M., Akinci, N., Becker, M., & Teschner, M.** (2011). *A Parallel SPH Implementation on Multi-Core CPUs*. **Computer Graphics Forum**, 30(1), 99–112.
   - DOI: [10.1111/j.1467-8659.2010.01832.x](https://doi.org/10.1111/j.1467-8659.2010.01832.x)
