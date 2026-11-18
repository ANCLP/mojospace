# mojospace API Reference

A pure Mojo library for spatial partitioning, neighbor searches, and particle physics acceleration.

---

## 1. Module `math`

### `struct Vec2`
2D Vector with high-performance arithmetic and vector geometry.

```mojo
from math.vec2 import Vec2
```

| Method | Signature | Description |
|---|---|---|
| `__init__` | `def __init__(out self, x: Float64 = 0.0, y: Float64 = 0.0)` | Initializes vector coordinates. |
| `+`, `-` | `def +(self, other: Vec2) -> Vec2` | Element-wise addition and subtraction. |
| `*`, `/` | `def *(self, scalar: Float64) -> Vec2` | Scalar multiplication and division. |
| `dot` | `def dot(self, other: Vec2) -> Float64` | Computes $x_1 x_2 + y_1 y_2$. |
| `length_squared` | `def length_squared(self) -> Float64` | Returns $\|\mathbf{v}\|^2 = x^2 + y^2$. |
| `length` | `def length(self) -> Float64` | Returns Euclidean length $\|\mathbf{v}\|$. |
| `normalized` | `def normalized(self) -> Vec2` | Returns unit vector; zero vector if length is zero. |
| `distance_to` | `def distance_to(self, other: Vec2) -> Float64` | Euclidean distance between two vectors. |
| `distance_squared_to` | `def distance_squared_to(self, other: Vec2) -> Float64` | Squared Euclidean distance. |

---

### `struct Vec3`
3D Vector with geometry and vector products.

```mojo
from math.vec3 import Vec3
```

| Method | Signature | Description |
|---|---|---|
| `__init__` | `def __init__(out self, x: Float64 = 0.0, y: Float64 = 0.0, z: Float64 = 0.0)` | Initializes 3D coordinates. |
| `+`, `-`, `*`, `/` | Operator overloads | Vector arithmetic and scalar scaling. |
| `dot` | `def dot(self, other: Vec3) -> Float64` | Computes $x_1 x_2 + y_1 y_2 + z_1 z_2$. |
| `cross` | `def cross(self, other: Vec3) -> Vec3` | Vector cross product $\mathbf{u} \times \mathbf{v}$. |
| `length_squared`, `length` | `def length_squared(self) -> Float64` | Squared norm and Euclidean norm. |
| `normalized` | `def normalized(self) -> Vec3` | Returns unit vector along $\mathbf{v}$. |
| `distance_to` | `def distance_to(self, other: Vec3) -> Float64` | Euclidean distance to target. |

---

### `struct AABB2D` and `struct AABB3D`
Axis-Aligned Bounding Box for spatial pruning.

| Method | Signature | Description |
|---|---|---|
| `contains` | `def contains(self, p: Vec) -> Bool` | Checks if point lies inside box bounds. |
| `intersects` | `def intersects(self, other: AABB) -> Bool` | Checks if two bounding boxes overlap. |
| `distance_squared_to` | `def distance_squared_to(self, p: Vec) -> Float64` | Squared minimum distance from point to box surface. |
| `expand` | `def expand(inout self, p: Vec)` | Encloses point into the bounding box. |

---

## 2. Module `particle`

### `struct Particle2D` and `struct Particle3D`
Standard particle state containers.

Fields:
- `id: Int`
- `pos: Vec` (position)
- `vel: Vec` (velocity)
- `acc: Vec` (acceleration)
- `mass: Float64` (inertial/gravitational mass)
- `radius: Float64` (collision / particle radius)
- `charge: Float64` (electrostatic charge)

---

## 3. Module `spatial_hash`

### `struct SpatialHashGrid2D` and `struct SpatialHashGrid3D`
Constant-time expected neighbor lookup grid using Teschner et al. (2003) prime hashing.

```mojo
from spatial_hash.hash_grid_3d import SpatialHashGrid3D

var grid = SpatialHashGrid3D(cell_size=0.5, table_size=10007)
```

| Method | Signature | Complexity | Description |
|---|---|---|---|
| `clear` | `def clear(inout self)` | $\mathcal{O}(\text{TableSize})$ | Clears buckets without reallocating arrays. |
| `insert` | `def insert(inout self, particle_id: Int, pos: Vec)` | $\mathcal{O}(1)$ | Hashes particle position and appends ID. |
| `build` | `def build(inout self, particles: List[Particle])` | $\mathcal{O}(N)$ | Inserts all particles from a list into grid. |
| `query_radius` | `def query_radius(self, center: Vec, radius: Float64, particles: List[Particle]) -> List[Int]` | $\mathcal{O}(1)$ exp. | Returns all particle IDs within radius $r$. |
| `query_pairs` | `def query_pairs(self, particles: List[Particle], radius: Float64) -> List[Tuple[Int, Int]]` | $\mathcal{O}(N)$ exp. | Returns unique collision pairs $(i, j)$ with $i < j$. |

---

## 4. Module `trees`

### `struct Quadtree` (2D) & `struct BarnesHutOctree` (3D)
Hierarchical spatial tree for $\mathcal{O}(N \log N)$ n-body force approximation and adaptive radius lookups.

```mojo
from trees.octree import BarnesHutOctree

var octree = BarnesHutOctree(bounds=box, max_depth=16)
octree.build(particles)
```

| Method | Signature | Complexity | Description |
|---|---|---|---|
| `build` | `def build(inout self, particles: List[Particle])` | $\mathcal{O}(N \log N)$ | Subdivides tree and aggregates center of mass. |
| `compute_force_barnes_hut` | `def compute_force_barnes_hut(self, target_idx: Int, particles: List[Particle], theta: Float64, G: Float64, softening: Float64) -> Vec` | $\mathcal{O}(\log N)$ | Computes gravitational force on particle using MAC $\theta$. |
| `query_radius` | `def query_radius(self, center: Vec, radius: Float64, particles: List[Particle]) -> List[Int]` | $\mathcal{O}(\log N + k)$ | Prunes subtrees using AABB distance metric. |

---

## 5. Module `physics`

### Forces & Kernels (`physics.force`)
- `def gravitational_force_direct(p1: Particle, p2: Particle, G: Float64, softening: Float64) -> Vec`: Direct pairwise Newtonian force with Plummer softening.
- `def coulomb_force(p1: Particle, p2: Particle, k_e: Float64, softening: Float64) -> Vec`: Coulomb electrostatic interaction.
- `def lennard_jones_force(r_vec: Vec, epsilon: Float64, sigma: Float64) -> Vec`: Lennard-Jones 12-6 force vector.
- `def sph_cubic_spline_kernel_3d(r: Float64, h: Float64) -> Float64`: Monaghan (1992) cubic spline kernel value.
- `def sph_spiky_gradient_3d(r_vec: Vec, h: Float64) -> Vec`: Müller et al. (2003) pressure gradient vector.

### Time Integrators (`physics.integrator`)
- `def integrate_velocity_verlet_step(inout particles: List[Particle], dt: Float64, forces: List[Vec])`: Second-order symplectic position and velocity advancement.
- `def integrate_symplectic_euler_step(inout particles: List[Particle], dt: Float64, forces: List[Vec])`: First-order symplectic momentum and coordinate update.
