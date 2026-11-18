# mojospace 🌌

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Pure Mojo](https://img.shields.io/badge/Language-Pure%20Mojo-ff69b4.svg)](https://modular.com/mojo)
[![Platform: Linux | macOS | WSL2](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS%20%7C%20WSL2-blue.svg)](https://modular.com)

A high-performance, lightweight 2D/3D spatial partitioning and neighbor search library written in **pure Mojo**. 

Designed for particle physics engines, Smoothed Particle Hydrodynamics (SPH), Discrete Element Method (DEM), collision detection, and astrophysical N-body gravitational simulations.

---

## ⚡ Features

- **🚀 Pure Mojo Implementation**: Zero external C/Python dependencies, leveraging Mojo's memory safety, compile-time performance, and SIMD hardware acceleration.
- **🔍 Spatial Hashing (2D & 3D)**: Constant-time $\mathcal{O}(1)$ expected radius neighbor lookup using the seminal [Teschner et al. (2003)](#references) XOR-prime hash formulation.
- **🌳 Barnes-Hut Trees (Quadtree & Octree)**: Hierarchical $\mathcal{O}(N \log N)$ n-body force calculation with Multipole Acceptance Criteria ($\theta$) and Plummer softening.
- **🌊 SPH & Physics Kernels**: Built-in Monaghan cubic spline smoothing kernels, Müller spiky pressure gradients, Lennard-Jones (12-6) molecular potentials, and symplectic time integrators (Velocity Verlet & Symplectic Euler).
- **📚 Rigorous Mathematical Documentation**: Complete mathematical formulas, derivations, complexity analysis, and academic citations in [`docs/`](docs/).

---

## 📐 Architecture & Math Overview

```
                      Physical Domain Space
                               │
         ┌─────────────────────┴─────────────────────┐
         ▼                                           ▼
  Spatial Hashing                             Barnes-Hut Trees
  (Fixed-cutoff, SPH, Contact)                (Long-range Gravity, Coulomb)
         │                                           │
  c = floor((x - min) / l)                    Quadtree (2D) / Octree (3D)
  H(c) = (c_x*p1 ^ c_y*p2 ^ c_z*p3) % N       Hierarchical Center of Mass (R_c, M_c)
  Query Stencil: 3^d neighboring cells        MAC Opening Angle: s / d < theta
         │                                           │
  Complexity: O(1) expected                   Complexity: O(N log N)
```

| Structure | Best Suited For | Target Interactions | Complexity |
|---|---|---|---|
| **Spatial Hash Grid** | Fluid (SPH), Contact, DEM, Collision | Short-range ($r \le r_{\text{cut}}$) | Construction: $\mathcal{O}(N)$<br>Query: $\mathcal{O}(1)$ exp. |
| **Barnes-Hut Tree** | Astrophysics, Celestial N-body, Coulomb | Long-range ($1/r^2$) | Construction: $\mathcal{O}(N \log N)$<br>Force: $\mathcal{O}(N \log N)$ |

---

## 🚀 Quickstart

### 1. 3D Spatial Hashing Neighbor Lookup

```mojo
from src.math.vec3 import Vec3
from src.particle.particle3d import Particle3D
from src.spatial_hash.hash_grid_3d import SpatialHashGrid3D
from collections import List

def main():
    # 1. Initialize grid with cell size matching cutoff radius
    var grid = SpatialHashGrid3D(cell_size=1.0, table_size=10007)

    # 2. Populate particles
    var particles = List[Particle3D]()
    particles.append(Particle3D(id=0, pos=Vec3(0.0, 0.0, 0.0)))
    particles.append(Particle3D(id=1, pos=Vec3(0.5, 0.2, 0.1)))
    particles.append(Particle3D(id=2, pos=Vec3(50.0, 50.0, 50.0)))

    # 3. Build spatial hash table
    grid.build(particles)

    # 4. Query all neighbors within radius 1.0 from particle 0
    var neighbors = grid.query_radius(particles[0].pos, 1.0, particles)
    print("Found neighbors count:", len(neighbors))
    for i in range(len(neighbors)):
        print("Neighbor particle ID:", neighbors[i])
```

### 2. Barnes-Hut 3D Gravitational Force Calculation

```mojo
from src.math.vec3 import Vec3
from src.math.bounds import AABB3D
from src.particle.particle3d import Particle3D
from src.trees.octree import BarnesHutOctree
from collections import List

def main():
    var bounds = AABB3D(Vec3(-100.0, -100.0, -100.0), Vec3(100.0, 100.0, 100.0))
    var octree = BarnesHutOctree(bounds=bounds, max_depth=10)

    var particles = List[Particle3D]()
    # ... populate particles ...
    octree.build(particles)

    # Calculate force on target particle using MAC theta = 0.5
    var force = octree.compute_force_barnes_hut(
        target_idx=0,
        particles=particles,
        theta=0.5,
        G=1.0,
        softening=0.01
    )
    print("Gravitational force:", force.x, force.y, force.z)
```

---

## 📖 In-Depth Mathematical Documentation

Check out the detailed scientific documentation in the [`docs/`](docs/) directory:
- [**Spatial Hashing Formulations**](docs/SPATIAL_HASHING.md): Coordinate mapping, Teschner prime hash constants, cell size tuning, collision behavior, and complexity proofs.
- [**Barnes-Hut & Hierarchical Trees**](docs/BARNES_HUT_AND_TREES.md): Quadtree & Octree subdivision, center of mass aggregation, Multipole Acceptance Criteria ($\theta$), Plummer softening, and tree pruning.
- [**Particle Physics & Smoothing Kernels**](docs/PARTICLE_PHYSICS_AND_KERNELS.md): Monaghan cubic spline SPH kernel, Müller spiky pressure gradient, Lennard-Jones 12-6 potential, and Velocity Verlet integration.
- [**API Reference**](docs/API_REFERENCE.md): Comprehensive signature and complexity reference for all structs and functions.

---

## 🛠️ Project Setup & Development

This project uses [Pixi](https://pixi.sh) to manage the Mojo environment across platforms.

### Installation & Running:
```bash
# Clone the repository
git clone https://github.com/ANCLP/mojospace.git
cd mojospace

# Run test suite
pixi run test

# Run neighbor search benchmark
pixi run bench

# Run N-body galactic simulation demo
pixi run demo-nbody
```

*(On Windows machines without native Mojo binaries, run the commands through WSL2).*

---

## 🔬 Citations & References

If you use `mojospace` in academic research, please cite:

```bibtex
@software{mojospace2026,
  author = {mojospace Contributors},
  title = {mojospace: Lightweight 2D/3D Spatial Hashing and Barnes-Hut Trees in Pure Mojo},
  year = {2026},
  url = {https://github.com/ANCLP/mojospace}
}
```

Key published references:
1. **Teschner, M., et al.** (2003). *Optimized Spatial Hashing for Collision Detection of Deformable Objects*. VMV 2003, pp. 47–54.
2. **Barnes, J., & Hut, P.** (1986). *A hierarchical $O(N \log N)$ force-calculation algorithm*. Nature, 324(6096), 446–449. [doi:10.1038/324446a0](https://doi.org/10.1038/324446a0)
3. **Monaghan, J. J.** (1992). *Smoothed Particle Hydrodynamics*. Ann. Rev. Astron. Astrophys., 30(1), 543–574. [doi:10.1146/annurev.aa.30.090192.002551](https://doi.org/10.1146/annurev.aa.30.090192.002551)
4. **Verlet, L.** (1967). *Computer "Experiments" on Classical Fluids. I.* Phys. Rev., 159(1), 98–103. [doi:10.1103/PhysRev.159.98](https://doi.org/10.1103/PhysRev.159.98)

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).
