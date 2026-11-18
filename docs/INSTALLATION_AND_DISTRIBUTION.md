# Installation & Distribution Guide

This document describes how to use and consume `mojospace` in your own projects without manually cloning the repository.

---

## 1. Zero-Clone Installation via Pixi (Recommended)

If you are developing a Mojo application or physics simulation using [Pixi](https://pixi.sh), you can add `mojospace` directly as a Git dependency in your project's `pixi.toml`:

```toml
[workspace]
channels = [
    "https://conda.modular.com/max",
    "https://repo.prefix.dev/modular-community",
    "conda-forge"
]
platforms = ["linux-64", "osx-arm64"]

[dependencies]
mojo = "*"
mojospace = { git = "https://github.com/anclp/mojospace" }
```

When you run `pixi install` or `pixi run`:
- Pixi automatically fetches and compiles `mojospace` in the background.
- Pixi places the compiled package into your environment's `$CONDA_PREFIX/lib/mojo/` directory.
- You can immediately import and use the library in your code:

```mojo
from mojospace.spatial_hash import SpatialHashGrid3D
from mojospace.trees import BarnesHutOctree
from mojospace.math import Vec3
```

---

## 2. Using Precompiled Binary Releases (`.mojoc`)

Every tagged GitHub release (e.g. `v0.1.0`) automatically compiles and publishes a standalone precompiled binary package:
**`mojospace.mojoc`** (~125 KB).

### How to use:
1. Go to the [Releases](https://github.com/anclp/mojospace/releases) page.
2. Download `mojospace.mojoc`.
3. Place `mojospace.mojoc` in your project root, or in your include path (`mojo -I /path/to/folder`).
4. Import directly in your Mojo scripts:

```mojo
import mojospace
from mojospace.spatial_hash.hash_grid_3d import SpatialHashGrid3D
from mojospace.math.vec3 import Vec3

def main():
    var grid = SpatialHashGrid3D(cell_size=1.0)
    print("mojospace initialized successfully.")
```

---

## 3. Building the Precompiled Binary Locally

To precompile `mojospace.mojoc` locally from source:

```bash
pixi run package
```

The output binary is written to `build/mojospace.mojoc`.

---

## 4. Future: Central Index Installation (`prefix.dev/modular-community`)

When `mojospace` is published to Prefix.dev, installation will be a single CLI command:

```bash
pixi add --channel https://prefix.dev/modular-community mojospace
```
