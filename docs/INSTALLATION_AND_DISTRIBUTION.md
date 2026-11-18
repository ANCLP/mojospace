# Installation & Distribution Guide

This document describes how to use and consume `mojospace` in your own projects without manually cloning the repository, as well as instructions for publishing the package to [prefix.dev](https://prefix.dev).

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

## 4. Publishing to `prefix.dev` Channel

`mojospace` is configured as a conda package using `pixi-build-mojo`.

### Automated Release via GitHub Actions (Secretless OIDC Workload Identity)
By linking your GitHub repository to your prefix.dev channel with Write access, releases publish automatically and securely without managing any static API tokens:

1. **Link Repository on prefix.dev**:
   - In your prefix.dev channel settings, connect the GitHub repository `anclp/mojospace` with **Write** permission.
2. **Automated Trigger**:
   - Creating a GitHub release or pushing a tag (`git tag v0.1.0 && git push origin v0.1.0`) triggers `.github/workflows/release.yml`.
   - The workflow uses `prefix-dev/setup-pixi` with `auth-workload-identity: true` to exchange a short-lived GitHub OIDC token with prefix.dev.
   - It precompiles `mojospace.mojoc`, generates Sigstore cryptographic attestations, and publishes `mojospace-*.conda` to `https://prefix.dev/anclp/mojospace`.

### Manual / Local Builds
> [!NOTE]
> If building packages locally on a Windows/WSL path containing spaces (e.g. `Startup Projects`), upstream `pixi-build-mojo` has an unquoted variable in its output path template. For local builds, run `pixi run package` (which outputs `build/mojospace.mojoc` directly), or test conda publishing inside a directory without spaces (e.g. `/tmp/mojospace`). GitHub Actions CI/CD runs in `/home/runner/work/...` which contains no spaces and is unaffected.

To test publishing locally:
```bash
pixi publish --dry-run --to https://prefix.dev/anclp/mojospace
```

### How Users Install from Your Channel:
Once published to prefix.dev, any user can install `mojospace` into their Pixi project:
```bash
pixi add --channel https://prefix.dev/anclp mojospace
```
Or by adding `"https://prefix.dev/anclp"` to their `channels` list in `pixi.toml`:
```toml
[workspace]
channels = [
    "https://prefix.dev/anclp",
    "https://conda.modular.com/max",
    "conda-forge"
]

[dependencies]
mojospace = "*"
```
