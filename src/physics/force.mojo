from math import sqrt
from src.math.vec2 import Vec2
from src.math.vec3 import Vec3
from src.particle.particle2d import Particle2D
from src.particle.particle3d import Particle3D

# ---------------------------------------------------------
# Gravitational Forces (Newtonian with Plummer Softening)
# ---------------------------------------------------------

fn gravitational_force_direct_2d(
    p1: Particle2D,
    p2: Particle2D,
    G: Float64 = 1.0,
    softening: Float64 = 0.01
) -> Vec2:
    """Computes exact direct gravitational force exerted by p2 on p1 in 2D."""
    let diff = p2.pos - p1.pos
    let r2 = diff.length_squared() + softening * softening
    let r = sqrt(r2)
    let mag = G * (p1.mass * p2.mass) / (r2 * r)
    return diff * (mag / r)

fn gravitational_force_direct_3d(
    p1: Particle3D,
    p2: Particle3D,
    G: Float64 = 1.0,
    softening: Float64 = 0.01
) -> Vec3:
    """Computes exact direct gravitational force exerted by p2 on p1 in 3D."""
    let diff = p2.pos - p1.pos
    let r2 = diff.length_squared() + softening * softening
    let r = sqrt(r2)
    let mag = G * (p1.mass * p2.mass) / (r2 * r)
    return diff * (mag / r)

# ---------------------------------------------------------
# Coulomb Electrostatic Force (with Softening)
# ---------------------------------------------------------

fn coulomb_force_3d(
    p1: Particle3D,
    p2: Particle3D,
    k_e: Float64 = 1.0,
    softening: Float64 = 0.01
) -> Vec3:
    """Computes electrostatic force exerted on p1 by p2 (repulsive for same signs)."""
    let diff = p1.pos - p2.pos
    let r2 = diff.length_squared() + softening * softening
    let r = sqrt(r2)
    let mag = k_e * (p1.charge * p2.charge) / (r2 * r)
    return diff * (mag / r)

# ---------------------------------------------------------
# Lennard-Jones (12-6) Potential & Force
# ---------------------------------------------------------

fn lennard_jones_force_3d(
    r_vec: Vec3,
    epsilon: Float64 = 1.0,
    sigma: Float64 = 1.0
) -> Vec3:
    """
    Computes Lennard-Jones (12-6) interatomic force:
    F(r) = 24 * eps / r^2 * (2 * (sig/r)^12 - (sig/r)^6) * r_vec
    """
    let r2 = r_vec.length_squared()
    if r2 == 0.0:
        return Vec3(0.0, 0.0, 0.0)

    let sig2 = sigma * sigma
    let inv_r2 = 1.0 / r2
    let s = (sig2 * inv_r2) * (sig2 * inv_r2) * (sig2 * inv_r2) # (sigma / r)^6
    let s2 = s * s                                              # (sigma / r)^12
    let factor = 24.0 * epsilon * inv_r2 * (2.0 * s2 - s)
    return r_vec * factor

# ---------------------------------------------------------
# Smoothed Particle Hydrodynamics (SPH) Smoothing Kernels
# ---------------------------------------------------------

fn sph_cubic_spline_kernel_2d(r: Float64, h: Float64) -> Float64:
    """
    Monaghan (1992) Cubic Spline Kernel in 2D.
    sigma_2 = 10 / (7 * pi * h^2)
    """
    let pi: Float64 = 3.141592653589793
    let sigma: Float64 = 10.0 / (7.0 * pi * h * h)
    let q: Float64 = r / h

    if q < 0.0 or q > 2.0:
        return 0.0
    elif q <= 1.0:
        return sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    else:
        let term = 2.0 - q
        return sigma * 0.25 * term * term * term

fn sph_cubic_spline_kernel_3d(r: Float64, h: Float64) -> Float64:
    """
    Monaghan (1992) Cubic Spline Kernel in 3D.
    sigma_3 = 1 / (pi * h^3)
    """
    let pi: Float64 = 3.141592653589793
    let sigma: Float64 = 1.0 / (pi * h * h * h)
    let q: Float64 = r / h

    if q < 0.0 or q > 2.0:
        return 0.0
    elif q <= 1.0:
        return sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    else:
        let term = 2.0 - q
        return sigma * 0.25 * term * term * term

fn sph_spiky_gradient_3d(r_vec: Vec3, h: Float64) -> Vec3:
    """
    Müller et al. (2003) Spiky Gradient Kernel for SPH fluid pressure forces.
    grad W = -45 / (pi * h^6) * (h - r)^2 * (r_vec / r)
    """
    let r = r_vec.length()
    if r <= 0.0 or r >= h:
        return Vec3(0.0, 0.0, 0.0)

    let pi: Float64 = 3.141592653589793
    let h6 = h * h * h * h * h * h
    let diff = h - r
    let factor = -45.0 / (pi * h6) * (diff * diff) / r
    return r_vec * factor
