from std.math import sqrt
from src.math.vec2 import Vec2
from src.math.vec3 import Vec3
from src.particle.particle2d import Particle2D
from src.particle.particle3d import Particle3D

# ---------------------------------------------------------
# Gravitational Forces (Newtonian with Plummer Softening)
# ---------------------------------------------------------

def gravitational_force_direct_2d(
    p1: Particle2D,
    p2: Particle2D,
    G: Float64 = 1.0,
    softening: Float64 = 0.01
) -> Vec2:
    """Computes exact direct gravitational force exerted by p2 on p1 in 2D."""
    var diff = p2.pos - p1.pos
    var r2 = diff.length_squared() + softening * softening
    var r = sqrt(r2)
    var mag = G * (p1.mass * p2.mass) / (r2 * r)
    return diff * (mag / r)

def gravitational_force_direct_3d(
    p1: Particle3D,
    p2: Particle3D,
    G: Float64 = 1.0,
    softening: Float64 = 0.01
) -> Vec3:
    """Computes exact direct gravitational force exerted by p2 on p1 in 3D."""
    var diff = p2.pos - p1.pos
    var r2 = diff.length_squared() + softening * softening
    var r = sqrt(r2)
    var mag = G * (p1.mass * p2.mass) / (r2 * r)
    return diff * (mag / r)

# ---------------------------------------------------------
# Coulomb Electrostatic Force (with Softening)
# ---------------------------------------------------------

def coulomb_force_3d(
    p1: Particle3D,
    p2: Particle3D,
    k_e: Float64 = 1.0,
    softening: Float64 = 0.01
) -> Vec3:
    """Computes electrostatic force exerted on p1 by p2 (repulsive for same signs)."""
    var diff = p1.pos - p2.pos
    var r2 = diff.length_squared() + softening * softening
    var r = sqrt(r2)
    var mag = k_e * (p1.charge * p2.charge) / (r2 * r)
    return diff * (mag / r)

# ---------------------------------------------------------
# Lennard-Jones (12-6) Potential & Force
# ---------------------------------------------------------

def lennard_jones_force_3d(
    r_vec: Vec3,
    epsilon: Float64 = 1.0,
    sigma: Float64 = 1.0
) -> Vec3:
    """
    Computes Lennard-Jones (12-6) interatomic force:
    F(r) = 24 * eps / r^2 * (2 * (sig/r)^12 - (sig/r)^6) * r_vec
    """
    var r2 = r_vec.length_squared()
    if r2 == 0.0:
        return Vec3(0.0, 0.0, 0.0)

    var sig2 = sigma * sigma
    var inv_r2 = 1.0 / r2
    var s = (sig2 * inv_r2) * (sig2 * inv_r2) * (sig2 * inv_r2) # (sigma / r)^6
    var s2 = s * s                                              # (sigma / r)^12
    var factor = 24.0 * epsilon * inv_r2 * (2.0 * s2 - s)
    return r_vec * factor

# ---------------------------------------------------------
# Smoothed Particle Hydrodynamics (SPH) Smoothing Kernels
# ---------------------------------------------------------

def sph_cubic_spline_kernel_2d(r: Float64, h: Float64) -> Float64:
    """
    Monaghan (1992) Cubic Spline Kernel in 2D.
    sigma_2 = 10 / (7 * pi * h^2)
    """
    var pi: Float64 = 3.141592653589793
    var sigma: Float64 = 10.0 / (7.0 * pi * h * h)
    var q: Float64 = r / h

    if q < 0.0 or q > 2.0:
        return 0.0
    elif q <= 1.0:
        return sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    else:
        var term = 2.0 - q
        return sigma * 0.25 * term * term * term

def sph_cubic_spline_kernel_3d(r: Float64, h: Float64) -> Float64:
    """
    Monaghan (1992) Cubic Spline Kernel in 3D.
    sigma_3 = 1 / (pi * h^3)
    """
    var pi: Float64 = 3.141592653589793
    var sigma: Float64 = 1.0 / (pi * h * h * h)
    var q: Float64 = r / h

    if q < 0.0 or q > 2.0:
        return 0.0
    elif q <= 1.0:
        return sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    else:
        var term = 2.0 - q
        return sigma * 0.25 * term * term * term

def sph_spiky_gradient_3d(r_vec: Vec3, h: Float64) -> Vec3:
    """
    Müller et al. (2003) Spiky Gradient Kernel for SPH fluid pressure forces.
    grad W = -45 / (pi * h^6) * (h - r)^2 * (r_vec / r)
    """
    var r = r_vec.length()
    if r <= 0.0 or r >= h:
        return Vec3(0.0, 0.0, 0.0)

    var pi: Float64 = 3.141592653589793
    var h6 = h * h * h * h * h * h
    var diff = h - r
    var factor = -45.0 / (pi * h6) * (diff * diff) / r
    return r_vec * factor
