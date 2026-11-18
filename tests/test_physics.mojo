from collections import List
from src.math.vec3 import Vec3
from src.particle.particle3d import Particle3D
from src.physics.force import (
    sph_cubic_spline_kernel_3d,
    sph_spiky_gradient_3d,
    lennard_jones_force_3d
)
from src.physics.integrator import integrate_velocity_verlet_step_3d

def assert_true(condition: Bool, msg: String):
    if not condition:
        print("[FAIL] " + msg)
    else:
        print("[PASS] " + msg)

def main():
    print("=== Testing Physics Module ===")

    # 1. SPH Cubic Spline Kernel tests
    var h: Float64 = 1.0
    var w_0 = sph_cubic_spline_kernel_3d(0.0, h)
    var w_1 = sph_cubic_spline_kernel_3d(0.5, h)
    var w_2 = sph_cubic_spline_kernel_3d(1.5, h)
    var w_cut = sph_cubic_spline_kernel_3d(2.1, h)

    assert_true(w_0 > w_1 and w_1 > w_2, "SPH Cubic spline kernel strictly decreases with distance")
    assert_true(w_cut == 0.0, "SPH Cubic spline kernel is strictly zero beyond 2*h")

    # 2. SPH Spiky Pressure Gradient tests
    var grad_cut = sph_spiky_gradient_3d(Vec3(1.2, 0.0, 0.0), h=1.0)
    assert_true(grad_cut == Vec3(0.0, 0.0, 0.0), "SPH Spiky gradient is zero beyond smoothing length h")

    var grad_inside = sph_spiky_gradient_3d(Vec3(0.5, 0.0, 0.0), h=1.0)
    assert_true(grad_inside.x < 0.0, "SPH Spiky gradient points in opposing direction to relative vector")

    # 3. Lennard-Jones (12-6) Force test
    # Equilibrium distance is r_eq = 2^(1/6) * sigma approx 1.122462
    var sigma: Float64 = 1.0
    var epsilon: Float64 = 1.0
    var r_eq = Vec3(1.122462048, 0.0, 0.0)
    var f_eq = lennard_jones_force_3d(r_eq, epsilon, sigma)
    var f_mag = f_eq.length()
    assert_true(f_mag < 1e-4, "Lennard-Jones force vanishes at equilibrium distance 2^(1/6)*sigma")

    # Strongly repulsive at r < r_eq
    var r_rep = Vec3(0.9, 0.0, 0.0)
    var f_rep = lennard_jones_force_3d(r_rep, epsilon, sigma)
    assert_true(f_rep.x > 0.0, "Lennard-Jones force is strongly repulsive at r < r_eq")

    # Attractive at r > r_eq
    var r_att = Vec3(1.5, 0.0, 0.0)
    var f_att = lennard_jones_force_3d(r_att, epsilon, sigma)
    assert_true(f_att.x < 0.0, "Lennard-Jones force is attractive at r > r_eq")

    # 4. Velocity Verlet Integrator test
    var particles = List[Particle3D]()
    particles.append(Particle3D(id=0, pos=Vec3(0.0, 0.0, 0.0), vel=Vec3(10.0, 0.0, 0.0), mass=2.0))
    var forces = List[Vec3]()
    # Constant force of 4.0 in X => acceleration a = F/m = 2.0
    forces.append(Vec3(4.0, 0.0, 0.0))

    var dt: Float64 = 0.1
    integrate_velocity_verlet_step_3d(particles, dt, forces)

    # Expected: x = 0 + 10*0.1 + 0 = 1.0 (with a_0 = 0 initially)
    # v = 10 + 0.5 * (0 + 2.0) * 0.1 = 10.1
    assert_true(particles[0].pos.x > 0.99 and particles[0].pos.x < 1.01, "Velocity Verlet step updates position correctly")
    assert_true(particles[0].vel.x > 10.09 and particles[0].vel.x < 10.11, "Velocity Verlet step updates velocity with new acceleration")

    print("=== All Physics Tests Completed ===")
