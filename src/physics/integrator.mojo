from std.collections import List
from src.math.vec2 import Vec2
from src.math.vec3 import Vec3
from src.particle.particle2d import Particle2D
from src.particle.particle3d import Particle3D

def integrate_symplectic_euler_step_2d(
    mut particles: List[Particle2D],
    dt: Float64,
    forces: List[Vec2]
):
    """
    First-order Symplectic Euler step:
    v(t + dt) = v(t) + (F / m) * dt
    x(t + dt) = x(t) + v(t + dt) * dt
    """
    for i in range(len(particles)):
        var inv_m = 1.0 / particles[i].mass
        var a = forces[i] * inv_m
        particles[i].vel = particles[i].vel + a * dt
        particles[i].pos = particles[i].pos + particles[i].vel * dt
        particles[i].acc = a

def integrate_symplectic_euler_step_3d(
    mut particles: List[Particle3D],
    dt: Float64,
    forces: List[Vec3]
):
    """3D First-order Symplectic Euler step."""
    for i in range(len(particles)):
        var inv_m = 1.0 / particles[i].mass
        var a = forces[i] * inv_m
        particles[i].vel = particles[i].vel + a * dt
        particles[i].pos = particles[i].pos + particles[i].vel * dt
        particles[i].acc = a

def integrate_velocity_verlet_step_3d(
    mut particles: List[Particle3D],
    dt: Float64,
    forces: List[Vec3]
):
    """
    Velocity Verlet numerical integration step (3D):
    x(t + dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt^2
    v(t + dt) = v(t) + 0.5 * (a(t) + a(t + dt)) * dt
    """
    var half_dt = 0.5 * dt
    for i in range(len(particles)):
        var inv_m = 1.0 / particles[i].mass
        var new_a = forces[i] * inv_m
        # Update position
        particles[i].pos = particles[i].pos + particles[i].vel * dt + particles[i].acc * (0.5 * dt * dt)
        # Update velocity with average acceleration
        particles[i].vel = particles[i].vel + (particles[i].acc + new_a) * half_dt
        particles[i].acc = new_a

def integrate_velocity_verlet_step_2d(
    mut particles: List[Particle2D],
    dt: Float64,
    forces: List[Vec2]
):
    """Velocity Verlet numerical integration step (2D)."""
    var half_dt = 0.5 * dt
    for i in range(len(particles)):
        var inv_m = 1.0 / particles[i].mass
        var new_a = forces[i] * inv_m
        particles[i].pos = particles[i].pos + particles[i].vel * dt + particles[i].acc * (0.5 * dt * dt)
        particles[i].vel = particles[i].vel + (particles[i].acc + new_a) * half_dt
        particles[i].acc = new_a
