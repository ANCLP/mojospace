from std.collections import List
from src.math.vec3 import Vec3
from src.math.bounds import AABB3D
from src.particle.particle3d import Particle3D
from src.trees.octree import BarnesHutOctree
from src.physics.integrator import integrate_velocity_verlet_step_3d

def main():
    print("=== mojospace: N-Body Gravitational Simulation (Barnes-Hut Octree) ===")

    # Initialize a system with a central massive body (star) and orbiting bodies (planets/dust)
    var particles = List[Particle3D]()
    
    # Central star (ID 0)
    particles.append(Particle3D(
        0,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        Vec3(),
        1000.0,
        1.0,
        0.0
    ))

    # Orbiting bodies in stable circular orbits
    particles.append(Particle3D(
        1,
        Vec3(10.0, 0.0, 0.0),
        Vec3(0.0, 10.0, 0.0),
        Vec3(),
        1.0,
        0.1,
        0.0
    ))
    particles.append(Particle3D(
        2,
        Vec3(0.0, 20.0, 0.0),
        Vec3(-7.07, 0.0, 0.0),
        Vec3(),
        1.0,
        0.1,
        0.0
    ))
    particles.append(Particle3D(
        3,
        Vec3(-30.0, 0.0, 0.0),
        Vec3(0.0, -5.77, 0.0),
        Vec3(),
        1.0,
        0.1,
        0.0
    ))

    var sim_bounds = AABB3D(Vec3(-100.0, -100.0, -100.0), Vec3(100.0, 100.0, 100.0))
    var octree = BarnesHutOctree(sim_bounds, max_depth=12)

    var dt: Float64 = 0.05
    var num_steps = 10

    print("Initial Positions:")
    for i in range(len(particles)):
        print("Body", particles[i].id, "pos: (", particles[i].pos.x, ",", particles[i].pos.y, ",", particles[i].pos.z, ")")

    print("\nAdvancing simulation steps...")
    for _ in range(num_steps):
        # 1. Build Barnes-Hut Octree for current particle configurations
        octree.build(particles)

        # 2. Compute gravitational forces using Barnes-Hut MAC (theta = 0.5)
        var forces = List[Vec3]()
        for i in range(len(particles)):
            var f = octree.compute_force_barnes_hut(
                target_idx=i,
                particles=particles,
                theta=0.5,
                G=1.0,
                softening=0.1
            )
            forces.append(f)

        # 3. Step forward with Velocity Verlet
        integrate_velocity_verlet_step_3d(particles, dt, forces)

    print("\nPositions after", num_steps, "steps:")
    for i in range(len(particles)):
        print("Body", particles[i].id, "pos: (", particles[i].pos.x, ",", particles[i].pos.y, ",", particles[i].pos.z, ")")

    print("Simulation finished successfully.")
