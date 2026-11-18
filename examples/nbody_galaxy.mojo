from collections import List
from src.math.vec3 import Vec3
from src.math.bounds import AABB3D
from src.particle.particle3d import Particle3D
from src.trees.octree import BarnesHutOctree
from src.physics.integrator import integrate_velocity_verlet_step_3d

fn main():
    print("=== mojospace: N-Body Gravitational Simulation (Barnes-Hut Octree) ===")

    # Initialize a system with a central massive body (star) and orbiting bodies (planets/dust)
    var particles = List[Particle3D]()
    
    # Central star (ID 0)
    particles.append(Particle3D(
        id=0,
        pos=Vec3(0.0, 0.0, 0.0),
        vel=Vec3(0.0, 0.0, 0.0),
        mass=1000.0,
        radius=1.0
    ))

    # Orbiting bodies in stable circular orbits
    particles.append(Particle3D(
        id=1,
        pos=Vec3(10.0, 0.0, 0.0),
        vel=Vec3(0.0, 10.0, 0.0),
        mass=1.0,
        radius=0.1
    ))
    particles.append(Particle3D(
        id=2,
        pos=Vec3(0.0, 20.0, 0.0),
        vel=Vec3(-7.07, 0.0, 0.0),
        mass=1.0,
        radius=0.1
    ))
    particles.append(Particle3D(
        id=3,
        pos=Vec3(-30.0, 0.0, 0.0),
        vel=Vec3(0.0, -5.77, 0.0),
        mass=1.0,
        radius=0.1
    ))

    let sim_bounds = AABB3D(Vec3(-100.0, -100.0, -100.0), Vec3(100.0, 100.0, 100.0))
    var octree = BarnesHutOctree(bounds=sim_bounds, max_depth=12)

    let dt: Float64 = 0.05
    let num_steps = 10

    print("Initial Positions:")
    for i in range(len(particles)):
        print("Body", particles[i].id, "pos: (", particles[i].pos.x, ",", particles[i].pos.y, ",", particles[i].pos.z, ")")

    print("\nAdvancing simulation steps...")
    for step in range(num_steps):
        # 1. Build Barnes-Hut Octree for current particle configurations
        octree.build(particles)

        # 2. Compute gravitational forces using Barnes-Hut MAC (theta = 0.5)
        var forces = List[Vec3]()
        for i in range(len(particles)):
            let f = octree.compute_force_barnes_hut(
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
