from collections import List
from src.math.vec3 import Vec3
from src.particle.particle3d import Particle3D
from src.spatial_hash.hash_grid_3d import SpatialHashGrid3D
from src.physics.force import sph_cubic_spline_kernel_3d

def main():
    print("=== mojospace: SPH Fluid Neighbor Search & Density Calculation ===")

    var smoothing_length: Float64 = 1.0
    var cutoff_radius: Float64 = 2.0 * smoothing_length

    # Initialize 3D Spatial Hash Grid with cell size equal to cutoff radius
    var grid = SpatialHashGrid3D(cell_size=cutoff_radius, table_size=10007)

    # Create a 3D block of fluid particles arranged on a lattice
    var particles = List[Particle3D]()
    var n_per_axis = 4
    var spacing: Float64 = 0.5
    var p_id = 0

    for i in range(n_per_axis):
        for j in range(n_per_axis):
            for k in range(n_per_axis):
                var pos = Vec3(
                    Float64(i) * spacing,
                    Float64(j) * spacing,
                    Float64(k) * spacing
                )
                particles.append(Particle3D(id=p_id, pos=pos, mass=0.02, radius=0.1))
                p_id += 1

    print("Created", len(particles), "SPH fluid particles.")

    # Build spatial hash table
    grid.build(particles)

    # Compute SPH density for each particle: rho_i = sum_j m_j * W(r_ij, h)
    var densities = List[Float64]()
    for i in range(len(particles)):
        var neighbors = grid.query_radius(particles[i].pos, cutoff_radius, particles)
        var rho: Float64 = 0.0
        for n_idx in range(len(neighbors)):
            var j = neighbors[n_idx]
            var r = particles[i].pos.distance_to(particles[j].pos)
            var w = sph_cubic_spline_kernel_3d(r, smoothing_length)
            rho += particles[j].mass * w
        densities.append(rho)

    print("SPH Density computation complete for all particles.")
    print("Center particle (ID 21) density:", densities[21])
    print("Corner particle (ID 0) density:", densities[0])
    print("Corner density is lower due to boundary vacuum, as physically expected.")
