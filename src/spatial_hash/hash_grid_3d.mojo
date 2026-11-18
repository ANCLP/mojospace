from math import floor
from collections import List
from src.math.vec3 import Vec3
from src.particle.particle3d import Particle3D

struct SpatialHashGrid3D:
    """
    3D Spatial Hash Grid based on Teschner et al. (2003).
    Provides O(1) expected neighbor search within interaction radius.
    """
    var cell_size: Float64
    var inv_cell_size: Float64
    var table_size: Int
    var buckets: List[List[Int]]

    fn __init__(out self, cell_size: Float64 = 1.0, table_size: Int = 10007):
        self.cell_size = cell_size
        self.inv_cell_size = 1.0 / cell_size
        self.table_size = table_size
        self.buckets = List[List[Int]]()
        for _ in range(table_size):
            self.buckets.append(List[Int]())

    fn hash_coord(self, cx: Int, cy: Int, cz: Int) -> Int:
        """
        Teschner et al. (2003) 3D prime hash:
        Hash(cx, cy, cz) = ((cx * p1) ^ (cy * p2) ^ (cz * p3)) mod N
        """
        let p1 = 73856093
        let p2 = 19349663
        let p3 = 83492791
        var h = (cx * p1) ^ (cy * p2) ^ (cz * p3)
        h = h % self.table_size
        if h < 0:
            h += self.table_size
        return h

    fn get_cell_coords(self, pos: Vec3) -> Tuple[Int, Int, Int]:
        """Maps continuous 3D position to discrete integer cell coordinates."""
        let cx = Int(floor(pos.x * self.inv_cell_size))
        let cy = Int(floor(pos.y * self.inv_cell_size))
        let cz = Int(floor(pos.z * self.inv_cell_size))
        return (cx, cy, cz)

    fn clear(inout self):
        """Resets all buckets without reallocating table memory."""
        for i in range(self.table_size):
            self.buckets[i].clear()

    fn insert(inout self, particle_idx: Int, pos: Vec3):
        """Inserts a particle index into its corresponding cell bucket."""
        let coords = self.get_cell_coords(pos)
        let h = self.hash_coord(coords[0], coords[1], coords[2])
        self.buckets[h].append(particle_idx)

    fn build(inout self, particles: List[Particle3D]):
        """Rebuilds the entire spatial grid from a list of particles."""
        self.clear()
        for i in range(len(particles)):
            self.insert(i, particles[i].pos)

    fn query_radius(
        self,
        center: Vec3,
        radius: Float64,
        particles: List[Particle3D]
    ) -> List[Int]:
        """
        Finds all particles within Euclidean distance radius from center.
        Inspects the 27-cell stencil (or larger if radius > cell_size).
        """
        var result = List[Int]()
        let r2 = radius * radius
        let span = Int(floor(radius * self.inv_cell_size)) + 1
        let coords = self.get_cell_coords(center)
        let cx = coords[0]
        let cy = coords[1]
        let cz = coords[2]

        for dx in range(-span, span + 1):
            for dy in range(-span, span + 1):
                for dz in range(-span, span + 1):
                    let h = self.hash_coord(cx + dx, cy + dy, cz + dz)
                    let bucket = self.buckets[h]
                    for k in range(len(bucket)):
                        let p_idx = bucket[k]
                        let dist2 = center.distance_squared_to(particles[p_idx].pos)
                        if dist2 <= r2:
                            result.append(p_idx)

        return result
