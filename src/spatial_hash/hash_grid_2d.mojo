from math import floor
from collections import List
from src.math.vec2 import Vec2
from src.particle.particle2d import Particle2D

struct SpatialHashGrid2D:
    """
    2D Spatial Hash Grid based on Teschner et al. (2003).
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

    fn hash_coord(self, cx: Int, cy: Int) -> Int:
        """
        Teschner et al. (2003) 2D prime hash:
        Hash(cx, cy) = ((cx * p1) ^ (cy * p2)) mod N
        """
        let p1 = 73856093
        let p2 = 19349663
        var h = (cx * p1) ^ (cy * p2)
        h = h % self.table_size
        if h < 0:
            h += self.table_size
        return h

    fn get_cell_coords(self, pos: Vec2) -> Tuple[Int, Int]:
        """Maps continuous position to discrete integer cell coordinates."""
        let cx = Int(floor(pos.x * self.inv_cell_size))
        let cy = Int(floor(pos.y * self.inv_cell_size))
        return (cx, cy)

    fn clear(inout self):
        """Resets all buckets without reallocating table memory."""
        for i in range(self.table_size):
            self.buckets[i].clear()

    fn insert(inout self, particle_idx: Int, pos: Vec2):
        """Inserts a particle index into its corresponding cell bucket."""
        let coords = self.get_cell_coords(pos)
        let h = self.hash_coord(coords[0], coords[1])
        self.buckets[h].append(particle_idx)

    fn build(inout self, particles: List[Particle2D]):
        """Rebuilds the entire spatial grid from a list of particles."""
        self.clear()
        for i in range(len(particles)):
            self.insert(i, particles[i].pos)

    fn query_radius(
        self,
        center: Vec2,
        radius: Float64,
        particles: List[Particle2D]
    ) -> List[Int]:
        """
        Finds all particles within Euclidean distance radius from center.
        Inspects the 9-cell stencil (or larger if radius > cell_size).
        """
        var result = List[Int]()
        let r2 = radius * radius
        let span = Int(floor(radius * self.inv_cell_size)) + 1
        let center_coords = self.get_cell_coords(center)
        let cx = center_coords[0]
        let cy = center_coords[1]

        for dx in range(-span, span + 1):
            for dy in range(-span, span + 1):
                let h = self.hash_coord(cx + dx, cy + dy)
                let bucket = self.buckets[h]
                for k in range(len(bucket)):
                    let p_idx = bucket[k]
                    let dist2 = center.distance_squared_to(particles[p_idx].pos)
                    if dist2 <= r2:
                        result.append(p_idx)

        return result
