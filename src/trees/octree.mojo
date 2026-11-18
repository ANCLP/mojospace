from math import max, sqrt
from collections import List
from src.math.vec3 import Vec3
from src.math.bounds import AABB3D
from src.particle.particle3d import Particle3D

@value
struct OctNode:
    var bounds: AABB3D
    var center_of_mass: Vec3
    var total_mass: Float64
    var particle_idx: Int
    var is_leaf: Bool
    var child_0: Int
    var child_1: Int
    var child_2: Int
    var child_3: Int
    var child_4: Int
    var child_5: Int
    var child_6: Int
    var child_7: Int

    fn __init__(out self, bounds: AABB3D):
        self.bounds = bounds
        self.center_of_mass = Vec3(0.0, 0.0, 0.0)
        self.total_mass = 0.0
        self.particle_idx = -1
        self.is_leaf = True
        self.child_0 = -1
        self.child_1 = -1
        self.child_2 = -1
        self.child_3 = -1
        self.child_4 = -1
        self.child_5 = -1
        self.child_6 = -1
        self.child_7 = -1

    fn get_child(self, index: Int) -> Int:
        if index == 0: return self.child_0
        if index == 1: return self.child_1
        if index == 2: return self.child_2
        if index == 3: return self.child_3
        if index == 4: return self.child_4
        if index == 5: return self.child_5
        if index == 6: return self.child_6
        return self.child_7

    fn set_child(inout self, index: Int, child_idx: Int):
        if index == 0: self.child_0 = child_idx
        elif index == 1: self.child_1 = child_idx
        elif index == 2: self.child_2 = child_idx
        elif index == 3: self.child_3 = child_idx
        elif index == 4: self.child_4 = child_idx
        elif index == 5: self.child_5 = child_idx
        elif index == 6: self.child_6 = child_idx
        else: self.child_7 = child_idx

struct BarnesHutOctree:
    """
    3D Barnes-Hut Octree for O(N log N) gravitational/electrostatic force
    computation and spatial radius neighbor search.
    """
    var nodes: List[OctNode]
    var root_bounds: AABB3D
    var max_depth: Int

    fn __init__(out self, bounds: AABB3D, max_depth: Int = 16):
        self.nodes = List[OctNode]()
        self.root_bounds = bounds
        self.max_depth = max_depth

    fn clear(inout self):
        self.nodes.clear()

    fn _get_octant(self, center: Vec3, pos: Vec3) -> Int:
        var oct = 0
        if pos.x >= center.x: oct |= 1
        if pos.y >= center.y: oct |= 2
        if pos.z >= center.z: oct |= 4
        return oct

    fn _subdivide(inout self, node_idx: Int):
        let b = self.nodes[node_idx].bounds
        let c = b.center()

        for k in range(8):
            let min_x = c.x if (k & 1) != 0 else b.min_pt.x
            let max_x = b.max_pt.x if (k & 1) != 0 else c.x

            let min_y = c.y if (k & 2) != 0 else b.min_pt.y
            let max_y = b.max_pt.y if (k & 2) != 0 else c.y

            let min_z = c.z if (k & 4) != 0 else b.min_pt.z
            let max_z = b.max_pt.z if (k & 4) != 0 else c.z

            let child_box = AABB3D(Vec3(min_x, min_y, min_z), Vec3(max_x, max_y, max_z))
            let new_idx = len(self.nodes)
            self.nodes.append(OctNode(child_box))
            self.nodes[node_idx].set_child(k, new_idx)

        self.nodes[node_idx].is_leaf = False

    fn _insert_particle(
        inout self,
        node_idx: Int,
        particle_idx: Int,
        particles: List[Particle3D],
        depth: Int
    ):
        let p = particles[particle_idx]

        # Update center of mass and total mass
        let new_mass = self.nodes[node_idx].total_mass + p.mass
        if new_mass > 0.0:
            let cm = (self.nodes[node_idx].center_of_mass * self.nodes[node_idx].total_mass + p.pos * p.mass) / new_mass
            self.nodes[node_idx].center_of_mass = cm
            self.nodes[node_idx].total_mass = new_mass

        if self.nodes[node_idx].is_leaf:
            if self.nodes[node_idx].particle_idx == -1:
                self.nodes[node_idx].particle_idx = particle_idx
                return
            elif depth < self.max_depth:
                let old_p = self.nodes[node_idx].particle_idx
                self.nodes[node_idx].particle_idx = -1
                self._subdivide(node_idx)

                let c = self.nodes[node_idx].bounds.center()
                let oct_old = self._get_octant(c, particles[old_p].pos)
                let child_old = self.nodes[node_idx].get_child(oct_old)
                self._insert_particle(child_old, old_p, particles, depth + 1)

                let oct_new = self._get_octant(c, p.pos)
                let child_new = self.nodes[node_idx].get_child(oct_new)
                self._insert_particle(child_new, particle_idx, particles, depth + 1)
                return
            else:
                return
        else:
            let c = self.nodes[node_idx].bounds.center()
            let oct = self._get_octant(c, p.pos)
            let child = self.nodes[node_idx].get_child(oct)
            self._insert_particle(child, particle_idx, particles, depth + 1)

    fn build(inout self, particles: List[Particle3D]):
        """Constructs octree from 3D particle array."""
        self.clear()
        if len(particles) == 0:
            return

        self.nodes.append(OctNode(self.root_bounds))
        for i in range(len(particles)):
            self._insert_particle(0, i, particles, 0)

    fn compute_force_barnes_hut(
        self,
        target_idx: Int,
        particles: List[Particle3D],
        theta: Float64 = 0.5,
        G: Float64 = 1.0,
        softening: Float64 = 0.01
    ) -> Vec3:
        """
        Computes gravitational force on particle target_idx using Barnes-Hut MAC.
        """
        if len(self.nodes) == 0:
            return Vec3(0.0, 0.0, 0.0)

        let target_p = particles[target_idx]
        var total_force = Vec3(0.0, 0.0, 0.0)
        let eps2 = softening * softening

        var stack = List[Int]()
        stack.append(0)

        while len(stack) > 0:
            let curr = stack.pop()
            let node = self.nodes[curr]

            if node.total_mass == 0.0:
                continue

            if node.is_leaf:
                if node.particle_idx >= 0 and node.particle_idx != target_idx:
                    let diff = node.center_of_mass - target_p.pos
                    let r2 = diff.length_squared() + eps2
                    let r = sqrt(r2)
                    let force_mag = G * (target_p.mass * node.total_mass) / (r2 * r)
                    total_force = total_force + diff * (force_mag / r)
            else:
                let diff = node.center_of_mass - target_p.pos
                let d2 = diff.length_squared()
                let d = sqrt(d2)
                let s_vec = node.bounds.size()
                let s = max(s_vec.x, max(s_vec.y, s_vec.z))

                if (s / d) < theta:
                    # Satisfies MAC opening angle criterion
                    let r2 = d2 + eps2
                    let r = sqrt(r2)
                    let force_mag = G * (target_p.mass * node.total_mass) / (r2 * r)
                    total_force = total_force + diff * (force_mag / r)
                else:
                    for k in range(8):
                        let child_idx = node.get_child(k)
                        if child_idx >= 0:
                            stack.append(child_idx)

        return total_force

    fn query_radius(
        self,
        center: Vec3,
        radius: Float64,
        particles: List[Particle3D]
    ) -> List[Int]:
        """Geometric radius query using AABB subtree distance pruning."""
        var results = List[Int]()
        if len(self.nodes) == 0:
            return results

        let r2 = radius * radius
        var stack = List[Int]()
        stack.append(0)

        while len(stack) > 0:
            let curr = stack.pop()
            let node = self.nodes[curr]

            if node.bounds.distance_squared_to(center) > r2:
                continue

            if node.is_leaf:
                if node.particle_idx >= 0:
                    let dist2 = center.distance_squared_to(particles[node.particle_idx].pos)
                    if dist2 <= r2:
                        results.append(node.particle_idx)
            else:
                for k in range(8):
                    let child_idx = node.get_child(k)
                    if child_idx >= 0:
                        stack.append(child_idx)

        return results
