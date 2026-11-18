from math import max, sqrt
from collections import List
from src.math.vec2 import Vec2
from src.math.bounds import AABB2D
from src.particle.particle2d import Particle2D

@value
struct QuadNode:
    var bounds: AABB2D
    var center_of_mass: Vec2
    var total_mass: Float64
    var particle_idx: Int
    var is_leaf: Bool
    var child_nw: Int
    var child_ne: Int
    var child_sw: Int
    var child_se: Int

    fn __init__(out self, bounds: AABB2D):
        self.bounds = bounds
        self.center_of_mass = Vec2(0.0, 0.0)
        self.total_mass = 0.0
        self.particle_idx = -1
        self.is_leaf = True
        self.child_nw = -1
        self.child_ne = -1
        self.child_sw = -1
        self.child_se = -1

struct Quadtree:
    """
    2D Quadtree supporting Barnes-Hut O(N log N) force approximation
    and geometric radius neighbor pruning.
    """
    var nodes: List[QuadNode]
    var root_bounds: AABB2D
    var max_depth: Int

    fn __init__(out self, bounds: AABB2D, max_depth: Int = 16):
        self.nodes = List[QuadNode]()
        self.root_bounds = bounds
        self.max_depth = max_depth

    fn clear(inout self):
        self.nodes.clear()

    fn _subdivide(inout self, node_idx: Int):
        let b = self.nodes[node_idx].bounds
        let c = b.center()

        let nw_box = AABB2D(Vec2(b.min_pt.x, c.y), Vec2(c.x, b.max_pt.y))
        let ne_box = AABB2D(c, b.max_pt)
        let sw_box = AABB2D(b.min_pt, c)
        let se_box = AABB2D(Vec2(c.x, b.min_pt.y), Vec2(b.max_pt.x, c.y))

        let idx_nw = len(self.nodes)
        self.nodes.append(QuadNode(nw_box))
        let idx_ne = len(self.nodes)
        self.nodes.append(QuadNode(ne_box))
        let idx_sw = len(self.nodes)
        self.nodes.append(QuadNode(sw_box))
        let idx_se = len(self.nodes)
        self.nodes.append(QuadNode(se_box))

        self.nodes[node_idx].child_nw = idx_nw
        self.nodes[node_idx].child_ne = idx_ne
        self.nodes[node_idx].child_sw = idx_sw
        self.nodes[node_idx].child_se = idx_se
        self.nodes[node_idx].is_leaf = False

    fn _get_quadrant(self, node_idx: Int, pos: Vec2) -> Int:
        let c = self.nodes[node_idx].bounds.center()
        let node = self.nodes[node_idx]
        if pos.x < c.x:
            if pos.y < c.y:
                return node.child_sw
            else:
                return node.child_nw
        else:
            if pos.y < c.y:
                return node.child_se
            else:
                return node.child_ne

    fn _insert_particle(
        inout self,
        node_idx: Int,
        particle_idx: Int,
        particles: List[Particle2D],
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
                # Empty leaf, insert particle here
                self.nodes[node_idx].particle_idx = particle_idx
                return
            elif depth < self.max_depth:
                # Leaf already occupied; subdivide and reinsert existing and new particle
                let old_p = self.nodes[node_idx].particle_idx
                self.nodes[node_idx].particle_idx = -1
                self._subdivide(node_idx)

                let quad_old = self._get_quadrant(node_idx, particles[old_p].pos)
                self._insert_particle(quad_old, old_p, particles, depth + 1)

                let quad_new = self._get_quadrant(node_idx, p.pos)
                self._insert_particle(quad_new, particle_idx, particles, depth + 1)
                return
            else:
                # Reached maximum depth, retain existing
                return
        else:
            let target_child = self._get_quadrant(node_idx, p.pos)
            self._insert_particle(target_child, particle_idx, particles, depth + 1)

    fn build(inout self, particles: List[Particle2D]):
        """Constructs quadtree from particle array."""
        self.clear()
        if len(particles) == 0:
            return

        self.nodes.append(QuadNode(self.root_bounds))
        for i in range(len(particles)):
            self._insert_particle(0, i, particles, 0)

    fn compute_force_barnes_hut(
        self,
        target_idx: Int,
        particles: List[Particle2D],
        theta: Float64 = 0.5,
        G: Float64 = 1.0,
        softening: Float64 = 0.01
    ) -> Vec2:
        """
        Computes gravitational force on particle target_idx using Barnes-Hut MAC.
        """
        if len(self.nodes) == 0:
            return Vec2(0.0, 0.0)

        let target_p = particles[target_idx]
        var total_force = Vec2(0.0, 0.0)
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
                let s = max(s_vec.x, s_vec.y)

                if (s / d) < theta:
                    # Satisfies MAC opening angle criterion
                    let r2 = d2 + eps2
                    let r = sqrt(r2)
                    let force_mag = G * (target_p.mass * node.total_mass) / (r2 * r)
                    total_force = total_force + diff * (force_mag / r)
                else:
                    if node.child_nw >= 0:
                        stack.append(node.child_nw)
                        stack.append(node.child_ne)
                        stack.append(node.child_sw)
                        stack.append(node.child_se)

        return total_force

    fn query_radius(
        self,
        center: Vec2,
        radius: Float64,
        particles: List[Particle2D]
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
                if node.child_nw >= 0:
                    stack.append(node.child_nw)
                    stack.append(node.child_ne)
                    stack.append(node.child_sw)
                    stack.append(node.child_se)

        return results
