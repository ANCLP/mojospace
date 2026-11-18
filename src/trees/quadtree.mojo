from std.math import max, sqrt
from std.collections import List
from src.math.vec2 import Vec2
from src.math.bounds import AABB2D
from src.particle.particle2d import Particle2D

@fieldwise_init
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

    def __init__(out self, bounds: AABB2D):
        self.bounds = bounds
        self.center_of_mass = Vec2()
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

    def __init__(out self, bounds: AABB2D, max_depth: Int = 16):
        self.nodes = List[QuadNode]()
        self.root_bounds = bounds
        self.max_depth = max_depth

    def clear(mut self):
        self.nodes.clear()

    def _subdivide(mut self, node_idx: Int):
        var b = self.nodes[node_idx].bounds
        var c = b.center()

        var nw_box = AABB2D(Vec2(b.min_pt.x, c.y), Vec2(c.x, b.max_pt.y))
        var ne_box = AABB2D(c, b.max_pt)
        var sw_box = AABB2D(b.min_pt, c)
        var se_box = AABB2D(Vec2(c.x, b.min_pt.y), Vec2(b.max_pt.x, c.y))

        var idx_nw = len(self.nodes)
        self.nodes.append(QuadNode(nw_box))
        var idx_ne = len(self.nodes)
        self.nodes.append(QuadNode(ne_box))
        var idx_sw = len(self.nodes)
        self.nodes.append(QuadNode(sw_box))
        var idx_se = len(self.nodes)
        self.nodes.append(QuadNode(se_box))

        self.nodes[node_idx].child_nw = idx_nw
        self.nodes[node_idx].child_ne = idx_ne
        self.nodes[node_idx].child_sw = idx_sw
        self.nodes[node_idx].child_se = idx_se
        self.nodes[node_idx].is_leaf = False

    def _get_quadrant(self, node_idx: Int, pos: Vec2) -> Int:
        var c = self.nodes[node_idx].bounds.center()
        var node = self.nodes[node_idx]
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

    def _insert_particle(
        mut self,
        node_idx: Int,
        particle_idx: Int,
        particles: List[Particle2D],
        depth: Int
    ):
        var p = particles[particle_idx]

        # Update center of mass and total mass
        var new_mass = self.nodes[node_idx].total_mass + p.mass
        if new_mass > 0.0:
            var cm = (self.nodes[node_idx].center_of_mass * self.nodes[node_idx].total_mass + p.pos * p.mass) / new_mass
            self.nodes[node_idx].center_of_mass = cm
            self.nodes[node_idx].total_mass = new_mass

        if self.nodes[node_idx].is_leaf:
            if self.nodes[node_idx].particle_idx == -1:
                # Empty leaf, insert particle here
                self.nodes[node_idx].particle_idx = particle_idx
                return
            elif depth < self.max_depth:
                # Leaf already occupied; subdivide and reinsert existing and new particle
                var old_p = self.nodes[node_idx].particle_idx
                self.nodes[node_idx].particle_idx = -1
                self._subdivide(node_idx)

                var quad_old = self._get_quadrant(node_idx, particles[old_p].pos)
                self._insert_particle(quad_old, old_p, particles, depth + 1)

                var quad_new = self._get_quadrant(node_idx, p.pos)
                self._insert_particle(quad_new, particle_idx, particles, depth + 1)
                return
            else:
                # Reached maximum depth, retain existing
                return
        else:
            var target_child = self._get_quadrant(node_idx, p.pos)
            self._insert_particle(target_child, particle_idx, particles, depth + 1)

    def build(mut self, particles: List[Particle2D]):
        """Constructs quadtree from particle array."""
        self.clear()
        if len(particles) == 0:
            return

        self.nodes.append(QuadNode(self.root_bounds))
        for i in range(len(particles)):
            self._insert_particle(0, i, particles, 0)

    def compute_force_barnes_hut(
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
            return Vec2()

        var target_p = particles[target_idx]
        var total_force = Vec2()
        var eps2 = softening * softening

        var stack = List[Int]()
        stack.append(0)

        while len(stack) > 0:
            var curr = stack.pop()
            var node = self.nodes[curr]

            if node.total_mass == 0.0:
                continue

            if node.is_leaf:
                if node.particle_idx >= 0 and node.particle_idx != target_idx:
                    var diff = node.center_of_mass - target_p.pos
                    var r2 = diff.length_squared() + eps2
                    var r = sqrt(r2)
                    var force_mag = G * (target_p.mass * node.total_mass) / (r2 * r)
                    total_force = total_force + diff * (force_mag / r)
            else:
                var diff = node.center_of_mass - target_p.pos
                var d2 = diff.length_squared()
                var d = sqrt(d2)
                var s_vec = node.bounds.size()
                var s = max(s_vec.x, s_vec.y)

                if (s / d) < theta:
                    # Satisfies MAC opening angle criterion
                    var r2 = d2 + eps2
                    var r = sqrt(r2)
                    var force_mag = G * (target_p.mass * node.total_mass) / (r2 * r)
                    total_force = total_force + diff * (force_mag / r)
                else:
                    if node.child_nw >= 0:
                        stack.append(node.child_nw)
                        stack.append(node.child_ne)
                        stack.append(node.child_sw)
                        stack.append(node.child_se)

        return total_force

    def query_radius(
        self,
        center: Vec2,
        radius: Float64,
        particles: List[Particle2D]
    ) -> List[Int]:
        """Geometric radius query using AABB subtree distance pruning."""
        var results = List[Int]()
        if len(self.nodes) == 0:
            return results^

        var r2 = radius * radius
        var stack = List[Int]()
        stack.append(0)

        while len(stack) > 0:
            var curr = stack.pop()
            var node = self.nodes[curr]

            if node.bounds.distance_squared_to(center) > r2:
                continue

            if node.is_leaf:
                if node.particle_idx >= 0:
                    var dist2 = center.distance_squared_to(particles[node.particle_idx].pos)
                    if dist2 <= r2:
                        results.append(node.particle_idx)
            else:
                if node.child_nw >= 0:
                    stack.append(node.child_nw)
                    stack.append(node.child_ne)
                    stack.append(node.child_sw)
                    stack.append(node.child_se)

        return results^
