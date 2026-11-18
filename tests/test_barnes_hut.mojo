from std.collections import List
from src.math.vec2 import Vec2
from src.math.vec3 import Vec3
from src.math.bounds import AABB2D, AABB3D
from src.particle.particle2d import Particle2D
from src.particle.particle3d import Particle3D
from src.trees.quadtree import Quadtree
from src.trees.octree import BarnesHutOctree
from src.physics.force import gravitational_force_direct_2d, gravitational_force_direct_3d

def assert_true(condition: Bool, msg: String):
    if not condition:
        print("[FAIL] " + msg)
    else:
        print("[PASS] " + msg)

def main():
    print("=== Testing Barnes-Hut Trees ===")

    # 1. Test 2D Quadtree Barnes-Hut
    var bounds_2d = AABB2D(Vec2(-50.0, -50.0), Vec2(50.0, 50.0))
    var quadtree = Quadtree(bounds_2d, max_depth=8)
    var p2d = List[Particle2D]()

    # Setup particles in 2D
    p2d.append(Particle2D(0, Vec2(0.0, 0.0), 10.0))
    p2d.append(Particle2D(1, Vec2(10.0, 0.0), 2.0))
    p2d.append(Particle2D(2, Vec2(10.2, 0.1), 2.0))
    p2d.append(Particle2D(3, Vec2(9.8, -0.1), 2.0))

    quadtree.build(p2d)

    # Compute exact direct force on particle 0
    var exact_force_2d = Vec2()
    for j in range(1, len(p2d)):
        exact_force_2d = exact_force_2d + gravitational_force_direct_2d(p2d[0], p2d[j], G=1.0, softening=0.01)

    # Compute Barnes-Hut force on particle 0 with theta = 0.5
    var bh_force_2d = quadtree.compute_force_barnes_hut(0, p2d, theta=0.5, G=1.0, softening=0.01)

    var rel_err_2d = (bh_force_2d - exact_force_2d).length() / exact_force_2d.length()
    assert_true(rel_err_2d < 0.05, "2D Quadtree Barnes-Hut force relative error < 5%")

    # 2. Test 3D Octree Barnes-Hut
    var bounds_3d = AABB3D(Vec3(-100.0, -100.0, -100.0), Vec3(100.0, 100.0, 100.0))
    var octree = BarnesHutOctree(bounds_3d, max_depth=10)
    var p3d = List[Particle3D]()

    # Target particle at origin
    p3d.append(Particle3D(0, Vec3(0.0, 0.0, 0.0), 5.0))
    # Distant cluster of 4 particles centered around (30, 30, 30)
    p3d.append(Particle3D(1, Vec3(30.0, 30.0, 30.0), 1.0))
    p3d.append(Particle3D(2, Vec3(30.5, 30.2, 29.8), 1.0))
    p3d.append(Particle3D(3, Vec3(29.8, 30.1, 30.3), 1.0))
    p3d.append(Particle3D(4, Vec3(30.2, 29.9, 30.1), 1.0))

    octree.build(p3d)

    var exact_force_3d = Vec3()
    for j in range(1, len(p3d)):
        exact_force_3d = exact_force_3d + gravitational_force_direct_3d(p3d[0], p3d[j], G=1.0, softening=0.01)

    var bh_force_3d = octree.compute_force_barnes_hut(0, p3d, theta=0.5, G=1.0, softening=0.01)

    var rel_err_3d = (bh_force_3d - exact_force_3d).length() / exact_force_3d.length()
    assert_true(rel_err_3d < 0.02, "3D Octree Barnes-Hut force relative error < 2%")

    # 3. Test Octree Radius Query
    var query_neighbors = octree.query_radius(Vec3(30.0, 30.0, 30.0), 2.0, p3d)
    assert_true(len(query_neighbors) == 4, "3D Octree radius query found all 4 cluster particles")

    print("=== All Barnes-Hut Tests Completed ===")
