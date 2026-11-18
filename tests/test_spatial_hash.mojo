from collections import List
from src.math.vec2 import Vec2
from src.math.vec3 import Vec3
from src.particle.particle2d import Particle2D
from src.particle.particle3d import Particle3D
from src.spatial_hash.hash_grid_2d import SpatialHashGrid2D
from src.spatial_hash.hash_grid_3d import SpatialHashGrid3D

fn assert_true(condition: Bool, msg: String):
    if not condition:
        print("[FAIL] " + msg)
    else:
        print("[PASS] " + msg)

fn main():
    print("=== Testing Spatial Hash Grids ===")

    # 1. Test 2D Spatial Hash Grid
    var grid2d = SpatialHashGrid2D(cell_size=1.0, table_size=1009)
    var p2d = List[Particle2D]()
    # Origin particle
    p2d.append(Particle2D(id=0, pos=Vec2(0.0, 0.0)))
    # Close neighbor within radius 0.8
    p2d.append(Particle2D(id=1, pos=Vec2(0.5, 0.5))) # dist ~ 0.707
    # Boundary neighbor at radius ~ 1.2
    p2d.append(Particle2D(id=2, pos=Vec2(0.8, 0.9))) # dist ~ 1.204
    # Distant particle
    p2d.append(Particle2D(id=3, pos=Vec2(10.0, 10.0)))

    grid2d.build(p2d)

    let neighbors_0_8 = grid2d.query_radius(Vec2(0.0, 0.0), 0.8, p2d)
    # Expect particles 0 and 1
    assert_true(len(neighbors_0_8) == 2, "2D query radius 0.8 found exactly 2 particles (self and id=1)")

    let neighbors_1_5 = grid2d.query_radius(Vec2(0.0, 0.0), 1.5, p2d)
    # Expect particles 0, 1, 2
    assert_true(len(neighbors_1_5) == 3, "2D query radius 1.5 found exactly 3 particles")

    # 2. Test 3D Spatial Hash Grid
    var grid3d = SpatialHashGrid3D(cell_size=2.0, table_size=5003)
    var p3d = List[Particle3D]()
    p3d.append(Particle3D(id=0, pos=Vec3(-1.0, -1.0, -1.0)))
    p3d.append(Particle3D(id=1, pos=Vec3(-0.5, -0.5, -0.5)))
    p3d.append(Particle3D(id=2, pos=Vec3(1.0, 1.0, 1.0)))
    p3d.append(Particle3D(id=3, pos=Vec3(100.0, 100.0, 100.0)))

    grid3d.build(p3d)

    let neighbors_3d = grid3d.query_radius(Vec3(-1.0, -1.0, -1.0), 1.0, p3d)
    # Expect particle 0 (dist 0) and particle 1 (dist ~ 0.866)
    assert_true(len(neighbors_3d) == 2, "3D query radius 1.0 found 2 negative-coordinate particles")

    let empty_neighbors = grid3d.query_radius(Vec3(50.0, 50.0, 50.0), 1.0, p3d)
    assert_true(len(empty_neighbors) == 0, "3D query in empty region returned 0 particles")

    print("=== All Spatial Hash Tests Completed ===")
