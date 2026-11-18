from src.math.vec2 import Vec2
from src.math.vec3 import Vec3
from src.math.bounds import AABB2D, AABB3D

fn assert_true(condition: Bool, msg: String):
    if not condition:
        print("[FAIL] Assertion failed: " + msg)
    else:
        print("[PASS] " + msg)

fn assert_approx(val: Float64, expected: Float64, tol: Float64, msg: String):
    let diff = val - expected
    let abs_diff = diff if diff >= 0.0 else -diff
    if abs_diff > tol:
        print("[FAIL] " + msg + " Expected: " + String(expected) + ", got: " + String(val))
    else:
        print("[PASS] " + msg)

fn main():
    print("=== Testing Math Module ===")

    # Test Vec2
    let v2_a = Vec2(3.0, 4.0)
    let v2_b = Vec2(1.0, 2.0)
    assert_approx(v2_a.length(), 5.0, 1e-6, "Vec2 length of (3,4) is 5")
    assert_approx(v2_a.dot(v2_b), 11.0, 1e-6, "Vec2 dot product (3*1 + 4*2 = 11)")
    let v2_norm = v2_a.normalized()
    assert_approx(v2_norm.length(), 1.0, 1e-6, "Vec2 normalized has length 1")

    # Test Vec3
    let v3_x = Vec3(1.0, 0.0, 0.0)
    let v3_y = Vec3(0.0, 1.0, 0.0)
    let v3_z = v3_x.cross(v3_y)
    assert_true(v3_z == Vec3(0.0, 0.0, 1.0), "Vec3 cross product (1,0,0) x (0,1,0) = (0,0,1)")
    assert_approx(v3_x.distance_to(v3_y), 1.41421356, 1e-5, "Vec3 distance between (1,0,0) and (0,1,0)")

    # Test AABB3D
    let box = AABB3D(Vec3(-10.0, -10.0, -10.0), Vec3(10.0, 10.0, 10.0))
    assert_true(box.contains(Vec3(0.0, 0.0, 0.0)), "AABB3D contains center point")
    assert_true(not box.contains(Vec3(15.0, 0.0, 0.0)), "AABB3D rejects point outside max bounds")
    assert_approx(box.distance_squared_to(Vec3(15.0, 0.0, 0.0)), 25.0, 1e-6, "AABB3D surface distance squared is 25")

    print("=== All Math Tests Completed ===")
