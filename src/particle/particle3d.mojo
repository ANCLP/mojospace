from src.math.vec3 import Vec3

@value
struct Particle3D:
    """3D Particle data container."""
    var id: Int
    var pos: Vec3
    var vel: Vec3
    var acc: Vec3
    var mass: Float64
    var radius: Float64
    var charge: Float64

    fn __init__(
        out self,
        id: Int = 0,
        pos: Vec3 = Vec3(0.0, 0.0, 0.0),
        vel: Vec3 = Vec3(0.0, 0.0, 0.0),
        acc: Vec3 = Vec3(0.0, 0.0, 0.0),
        mass: Float64 = 1.0,
        radius: Float64 = 0.1,
        charge: Float64 = 0.0
    ):
        self.id = id
        self.pos = pos
        self.vel = vel
        self.acc = acc
        self.mass = mass
        self.radius = radius
        self.charge = charge
