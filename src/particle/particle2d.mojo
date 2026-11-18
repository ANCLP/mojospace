from src.math.vec2 import Vec2

@value
struct Particle2D:
    """2D Particle data container."""
    var id: Int
    var pos: Vec2
    var vel: Vec2
    var acc: Vec2
    var mass: Float64
    var radius: Float64
    var charge: Float64

    fn __init__(
        out self,
        id: Int = 0,
        pos: Vec2 = Vec2(0.0, 0.0),
        vel: Vec2 = Vec2(0.0, 0.0),
        acc: Vec2 = Vec2(0.0, 0.0),
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
