from src.math.vec2 import Vec2

@fieldwise_init
struct Particle2D(ImplicitlyCopyable):
    """2D Particle data container."""
    var id: Int
    var pos: Vec2
    var vel: Vec2
    var acc: Vec2
    var mass: Float64
    var radius: Float64
    var charge: Float64

    def __init__(out self, id: Int, pos: Vec2):
        self.id = id
        self.pos = pos
        self.vel = Vec2()
        self.acc = Vec2()
        self.mass = 1.0
        self.radius = 0.1
        self.charge = 0.0

    def __init__(out self, id: Int, pos: Vec2, mass: Float64):
        self.id = id
        self.pos = pos
        self.vel = Vec2()
        self.acc = Vec2()
        self.mass = mass
        self.radius = 0.1
        self.charge = 0.0
