from src.math.vec3 import Vec3

@fieldwise_init
struct Particle3D:
    """3D Particle data container."""
    var id: Int
    var pos: Vec3
    var vel: Vec3
    var acc: Vec3
    var mass: Float64
    var radius: Float64
    var charge: Float64

    def __init__(out self, id: Int, pos: Vec3):
        self.id = id
        self.pos = pos
        self.vel = Vec3()
        self.acc = Vec3()
        self.mass = 1.0
        self.radius = 0.1
        self.charge = 0.0

    def __init__(out self, id: Int, pos: Vec3, mass: Float64):
        self.id = id
        self.pos = pos
        self.vel = Vec3()
        self.acc = Vec3()
        self.mass = mass
        self.radius = 0.1
        self.charge = 0.0
