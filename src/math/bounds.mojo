from src.math.vec2 import Vec2
from src.math.vec3 import Vec3

@fieldwise_init
struct AABB2D(ImplicitlyCopyable):
    """2D Axis-Aligned Bounding Box."""
    var min_pt: Vec2
    var max_pt: Vec2

    def __init__(out self):
        self.min_pt = Vec2()
        self.max_pt = Vec2()

    def size(self) -> Vec2:
        return self.max_pt - self.min_pt

    def center(self) -> Vec2:
        return (self.min_pt + self.max_pt) * 0.5

    def contains(self, p: Vec2) -> Bool:
        return (p.x >= self.min_pt.x and p.x <= self.max_pt.x and
                p.y >= self.min_pt.y and p.y <= self.max_pt.y)

    def intersects(self, other: AABB2D) -> Bool:
        return not (self.max_pt.x < other.min_pt.x or self.min_pt.x > other.max_pt.x or
                    self.max_pt.y < other.min_pt.y or self.min_pt.y > other.max_pt.y)

    def distance_squared_to(self, p: Vec2) -> Float64:
        """Minimum squared distance from point p to the bounding box surface."""
        var d2: Float64 = 0.0
        if p.x < self.min_pt.x:
            var d = self.min_pt.x - p.x
            d2 += d * d
        elif p.x > self.max_pt.x:
            var d = p.x - self.max_pt.x
            d2 += d * d

        if p.y < self.min_pt.y:
            var d = self.min_pt.y - p.y
            d2 += d * d
        elif p.y > self.max_pt.y:
            var d = p.y - self.max_pt.y
            d2 += d * d

        return d2

    def expand(mut self, p: Vec2):
        if p.x < self.min_pt.x:
            self.min_pt.x = p.x
        if p.x > self.max_pt.x:
            self.max_pt.x = p.x
        if p.y < self.min_pt.y:
            self.min_pt.y = p.y
        if p.y > self.max_pt.y:
            self.max_pt.y = p.y


@fieldwise_init
struct AABB3D(ImplicitlyCopyable):
    """3D Axis-Aligned Bounding Box."""
    var min_pt: Vec3
    var max_pt: Vec3

    def __init__(out self):
        self.min_pt = Vec3()
        self.max_pt = Vec3()

    def size(self) -> Vec3:
        return self.max_pt - self.min_pt

    def center(self) -> Vec3:
        return (self.min_pt + self.max_pt) * 0.5

    def contains(self, p: Vec3) -> Bool:
        return (p.x >= self.min_pt.x and p.x <= self.max_pt.x and
                p.y >= self.min_pt.y and p.y <= self.max_pt.y and
                p.z >= self.min_pt.z and p.z <= self.max_pt.z)

    def intersects(self, other: AABB3D) -> Bool:
        return not (self.max_pt.x < other.min_pt.x or self.min_pt.x > other.max_pt.x or
                    self.max_pt.y < other.min_pt.y or self.min_pt.y > other.max_pt.y or
                    self.max_pt.z < other.min_pt.z or self.min_pt.z > other.max_pt.z)

    def distance_squared_to(self, p: Vec3) -> Float64:
        """Minimum squared distance from point p to the bounding box surface."""
        var d2: Float64 = 0.0
        if p.x < self.min_pt.x:
            var d = self.min_pt.x - p.x
            d2 += d * d
        elif p.x > self.max_pt.x:
            var d = p.x - self.max_pt.x
            d2 += d * d

        if p.y < self.min_pt.y:
            var d = self.min_pt.y - p.y
            d2 += d * d
        elif p.y > self.max_pt.y:
            var d = p.y - self.max_pt.y
            d2 += d * d

        if p.z < self.min_pt.z:
            var d = self.min_pt.z - p.z
            d2 += d * d
        elif p.z > self.max_pt.z:
            var d = p.z - self.max_pt.z
            d2 += d * d

        return d2

    def expand(mut self, p: Vec3):
        if p.x < self.min_pt.x:
            self.min_pt.x = p.x
        if p.x > self.max_pt.x:
            self.max_pt.x = p.x
        if p.y < self.min_pt.y:
            self.min_pt.y = p.y
        if p.y > self.max_pt.y:
            self.max_pt.y = p.y
        if p.z < self.min_pt.z:
            self.min_pt.z = p.z
        if p.z > self.max_pt.z:
            self.max_pt.z = p.z
