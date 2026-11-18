from std.math import sqrt

@fieldwise_init
struct Vec2(ImplicitlyCopyable):
    """2D Vector with standard vector arithmetic and Euclidean metrics."""
    var x: Float64
    var y: Float64

    def __init__(out self):
        self.x = 0.0
        self.y = 0.0

    def __add__(self, other: Vec2) -> Vec2:
        return Vec2(self.x + other.x, self.y + other.y)

    def __sub__(self, other: Vec2) -> Vec2:
        return Vec2(self.x - other.x, self.y - other.y)

    def __mul__(self, scalar: Float64) -> Vec2:
        return Vec2(self.x * scalar, self.y * scalar)

    def __truediv__(self, scalar: Float64) -> Vec2:
        return Vec2(self.x / scalar, self.y / scalar)

    def __neg__(self) -> Vec2:
        return Vec2(-self.x, -self.y)

    def __eq__(self, other: Vec2) -> Bool:
        return self.x == other.x and self.y == other.y

    def __ne__(self, other: Vec2) -> Bool:
        return not (self == other)

    def dot(self, other: Vec2) -> Float64:
        """Computes the dot product of two vectors."""
        return self.x * other.x + self.y * other.y

    def length_squared(self) -> Float64:
        """Returns the squared Euclidean length."""
        return self.dot(self)

    def length(self) -> Float64:
        """Returns the Euclidean length (magnitude)."""
        return sqrt(self.length_squared())

    def normalized(self) -> Vec2:
        """Returns the unit vector in the direction of self."""
        var l = self.length()
        if l == 0.0:
            return Vec2(0.0, 0.0)
        return self / l

    def distance_squared_to(self, other: Vec2) -> Float64:
        """Computes squared Euclidean distance to target vector."""
        var dx = self.x - other.x
        var dy = self.y - other.y
        return dx * dx + dy * dy

    def distance_to(self, other: Vec2) -> Float64:
        """Computes Euclidean distance to target vector."""
        return sqrt(self.distance_squared_to(other))
