from math import sqrt

@value
struct Vec2:
    """2D Vector with standard vector arithmetic and Euclidean metrics."""
    var x: Float64
    var y: Float64

    fn __init__(out self, x: Float64 = 0.0, y: Float64 = 0.0):
        self.x = x
        self.y = y

    fn __add__(self, other: Vec2) -> Vec2:
        return Vec2(self.x + other.x, self.y + other.y)

    fn __sub__(self, other: Vec2) -> Vec2:
        return Vec2(self.x - other.x, self.y - other.y)

    fn __mul__(self, scalar: Float64) -> Vec2:
        return Vec2(self.x * scalar, self.y * scalar)

    fn __truediv__(self, scalar: Float64) -> Vec2:
        return Vec2(self.x / scalar, self.y / scalar)

    fn __neg__(self) -> Vec2:
        return Vec2(-self.x, -self.y)

    fn __eq__(self, other: Vec2) -> Bool:
        return self.x == other.x and self.y == other.y

    fn __ne__(self, other: Vec2) -> Bool:
        return not (self == other)

    fn dot(self, other: Vec2) -> Float64:
        """Computes the dot product of two vectors."""
        return self.x * other.x + self.y * other.y

    fn length_squared(self) -> Float64:
        """Returns the squared Euclidean length."""
        return self.dot(self)

    fn length(self) -> Float64:
        """Returns the Euclidean length (magnitude)."""
        return sqrt(self.length_squared())

    fn normalized(self) -> Vec2:
        """Returns the unit vector in the direction of self."""
        let l = self.length()
        if l == 0.0:
            return Vec2(0.0, 0.0)
        return self / l

    fn distance_squared_to(self, other: Vec2) -> Float64:
        """Computes squared Euclidean distance to target vector."""
        let dx = self.x - other.x
        let dy = self.y - other.y
        return dx * dx + dy * dy

    fn distance_to(self, other: Vec2) -> Float64:
        """Computes Euclidean distance to target vector."""
        return sqrt(self.distance_squared_to(other))
