from std.math import sqrt

@fieldwise_init
struct Vec3(ImplicitlyCopyable):
    """3D Vector with standard vector arithmetic, cross product, and metrics."""
    var x: Float64
    var y: Float64
    var z: Float64

    def __init__(out self):
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0

    def __add__(self, other: Vec3) -> Vec3:
        return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: Vec3) -> Vec3:
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, scalar: Float64) -> Vec3:
        return Vec3(self.x * scalar, self.y * scalar, self.z * scalar)

    def __truediv__(self, scalar: Float64) -> Vec3:
        return Vec3(self.x / scalar, self.y / scalar, self.z / scalar)

    def __neg__(self) -> Vec3:
        return Vec3(-self.x, -self.y, -self.z)

    def __eq__(self, other: Vec3) -> Bool:
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __ne__(self, other: Vec3) -> Bool:
        return not (self == other)

    def dot(self, other: Vec3) -> Float64:
        """Computes the dot product of two vectors."""
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other: Vec3) -> Vec3:
        """Computes the 3D cross product self x other."""
        return Vec3(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x
        )

    def length_squared(self) -> Float64:
        """Returns the squared Euclidean length."""
        return self.dot(self)

    def length(self) -> Float64:
        """Returns the Euclidean length (magnitude)."""
        return sqrt(self.length_squared())

    def normalized(self) -> Vec3:
        """Returns the unit vector in the direction of self."""
        var l = self.length()
        if l == 0.0:
            return Vec3(0.0, 0.0, 0.0)
        return self / l

    def distance_squared_to(self, other: Vec3) -> Float64:
        """Computes squared Euclidean distance to target vector."""
        var dx = self.x - other.x
        var dy = self.y - other.y
        var dz = self.z - other.z
        return dx * dx + dy * dy + dz * dz

    def distance_to(self, other: Vec3) -> Float64:
        """Computes Euclidean distance to target vector."""
        return sqrt(self.distance_squared_to(other))
