from math import sqrt

@value
struct Vec3:
    """3D Vector with standard vector arithmetic, cross product, and metrics."""
    var x: Float64
    var y: Float64
    var z: Float64

    fn __init__(out self, x: Float64 = 0.0, y: Float64 = 0.0, z: Float64 = 0.0):
        self.x = x
        self.y = y
        self.z = z

    fn __add__(self, other: Vec3) -> Vec3:
        return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)

    fn __sub__(self, other: Vec3) -> Vec3:
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

    fn __mul__(self, scalar: Float64) -> Vec3:
        return Vec3(self.x * scalar, self.y * scalar, self.z * scalar)

    fn __truediv__(self, scalar: Float64) -> Vec3:
        return Vec3(self.x / scalar, self.y / scalar, self.z / scalar)

    fn __neg__(self) -> Vec3:
        return Vec3(-self.x, -self.y, -self.z)

    fn __eq__(self, other: Vec3) -> Bool:
        return self.x == other.x and self.y == other.y and self.z == other.z

    fn __ne__(self, other: Vec3) -> Bool:
        return not (self == other)

    fn dot(self, other: Vec3) -> Float64:
        """Computes the dot product of two vectors."""
        return self.x * other.x + self.y * other.y + self.z * other.z

    fn cross(self, other: Vec3) -> Vec3:
        """Computes the 3D cross product self x other."""
        return Vec3(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x
        )

    fn length_squared(self) -> Float64:
        """Returns the squared Euclidean length."""
        return self.dot(self)

    fn length(self) -> Float64:
        """Returns the Euclidean length (magnitude)."""
        return sqrt(self.length_squared())

    fn normalized(self) -> Vec3:
        """Returns the unit vector in the direction of self."""
        let l = self.length()
        if l == 0.0:
            return Vec3(0.0, 0.0, 0.0)
        return self / l

    fn distance_squared_to(self, other: Vec3) -> Float64:
        """Computes squared Euclidean distance to target vector."""
        let dx = self.x - other.x
        let dy = self.y - other.y
        let dz = self.z - other.z
        return dx * dx + dy * dy + dz * dz

    fn distance_to(self, other: Vec3) -> Float64:
        """Computes Euclidean distance to target vector."""
        return sqrt(self.distance_squared_to(other))
