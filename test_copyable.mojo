@fieldwise_init
struct Point(ImplicitlyCopyable):
    var x: Float64
    var y: Float64

    def __init__(out self):
        self.x = 0.0
        self.y = 0.0

def main():
    var p1 = Point()
    var p2 = p1
    print("p2 x:", p2.x)
