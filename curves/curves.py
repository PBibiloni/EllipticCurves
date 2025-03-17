from .finitegroup import FiniteGroup


class WeierstrassCurve:
    """ Curve in normal form y^2 = x^3 + ax + b. """
    def __init__(self, a, b, field):
        self.a = field.from_coefficients(a)
        self.b = field.from_coefficients(b)
        self.field = field
        assert field.prime > 3, 'Addition operation not implemented for F_2 (specifically slope between points)'

    def point(self, x, y, check_belongs=True):
        P = WeierstrassCurve.Point(self, x, y)
        if check_belongs and not self.check_belongs(P):
            raise ValueError(f'Point {P} does not belong to curve {self}')
        return P

    def neutral_element(self):
        return WeierstrassCurve.Point(self, None, None)

    def check_belongs(self, point):
        # Check point at infinity first
        if point.x is None and point.y is None:
            return True
        elif point.x is None or point.y is None:
            return False
        else:
            return point.y ** 2 == point.x ** 3 + self.a * point.x + self.b

    def slope(self, point_A, point_B):
        """ Slope of line containing both points. """
        if point_A.x is None or point_B is None:    # Either of them is the point at infinity
            return None
        if point_A.x == point_B.x and point_A.y == point_B.y == 0:  # Doubling point with y = 0 -> infty
            return None
        if point_A.x == point_B.x and point_A.y == point_B.y != 0:  # Doubling point with y != 0 -> tangent to E
            return (3 * point_A.x ** 2 + self.a) / (2 * point_A.y)
        if point_A.x == point_B.x and point_A.y != point_B.y:       # Pair of points (x, y) and (x, -y).
            return None
        return (point_B.y - point_A.y) / (point_B.x - point_A.x)    # General case for x1 != x2

    def point_generator(self):
        yield WeierstrassCurve.Point(self, None, None)
        for x in self.field.field_elements():
            y_squared = x ** 3 + self.a * x + self.b
            for y in self.field.square_root(y_squared):
                yield WeierstrassCurve.Point(self, x, y)

    def get_all_points(self):
        return list(self.point_generator())

    def as_group(self):
        return FiniteGroup(
            group_elements=self.get_all_points(),
            identity_element=WeierstrassCurve.Point(self, None, None),
            operation=lambda P, Q: P + Q,
            inverse=lambda P: -P,
        )

    def __repr__(self):
        return f'y^2 = x^3 + ({self.a})·x + ({self.b}) (for x, y in {self.field})'

    class Point:
        def __init__(self, curve, x, y):
            if x is None and y is None:
                # Point at infinity
                self.x = None
                self.y = None
            elif x is None or y is None:
                raise ValueError('Both x and y must be None or not None')
            else:
                self.x = curve.field.from_coefficients(x)
                self.y = curve.field.from_coefficients(y)
            self.curve = curve

        def __add__(self, other):
            assert isinstance(other, self.__class__)
            if self.x is None:
                return other
            if other.x is None:
                return self
            if self.x == other.x and self.y != other.y:
                return self.__class__(self.curve, None, None)
            lmbda = self.curve.slope(self, other)
            if lmbda is None:   # Doubling a point with lmbd = infty
                return self.__class__(self.curve, None, None)
            x = lmbda ** 2 - self.x - other.x
            y = lmbda * (self.x - x) - self.y
            return self.__class__(self.curve, x, y)

        def __neg__(self):
            return self.__class__(self.curve, self.x, -self.y)

        def __sub__(self, other):
            assert isinstance(other, self.__class__)
            return self + (-other)

        def __mul__(self, other):
            assert isinstance(other, int)
            base = self if other >= 0 else -self
            times = abs(other)
            # Double-and-add algorithm
            value = self.curve.neutral_element()
            for bit in bin(times)[2:]:    # Skip '0b' prefix
                value = value + value
                if bit == '1':
                    value = value + base
            return value

        def __rmul__(self, other):
            return self.__mul__(other)

        def __repr__(self):
            return f'({self.x}, {self.y})'

        def __eq__(self, other):
            assert isinstance(other, self.__class__)
            if self.curve != other.curve:
                return False
            if self.x is None and other.x is None and self.y is None and other.y is None:
                return True
            return self.x == other.x and self.y == other.y

        def __hash__(self):
            return hash((self.x, self.y, self.curve))
