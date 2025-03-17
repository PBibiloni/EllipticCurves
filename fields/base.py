from curves.finitegroup import FiniteGroup


class BaseFiniteField:

    def __init__(self):
        self._square_roots = None

    def from_coefficients(self, *coefficients):
        raise NotImplementedError

    def field_elements(self):
        raise NotImplementedError

    def addition(self, element_A, element_B):
        raise NotImplementedError

    def negation(self, element):
        raise NotImplementedError

    def multiplication(self, element_A, element_B):
        raise NotImplementedError

    def inverse(self, element):
        raise NotImplementedError

    def square_root(self, n):
        if self._square_roots is None:
            # Precompute all square roots in field
            self._square_roots = {}
            for element in self.field_elements():
                squared = element * element
                self._square_roots[squared] = self._square_roots.get(squared, []) + [element]
        return self._square_roots.get(self.from_coefficients(n), [])

    def nth_roots(self, n):
        return [
            element
            for element in self.field_elements()
            if (element ** n) == self.from_coefficients(1)
        ]

    def primitive_nth_roots(self, n):
        primitive_roots = self.nth_roots(n)
        orders = {}
        for root in primitive_roots:
            orders[root] = 1
            result = root
            while result != self.from_coefficients(1):
                result = result * root
                orders[root] += 1
        return [root for root, order in orders.items() if order == n]

    def group_nth_roots(self, n):
        return FiniteGroup(
            group_elements=self.nth_roots(n),
            identity_element=self.from_coefficients(1),
            operation=lambda a, b: a * b,
            inverse=lambda a: 1 / a,
        )


class FieldElement:
    """ Implement pair a + b · alpha, where alpha is a 3th primitive root of unity. """
    def __init__(self, field, *coefficients):
        if any(not isinstance(value, int) for value in coefficients):
            raise ValueError(f'All coefficients must be integer (but {coefficients=}).')
        self.coefficients = coefficients
        self.field = field

    def _sanitize_other(self, other):
        if isinstance(other, self.__class__) and self.field == other.field:
            return other
        elif isinstance(other, self.__class__) and self.field != other.field:
            raise ValueError('Fields do not match')
        return self.field.from_coefficients(other)

    def __add__(self, other):
        return self.field.addition(self, self._sanitize_other(other))

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self.field.addition(self, self.field.negation(self._sanitize_other(other)))

    def __neg__(self):
        return self.field.negation(self)

    def __rsub__(self, other):
        return self._sanitize_other(other) - self

    def __mul__(self, other):
        return self.field.multiplication(self, self._sanitize_other(other))

    def __rmul__(self, other):
        return self * other

    def __pow__(self, other):
        assert isinstance(other, int)
        base = self if other >= 0 else self.field.inverse(self)
        exp = abs(other)
        # double-and-add exponentiation:
        result = 1
        for bit in bin(exp)[2:]:    # Avoid '0b' prefix
            result = result * result
            if bit == '1':
                result = result * base
        return result

    def __truediv__(self, other):
        return self.field.multiplication(self, self.field.inverse(self._sanitize_other(other)))

    def __rtruediv__(self, other):
        return self._sanitize_other(other) / self

    def __eq__(self, other):
        try:
            other = self._sanitize_other(other)
        except ValueError:
            return False
        return self.coefficients == other.coefficients and self.field == other.field

    def __repr__(self):
        n = len(self.coefficients)
        if n == 1:
            return f'{self.coefficients[0]}'
        elif n == 2:
            return f'({self.coefficients[0]}+{self.coefficients[1]}·α)'
        else:
            return f'({self.coefficients[0]}+' + '+'.join(f'{self.coefficients[1]}·α_{i})' for i in range(1, n)) + ")"

    def __hash__(self):
        return hash((self.coefficients, self.field))
