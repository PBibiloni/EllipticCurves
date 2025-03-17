from fields.base import BaseFiniteField
from fields.utils import bezout_identity_Z
from fields.primeorder import FieldElement


class FiniteFieldExtension3thPrimitiveRoot(BaseFiniteField):
    def __init__(self, prime):
        super().__init__()
        self.prime = prime
        assert prime % 3 == 2, 'This implementation is only defined for primes p such that p % 3 == 2'

    def from_coefficients(self, *coefficients):
        if len(coefficients) == 1 and isinstance(coefficients[0], FieldElement):
            if coefficients[0].field != self:
                raise ValueError(f'Original field ({coefficients[0].field}) do not match this field ({self})')
            return coefficients[0]
        elif len(coefficients) == 1 and isinstance(coefficients[0], FieldElement):
            if coefficients[0].field.prime != self.prime:
                raise ValueError(f'Original field ({coefficients[0].field}) do not match this field ({self})')
            return FieldElement(self, coefficients[0].value, 0)
        elif len(coefficients) == 1 and isinstance(coefficients[0], int):
            return FieldElement(self, coefficients[0], 0)
        elif len(coefficients) == 2 and all(isinstance(value, int) for value in coefficients):
            return FieldElement(self, coefficients[0], coefficients[1])
        raise ValueError(f'Not implemented for {coefficients=}')

    def field_elements(self):
        return [FieldElement(self, i, j)
                for i in range(self.prime)
                for j in range(self.prime)]

    def addition(self, element_A, element_B):
        return self.from_coefficients(
            (element_A.coefficients[0] + element_B.coefficients[0]) % self.prime,
            (element_A.coefficients[1] + element_B.coefficients[1]) % self.prime,
        )

    def negation(self, element):
        return self.from_coefficients(
            (-element.coefficients[0]) % self.prime,
            (-element.coefficients[1]) % self.prime,
        )

    def multiplication(self, element_A, element_B):
        # Take into account that alpha^2 = - alpha - 1
        return self.from_coefficients(
            (element_A.coefficients[0] * element_B.coefficients[0] - element_A.coefficients[1] * element_B.coefficients[1]) % self.prime,
            (element_A.coefficients[0] * element_B.coefficients[1] + element_A.coefficients[1] * element_B.coefficients[0] - element_A.coefficients[1] * element_B.coefficients[1]) % self.prime,
        )

    def inverse(self, element):
        # Compute 1/element as (a+b*alpha^2)/ (a^2 - ab + b^2), where alpha^2 = -alpha - 1
        norm = element.coefficients[0] ** 2 - element.coefficients[0] * element.coefficients[1] + element.coefficients[1] ** 2
        inv_norm, _, gcd = bezout_identity_Z(norm, self.prime)
        assert gcd == 1, f'{norm} is not invertible: could not invert {element}.'
        return self.from_coefficients(
            (element.coefficients[0] - element.coefficients[1]) * inv_norm % self.prime,
            (-element.coefficients[1]) * inv_norm % self.prime,
        )

    def __eq__(self, other):
        return self.prime == other.prime

    def __hash__(self):
        return hash(self.prime)

    def __repr__(self):
        return f'F_{self.prime}(α), st. α^3 = 1 and α != 1'
