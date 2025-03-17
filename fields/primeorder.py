from fields.base import BaseFiniteField, FieldElement
from fields.utils import bezout_identity_Z


class FiniteFieldPrimeOrder(BaseFiniteField):
    def __init__(self, prime):
        super().__init__()
        self.prime = prime

    def from_coefficients(self, *coefficients):
        assert len(coefficients) == 1, 'Only one value can be converted to a field element.'
        value = coefficients[0]
        if isinstance(value, FieldElement) and value.field == self:
            return value
        elif isinstance(value, FieldElement):
            raise ValueError(f'Original field ({value.field}) do not match this field ({self})')
        return FieldElement(self, value)

    def field_elements(self):
        return [FieldElement(self, i) for i in range(self.prime)]

    def addition(self, element_A, element_B):
        return self.from_coefficients((element_A.coefficients[0] + element_B.coefficients[0]) % self.prime)
    
    def negation(self, element):
        return self.from_coefficients((-element.coefficients[0]) % self.prime)
    
    def multiplication(self, element_A, element_B):
        return self.from_coefficients((element_A.coefficients[0] * element_B.coefficients[0]) % self.prime)
    
    def inverse(self, element):
        alpha, beta, gcd = bezout_identity_Z(element.coefficients[0], self.prime)
        return self.from_coefficients(alpha % self.prime)

    def __eq__(self, other):
        return self.prime == other.prime

    def __hash__(self):
        return hash(self.prime)

    def __repr__(self):
        return f'F_{self.prime}'
