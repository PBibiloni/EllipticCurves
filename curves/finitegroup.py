import sympy
from sympy import divisors


class FiniteGroup:
    def __init__(self, group_elements, identity_element, operation, inverse, order_by_element=None):
        self.group_elements = group_elements
        self.identity_element = identity_element
        self.operation = operation
        self.inverse = inverse
        self._order_of_group = len(group_elements)
        self._order_of_group_divisors = divisors(self._order_of_group)  # [d1, d2, ...]
        self._order_of_group_divisors.sort()
        self._order_by_element = {} if order_by_element is None else order_by_element

    def order_of_group(self):
        if self._order_of_group is None:
            self._order_of_group = len(self.group_elements())
        return self._order_of_group

    def get_nontrivial_element(self):
        if self.order_of_group() == 1:
            raise ValueError('Group has only one element')
        elements = self.group_elements()
        elements.sort(key=lambda x: self.order(x))
        return elements[1]

    def order(self, element):
        if element not in self._order_by_element:
            for d in self._order_of_group_divisors:
                if self.exponentiation(element, d) == self.identity_element:
                    self._order_by_element[element] = d
                    return d
            raise ValueError(f'Element {element} has order not divided by |G| = {self.order_of_group()}.')
        return self._order_by_element[element]

    def exponentiation(self, element, exp):
        """ Repeatedly apply the group operation of one element with itself. """
        assert isinstance(exp, int)
        base = element if exp >= 0 else self.inverse(element)
        exp = abs(exp)
        # Double-and-add algorithm
        result = self.identity_element
        for bit in bin(exp)[2:]:    # Avoid '0b' prefix and left-most bit (always set to 1).
            result = self.operation(result, result)
            if bit == '1':
                result = self.operation(result, base)
        return result

    def cyclic_subgroup(self, generator):
        elements = [generator]
        while elements[-1] != self.identity_element:
            elements.append(self.operation(elements[-1], generator))
        return FiniteGroup([elements[-1]] + elements[:-1], self.identity_element, self.operation, self.inverse)

    def subgroup_generated_by(self, elements):
        subgroup_elements = {self.identity_element}

        for e in elements:
            factor = self.identity_element
            cosets = []
            for _ in range(self.order(e)):
                factor = self.operation(factor, e)
                cosets.append(set(self.operation(factor, element) for element in subgroup_elements))
            for _ in range(self.order(e)):
                subgroup_elements = subgroup_elements | cosets.pop(0)

        return FiniteGroup(
            group_elements=subgroup_elements,
            identity_element=self.identity_element,
            operation=self.operation,
            inverse=self.inverse
        )

    def n_torsion_subgroup(self, n):
        """ Returns the elements of order n. """
        ntorsion_elements = [element for element in self.group_elements if n % self.order(element) == 0]
        return FiniteGroup(
            group_elements=ntorsion_elements,
            identity_element=self.identity_element,
            operation=self.operation,
            inverse=self.inverse,
            order_by_element={element: self._order_by_element[element] for element in ntorsion_elements}
        )

    def quotient_group(self, subgroup):
        remaining_elements = self.group_elements
        cosets = {}     # Dictionary {representative: coset}
        element_to_representative = {}
        while len(remaining_elements) > 0:
            coset_representative = remaining_elements.pop(0)
            coset = tuple(
                self.operation(coset_representative, element)
                for element in subgroup.group_elements
            )
            element_to_representative = element_to_representative | {   # Add to the dictionary
                element: coset_representative
                for element in coset
            }
            cosets[coset_representative] = coset

            remaining_elements = [
                element
                for element in remaining_elements
                if element not in coset
            ]

        return FiniteGroup(
            group_elements=list(cosets.values()),
            identity_element=cosets[element_to_representative[self.identity_element]],
            operation=lambda c1, c2: cosets[element_to_representative[self.operation(c1[0], c2[0])]],
            inverse=lambda c1: cosets[element_to_representative[self.inverse(c1[0])]]
        )

    def classify_finite_abelian_group(self):
        subgroups = []
        for prime, exp in sympy.factorint(self.order_of_group()).items():
            orders = set(prime ** i for i in range(exp + 1))
            subgroups += [
                FiniteGroup(
                    group_elements=[element for element in self.group_elements if self.order(element) in orders],
                    identity_element=self.identity_element,
                    operation=self.operation,
                    inverse=self.inverse,
                    order_by_element=self._order_by_element
                )
            ]
        idcs = []
        for subgroup in subgroups:
            if subgroup.order_of_group() == 1:
                pass
            elif sympy.isprime(subgroup.order_of_group()):
                idcs += [subgroup.order_of_group()]
            else:
                current_group = subgroup
                while current_group.order_of_group() > 1:
                    order_by_element = {e: current_group.order(e) for e in current_group.group_elements}
                    orders = set(order_by_element.values())
                    orders.remove(1)
                    p = min(orders)
                    max_power_p = max(orders)
                    if max_power_p == current_group.order_of_group():
                        idcs += [max_power_p]
                        break
                    elif max_power_p * p == current_group.order_of_group():
                        idcs += [max_power_p, p]
                        break
                    else:
                        generator = (e for e, o in order_by_element.items() if o == max_power_p).__next__()
                        subgroup = current_group.cyclic_subgroup(generator)
                        current_group = current_group.quotient_group(subgroup)
                        idcs += [subgroup.order_of_group()]
        return idcs or [1]

    def __repr__(self):
        return f'Group (order {self.order_of_group()})'
