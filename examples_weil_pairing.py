import math
import random

import sympy

from curves.curves import WeierstrassCurve
from curves.weil import weil_pairing
from fields.extension import FiniteFieldExtension3thPrimitiveRoot
from fields.primeorder import FiniteFieldPrimeOrder


if __name__ == '__main__':
    p = 3
    while True:
        p = sympy.nextprime(p)
        print(f'\n\n\nChecking {p=}.')
        if not sympy.isprime(p):
            print(f'   {p=} is not prime.')
            continue
        if p < 5:
            print(f'   {p=} is too small.')
            continue
        if not p % 3 == 2:
            print(f'   {p=} = 2 (mod 3), therefore F_{p=} already has primitive 3-roots: we will not find E[n] in E(F_p(alpha)).')
            continue

        # Define n at "random" so that it is square free
        n = math.prod(factor for factor in sympy.primefactors(p+1) if random.random() < 0.5)
        print(f'    Let\'s consider {n=}.')

        field_Fp = FiniteFieldPrimeOrder(prime=p)
        field_Fp2 = FiniteFieldExtension3thPrimitiveRoot(prime=p)
        print(f'    {n}-th roots of {field_Fp} --> {field_Fp.nth_roots(n)}')
        print(f'    {n}-th roots of {field_Fp2} --> {field_Fp2.nth_roots(n)}')

        print()
        # Study E(F_p)
        curve_Fp = WeierstrassCurve(a=0, b=1, field=field_Fp)
        E_Fp = curve_Fp.as_group()
        E_Fp_ntorsion = E_Fp.n_torsion_subgroup(n)
        print(f'    Order(E(F_p)) = {E_Fp.order_of_group()}.')
        print(f'    E(F_p) <-> ' + " x ".join(f'Z_{{{idx}}}' for idx in E_Fp.classify_finite_abelian_group()))
        print(f'    -> {n}-torsion subgroup: {E_Fp_ntorsion.group_elements}.')

        # Compute Weil pairings on E(F_p)
        print(f'    Weil pairings in E(F_p) {n}-torsion group:')
        for _ in range(6):
            P = random.choice(E_Fp_ntorsion.group_elements)
            Q = random.choice(E_Fp_ntorsion.group_elements)
            w = weil_pairing(P, Q, n)
            print(f'     -> e_{n}({P=}, {Q=}) = {w}.')
            print(f'        -> {w}^{n} = {w ** n}.')

        print()
        # Study E(F_p(alpha)), for alpha as 3-th primitive root unit for F_p
        curve_Fp2 = WeierstrassCurve(a=0, b=1, field=field_Fp2)
        E_Fp2 = curve_Fp2.as_group()
        E_Fp2_ntorsion = E_Fp2.n_torsion_subgroup(n)
        print(f'    Order(E(F_p(alpha))) = {E_Fp2.order_of_group()}.')
        print(f'    E(F_p(alpha)) <-> ' + " x ".join(f'Z_{{{idx}}}' for idx in E_Fp2.classify_finite_abelian_group()))
        print(f'    -> {n}-torsion subgroup: {E_Fp2_ntorsion.group_elements}.')

        # Compute Weil pairings on E(F_p)
        print(f'    Weil pairings in E(F_p(alpha)) {n}-torsion group:')
        for _ in range(6):
            P = random.choice(E_Fp2_ntorsion.group_elements)
            Q = random.choice(E_Fp2_ntorsion.group_elements)
            w = weil_pairing(P, Q, n)
            print(f'     -> e_{n}({P=}, {Q=}) = {w}.')
            print(f'        -> ({w})^{n} = {w ** n}.')
