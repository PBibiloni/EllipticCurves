import random

import sympy

from curves.curves import WeierstrassCurve
from curves.finitegroup import FiniteGroup
from curves.weil import weil_pairing
from fields.utils import bezout_identity_Z
from fields.extension import FiniteFieldExtension3thPrimitiveRoot
from fields.primeorder import FiniteFieldPrimeOrder


def keygen(bits_of_security):
    # Generate q1, q2.
    q1 = sympy.nextprime(2 ** (bits_of_security + random.random()))
    q2 = sympy.nextprime(2 ** (bits_of_security + random.random()))
    if q1 == q2:
        q2 = sympy.nextprime(q2)
    n = q1 * q2
    # Find p, compute G = curve_Fp
    p = find_smallest_p(q1, q2)
    field_Fp = FiniteFieldPrimeOrder(prime=p)
    # Find g, u random generators in G.
    curve_Fp = WeierstrassCurve(a=0, b=1, field=field_Fp)
    g = find_point_order_n(curve_Fp, q1, q2, p)
    u = g * random.randint(1, q1*q2-1)
    netrual_element = curve_Fp.neutral_element()
    while q1 * u == netrual_element or q2 * u == netrual_element:
        u = random.randint(1, p-1) * g
    # Find h
    h = q2 * u
    # Define modified weil pairing
    field_Fp2 = FiniteFieldExtension3thPrimitiveRoot(prime=p)
    curve_Fp2 = WeierstrassCurve(a=0, b=1, field=field_Fp2)
    S = curve_Fp2.point(x=field_Fp2.from_coefficients(0), y=field_Fp2.from_coefficients(1))
    def e(P, Q):
        return modified_weil_pairing(curve_Fp2, P, Q, n, S=S)
    # Compute G1 (from weil pairing)
    G1_gen = e(g, u)
    G1_elements = [field_Fp2.from_coefficients(1)]
    for i in range(1, n):
        G1_elements.append(G1_gen ** i)
    G1 = FiniteGroup(
        group_elements=G1_elements,
        identity_element=1,
        operation=lambda a, b: a * b,
        inverse=lambda a: 1 / a,
    )
    # Return private key and secret key
    pk = (n, curve_Fp, G1, e, g, h)
    sk = (q1, )
    return pk, sk


def modified_weil_pairing(curve_Fp2, P, Q, n, S=None):
    assert curve_Fp2.a == 0 and curve_Fp2.b == 1, 'Implementation only for y^2 = x^3 + 1 due to distortion map.'
    Fp2 = curve_Fp2.field
    # Cast points into Fp2
    P_Fp2 = curve_Fp2.point(
        x=Fp2.from_coefficients(P.x.coefficients[0]) if P.x is not None else None,
        y=Fp2.from_coefficients(P.y.coefficients[0]) if P.y is not None else None
    )
    Q_Fp2_distorted = curve_Fp2.point(
        x=Fp2.from_coefficients(0, Q.x.coefficients[0]) if Q.x is not None else None,    # Applying distortion map: (x, y) -> (x · alpha, y)
        y=Fp2.from_coefficients(Q.y.coefficients[0]) if Q.y is not None else None
    )
    # Compute
    return weil_pairing(P_Fp2, Q_Fp2_distorted, n, S)


def find_smallest_p(q_1, q_2):
    n = q_1 * q_2
    k = 0
    while True:
        k += 1
        p = n*k - 1
        if p < 5:
            pass
        elif not sympy.isprime(p):
            pass
        elif p % 3 != 2:
            pass
        else:
            return p


def find_point_order_n(curve, q1, q2, p):
    assert curve.field.prime % 3 == 2, 'Implementation only for p = 2 (mod 3).'
    assert curve.a == 0 and curve.b == 1, 'Implementation only for y^2 = x^3 + 1.'
    neutral_element = curve.neutral_element()
    for y_int in range(p):
        y = curve.field.from_coefficients(y_int)
        x_cube = y ** 2 - 1
        # x^3 can be computed as (x^3)^e, where e = 1/3 (mod p-1)
        alpha, _, _ = bezout_identity_Z(3, p - 1)
        x = x_cube ** (alpha % (p-1))
        P = curve.point(x=x, y=y)
        if (q1*q2) * P == neutral_element and q1 * P != neutral_element and q2 * P != neutral_element:
            return P
    raise ValueError(f'Could not find a point of order {q1*q2} in curve {curve}.')
