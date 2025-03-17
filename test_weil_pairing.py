from curves.curves import WeierstrassCurve
from curves.weil import weil_pairing, f
from fields.primeorder import FiniteFieldPrimeOrder

if __name__ == '__main__':
    # Example in Silverman
    field = FiniteFieldPrimeOrder(prime=631)
    curve = WeierstrassCurve(a=30, b=34, field=field)
    P = curve.point(36, 60)
    Q = curve.point(121, 387)
    S = curve.point(0, 36)
    assert f(curve, P, Q+S, n=5) == curve.field.from_coefficients(103)
    w = weil_pairing(P, Q, n=5, S=S)
    assert w == curve.field.from_coefficients(242), f'weil_pairing(curve, P, Q, n=7, S=S) = {w} (but expected to be 242)'

    # Example in 2011 - Aftuck - The Weil pairing on elliptic curves and its cryptographic applications
    field = FiniteFieldPrimeOrder(prime=1009)
    curve = WeierstrassCurve(a=37, b=0, field=field)
    P = curve.point(8, 703)
    Q = curve.point(49, 20)
    S = curve.point(0, 0)
    w = weil_pairing(P, Q, n=7, S=S)
    assert w == curve.field.from_coefficients(105), f'weil_pairing(curve, P, Q, n=7, S=S) = {w} (but expected to be 105)'
