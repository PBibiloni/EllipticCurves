def bezout_identity_Z(a, b):
    # Returns alpha, beta, gcd for a * alpha + b * beta = gcd(a, b)
    if b == 0:
        return 1, 0, a
    else:
        alpha, beta, gcd = bezout_identity_Z(b, a % b)
        return beta, alpha - (a // b) * beta, gcd
