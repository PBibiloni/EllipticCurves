from bgn.keygen import keygen

if __name__ == '__main__':
    for b in [2, 3, 4, 5, 6, 7, 8]:
        pk, sk = keygen(bits_of_security=b)
        print(f'For {b=}:')
        print(f'    {pk=}')
        print(f'    {sk=}')
