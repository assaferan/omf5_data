hecke_eran = eval(open("hecke_all.dat").read())
hecke_gustavo = eval(open("hecke_3_0.dat").read())

def first_prime(N):
    p = 2
    while (N % p == 0):
        p = next_prime(p)
    return p

def check(N,d):
    if (hecke_eran[N][first_prime(N)][d] == []):
        if (hecke_gustavo[N] == {}):
            return True
        return Matrix(hecke_eran[N][first_prime(N)][d]).charpoly() == Matrix(hecke_gustavo[N][first_prime(N)][d]).charpoly()

def verify(N):
    if (N == 0):
        return True
    for d in divisors(N):
        if is_squarefree(d):
            if (not check(N,d)):
                return False
    return True

all([verify(N) for N in range(1000)])
