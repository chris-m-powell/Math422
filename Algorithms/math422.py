import random


def naturals(n):
    """Print first n natural numbers"""
    i = 0
    while i < n+1:
        print(i)
        i += 1


def triangular(n):
    """Compute nth triangular number using while"""
    sum, i = 0, 0
    while i < n:
        i += 1
        sum = sum + i
    return(sum)


def special_sum(n):
    """Compute sum from 1 to n, not multiples of 3, using while."""
    sum, i = 0, 0
    while (i < n):
        i += 1
        if i % 3 != 0:
            sum = sum + i
    return(sum)


def dumber_multiply(a, b):
    """Multiply positive integers a and b using for."""
    p = 0
    for i in range(1, b + 1):
        if b != 0:
            p, b = p + a, b - 1
    return p


def factorial(n):
    """Compute n factorial using while."""
    p, i = 1, 0
    while i < n:
        i += 1
        p = p * i
    return(p)


def fibonacci(n):
    """Compute nth fibonacci number using while."""
    if n == 1:
        return 0
    a, b = 0, 1
    i = 0
    while i < n:
        i += 1
        a, b = b, a + b
    return(a)


def pisano(n):
    """Compute the sum of the nth and 2nth fibonacci numbers."""
    return fibonacci(n) + fibonacci(2*n)


def sod(x):
    """Returns sum of digits"""
    s = 0
    while x > 0:
        s = s + (x % 10)
        x = x // 10
    return s


def collatz(n):
    """Compute Length of Termination and Terminating value"""
    A, i = [], 0
    while n not in A:
        A, i = A + [n], i + 1
        if n % 2 == 0:
            n = n // 2
        else:
            n = (3 * n) + 1
    return A, i, A[i-1]


def collatz_table(k):
    """Print table for Length of Termation and Terminating values"""
    for n in range(1, k+1):
        print(n, collatz(n))


def prime_factors(n):
    """Return prime factors of positive integer n"""
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    return factors


def gcd(a, b):
    """Return gcd(a,b)"""
    while b:
        a, b = b, a % b
    return a


def xgcd(a, b):
    """Return (g, x, y) such that a*x + b*y = g = gcd(a, b)"""
    if b == 0:
        return a, 1, 0
    x, g, v, w = 1, a, 0, b
    while w != 0:
        x, g, v, w = v, w, x - (g // w) * v, g % w
    x = x % (b // g)
    return g, x, (g - (a * x)) // b


def linear_cong(a, c, m):
    """Return {x in Z | ax = c (mod m)}"""
    g = xgcd(a, m)[0]
    S = []
    if c % g != 0:
        return S  # solution-set S is empty
    else:
        x = (c * xgcd(a, m)[1]) // g
        for i in range(0, g):
            S = S + [(x + i * (m // g)) % m]
    return S  # return g incongruent solutions


def feb21(n):
    """Return day of week of Feb 21 for year 2000 < n < 2100"""
    days = ['Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat']
    return days[((n % 2000) % 7 + (1 + ((n - 1) % 2000) // 4) % 7) % 7]


def baseb(n, b):
    """Given positive integers n and b, return n in base b"""
    d = []
    while n:
        d += [n % b]
        n //= b
    return d[::-1]  # Return reversed list


def totient(n):
    """Given positive integer n, return totient(n)"""
    t = 0
    for i in range(1, n):
        if gcd(i, n) == 1:
            t += 1
    return t


def crt(a, b, m, n):
    """Given integers a,b w/ gcd(m,n)=1, return unique x cong a,b (mod m,n)"""
    (g, r, s) = xgcd(m, n)
    return ((a * s * n) + (b * r * m)) % (m * n)


def expmod(a, k, m):
    """compute a^k mod m"""
    b = 1
    while k:
        if k % 2 == 1:
            b = (b * a) % m
        a, k = (a ** 2) % m, k // 2
    return b


def modroot(k, b, m):
    """return x such that x^k cong b mod m if gcd(b,m)=1, gcd(k, totient(m))"""
    if gcd(b, m) == 1:
        (g, u, v) = xgcd(k, totient(m))
        if g == 1:
            return expmod(b, u, m)


def modrootpq(k, b, m, p, q):
    """given p,q, return x such that x^k cong b mod m=pq if gcd(b,m)=1, gcd(k, totient(m))"""
    if gcd(b, m) == 1:
        (g, u, v) = xgcd(k, totient(p)*totient(q))
        if g == 1:
            return expmod(b, u, m)


def probablyprime(n):
    """Returns whether integer n is likely prime or composite"""
    A = []
    for i in range(10):
        A += [random.randint(2, n-1)]
    for i in A:
        if expmod(A[i], n-1, n) % n != 1:
            return str(n) + " is composite"
    return str(n) + " is prime"


def rsa_decrypt(k, B, m, p, q):
    "Decrypt RSA ciphertext B given k, m, p, q"
    a = ""
    for b in B:
        a += str(modrootpq(k, b, m, p, q))
    A, a = [chr(int(a)+54) for a in [a[i:i + 2] for i in range(0, len(a), 2)]], ""
    for i in A:
        a += i
    return a


def residue(n, p):
    """Return each a ≡ b^n mod p for some integer b and modulus p; n=2,3"""
    R = []
    for i in range(p):
        R += [(i**n) % p]
    return set(R)


def invmod(a, m):
    """ """
    g, x, y = xgcd(a, m)
    if g != 1:
        return "Modular inverse does not exit"
    else:
        return x % m


def CRT(L):
    """Given (a1,m1),...,(an,mn) with gcd(mi,mj)=1, return x such that x = ai (mod mi) for all i"""
    x, M = 0, 1
    for i in range(len(L)):
        M *= L[i][1]
    for i in range(len(L)):
        x += L[i][0] * (M // L[i][1]) * (xgcd(M//L[i][1], L[i][1])[1] % L[i][1])
    return x % M

L = [[4, 5], [2, 7]]
print(CRT(L))

def hw11(a, p):
    """Return solution pair to x^2=a (mod p) for prime p=5 (mod 8), where a is a QR (mod p)"""
    if a not in residue(2, p):
        return str(a)+" is not a quadratic residue modulo "+str(p)+"." 
    if expmod(a, ((p - 1) // 4), p) == 1:
        x = expmod(a, ((p + 3) // 8), p)
        return x, -x % p
    else:
        x = ((2 * a) * expmod(4 * a, ((p - 5) // 8), p)) % p
        return x, -x % p


def descent(A, B, p):
    "return integers (A,B) s.t. A^2+B^2=p; p prime congruent to 1 mod 4, via Fermat's descent method"
    M = ((A ** 2) + (B ** 2)) // p
    if ((A ** 2) + (B ** 2)) % p != 0:
        return "No solution exists."
    else:
        while M > 1:
            u = A % M
            while u > (M // 2):
                u = u - M
            v = B % M
            while v > (M // 2):
                v = v - M
            print(u,v)
            print(A,B)
            A, B = ((u * A) + (v * B)) // M, ((v * A) - (u * B)) // M
            M = ((A ** 2) + (B ** 2)) // p
        return A, B


def ord(a, p):
    for e in range(1, p):
        if (a ** e) % p == 1:
            return e


# print(descent(259, 1, 1973))
