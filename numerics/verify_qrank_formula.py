# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Closed form for the Q-rank of the ring frequencies: rank(N) = phi(2N)/2.

The subtorus dimension d_Q(N) = dim_Q span{2 sin(pi r/N): r=1..N-1} equals phi(2N)/2, the dimension
of the "imaginary" (minus-eigenspace of complex conjugation) part of Q(zeta_{2N}), since
2 i sin(pi r/N) = zeta_{2N}^r - zeta_{2N}^{-r}. Consequently the full-rank case --- where Bohr gives
A_N = U_N exactly --- is rank(N) = floor(N/2), i.e.

    phi(2N)/2 = floor(N/2)   <=>   N is prime or a power of two,

exactly the classes of Theorems 1' and 1''. This pins the whole composite map with one formula.
Verified by an exact (mod large prime) rank against phi(2N)/2 for N = 3..64.
"""

from __future__ import annotations


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def totient(n: int) -> int:
    result, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            result -= result // p
        p += 1
    if m > 1:
        result -= result // m
    return result


def poly_div_exact(a: list[int], b: list[int]) -> list[int]:
    a = a[:]
    db = len(b) - 1
    q = [0] * (len(a) - db)
    for i in range(len(a) - 1, db - 1, -1):
        c = a[i]
        if c:
            q[i - db] = c
            for j in range(len(b)):
                a[i - db + j] -= c * b[j]
    return q


def cyclotomic(n: int, cache: dict[int, list[int]]) -> list[int]:
    if n in cache:
        return cache[n]
    num = [-1] + [0] * (n - 1) + [1]
    for d in divisors(n):
        if d < n:
            num = poly_div_exact(num, cyclotomic(d, cache))
    cache[n] = num
    return num


def power_mod(e: int, phi: list[int]) -> list[int]:
    deg = len(phi) - 1
    p = [0] * max(e + 1, deg)
    p[e] = 1
    for i in range(len(p) - 1, deg - 1, -1):
        c = p[i]
        if c:
            for j in range(deg + 1):
                p[i - deg + j] -= c * phi[j]
    return p[:deg]


def ring_rank(n: int, q: int = (1 << 61) - 1) -> int:
    cache: dict[int, list[int]] = {}
    phi = cyclotomic(2 * n, cache)
    rows = []
    for r in range(1, n):
        hi = power_mod(r % (2 * n), phi)
        lo = power_mod((2 * n - r) % (2 * n), phi)
        rows.append([(h - lo_) % q for h, lo_ in zip(hi, lo, strict=True)])
    mat = rows
    pr = 0
    for col in range(len(mat[0])):
        piv = next((rr for rr in range(pr, len(mat)) if mat[rr][col] % q), None)
        if piv is None:
            continue
        mat[pr], mat[piv] = mat[piv], mat[pr]
        inv = pow(mat[pr][col], q - 2, q)
        mat[pr] = [(x * inv) % q for x in mat[pr]]
        for rr in range(len(mat)):
            if rr != pr and mat[rr][col]:
                f = mat[rr][col]
                mat[rr] = [(mat[rr][k] - f * mat[pr][k]) % q for k in range(len(mat[0]))]
        pr += 1
    return pr


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def main() -> int:
    print("=" * 62)
    print("Q-rank of {2 sin(pi r/N)} = phi(2N)/2  (subtorus dimension)")
    print("=" * 62)
    ok = True
    full_iff = True
    mism = []
    for n in range(3, 65):
        rk = ring_rank(n)
        formula = totient(2 * n) // 2
        if rk != formula:
            mism.append((n, rk, formula))
            ok = False
        full = rk == n // 2
        expect_full = is_prime(n) or (n & (n - 1) == 0)
        if full != expect_full:
            full_iff = False
    print(f"\n  rank(N) == phi(2N)/2 for all N in 3..64 : {'PASS' if ok else f'FAIL {mism}'}")
    print(f"  rank(N) == floor(N/2)  <=>  N prime or 2^m : {'PASS' if full_iff else 'FAIL'}")
    print("\n  sample:")
    print(f"  {'N':>4} {'rank':>5} {'phi(2N)/2':>10} {'floor(N/2)':>11} {'full?':>6} {'class':>10}")
    for n in (5, 7, 8, 12, 15, 16, 24, 30, 36):
        rk = ring_rank(n)
        cls = "prime" if is_prime(n) else ("2^m" if (n & (n - 1)) == 0 else "composite")
        print(f"  {n:>4} {rk:>5} {totient(2 * n) // 2:>10} {n // 2:>11} "
              f"{rk == n // 2!s:>6} {cls:>10}")

    print("\n" + "=" * 62)
    print("RESULT:", "Q-RANK FORMULA VERIFIED" if ok and full_iff else "CHECK FAILED")
    print("d_Q(N)=phi(2N)/2; full rank (A_N=U_N) exactly for N prime or 2^m. One formula closes the map.")
    print("=" * 62)
    return 0 if ok and full_iff else 1


if __name__ == "__main__":
    raise SystemExit(main())
