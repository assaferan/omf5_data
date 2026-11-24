# Hecke eigenvalues beyond p < 100

Raw omf5 output for Hecke eigenvalues at primes beyond those in **hecke_evs_3_0**, which stops at p < 100 for T(p) and at p <= 7 for T(p^2). These runs were made to compute spinor L-functions with the certified L-function calculator (edgarcosta/lfunctions). Each file is the unedited standard output of one omf5 run. The matching `.err` file is its standard error, ending with the wall time.

Each output is a dictionary keyed by the conductor d of the spinor character, as in **hecke_rows**. The value is a list of the eigenforms in that space. For each form:
- `field_poly` gives the number of coefficients, then the coefficients from the constant term up. For example, `'2  7 1'` is x + 7.
- `lambda_p` is the T(p) eigenvalue, as a polynomial in a root `a` of `field_poly`.
- `lambda_p_square` is the T(p,2) eigenvalue. It is printed only when the run used `-square`.

Types print as U because the runs do not use `-prec`. At p | N the T(p) value is omf5's eigenvalue lambda. For p exactly dividing N, the coefficient a_p of the L-function is c = lambda + eps_p p^a, where eps_p is the Atkin-Lehner sign and a = k - 3 + j/2 (a = 0 for weight (3,0)).

## 3_0/N61: weight (3,0), level 61

omf5 commit 0d2c6a1, lattice `-quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8` (disc 61), `-format=GG`.

- **raw_tp/q{p}.txt**: T(p) for all 180 primes p <= 1069, including p = 61. T(p,2) is included for p <= 31. The command is in `raw_tp/one.sh`.
- **raw_sq/s{p}.txt**: T(p) and T(p,2) for the 13 primes 37 <= p <= 97, p != 61. The command is in `raw_sq/one.sh`.

The space at d = 1 has dimension 8:

| field_poly | orbit |
|---|---|
| `'2  7 1'` (lambda_2 = -7) | 2.K.61.3.0.a.a (type G) |
| `'7  2026 -5205 4471 -1714 322 -29 1'` (degree 6) | 2.K.61.3.0.a.b (Saito-Kurokawa lifts) |
| `'2  -15 1'` (lambda_2 = 15) | the Eisenstein series |

At p = 61, 2.K.61.3.0.a.a has lambda = 145 and eps_61 = -1, so c = 144. The bad Euler factor is (1 - 83x + 61^3 x^2)(1 - 61x).
