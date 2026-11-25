# Hecke eigenvalues beyond p < 100

Raw omf5 output for Hecke eigenvalues at primes beyond those in **hecke_evs_3_0**, which stops at p < 100 for T(p) and at p <= 7 for T(p^2). These runs were made to compute spinor L-functions with the certified L-function calculator (edgarcosta/lfunctions). Each file is the unedited standard output of one omf5 run. The matching `.err` file is its standard error, ending with the wall time.

Each output is a dictionary keyed by the conductor d of the spinor character, as in **hecke_rows**. The value is a list of the eigenforms in that space. For each form:
- `field_poly` gives the number of coefficients, then the coefficients from the constant term up. For example, `'2  7 1'` is x + 7.
- `lambda_p` is the T(p) eigenvalue, as a polynomial in a root `a` of `field_poly`.
- `lambda_p_square` is the T(p,2) eigenvalue. It is printed only when the run used `-square`.

Types print as U because the runs do not use `-prec`. For p exactly dividing N, the printed T(p) value is the coefficient a_p = c of the L-function. The Euler factor is (1 - (c + eps_p p^(a+1)) x + p^w x^2)(1 + eps_p p^(a+1) x), where eps_p is the Atkin-Lehner sign, a = k - 3 + j/2 and w = 2k + j - 3 (a = 0 and w = 3 for weight (3,0)).

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

At p = 61, 2.K.61.3.0.a.a has c = 145 and eps_61 = -1, so the bad Euler factor is (1 - 84x + 61^3 x^2)(1 - 61x). With it, PARI's lfuncheckfeq on 1084 coefficients holds with root number +1 to 2^-191, and the L-function calculator certifies the L-function (rank 0, root number +1). The factor (1 - 83x + 61^3 x^2)(1 - 61x), from c = 145 + eps_61, holds only to 2^-47 and is not certified.

## Reduced files

Where noted below, a file is reduced. Orbits of degree at most 20 are kept as omf5 printed them. For a larger orbit, `lambda_p` and `lambda_p_square` are replaced by `trace_lambda_p` and `trace_lambda_p_square`, their traces from the Hecke field to Q, and `field_poly` is kept. This follows **hecke_evs_3_0**, which stores only traces for orbits of degree greater than 20. The raw output would be 4 to 500 times larger. Almost all of it is the eigenvalues of these orbits, printed as polynomials of degree up to 265 with large rational coefficients.

## 3_0/N499, 3_0/N997: weight (3,0), levels 499 and 997 (reduced)

T(p) at p = 53, 101 and 199 for every orbit, from timing runs. omf5 commit 0d2c6a1:

    omf5 -genus=qf5db.sage -isom -disc=N -format=GG -p=p

`p{p}.dat` is the reduced output. `p{p}.err` is the standard error, ending with the wall time. At p = 53 the traces agree with `trace_lambda_p` in **hecke_evs_3_0**.

## 4_0/N31: weight (4,0), level 31

omf5 weight (1,0), which is the paramodular weight (4,0). omf5 commit 0d2c6a1:

    omf5 -genus=qf5plain.txt -disc=31 -weight=1 -prec=200     (prec.txt: T(p) for p < 200, p != 31)
    omf5 -genus=qf5plain.txt -disc=31 -weight=1 -p=31         (p31.txt)

`qf5plain.txt` is `qf5.db` of this repository. Both forms have spinor conductor 31:

| field_poly | orbit |
|---|---|
| `'2  15 1'` (lambda_2 = -15) | 2.K.31.4.0.a.a (type G) |
| degree 5, x^5 - 51x^4 + ... | 2.K.31.4.0.a.b (Saito-Kurokawa lifts of 31.6.a.a) |

The Atkin-Lehner sign of 2.K.31.4.0.a.a is eps_31 = +1. At p = 31, a = 1 and c = 7527, and the bad Euler factor is (1 - (c + 31^2) x + 31^5 x^2)(1 + 31^2 x). With this factor, PARI's lfuncheckfeq on the first 199 coefficients holds with root number +1, to 2^-57.

## 3_0/bad_primes: weight (3,0), the primes dividing the level (reduced)

T(p) and T(p,2) at the primes p dividing N, at 882 levels N <= 1000. **hecke_evs_3_0** stores `NULL` there. At the levels without a full run, the sweep also covers three good primes, which are used to match the omf5 forms with the stored orbits. There are 3,094 runs. omf5 commit f3e0c64, with the default output format:

    omf5 -genus=qf5db.sage -isom -disc=N -p=p -square

`N_p.dat` is the reduced output of one run. `times.txt` has one line per run: N, p, the exit code, the wall time and the peak memory. `one.sh` and `one2.sh` are the scripts of the two sweeps.

- The runs at p = N for the prime levels 463, 491, 503, 509, 547, 587 and 631 took 16 to 30 hours each. The runs at p = N for 797, 811 and 997 are still running and are not included yet.
- For p not dividing N and p < 100, the degrees and T(p) traces of every run agree with **hecke_evs_3_0**, and every T(p,2) trace stored there also occurs in the run, with one exception. At level 488 (p = 3, 5, 7), commit f3e0c64 prints one old form with a rational field. **hecke_evs_3_0** gives it the field x^2 - 20x + 88, after the fix of old-form fields in omf5 PR #75, which f3e0c64 predates. The newforms are not affected.
