"""
omf5_data rows for the paramodular forms of weight (k, j) != (3, 0), from runs of omf5 with a
weight at precision 200:

    src/omf5 -genus=qf5plain.txt -disc=N -weight=a,b -prec=200 > RUNDIR/N_W.out

omf5's weight (a, b) gives the paramodular weight (k, j) = (a - b + 3, 2b), and W is a for
b = 0 and a,b otherwise: W = 1 is (4, 0), W = 1,1 is (3, 2), W = 2,1 is (4, 2).  omf5
normalizes T(p) and T(p,2) by p^a and p^(2a), so that the eigenvalues are algebraic integers
and the spinor Euler factor at a good prime is

    1 - T x + p (T2 + p^(2a) (p^2 + 1)) x^2 - p^(2a+3) T x^3 + p^(4a+6) x^4.

For each level N it writes OUTDIR/hecke_ev_K_J_N.dat, in the format of hecke_evs_3_0/data
(see the README), with one record per form of the run:

  - 'aut_rep_type' is omf5's type, except that a form of type G whose eigenvalues at the
    primes p not dividing N (T(p) for p < 200, T(p,2) for p^2 < 200) are those of a form
    of type G of a run of level M, M | N, M < N, gets 'O': it is an old form.  (They are
    compared by their characteristic polynomials: see signature().)
    So the forms of type G are the new non-lifts.
  - The Atkin-Lehner fields are the signs of the paramodular form: for p | N, the sign is
    e_p if p does not divide the conductor d of the spinor character of the form in the
    run, and -e_p if it does, e_p the Eichler invariant at p of the lattice (as for weight
    0, and for weights (a, b) by the functional equation).  e_p is read off the weight-0
    data: it is the sign at p of the constant form (T(p) = (p + 1)(p^2 + 1)), whose
    spinor character is trivial, in WT0DIR/hecke_ev_3_0_N.dat.
  - 'trace_lambda_p' and 'trace_lambda_p_square' are the traces of T(p), p < 100, and of
    T(p,2), p^2 < 200, with 'NULL' at p | N.
  - The new forms of type G of degree at most 20 get the Hecke ring (the order generated
    by these eigenvalues), its index in the ring of integers, nbound and the eigenvalue
    coordinates, as add_coordinates_14_20.orbit_data computes them; their field records
    (field_disc, its factorization, polredabs) go to OUTDIR/fields_K_J_N.json.
    The new forms of type G of larger degree get 'trace_lambda_p_squared' instead: the
    traces of T(p)^2 for p^2 < 200 ('NULL' at p | N), from which the loader gets the
    traces of T(p^2).

There is no run of level 1 (no lattice of discriminant 1), so an old form coming from
level 1 is not recognized: this needs S_{k,j}(Sp(4, Z)) = 0, which holds for (4, 0),
(3, 2) and (4, 2), and is checked against the database by the loader.

OUTDIR/summary_K_J.json has, for each level, the dimensions of the new forms of type G for
each Atkin-Lehner string, and the number of forms of each type.

    sage -python convert_weight_runs.py --weight=W RUNDIR WT0DIR OUTDIR N [N ...]

The runs of hecke_evs_4_0 and hecke_evs_3_2 used omf5 855b651.
"""
import json
import os
import re
import sys
import time

from sage.all import PolynomialRing, QQ, NumberField, prime_divisors, prime_range, sage_eval

from add_coordinates_14_20 import orbit_data, KEYS

PREC = 200
B = 100           # T(p) stored for p < B, as in hecke_evs_3_0
MAX_RING_DEGREE = 20


def load(fname):
    with open(fname) as f:
        text = f.read()
    return eval(re.sub(r'(\d)L', r'\1', text[text.index('{'):text.rindex('}') + 1]))


def weight_of(W):
    a, b = (int(x) for x in W.split(',')) if ',' in W else (int(W), 0)
    return a, b, a - b + 3, 2 * b


def run_forms(rundir, N, W):
    """[(d, type, field_poly coefficients, F, {p: T(p)}, {p: T(p,2)})] of the run of level N."""
    fname = os.path.join(rundir, '%d_%s.out' % (N, W))
    if not os.path.exists(fname) or os.path.getsize(fname) == 0:
        return None
    Qx = PolynomialRing(QQ, 'x')
    good = [p for p in prime_range(PREC) if N % p]
    good_sq = [p for p in prime_range(PREC) if N % p and p * p < PREC]
    out = []
    for d, fs in load(fname).items():
        for f in fs:
            c = f['field_poly'].split()
            assert int(c[0]) == len(c) - 1, (N, d, f['field_poly'])
            pol = Qx([int(x) for x in c[1:]])
            F = NumberField(pol, 'a')
            a = F.gen()
            assert len(f['lambda_p']) == len(good) and len(f['lambda_p_square']) == len(good_sq), (N, d)
            lp = {p: F(sage_eval(s, locals={'a': a})) for p, s in zip(good, f['lambda_p'])}
            lp2 = {p: F(sage_eval(s, locals={'a': a})) for p, s in zip(good_sq, f['lambda_p_square'])}
            out.append((int(d), f['aut_rep_type'], [int(x) for x in c[1:]], F, lp, lp2))
    return out


def signature(form, N):
    """
    The characteristic polynomials of the eigenvalues at the primes not dividing N, and of
    T(q1) + 2 T(q2) + 3 T(q3) for the first three of them: equal for an old form and its
    source.  (Not the field_poly: omf5 takes the minimal polynomial of the first eigenvalue
    that generates the field, so it depends on which primes divide the level.)
    """
    _, _, _, F, lp, lp2 = form
    good = [p for p in sorted(lp) if N % p]
    comb = sum((i + 1) * lp[p] for i, p in enumerate(good[:3]))
    return (tuple(tuple(lp[p].charpoly().list()) for p in good),
            tuple(tuple(lp2[p].charpoly().list()) for p in sorted(lp2) if N % p),
            tuple(comb.charpoly().list()))


def eichler_invariants(wt0dir, N):
    """{p: e_p}, from the sign of the constant form in the weight-0 data."""
    with open(os.path.join(wt0dir, 'hecke_ev_3_0_%d.dat' % N)) as f:
        rows = eval(f.read())
    ps = prime_range(B)
    const = [r for r in rows if len(r['field_poly']) == 2 and
             all(t == 'NULL' or t == (p + 1) * (p * p + 1) for p, t in zip(ps, r['trace_lambda_p']))]
    assert len(const) == 1, "level %d: %d constant forms in the weight-0 data" % (N, len(const))
    return {p: (1 if s == '+' else -1) for p, s in zip(prime_divisors(N), const[0]['atkin_lehner_string'])}


def process(N, W, runs, wt0dir, outdir, summary):
    a, b, k, j = weight_of(W)
    forms = runs[N]
    e = eichler_invariants(wt0dir, N) if forms else {}
    ps, sq_ps = prime_range(B), [p for p in prime_range(PREC) if p * p < PREC]
    rows, records, dims, types = [], [], {}, {}
    for form in forms:
        d, typ, coeffs, F, lp, lp2 = form
        if typ == 'G':
            for M in sorted(runs):
                if M < N and N % M == 0 and runs[M] is not None:
                    s = signature(form, N)
                    if any(g[1] == 'G' and signature(g, N) == s for g in runs[M]):
                        typ = 'O'
                        break
        al = [[int(p), e[p] * (-1 if d % p == 0 else 1)] for p in prime_divisors(N)]
        al_string = ''.join('+' if s == 1 else '-' for _, s in al)
        lambda_p = ['NULL' if N % p == 0 else lp[p] for p in ps]
        lambda_p_square = ['NULL' if N % p == 0 else lp2[p] for p in sq_ps]
        row = {'field_poly': coeffs, 'aut_rep_type': typ, 'atkin_lehner_eigenvals': al,
               'atkin_lehner_string': al_string,
               'trace_lambda_p': ['NULL' if x == 'NULL' else int(x.trace()) for x in lambda_p],
               'trace_lambda_p_square': ['NULL' if x == 'NULL' else int(x.trace()) for x in lambda_p_square]}
        if typ == 'G' and F.degree() > MAX_RING_DEGREE:
            # no coordinates: the traces of T(p)^2, which the trace of T(p^2) needs
            row['trace_lambda_p_squared'] = ['NULL' if x == 'NULL' else int((lp[p]**2).trace())
                                             for p, x in zip(sq_ps, lambda_p_square)]
        if typ == 'G' and F.degree() <= MAX_RING_DEGREE:
            t0 = time.time()
            new, record = orbit_data(N, F, lambda_p, lambda_p_square)
            row.update(new)
            row = {key: row[key] for key in KEYS}
            records.append(dict({'level': N, 'field_poly': coeffs}, **record))
            print("level %d: G orbit of dimension %d, index %d, nbound %d (%.1fs)"
                  % (N, F.degree(), new['hecke_ring_index'], new['hecke_ring_generator_nbound'], time.time() - t0),
                  flush=True)
        if typ == 'G':
            dims[al_string] = dims.get(al_string, 0) + int(F.degree())
        types[typ] = types.get(typ, 0) + 1
        rows.append(row)
    with open(os.path.join(outdir, 'hecke_ev_%d_%d_%d.dat' % (k, j, N)), 'w') as f:
        f.write(str(rows) + '\n')
    with open(os.path.join(outdir, 'fields_%d_%d_%d.json' % (k, j, N)), 'w') as f:
        json.dump(records, f)
    summary[int(N)] = {'new_G_dims_by_AL': dims, 'types': types}


def main(args):
    assert args[0].startswith('--weight=')
    W = args[0][len('--weight='):]
    rundir, wt0dir, outdir = args[1:4]
    levels = [int(N) for N in args[4:]]
    a, b, k, j = weight_of(W)
    os.makedirs(outdir, exist_ok=True)
    # the runs of every divisor of every level, for the old forms
    need = sorted({M for N in levels for M in range(1, N + 1) if N % M == 0})
    runs = {M: run_forms(rundir, M, W) for M in need}
    summary_file = os.path.join(outdir, 'summary_%d_%d.json' % (k, j))
    summary = {}
    if os.path.exists(summary_file):
        with open(summary_file) as f:
            summary = {int(N): v for N, v in json.load(f).items()}
    for N in levels:
        if runs[N] is None:
            print("level %d: no run" % N, flush=True)
            continue
        missing = [M for M in need if M < N and N % M == 0 and runs[M] is None and M > 1]
        if missing:
            # an old form could come from a level without a run
            print("level %d: no run at the divisors %s, skipped" % (N, missing), flush=True)
            continue
        process(N, W, runs, wt0dir, outdir, summary)
        print("level %d:" % N, summary[N], flush=True)
    with open(summary_file, 'w') as f:
        json.dump(summary, f, sort_keys=True)


if __name__ == '__main__':
    main(sys.argv[1:])
