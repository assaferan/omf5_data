"""
Replace the old forms (type O) whose field was too small (omf5 PR #75).

When the Hecke operators at the primes up to 19 did not split a space, omf5
gave the old form the field of the last operator it tried, and read its other
eigenvalues off an eigenvector of that operator alone.  When that field is
smaller than the field of the forms in the space, those eigenvalues are not
those of any eigenvector: witness level 460, where the stored orbit x - 62 has
T_3 = -7, while T_3 has the characteristic polynomial (x^2 - 4x - 41)^3 on its
space.  omf5 PR #75 gives such a form the field of its forms (from the
operator with the largest irreducible factor, or from a combination of the
operators).  Every old form of omf5_data up to level 999 was checked against
the Hecke matrices at the good primes up to 19: these 25 orbits are the ones
which fail.

For each orbit (N, c, i) of ORBITS (the level, the conductor in the output of
omf5, and the index of the form in its list):

  - the stored orbit is the one with the field_poly and the traces at the good
    p < 20 of form i of conductor c in BASEDIR/N.base (where these are integers,
    or rounded down in degree 1, as parse_omf5 stored them), the output of omf5
    before PR #75 at precision 20 (each matches exactly one);
  - the new orbit is form i of conductor c in RUNDIR/N.out, the output of omf5
    after PR #75 at precision 100,

        src/omf5 -genus=qf5db.sage -isom -disc=N -prec=100

    whose eigenvalues at p < 20 must be those of form i of conductor c in
    NEWDIR/N.new (after PR #75, precision 20);
  - field_poly, trace_lambda_p and trace_lambda_p_square are replaced, and the
    orbit is given eigenvalue coordinates and its Hecke ring, as in
    add_coordinates_14_20.py (orbit_data), as every orbit of type O of degree
    at most 13 should have (none of the 24 stored orbits had them).
    aut_rep_type stays 'O' and the Atkin-Lehner fields are kept.  (At levels
    768 and 882 the new form has multiplicity one, and omf5 now gives it the
    type P, but it is the Saito-Kurokawa lift of a form of lower level, 256.4.a.n
    and 441.4.a.x, so it is not in the paramodular newspace.)

    sage -python fix_old_form_fields.py BASEDIR NEWDIR RUNDIR DATADIR OUTDIR

Writes OUTDIR/hecke_ev_3_0_N.dat for each level of ORBITS, and prints one line
per orbit.  One of these orbits, at the levels 488 and 976, is also in
data_nl_200 (this script does not change it; fix_nl_200_old_fields.py does).
None has degree 14 to 20, so fields.json is unchanged.
"""
import os
import re
import sys

from sage.all import PolynomialRing, QQ, NumberField, floor, prime_range, sage_eval, sqrt

from add_coordinates_14_20 import orbit_data, KEYS

B = 100

ORBITS = [(460, 1, 12), (468, 1, 13), (486, 1, 26), (488, 122, 10), (540, 1, 21), (552, 1, 16),
          (580, 1, 10), (644, 1, 13), (696, 1, 21), (736, 1, 36), (748, 187, 5), (756, 1, 20),
          (756, 1, 21), (756, 21, 15), (768, 1, 41), (780, 1, 14), (812, 1, 11), (828, 1, 13),
          (882, 1, 26), (920, 1, 22), (920, 2, 8), (920, 10, 22), (928, 1, 28), (936, 1, 22),
          (976, 122, 22)]


def load(fname):
    with open(fname) as f:
        return eval(re.sub(r'(\d)L', r'\1', f.read()))


def form(f):
    """The field and the eigenvalues of a form of the output of omf5."""
    Qx = PolynomialRing(QQ, 'x')
    pol = Qx([QQ(c) for c in f['field_poly'].split()[1:]])
    F = NumberField(pol, 'a')
    a = F.gen()
    ev = lambda s: 'NULL' if s == 'NULL' else F(sage_eval(s, locals={'a': a}))
    return pol, F, [ev(s) for s in f['lambda_p']], [ev(s) for s in f['lambda_p_square']]


def coeffs(pol):
    return [int(c) for c in pol.coefficients(sparse=False)]


def process(N, orbits, basedir, newdir, rundir, datadir, outdir):
    ps = prime_range(B)
    sq_ps = prime_range(sqrt(B))
    good20 = [p for p in prime_range(20) if N % p]
    fname = os.path.join(datadir, "hecke_ev_3_0_%d.dat" % N)
    with open(fname) as f:
        newline = f.read().endswith("\n")
    entries = load(fname)
    base, new20, run = (load(os.path.join(d, "%d.%s" % (N, s)))
                        for d, s in [(basedir, 'base'), (newdir, 'new'), (rundir, 'out')])
    done = set()
    for c, i in orbits:
        pol0, F0, lp0, _ = form(base[c][i])
        # 10 of these forms had non-integral eigenvalues, which parse_omf5 stored
        # rounded down (in degree 1) or with their coordinates rounded down; so
        # compare the traces which are integers, and the rounded down ones in degree 1
        tr0 = [x.trace() for x in lp0]
        same = lambda t, t0: t == t0 if t0.is_integer() else (pol0.degree() > 1 or t == floor(t0))
        m = [j for j, e in enumerate(entries) if j not in done and e['field_poly'] == coeffs(pol0)
             and all(same(e['trace_lambda_p'][ps.index(p)], t0) for p, t0 in zip(good20, tr0))]
        assert len(m) == 1, "level %d: stored orbits %s for form %d of conductor %d" % (N, m, i, c)
        j = m[0]
        done.add(j)
        e = entries[j]
        assert e['aut_rep_type'] == 'O', (N, j)
        pol, F, lp, lp2 = form(run[c][i])
        assert run[c][i]['field_poly'] == new20[c][i]['field_poly'], (N, c, i)
        _, _, lp20, _ = form(new20[c][i])
        assert lp[:len(lp20)] == lp20, (N, c, i)
        good = [p for p in ps if N % p]
        good_sq = [p for p in sq_ps if N % p]
        assert len(lp) == len(good) and len(lp2) == len(good_sq) and 'NULL' not in lp + lp2, (N, c, i)
        lambda_p = ['NULL' if N % p == 0 else lp[good.index(p)] for p in ps]
        lambda_p_square = ['NULL' if N % p == 0 else lp2[good_sq.index(p)] for p in sq_ps]
        tr = lambda xs: ['NULL' if x == 'NULL' else int(x.trace()) for x in xs]
        new = dict(e)
        new['field_poly'] = coeffs(pol)
        new['trace_lambda_p'] = tr(lambda_p)
        new['trace_lambda_p_square'] = tr(lambda_p_square)
        if pol.degree() <= 13:
            ring, _ = orbit_data(N, F, lambda_p, lambda_p_square)
            for k in KEYS:
                if k in ring:
                    new[k] = ring[k]
        entries[j] = new
        print("level %d: orbit %d (form %d of conductor %d): field %s -> %s, coordinates %s"
              % (N, j, i, c, coeffs(pol0), coeffs(pol), 'yes' if pol.degree() <= 13 else 'no'), flush=True)
    with open(os.path.join(outdir, "hecke_ev_3_0_%d.dat" % N), "w") as f:
        f.write(str(entries) + ("\n" if newline else ""))


if __name__ == '__main__':
    basedir, newdir, rundir, datadir, outdir = sys.argv[1:6]
    levels = sorted(set(N for N, _, _ in ORBITS))
    if len(sys.argv) > 6:
        levels = [int(N) for N in sys.argv[6:]]
    for N in levels:
        process(N, [(c, i) for M, c, i in ORBITS if M == N], basedir, newdir, rundir, datadir, outdir)
