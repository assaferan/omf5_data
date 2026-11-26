"""
Replace the two orbits of data_nl_200 which fix_old_form_fields.py missed.

fix_old_form_fields.py (PR #9) replaced 25 old forms whose field was too small
(omf5 PR #75) in hecke_evs_3_0/data, but its docstring wrongly says that none
of them is in data_nl_200.  One of them is: the form of conductor 122 with the
field x - 434 at the levels 488 and 976, which is orbit 7 of
hecke_ev_3_0_nl_200_488.dat and orbit 20 of hecke_ev_3_0_nl_200_976.dat.  Its
field is x^2 - 20x + 88, as in hecke_evs_3_0/data since PR #9.

For each orbit (N, c, i, j) of ORBITS (the level, the conductor and the index
of the form in the output of omf5, and the index of the orbit in
hecke_evs_3_0/data):

  - the stored orbit of data_nl_200 is the one with the old field x - 434 and
    the Atkin-Lehner signs of orbit j of hecke_evs_3_0/data (exactly one);
  - the new orbit is form i of conductor c in RUNDIR/N.out, the output of omf5
    after PR #75 at precision 200,

        src/omf5 -genus=qf5db.sage -isom -disc=N -prec=200

    whose traces at p < 100, and of lambda_p_square, must be those of orbit j
    of hecke_evs_3_0/data;
  - field_poly, lambda_p, lambda_p_square and the Hecke ring are replaced, as
    in add_coordinates_14_20.py (orbit_data) with the primes up to 200.
    aut_rep_type stays 'O' and the Atkin-Lehner fields are kept.

    sage -python fix_nl_200_old_fields.py RUNDIR DATADIR OUTDIR

DATADIR is hecke_evs_3_0.  Writes OUTDIR/hecke_ev_3_0_nl_200_N.dat for both
levels, and prints one line per orbit.  The field has degree 2, so
fields.json is unchanged.
"""
import os
import re
import sys

from sage.all import PolynomialRing, QQ, NumberField, prime_range, sage_eval

import add_coordinates_14_20
from add_coordinates_14_20 import orbit_data, KEYS

B = 200
OLD_FIELD = [-434, 1]
ORBITS = [(488, 122, 10, 25), (976, 122, 22, 59)]


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


def process(N, c, i, j, rundir, datadir, outdir):
    # orbit_data takes the primes up to add_coordinates_14_20.B
    add_coordinates_14_20.B = B
    ps = prime_range(B)
    sq_ps = [p for p in ps if p * p < B]
    fname = os.path.join(datadir, "data_nl_200", "hecke_ev_3_0_nl_200_%d.dat" % N)
    with open(fname) as f:
        newline = f.read().endswith("\n")
    entries = load(fname)
    data = load(os.path.join(datadir, "data", "hecke_ev_3_0_%d.dat" % N))[j]
    m = [k for k, e in enumerate(entries) if e['field_poly'] == OLD_FIELD
         and e['atkin_lehner_string'] == data['atkin_lehner_string']]
    assert len(m) == 1, "level %d: stored orbits %s" % (N, m)
    k = m[0]
    e = entries[k]
    assert e['aut_rep_type'] == 'O', (N, k)
    assert e['atkin_lehner_eigenvals'] == data['atkin_lehner_eigenvals'], (N, k)
    run = load(os.path.join(rundir, "%d.out" % N))
    pol, F, lp, lp2 = form(run[c][i])
    good = [p for p in ps if N % p]
    good_sq = [p for p in sq_ps if N % p]
    assert len(lp) == len(good) and len(lp2) == len(good_sq) and 'NULL' not in lp + lp2, (N, c, i)
    lambda_p = ['NULL' if N % p == 0 else lp[good.index(p)] for p in ps]
    lambda_p_square = ['NULL' if N % p == 0 else lp2[good_sq.index(p)] for p in sq_ps]
    tr = lambda xs: ['NULL' if x == 'NULL' else int(x.trace()) for x in xs]
    assert [int(a) for a in pol.coefficients(sparse=False)] == data['field_poly'], (N, j)
    n = len(data['trace_lambda_p'])
    assert tr(lambda_p)[:n] == data['trace_lambda_p'], (N, j)
    n2 = len(data['trace_lambda_p_square'])
    assert tr(lambda_p_square)[:n2] == data['trace_lambda_p_square'], (N, j)
    assert len(lambda_p) == len(e['lambda_p']) and len(lambda_p_square) == len(e['lambda_p_square']), (N, k)
    ring, _ = orbit_data(N, F, lambda_p, lambda_p_square)
    new = dict(e)
    new['field_poly'] = data['field_poly']
    for key in KEYS:
        if key in ring:
            new[key] = ring[key]
    assert set(new) == set(e), (N, k, set(new) ^ set(e))
    entries[k] = new
    print("level %d: orbit %d (form %d of conductor %d): field %s -> %s"
          % (N, k, i, c, OLD_FIELD, new['field_poly']), flush=True)
    with open(os.path.join(outdir, "hecke_ev_3_0_nl_200_%d.dat" % N), "w") as f:
        f.write(str(entries) + ("\n" if newline else ""))


if __name__ == '__main__':
    rundir, datadir, outdir = sys.argv[1:4]
    for N, c, i, j in ORBITS:
        process(N, c, i, j, rundir, datadir, outdir)
