"""
Correct the eigenvalues that parse_omf5 stored rounded down.

parse_omf5 (omf5/scripts/parse_omf.py, before omf5 PR #52) wrote each eigenvalue
in the basis of the order spanned by the first few lambda_p, rounding its
coordinates down when the eigenvalue was not in that order.  So in some orbits
whose stored hecke_ring_index was > 1, some stored eigenvalues (and their traces)
are wrong: each is the true eigenvalue with its coordinates in the stored basis
rounded down.  Witness: level 992, field x^2 - 8x - 4, where lambda_23 has trace
6 stored and 10 true.

This script compares every orbit with eigenvalue coordinates against a full run
of omf5 at precision 100,

    src/omf5 -genus=qf5db.sage -isom -disc=N -prec=100

(all forms, lifts and oldforms included).  An orbit is correct if some form of the
run with the same field_poly has exactly its eigenvalues.  Otherwise the run must
have exactly one eigenvalue system with the same field_poly whose eigenvalues, with their
coordinates in the basis of OLDDATADIR (the omf5_data of commit 5693c84a, as
parse_omf5 wrote it) rounded down, are the stored ones, and the orbit is
rewritten from that form: the eigenvalues, their traces, and the Hecke ring with
its index and nbound, as in add_coordinates_14_20.py (orbit_data).  The other
fields, and every other orbit, are unchanged.

    sage -python fix_floored_eigenvalues.py RUNDIR OLDDATADIR DATADIR OUTDIR N [N ...]

Writes OUTDIR/hecke_ev_3_0_N.dat for each level with a corrected orbit, and
prints one line per corrected orbit and a summary line per level.
"""
import os
import re
import sys

from sage.all import Matrix, PolynomialRing, QQ, ZZ, NumberField, floor, prime_range, sage_eval, sqrt, vector

from add_coordinates_14_20 import orbit_data

B = 100


def load(fname):
    with open(fname) as f:
        return eval(re.sub(r'(\d)L', r'\1', f.read()))


def basis_of(e, F):
    a = F.gen()
    return [sum(n * a**i for i, n in enumerate(nums)) / den
            for nums, den in zip(e['hecke_ring_numerators'], e['hecke_ring_denominators'])]


def process(N, rundir, olddir, datadir, outdir):
    ps = prime_range(B)
    sq = prime_range(sqrt(B))
    good = [p for p in ps if N % p]
    good_sq = [p for p in sq if N % p]
    Qx = PolynomialRing(QQ, 'x')
    run = []
    for fs in load(os.path.join(rundir, '%d.dat' % N)).values():
        for f in fs:
            run.append(([int(c) for c in f['field_poly'].split()[1:]], f))
    fname = os.path.join(datadir, 'hecke_ev_3_0_%d.dat' % N)
    with open(fname) as f:
        newline = f.read().endswith('\n')
    entries = load(fname)
    old_entries = load(os.path.join(olddir, 'hecke_ev_3_0_%d.dat' % N))
    assert len(old_entries) == len(entries)
    fixed = 0
    for idx, (e, old) in enumerate(zip(entries, old_entries)):
        if not old.get('lambda_p'):
            continue
        assert e['field_poly'] == old['field_poly'] and e['atkin_lehner_string'] == old['atkin_lehner_string']
        F = NumberField(Qx(e['field_poly']), 'a')
        a = F.gen()
        Bas = basis_of(old, F)
        Binv = Matrix(QQ, [list(b) for b in Bas])**(-1)
        stored = {p: sum(c * b for c, b in zip(v, Bas)) for p, v in zip(ps, old['lambda_p']) if v != 'NULL'}
        stored.update({p * p: sum(c * b for c, b in zip(v, Bas))
                       for p, v in zip(sq, old['lambda_p_square']) if v != 'NULL'})
        exact, floored = [], []
        for fp, f in run:
            if fp != e['field_poly']:
                continue
            true = {p: F(sage_eval(s, locals={'a': a})) for p, s in zip(good, f['lambda_p'])}
            true.update({p * p: F(sage_eval(s, locals={'a': a})) for p, s in zip(good_sq, f['lambda_p_square'])})
            assert set(true) == set(stored), (N, idx)
            if all(true[j] == x for j, x in stored.items()):
                exact.append(true)
            elif all(vector(QQ, list(x)) * Binv == vector(ZZ, [floor(c) for c in vector(QQ, list(true[j])) * Binv])
                     for j, x in stored.items()):
                floored.append(true)
        if exact:
            continue
        # (an eigenvalue system can occur in the run more than once, in the spaces of different spinor characters)
        floored = [t for i, t in enumerate(floored) if t not in floored[:i]]
        assert len(floored) == 1, "level %d orbit %d: %d forms of the run round down to it" % (N, idx, len(floored))
        true = floored[0]
        # the current entry has the same eigenvalues as the old one (recompute_hecke_rings.py changes only the basis)
        cur = basis_of(e, F)
        assert all(sum(c * b for c, b in zip(v, cur)) == stored[p]
                   for p, v in zip(ps, e['lambda_p']) if v != 'NULL'), (N, idx)
        lambda_p = ['NULL' if N % p == 0 else true[p] for p in ps]
        lambda_p_square = ['NULL' if N % p == 0 else true[p * p] for p in sq]
        new, _ = orbit_data(N, F, lambda_p, lambda_p_square)
        new['trace_lambda_p'] = ['NULL' if x == 'NULL' else int(x.trace()) for x in lambda_p]
        new['trace_lambda_p_square'] = ['NULL' if x == 'NULL' else int(x.trace()) for x in lambda_p_square]
        wrong = [('T(%d)' % j if j in ps else 'T(%d,2)' % ZZ(j).sqrt()) for j in sorted(stored) if stored[j] != true[j]]
        out = dict(e)     # keeps the key order
        out.update({k: new[k] for k in new if k in e})
        print("level %d orbit %d type %s deg %d: index %d -> %d, nbound %d -> %d, wrong %s"
              % (N, idx, e['aut_rep_type'], F.degree(), e['hecke_ring_index'], out['hecke_ring_index'],
                 e['hecke_ring_generator_nbound'], out['hecke_ring_generator_nbound'], ' '.join(wrong)), flush=True)
        entries[idx] = out
        fixed += 1
    if fixed:
        with open(os.path.join(outdir, 'hecke_ev_3_0_%d.dat' % N), 'w') as f:
            f.write(str(entries) + ('\n' if newline else ''))
    print('level %d: %d corrected' % (N, fixed), flush=True)
    return fixed


if __name__ == '__main__':
    rundir, olddir, datadir, outdir = sys.argv[1:5]
    total = 0
    for N in map(int, sys.argv[5:]):
        total += process(N, rundir, olddir, datadir, outdir)
    print('orbits corrected:', total)
