"""
Correct the automorphic types and Atkin-Lehner signs of hecke_evs_3_0/data_nl_200.

Every record there was labelled 'G', though many are lifts or oldforms, and its
Atkin-Lehner fields held the spinor character of the orthogonal modular form, not the signs
of the paramodular form (which hecke_evs_3_0/data has had since "Store paramodular
Atkin-Lehner signs").  The two differ exactly at the primes where the quaternion algebra of
the genus ramifies, an odd number of primes of the level.

Each record is matched to the form of hecke_evs_3_0/data at the same level with the same
field polynomial, the same traces of lambda_p for p < 100, and the converted signs; its
'aut_rep_type', 'atkin_lehner_string' and 'atkin_lehner_eigenvals' are copied from it.
The set of primes where the signs change is not read from the lattice: for each level it is
the only set of odd size for which every record matches exactly one form of data, and
different records match different forms.  The script fails if there is not exactly one.
A level whose records already agree with their forms of data is skipped.

    sage -python fix_nl_200_types.py [N ...]      # all levels by default; run in the repository root
"""
import glob
import re
import sys
from collections import Counter

from sage.all import Matrix, QQ, NumberField, polygen, vector
from rebasis_hecke_rings import dump, load, rows_of


def traces(e):
    """The traces of lambda_p (stored only for the records without coordinates)."""
    if not e.get('lambda_p'):
        return e['trace_lambda_p']
    d = len(e['field_poly']) - 1
    if d == 1:
        tv = vector(QQ, [1])
    else:
        x = polygen(QQ)
        F = NumberField(sum(a * x**i for i, a in enumerate(e['field_poly'])), 'a')
        tv = vector(QQ, [(F.gen()**i).trace() for i in range(d)])
    B = rows_of(e['hecke_ring_numerators'], e['hecke_ring_denominators']) if e.get('hecke_ring_numerators') \
        else Matrix.identity(QQ, d)
    return ['NULL' if v == 'NULL' else int(vector(QQ, [QQ(t) for t in v]) * B * tv) for v in e['lambda_p']]


def convert(s, flips):
    return ''.join(('-' if c == '+' else '+') if i in flips else c for i, c in enumerate(s))


def process(N, stats):
    fname = 'hecke_evs_3_0/data_nl_200/hecke_ev_3_0_nl_200_%d.dat' % N
    with open(fname) as f:
        text = f.read()
    nl = "\n" if text.endswith("\n") else ""
    entries = load(text)
    if not entries:
        return
    with open('hecke_evs_3_0/data/hecke_ev_3_0_%d.dat' % N) as f:
        data = load(f.read())
    cands = []
    for e in entries:
        tr = traces(e)
        cands.append([j for j, f in enumerate(data) if f['field_poly'] == e['field_poly']
                      and f['trace_lambda_p'] == tr[:len(f['trace_lambda_p'])]])
    patterns = set(frozenset(i for i, (s, t) in enumerate(zip(e['atkin_lehner_string'], data[j]['atkin_lehner_string']))
                             if s != t) for e, c in zip(entries, cands) for j in c)
    def fit(fl):
        m = [[j for j in c if data[j]['atkin_lehner_string'] == convert(e['atkin_lehner_string'], fl)]
             for e, c in zip(entries, cands)]
        if all(len(x) == 1 for x in m) and len(set(x[0] for x in m)) == len(m):
            return [x[0] for x in m]

    done = fit(frozenset())
    if done is not None and all(e[k] == data[j][k] for e, j in zip(entries, done)
                                for k in ('aut_rep_type', 'atkin_lehner_string', 'atkin_lehner_eigenvals')):
        stats['records already corrected'] += len(entries)
        return
    good = {fl: fit(fl) for fl in patterns if len(fl) % 2 and fit(fl) is not None}
    assert len(good) == 1, "level %d: %d sets of primes fit" % (N, len(good))
    (fl, match), = good.items()
    changed = 0
    for i, j in enumerate(match):
        e, f = entries[i], data[j]
        assert [p for p, _ in e['atkin_lehner_eigenvals']] == [p for p, _ in f['atkin_lehner_eigenvals']]
        new = dict(e)     # keeps the key order
        for k in ('aut_rep_type', 'atkin_lehner_string', 'atkin_lehner_eigenvals'):
            new[k] = f[k]
        stats[(e['aut_rep_type'], f['aut_rep_type'])] += 1
        if new != e:
            entries[i] = new
            changed += 1
    stats['records'] += len(entries)
    stats['changed'] += changed
    if changed:
        assert dump(load(text)) + nl == text, "%s is not reproduced from what was read" % fname
        with open(fname, 'w') as f:
            f.write(dump(entries) + nl)
    primes = [p for p, _ in entries[0]['atkin_lehner_eigenvals']]
    print("level %d: signs change at %s; %d of %d records changed" % (
        N, sorted(primes[i] for i in fl), changed, len(entries)), flush=True)


if __name__ == '__main__':
    levels = sorted(map(int, sys.argv[1:])) or sorted(
        int(re.search(r'_(\d+)\.dat$', f).group(1)) for f in glob.glob('hecke_evs_3_0/data_nl_200/*.dat'))
    stats = Counter()
    for N in levels:
        process(N, stats)
    print({str(k): v for k, v in stats.items()}, flush=True)
