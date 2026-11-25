"""
Type as oldforms ('O') the Saito-Kurokawa lifts that are new at a lower level.

The Saito-Kurokawa lift of a newform f in S_4(Gamma_0(M)) is a paramodular newform of level M,
so at a level N > M it is an oldform.  Many such lifts in hecke_evs_3_0/data were typed 'O',
as the README defines, but others were typed 'P' (for example the lifts of the level-256
forms at level 768, where there is no space at level 256 to see them in, but also the lift of
16.4.a.a at level 976).  The LMFDB has none of them in its newspaces: every one of its P
newforms is the lift of a form of its own level.

Each record of type 'P' is matched to the newforms f of weight 4, trivial character, level
dividing N and dimension the degree of its field, with
    tr lambda_p = tr a_p(f) + d (p + p^2)
at every good p < 100.  It must match exactly one; if f has level < N the record becomes 'O'.
Nothing else in the record changes.  The newforms come from the LMFDB's mf_newforms (read
only), so run this from an LMFDB checkout:

    PYTHONPATH=/path/to/lmfdb sage -python fix_old_lifts.py [N ...]   # all levels by default; run in the repository root

The same records in hecke_evs_3_0/data_nl_200 (matched as in fix_nl_200_types.py, by field
polynomial, traces of lambda_p for p < 100 and Atkin-Lehner signs, which that script made the
same) get the same type.
"""
import glob
import re
import sys
from collections import Counter

from sage.all import divisors, prime_range
from lmfdb import db
from rebasis_hecke_rings import dump, load
from fix_nl_200_types import traces

PS = prime_range(100)
_classical = {}


def classical(M, d):
    if (M, d) not in _classical:
        _classical[(M, d)] = list(db.mf_newforms.search(
            {'level': M, 'weight': 4, 'char_orbit_index': 1, 'dim': d}, ['label', 'level', 'traces']))
    return _classical[(M, d)]


def lifted_form(N, e):
    d = len(e['field_poly']) - 1
    tr = e['trace_lambda_p']
    good = [(i, p) for i, p in enumerate(PS) if N % p and i < len(tr) and tr[i] != 'NULL']
    assert len(good) >= 10, (N, e['field_poly'])
    hits = [f for M in divisors(N) for f in classical(M, d)
            if all(int(tr[i]) == f['traces'][p - 1] + d * (p + p * p) for i, p in good)]
    assert len(hits) == 1, "level %d, field %s: %d classical forms match" % (N, e['field_poly'], len(hits))
    return hits[0]


def process(N, stats):
    fname = 'hecke_evs_3_0/data/hecke_ev_3_0_%d.dat' % N
    with open(fname) as f:
        text = f.read()
    nl = "\n" if text.endswith("\n") else ""
    entries = load(text)
    changed = []
    for i, e in enumerate(entries):
        if e['aut_rep_type'] != 'P':
            continue
        f = lifted_form(N, e)
        if f['level'] == N:
            stats['P, lift of a form of level N'] += 1
            continue
        new = dict(e)     # keeps the key order
        new['aut_rep_type'] = 'O'
        entries[i] = new
        changed.append((i, f['label']))
    stats['P -> O'] += len(changed)
    if changed:
        assert dump(load(text)) + nl == text, "%s is not reproduced from what was read" % fname
        with open(fname, 'w') as f:
            f.write(dump(entries) + nl)
        print("level %d: %s" % (N, ", ".join("#%d (lift of %s)" % c for c in changed)), flush=True)
    process_nl_200(N, entries, stats)


def process_nl_200(N, data, stats):
    fname = 'hecke_evs_3_0/data_nl_200/hecke_ev_3_0_nl_200_%d.dat' % N
    try:
        with open(fname) as f:
            text = f.read()
    except FileNotFoundError:
        return
    nl = "\n" if text.endswith("\n") else ""
    entries = load(text)
    changed = 0
    for i, e in enumerate(entries):
        if e['aut_rep_type'] != 'P':      # only a P record can change
            continue
        tr = traces(e)
        m = [f for f in data if f['field_poly'] == e['field_poly'] and f['trace_lambda_p'] == tr[:len(f['trace_lambda_p'])]
             and f['atkin_lehner_string'] == e['atkin_lehner_string']]
        assert len(m) == 1, "%s, record %d: %d records of data match" % (fname, i, len(m))
        if e['aut_rep_type'] != m[0]['aut_rep_type']:
            assert (e['aut_rep_type'], m[0]['aut_rep_type']) == ('P', 'O'), (fname, i)
            new = dict(e)
            new['aut_rep_type'] = 'O'
            entries[i] = new
            changed += 1
    stats['data_nl_200: P -> O'] += changed
    if changed:
        assert dump(load(text)) + nl == text, "%s is not reproduced from what was read" % fname
        with open(fname, 'w') as f:
            f.write(dump(entries) + nl)


if __name__ == '__main__':
    levels = sorted(map(int, sys.argv[1:])) or sorted(
        int(re.search(r'_(\d+)\.dat$', f).group(1)) for f in glob.glob('hecke_evs_3_0/data/hecke_ev_3_0_*.dat'))
    stats = Counter()
    for N in levels:
        process(N, stats)
    print(dict(stats), flush=True)
