"""Command line: depclaims build | score | lineage | published."""
from __future__ import annotations

import argparse, json, sys

from .core import Atlas, build


def main(argv=None):
    ap = argparse.ArgumentParser(prog='depclaims', description='Score published gene-line growth claims against DepMap screens.')
    sub = ap.add_subparsers(dest='cmd', required=True)
    b = sub.add_parser('build', help='build the local cache from DepMap downloads')
    for k in ('chronos', 'demeter2', 'expression', 'models', 'ceg', 'neg'):
        b.add_argument(f'--{k}', required=True)
    b.add_argument('--out', required=True)
    for name, hlp in (('score', 'screen evidence for one gene-line claim'), ('lineage', 'gene effects across the lines of a lineage'),
                      ('published', 'published claims in the atlas')):
        p = sub.add_parser(name, help=hlp)
        p.add_argument('--cache', required=True); p.add_argument('--claims')
        p.add_argument('gene'); p.add_argument('line' if name != 'lineage' else 'lineage', nargs='?' if name == 'published' else None)
    a = ap.parse_args(argv)
    if a.cmd == 'build':
        print(build(a.chronos, a.demeter2, a.expression, a.models, a.ceg, a.neg, a.out)); return 0
    A = Atlas(a.cache, a.claims)
    if a.cmd == 'score':
        json.dump(A.score(a.gene, a.line).as_dict(), sys.stdout, indent=2); print()
    elif a.cmd == 'lineage':
        print(A.lineage_profile(a.gene, a.lineage).to_string())
    else:
        print(A.published(a.gene, a.line).to_csv(sep='\t', index=False), end='')
    return 0


if __name__ == '__main__':
    sys.exit(main())
