#!/bin/bash
# M5: rule-based claim tagging over candidate abstracts, PubTator3 gene + cell-line annotations, Cellosaurus->DepMap map
set -euo pipefail
cd ~/depclaims
~/Documents/Bioinformatics/bin/python3 scripts/02_tag_claims.py work/cand ~/echo/raw/pubtator3/gene2pubtator3.gz cellline2pubtator3.gz cvcl2ach.tsv work/claims 16
