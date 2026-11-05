#!/bin/bash
# M5: candidate abstracts from the full PubMed 2026 baseline (echo's verified copy), 16 processes
set -euo pipefail
~/Documents/Bioinformatics/bin/python3 ~/depclaims/scripts/01_extract_candidates.py ~/echo/raw/baseline ~/depclaims/work/cand 16
