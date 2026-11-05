#!/bin/zsh
# Re-run the unchanged pipeline on the Claude Sonnet 5.5 reader output (notes/analysis_plan_review2_2026-11-04.md, section E).
# The Qwen-based outputs are kept under work/qwen/ for the between-reader comparison. Steps that need the M5 (figure data) are started
# with mjob and waited for. Usage: scripts/33_rerun_reader.sh
set -euo pipefail
cd /Volumes/ResearchData/Claude/depclaims
PY=/opt/homebrew/Cellar/micromamba/2.6.2/envs/admg/bin/python
MJ=/Volumes/ResearchData/Claude/remote/m5/mjob
S=work/creader/all_sonnet.jsonl
[[ $(wc -l < $S) -eq 32939 ]] || { echo "reader output incomplete"; exit 1; }
if [[ ! -d work/qwen ]]; then
  mkdir -p work/qwen/db_tables
  cp work/reader_all.jsonl work/reader_mapped_explore.parquet work/qwen/
  cp -R release work/qwen/release; cp -R db/parquet work/qwen/db_parquet; cp -R figures/main work/qwen/figures_main
  cp -R work/fig17/out work/qwen/fig17_out; cp -R work/fig work/qwen/fig
  for t in confirm_results independent_repeats robustness_additional readout_platform_additional gdsc_drug_additional; do
    duckdb db/depclaims.duckdb -c "copy (select * from $t) to 'work/qwen/db_tables/$t.parquet' (format parquet)" || true
  done
fi
cp $S work/reader_all.jsonl
$PY -I scripts/07_map_reader.py work/reader_all.jsonl work/reader_mapped_explore.parquet
DEPCLAIMS_READER='Claude Sonnet 5.5 reader (prompt of reader v2)' $PY scripts/12_confirm.py > work/confirm_run_sonnet.log 2>&1
duckdb db/depclaims.duckdb -c "copy (select * from confirm_results) to 'work/confirm_results_sonnet.parquet' (format parquet)"
# author lists for all candidate papers (superset of every claim set)
[[ -f work/dcauth_all/out/authors.parquet ]] && { cp work/dcauth/out/authors.parquet work/qwen/authors_claims.parquet 2>/dev/null || true; cp work/dcauth_all/out/authors.parquet work/dcauth/out/authors.parquet; }
$PY scripts/15_independent_repeats.py > work/repeats_sonnet.log 2>&1
$PY scripts/16_figdata_line.py > work/figdata_line_sonnet.log 2>&1
# figure data on the M5 (same inputs as the first run, new claim tables)
IN=(); for f in raw/gene_effect_26q1_chronos.csv raw/CEGv2.txt raw/NEGv1.txt raw/d2_gene_effect.csv raw/Model.csv \
   raw/OmicsExpressionProteinCodingGenesTPMLogp1.csv db/parquet/claims_confirmatory.parquet db/parquet/claims_development.parquet; do IN+=(--in $f); done
$MJ run fig17s scripts/m5/17_figdata_calib.py --conda admg "${IN[@]}" --fetch work/fig17s > work/fig17s.log 2>&1
rm -rf work/fig17/out && mkdir -p work/fig17 && cp -R work/fig17s/out work/fig17/out
$PY scripts/25_atlas_release.py > work/release_sonnet.log 2>&1
$PY scripts/28_robustness.py > work/robustness_sonnet.log 2>&1
$PY scripts/30_readout_platform.py > work/readout_platform_sonnet.log 2>&1
$PY scripts/31_gdsc_drug.py > work/gdsc_sonnet.log 2>&1
$PY scripts/34_reader_agreement.py > work/reader_agreement.log 2>&1
echo "pipeline re-run done"
