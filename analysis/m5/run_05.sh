#!/bin/bash
# M5: frozen reader v2 over every paper with a named DepMap line (resumable)
set -euo pipefail
cd ~/depclaims
venv_mlx/bin/python 05_llm_read.py ~/models/Qwen3.5-9B-MLX-4bit llm_items_all.parquet reader_all.jsonl 8
