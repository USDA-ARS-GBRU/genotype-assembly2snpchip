# Common Workflows

## Update a workflow flag or command

Check all of these:

- `sbatch/call_panel_variants_and_gtcheck.sbatch`
- `sbatch/call_panel_variants_and_gtcheck_array.sbatch`
- `docs/step-2-prepare-panel-and-call-sites.md`
- `docs/step-3-compare-and-summarize.md`
- `examples/soy50k_gtcheck_from_asm20.sbatch`
- `examples/README*.md`

## Update a plotting change

Check:

- `scripts/plot_gtcheck_summary.py`
- `docs/visualization.md`
- `README.md` if quick-start figures or captions change
- `examples/figures/` if example outputs should be regenerated

## Update environment/tool requirements

Check:

- `environment.yml`
- `pixi.toml`
- `requirements.txt` if Python-only dependency changes
- `README.md`
- `docs/setup.md`
- `docs/hpc_notes.md`

## Review a pull request or bug report

Start with:

- `rg` for the exact flag, script, or output being discussed
- `git log -- <file>` if behavior may have changed recently
- repo docs plus `examples/` if the issue smells like legacy soybean carryover

## Add new post-processing output

Update:

- canonical script in `scripts/`
- example outputs if bundled examples should illustrate it
- `README.md` quick-start expected outputs
- deeper `docs/` page explaining interpretation

