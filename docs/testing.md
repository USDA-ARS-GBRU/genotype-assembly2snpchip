# Testing

## Tiny Regression Checks

The repository includes small fixture-based tests for the main analysis helpers.

### Parser

```bash
bash tests/run_tiny_test.sh
```

### Summary plots

```bash
bash tests/run_tiny_plot_test.sh
```

### PCA / MDS

```bash
bash tests/run_tiny_pca_test.sh
```

### GRIN enrichment

```bash
bash tests/run_tiny_grin_test.sh
```

## Pixi Shortcuts

```bash
pixi run test-parser
pixi run test-plots
pixi run test-pca
pixi run test-grin
```

## Matplotlib Cache Note

On some HPC systems, set `MPLCONFIGDIR` to a writable temporary directory before plotting:

```bash
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-$USER"
mkdir -p "$MPLCONFIGDIR"
```

This prevents matplotlib from trying to write cache files under a restricted home directory.
