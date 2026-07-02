#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

require_fixed() {
    local needle="$1"
    local file="$2"
    if ! grep -Fq -- "$needle" "$file"; then
        fail "Missing expected text '$needle' in $file"
    fi
}

require_regex() {
    local pattern="$1"
    local file="$2"
    if ! grep -Eq -- "$pattern" "$file"; then
        fail "Missing expected pattern '$pattern' in $file"
    fi
}

for file in \
    sbatch/call_panel_variants_and_gtcheck.sbatch \
    sbatch/call_panel_variants_and_gtcheck_array.sbatch
do
    require_fixed "--keep-refs" "$file"
    require_fixed "+fixploidy" "$file"
    require_regex '^[[:space:]]+-i[[:space:]]*(\\)?[[:space:]]*$' "$file"
    if grep -Eq -- '^[[:space:]]+-i[[:space:]]+1([[:space:]]|\\|$)' "$file"; then
        fail "Found forbidden '-i 1' in $file"
    fi
done

for file in README.md docs/setup.md docs/hpc_notes.md
do
    require_fixed "bcftools >= 1.23" "$file"
done

require_fixed "--insert-missed" docs/step-2-prepare-panel-and-call-sites.md
require_fixed "--keep-refs" docs/step-3-compare-and-summarize.md

echo "Workflow invariants check passed."
