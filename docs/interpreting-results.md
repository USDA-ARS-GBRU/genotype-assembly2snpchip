# Interpreting Results

## First Look

Three values deserve immediate attention:

1. `sites_compared`
2. `match_fraction`
3. the gap between rank 1 and rank 2

A very high match fraction based on a handful of sites is weak evidence. A slightly lower match fraction based on thousands of sites may be much stronger.

## Common Patterns

### Strong Match

```text
rank 1   match_fraction 0.995   sites_compared 18000
rank 2   match_fraction 0.962   sites_compared 17950
```

This is usually reassuring.

### Ambiguous Match

```text
rank 1   match_fraction 0.991   sites_compared 16000
rank 2   match_fraction 0.989   sites_compared 16010
rank 3   match_fraction 0.988   sites_compared 15980
```

This often means related lines, duplicates, or near-duplicates in the panel.

### Weak Match

```text
rank 1   match_fraction 0.780   sites_compared 300
rank 2   match_fraction 0.775   sites_compared 290
```

This is not enough evidence for a confident identity assignment.

## Why Site Counts Differ

`sites_compared` usually does not equal the total panel marker count.

Sites can be lost because:

- the assembly does not align over the marker
- mapping quality is low
- the genotype was set to missing
- chromosome names are incompatible
- REF or ALT representation does not match cleanly
- the locus is structurally different in the assembly

## Common Problems and Fixes

### Very few compared sites

Check:

- FASTA and VCF coordinate compatibility
- BAM indexing
- panel indexing
- `MIN_MQ` and `MIN_GQ`
- whether the query VCF contains many `./.` calls

### Top two hits are nearly tied

That may be biologically real. Treat it as a cluster-level clue unless rank 1 clearly separates from rank 2.

### Expected sample is absent from the panel

The workflow can still find the closest available genotype, but it cannot prove identity to a sample missing from the database.
