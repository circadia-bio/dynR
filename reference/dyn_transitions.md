# State transition probabilities (Markov analysis)

From a long-format data frame of cluster labels ordered in time, compute
first-order Markov transition probabilities between brain states: for
each source state, the fraction of transitions that lead to each target
state.

## Usage

``` r
dyn_transitions(tbl, vars, cVar, sortBy, groupBy, remIntra = FALSE)
```

## Arguments

- tbl:

  A data frame in long format (one row per timepoint).

- vars:

  Character vector of grouping / covariate column names to preserve in
  the output (e.g. `c("sub", "ses", "age", "sex")`).

- cVar:

  Character. Name of the column holding integer cluster labels.

- sortBy:

  Character vector of column names used to sort rows within each group
  before transitions are computed (typically
  `c("sub", "ses", "ttime")`).

- groupBy:

  Character vector of column names that define independent sequences
  (typically `c("sub", "ses")`). Transitions are never computed across
  group boundaries.

- remIntra:

  Logical. If `TRUE`, self-transitions (state → same state) are removed
  before normalising probabilities. Default `FALSE`.

## Value

A nested tibble with columns `tag`, `source`, `target`, and `data`.
`tag` encodes the transition as `"<source>_<target>"`. Each `data`
element is a per-group tibble with columns inherited from `vars` plus:

- n:

  Raw transition count.

- tot:

  Total transitions out of `source` for that group.

- nCount:

  Transition probability (`n / tot`).

## Details

Ported from `clusters_markov()` in the `neonatal_dfc` analysis pipeline
(França et al., *Nat Commun*).

## Examples

``` r
set.seed(1)
df <- data.frame(
  sub   = rep(c("A", "B"), each = 50),
  ses   = 1L,
  ttime = rep(seq_len(50), 2),
  clus4 = sample(1:4, 100, replace = TRUE)
)
tr <- dyn_transitions(
  df,
  vars    = c("sub", "ses"),
  cVar    = "clus4",
  sortBy  = c("sub", "ses", "ttime"),
  groupBy = c("sub", "ses")
)
#> Error in dyn_transitions(df, vars = c("sub", "ses"), cVar = "clus4", sortBy = c("sub",     "ses", "ttime"), groupBy = c("sub", "ses")): could not find function "dyn_transitions"
tr
#> Error: object 'tr' not found
```
