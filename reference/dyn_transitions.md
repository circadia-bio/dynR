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
tr
#> # A tibble: 16 × 4
#> # Groups:   tag, source, target [16]
#>    source target tag   data            
#>     <int>  <int> <chr> <list>          
#>  1      1      1 1_1   <tibble [2 × 5]>
#>  2      1      2 1_2   <tibble [2 × 5]>
#>  3      1      3 1_3   <tibble [2 × 5]>
#>  4      1      4 1_4   <tibble [2 × 5]>
#>  5      2      1 2_1   <tibble [2 × 5]>
#>  6      2      2 2_2   <tibble [2 × 5]>
#>  7      2      3 2_3   <tibble [2 × 5]>
#>  8      2      4 2_4   <tibble [2 × 5]>
#>  9      3      1 3_1   <tibble [2 × 5]>
#> 10      3      2 3_2   <tibble [2 × 5]>
#> 11      3      3 3_3   <tibble [2 × 5]>
#> 12      3      4 3_4   <tibble [2 × 5]>
#> 13      4      1 4_1   <tibble [2 × 5]>
#> 14      4      2 4_2   <tibble [2 × 5]>
#> 15      4      3 4_3   <tibble [2 × 5]>
#> 16      4      4 4_4   <tibble [2 × 5]>
```
