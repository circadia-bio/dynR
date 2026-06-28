# Dynamic phase-locking matrix (dPL)

Compute instantaneous phase-locking matrices from parcel-level phase
time series and extract leading eigenvectors via
[`get_leida()`](https://dynr.circadia-lab.uk/reference/get_leida.md).
The first and last 10 timepoints are discarded to avoid edge effects
from the Hilbert transform.

## Usage

``` r
dyn_phase_lock(phases)
```

## Arguments

- phases:

  Numeric matrix \[N × Tmax\]. Instantaneous phases in radians, as
  returned by
  [`hilbert_phases()`](https://dynr.circadia-lab.uk/reference/hilbert_phases.md).

## Value

A list with:

- sync_conn:

  Array \[N, N, Tmax-20\]. Instantaneous phase-locking matrices, one per
  (trimmed) timepoint.

- leida:

  Matrix \[Tmax-20, N\]. Leading eigenvectors from
  [`get_leida()`](https://dynr.circadia-lab.uk/reference/get_leida.md).

## References

Cabral, J. et al. (2017). Cognitive performance in healthy older adults
relates to spontaneous switching between states of functional
connectivity during rest. *Scientific Reports*, 7(1), 5135.
[doi:10.1038/s41598-017-05425-7](https://doi.org/10.1038/s41598-017-05425-7)

Lord, L.-D. et al. (2019). Dynamical exploration of the repertoire of
brain networks at rest is modulated by psilocybin. *NeuroImage*, 199,
127–142.
[doi:10.1016/j.neuroimage.2019.05.060](https://doi.org/10.1016/j.neuroimage.2019.05.060)

## Examples

``` r
set.seed(1)
ts <- matrix(rnorm(10 * 200), nrow = 10, ncol = 200)
phases <- hilbert_phases(ts)
res <- dyn_phase_lock(phases)
dim(res$sync_conn)  # 10 x 10 x 180
#> [1]  10  10 180
dim(res$leida)      # 180 x 10
#> [1] 180  10
```
