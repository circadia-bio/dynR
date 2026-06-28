# Leading eigenvector decomposition (LEiDA)

Extract the leading eigenvector from a series of instantaneous
phase-locking matrices. The sign of each eigenvector is normalised so
that its sum is non-positive, following the LEiDA convention.

## Usage

``` r
get_leida(sync_conn)
```

## Arguments

- sync_conn:

  Numeric array \[N, N, Tmax\]. Instantaneous phase-locking matrices, as
  returned by
  [`dyn_phase_lock()`](https://dynr.circadia-lab.uk/reference/dyn_phase_lock.md).

## Value

Numeric matrix \[Tmax × N\]. Leading eigenvectors, one per timepoint.

## References

Cabral, J. et al. (2017). Cognitive performance in healthy older adults
relates to spontaneous switching between states of functional
connectivity during rest. *Scientific Reports*, 7(1), 5135.
[doi:10.1038/s41598-017-05425-7](https://doi.org/10.1038/s41598-017-05425-7)

Lord, L.-D. et al. (2019). Dynamical exploration of the repertoire of
brain networks at rest is modulated by psilocybin. *NeuroImage*, 199,
127–142.
[doi:10.1016/j.neuroimage.2019.05.060](https://doi.org/10.1016/j.neuroimage.2019.05.060)
