# Functional connectivity matrix (200 parcels)

Static functional connectivity matrix derived from
[ts](https://dynr.circadia-lab.uk/reference/ts.md) by computing the
full-timeseries Pearson correlation across all 200 parcel pairs. Serves
as a ground-truth reference for validating sliding-window and
phase-based dynFC methods when the window spans the full timeseries.

## Usage

``` r
fc
```

## Format

A numeric matrix with 200 rows and 200 columns. Diagonal is 1;
off-diagonal values are Pearson correlations in \[-1, 1\].

## Source

Derived from [ts](https://dynr.circadia-lab.uk/reference/ts.md).
