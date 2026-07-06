# Plot a brain state sequence

Tile plot of brain state labels over time (or windows), coloured with
the dynR state palette (periwinkle, brick red, mauve, deep indigo, dusty
rose).

## Usage

``` r
plot_state_sequence(
  states,
  x_label = "Timepoint",
  title = "Brain state sequence",
  palette = NULL
)
```

## Arguments

- states:

  Integer or factor vector of state labels.

- x_label:

  Character. x-axis label. Default `"Timepoint"`.

- title:

  Character. Plot title. Default `"Brain state sequence"`.

- palette:

  Character vector of colours. Defaults to the dynR 5-colour state
  palette, recycled if `K > 5`.

## Value

A `ggplot` object.

## Examples

``` r
set.seed(1)
states <- sample(1:4, 100, replace = TRUE)
plot_state_sequence(states)
```
