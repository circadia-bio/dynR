#' State transition probabilities (Markov analysis)
#'
#' From a long-format data frame of cluster labels ordered in time, compute
#' first-order Markov transition probabilities between brain states: for each
#' source state, the fraction of transitions that lead to each target state.
#'
#' Ported from `clusters_markov()` in the `neonatal_dfc` analysis pipeline
#' (França et al., *Nat Commun*).
#'
#' @param tbl A data frame in long format (one row per timepoint).
#' @param vars Character vector of grouping / covariate column names to
#'   preserve in the output (e.g. `c("sub", "ses", "age", "sex")`).
#' @param cVar Character. Name of the column holding integer cluster labels.
#' @param sortBy Character vector of column names used to sort rows within each
#'   group before transitions are computed (typically `c("sub", "ses",
#'   "ttime")`).
#' @param groupBy Character vector of column names that define independent
#'   sequences (typically `c("sub", "ses")`). Transitions are never computed
#'   across group boundaries.
#' @param remIntra Logical. If `TRUE`, self-transitions (state → same state)
#'   are removed before normalising probabilities. Default `FALSE`.
#'
#' @return A nested tibble with columns `tag`, `source`, `target`, and `data`.
#'   `tag` encodes the transition as `"<source>_<target>"`. Each `data`
#'   element is a per-group tibble with columns inherited from `vars` plus:
#'   \item{n}{Raw transition count.}
#'   \item{tot}{Total transitions out of `source` for that group.}
#'   \item{nCount}{Transition probability (`n / tot`).}
#'
#' @importFrom dplyr rename all_of group_by arrange mutate ungroup filter summarise n
#' @importFrom tidyr nest
#' @importFrom rlang syms
#' @export
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   sub   = rep(c("A", "B"), each = 50),
#'   ses   = 1L,
#'   ttime = rep(seq_len(50), 2),
#'   clus4 = sample(1:4, 100, replace = TRUE)
#' )
#' tr <- dyn_transitions(
#'   df,
#'   vars    = c("sub", "ses"),
#'   cVar    = "clus4",
#'   sortBy  = c("sub", "ses", "ttime"),
#'   groupBy = c("sub", "ses")
#' )
#' tr
utils::globalVariables(c("clus", "tag", "target", "tot"))

dyn_transitions <- function(tbl, vars, cVar, sortBy, groupBy,
                            remIntra = FALSE) {
  vars_ext_ <- rlang::syms(c(vars, "source", "target", "tag"))
  vars_src_ <- rlang::syms(c(vars, "source"))
  grp_      <- rlang::syms(groupBy)
  srt_      <- rlang::syms(sortBy)

  mk <- tbl |>
    dplyr::rename(clus = dplyr::all_of(cVar)) |>
    dplyr::group_by(!!!grp_) |>
    dplyr::arrange(!!!srt_) |>
    dplyr::mutate(
      source = dplyr::lag(clus, n = 1L),
      target = clus,
      tag    = paste0(source, "_", target)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(source))

  if (remIntra) {
    n_clus <- max(as.integer(mk$clus), na.rm = TRUE)
    inward <- paste0(seq(0L, n_clus), "_", seq(0L, n_clus))
    mk     <- dplyr::filter(mk, !tag %in% inward)
  }

  mk |>
    dplyr::group_by(!!!vars_ext_) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::group_by(!!!vars_src_) |>
    dplyr::mutate(tot = sum(n), nCount = n / tot) |>
    dplyr::ungroup() |>
    dplyr::group_by(tag, source, target) |>
    tidyr::nest()
}
