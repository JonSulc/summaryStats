#' @import data.table

assign_standardized_names <- function(
    genomic_ranges,
    all_colnames = summary_stats_column_names,
    ...
) {
  old_names <- get_standard_colnames(
    names(genomic_ranges),
    all_colnames = all_colnames,
    ...
  )

  if (length(unlist(old_names)) == 0) return(invisible(genomic_ranges))

  data.table::setnames(
    genomic_ranges,
    old = unlist(old_names),
    new = names(unlist(old_names)),
    skip_absent = TRUE
  )
}
