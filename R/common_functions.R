#' @import data.table

#' @export
validate_chr <- function(x) UseMethod("validate_chr")
#' @export
validate_pos <- function(x) UseMethod("validate_pos")
#' @export
validate_start_end <- function(x) UseMethod("validate_start_end")
is_positive_integer <- function(x, cnames) {
  x_is_positive <- lapply(cnames, \(cname) {
    x[
      !is.na(get(cname)),
      all(0 < get(cname))
    ]
  }) |>
    unlist()
  x_is_integer <- sapply(cnames, \(cname) is.integer(x[[cname]])) |>
    unlist()
  if (!(all(c(x_is_integer & x_is_positive)))) {
    stop("All position values (",
         sprintf("'%s'", cnames) |> paste(collapse = ", "),
         ") must be positive integers")
  }
  invisible(x)
}

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

assign_correct_chr <- function(summary_stats) {
  if (!is.character(summary_stats$chr)) {
    summary_stats[
      ,
      chr := as.character(chr) # Force conversion in case DT is empty or NA
    ][]
  }
  summary_stats[
    substr(chr, 1, 3) != "chr",
    chr := paste0("chr", chr)
  ][]
  invisible(summary_stats)
}
