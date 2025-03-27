#' @import data.table

#' @export
validate_chr <- function(x) UseMethod("validate_chr")
#' @export
validate_pos <- function(x) UseMethod("validate_pos")
#' @export
validate_start_end <- function(x) UseMethod("validate_start_end")
are_only_positive_integers <- function(x, cnames) {
  stopifnot(data.table::is.data.table(x))
  stopifnot(all(cnames %in% names(x)))
  x_is_positive <- x[
    ,
    !sapply(mget(cnames), is_any_non_positive)
  ]
  x_is_integer <- sapply(cnames, \(cname) is.integer(x[[cname]])) |>
    unlist()
  if (!(all(c(x_is_integer & x_is_positive)))) {
    stop("All position values must be positive integers",
         " (",
         sprintf("'%s'", cnames[!x_is_integer | !x_is_positive]) |> paste(collapse = ", "),
         " contain non-conform values)")
  }
  invisible(x)
}

is_any_non_positive <- function(x, default = FALSE) {
  if (all(is.na(x))) return(default)
  x[!is.na(x)] <= 0
}

convert_to_dt <- function(dt) {
  if (data.table::is.data.table(dt)) return(data.table::copy(dt))
  data.table::as.data.table(dt)
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
  stopifnot("chr" %in% colnames(summary_stats))
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

assign_correct_pos <- function(genomic_positions, position_columns) {
  stopifnot(all(position_columns %in% colnames(genomic_positions)))
  if (any(genomic_positions[,
                            sapply(.SD, is_any_non_positive),
                            .SDcols = position_columns])) {
    warning("Some position values were not positive, converted to NA")
  }
  genomic_positions[
    ,
    (position_columns) := lapply(.SD, make_positive_integer),
    .SDcols = position_columns
  ][]
  data.table::setattr(genomic_positions, "pos", position_columns)
}

make_positive_integer <- function(x) {
  x <- as.integer(x)
  x[x <= 0L] <- NA_integer_
  x
}
