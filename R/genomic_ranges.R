#' @import data.table

#' @export
is_genomic_ranges <- function(object) inherits(object, "genomic_ranges")

#' @export
new_genomic_ranges <- function(
    chromosome,
    start = NULL,
    end   = start,
    build = c("b38", "b37", "b36"),
    required_columns = c("chr", "start", "end")
) {
  build <- match.arg(tolower(build), build)

  if (is.data.frame(chromosome)) {
    if (is.data.table(chromosome)) {
      genomic_ranges <- data.table::copy(chromosome)
    } else {
      genomic_ranges <- data.table::as.data.table(chromosome)
    }
    assign_standardized_names(genomic_ranges)
  } else {
    genomic_ranges <- data.table::data.table(
      chr   = chromosome,
      start = as.integer(start),
      end   = as.integer(end)
    )
  }

  if (!all(required_columns %chin% names(genomic_ranges))) {
    stop("Missing required columns:",
         required_columns[!required_columns %chin% names(genomic_ranges)] |>
           paste(collapse = ", "))
  }

  if (nrow(genomic_ranges) == 0) {
    if (!ignore_warning)
      warning("No variants in data.")
  }
  genomic_ranges
}
