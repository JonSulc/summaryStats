#' @import data.table

#' @export
is_genomic_ranges <- function(object) inherits(object, "genomic_ranges")

as_genomic_ranges <- function(object) {
  data.table::setattr(object,
                      "class",
                      unique(c("genomic_ranges", class(object))))
}

#' @export
new_genomic_ranges <- function(
    chr,
    start = NULL,
    end   = start,
    build = c("b38", "b37", "b36"),
    required_columns = c("chr", "start", "end")
) {
  build <- match.arg(tolower(build), build)

  if (is.data.frame(chr)) {
    if (is.data.table(chr)) {
      genomic_ranges <- data.table::copy(chr)
    } else {
      genomic_ranges <- data.table::as.data.table(chr)
    }
    assign_standardized_names(genomic_ranges)
  } else {
    genomic_ranges <- data.table::data.table(
      chr   = chr,
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
  data.table::setattr(genomic_ranges,
                      "class",
                      unique(c("genomic_ranges", class(genomic_ranges))))
  genomic_ranges
}
