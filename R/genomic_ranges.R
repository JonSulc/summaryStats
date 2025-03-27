#' @import data.table

validate_start_less_than_end <- function(object) {
  stopifnot(all(c("start", "end") %in% colnames(object)))
  object[
    !is.na(start) & !is.na(end),
    {
      if (0 < .N) {
        if (any(end < start)) {
          stop("All 'start' values must be lower than (or equal to) 'end' values.")
        }
      }
    }
  ]
  invisible(object)
}

#' @export
validate_genomic_ranges <- function(object) {
  if (!all(c("start", "end") %in% colnames(object))) {
    stop("Genomic ranges must include columns 'start' and 'end'.",
         " Current columns: ",
         paste(colnames(object), collapse = ", "))
  }

  validate_start_less_than_end(object)
  validate_genomic_positions(object)
}


#' @export
is_genomic_ranges <- function(object) inherits(object, "genomic_ranges")

as_genomic_ranges <- function(object) {
  as_genomic_positions(object)
  data.table::setattr(object,
                      "class",
                      unique(c("genomic_ranges", class(object))))
}

#' @export
new_genomic_ranges <- function(
    dt,
    build = c("b38", "b37", "b36"),
    ...
) {
  new_genomic_positions(
    dt,
    build = build,
    position_columns = c("start", "end"),
    ...,
    class = "genomic_ranges"
  )
}
