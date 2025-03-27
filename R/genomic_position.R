#' @import data.table

#' @export
validate_chr.genomic_positions <- function(genomic_positions) {
  stopifnot("chr" %in% colnames(genomic_positions))
  if (!is.character(genomic_positions$chr)) {
    stop("Chromosome ('chr') must be in character format starting with 'chr',",
         " e.g., 'chr1'")
  }
  if (
    genomic_positions[
      !is.na(chr),
      any(substr(chr, 1, 3) != "chr")
    ]
  ) {
    stop("All chromosome values ('chr') must start with 'chr', e.g., 'chr1'")
  }
  invisible(genomic_positions)
}

#' @export
validate_pos.genomic_positions <- function(genomic_positions) {
  if (is.null(attr(genomic_positions, "pos"))) {
    stop("Object does not appear to have a valid position column",
         " (attr(object, 'pos') is NULL).")
  }
  pos_columns <- attr(genomic_positions, "pos")
  if (!all(pos_columns %in% colnames(genomic_positions))) {
    stop("Columns ",
         sprintf("'%s'",
                 pos_columns[!pos_columns %in% colnames(genomic_positions)]) |>
                               paste(collapse = ", "),
         " are indicated as position columns but are absent.")
  }
  are_only_positive_integers(genomic_positions, pos_columns)
  invisible(genomic_positions)
}

#' @export
validate_genomic_positions <- function(object) {
  if (!"chr" %in% colnames(object)) {
    stop("Genomic ranges must include a column for 'chr'. Current columns: ",
         paste(colnames(object), collapse = ", "))
  }

  validate_chr(object)
  validate_pos(object)
}


#' @export
is_genomic_positions <- function(object) inherits(object, "genomic_positions")

as_genomic_positions <- function(object) {
  data.table::setattr(object,
                      "class",
                      unique(c("genomic_positions", class(object))))
}

#' @export
new_genomic_positions <- function(
  data,
  build = c("b38", "b37", "b36"),
  position_columns = intersect(c("pos", "start", "end"),
                               colnames(data)),
  ...,
  class = character()
) {
  build <- match.arg(tolower(build), build)

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }

  if (is.data.table(data)) {
    genomic_positions <- data.table::copy(data)
  } else {
    genomic_positions <- data.table::as.data.table(data)
  }
  assign_standardized_names(genomic_positions)

  assign_correct_chr(genomic_positions)
  assign_correct_pos(genomic_positions, position_columns)

  if (nrow(genomic_positions) == 0) {
    warning("No observations in data.")
  }

  data.table::setattr(genomic_positions,
                      "build",
                      build)

  data.table::setattr(genomic_positions,
                      "class",
                      unique(c(class, "genomic_positions", class(genomic_positions))))

  genomic_positions
}
