#' @import data.table

#' @export
validate_chr.genomic_ranges <- function(genomic_ranges) {
  stopifnot("chr" %in% colnames(genomic_ranges))
  if (!is.character(genomic_ranges$chr)) {
    stop("Chromosome ('chr') must be in character format starting with 'chr',",
         " e.g., 'chr1'")
  }
  if (
    genomic_ranges[
      !is.na(chr),
      any(substr(chr, 1, 3) != "chr")
    ]
  ) {
    stop("All chromosome values ('chr') must start with 'chr', e.g., 'chr1'")
  }
  invisible(genomic_ranges)
}
#' @export
validate_pos.genomic_ranges <- function(
  genomic_ranges,
  cnames = intersect(c("pos", "start", "end"), colnames(genomic_ranges))
) {
  if (!all(c("pos") %in% colnames(genomic_ranges))) {
    return(invisible(genomic_ranges))
  }
  is_positive_integer(genomic_ranges, "pos")
  invisible(genomic_ranges)
}
#' @export
validate_start_end.genomic_ranges <- function(
  genomic_ranges
) {
  if (!all(c("start", "end") %in% colnames(genomic_ranges))) {
    return(invisible(genomic_ranges))
  }
  is_positive_integer(genomic_ranges, c("start", "end"))
  if (
    genomic_ranges[
      !is.na(start) & !is.na(end),
      any(end < start)
    ]
  ) {
    stop("For all entries, the start must be lower than the end.")
  }
}

#' @export
validate_genomic_ranges <- function(object) {
  if (!"chr" %in% colnames(object)
      | (!all(c("start", "end") %in% colnames(object))
         & !"pos" %in% colnames(object))) {
    stop("Genomic ranges must include columns for 'chr' and either 'pos' or",
         " both 'start' and 'end'. Current columns: ",
         paste(colnames(object), collapse = ", "))
  }

  validate_chr(object)
  validate_pos(object)
  validate_start_end(object)

  invisible(object)
}


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
      chr   = as.character(chr),
      start = as.integer(start),
      end   = as.integer(end)
    )
  }

  if (!all(required_columns %chin% names(genomic_ranges))) {
    stop("Missing required columns:",
         required_columns[!required_columns %chin% names(genomic_ranges)] |>
           paste(collapse = ", "))
  }

  assign_correct_chr(genomic_ranges)
  assign_correct_start_end(genomic_ranges)

  if (nrow(genomic_ranges) == 0) {
    if (!ignore_warning)
      warning("No variants in data.")
  }
  data.table::setattr(genomic_ranges,
                      "class",
                      unique(c("genomic_ranges", class(genomic_ranges))))

  genomic_ranges
}

assign_correct_start_end <- function(
  genomic_ranges
) {
  stopifnot(data.table::is.data.table(genomic_ranges))
  stopifnot(all(c("start", "end") %in% colnames(genomic_ranges)))

  genomic_ranges[
    ,
    c("start", "end") := lapply(.(start, end), as.integer)
  ][]

  swapped <- genomic_ranges[, end < start]
  if (any(swapped)) {
    warning("Some entries have 'end' values smaller than 'start', e.g.:\n",
            genomic_ranges[swapped][1, ], "\n",
            "Swapping 'start' and 'end' in those cases...")
    genomic_ranges[swapped, c("start", "end") := .(end, start)]
  }

  invisible(genomic_ranges)
}
