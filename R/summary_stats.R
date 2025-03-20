#' @import data.table
#' @importFrom stringr str_match str_to_upper str_detect

#' @export
validate_summary_stats <- function(
  summary_stats
) {
  # stopifnot(inherits(summary_stats, "summary_stats"))
  stopifnot(inherits(summary_stats, "data.table"))

  # Check missing columns
  missing_columns <- setdiff(
    c("variant_id", "chr", "pos", "ref", "alt", "effect", "effect_se", "pval"),
    colnames(summary_stats)
  )
  if (length(missing_columns) != 0) {
    stop("Required columns are missing: ",
         paste(missing_columns, collapse = ", "))
  }

  # Check variant IDs
  if (any(duplicated(summary_stats$variant_id)))
    stop("Some variant_id values are not unique")

  # Check chromosome
  summary_stats[
    !stringr::str_detect(chr,
                         "^chr([0-9]{1,2}|[XY]|[mM][tT])$"),
    {
      if (.N != 0) {
        stop("Chromosome format does not match expect (e.g., 'chr1'): ",
             head(unique(chr)) |>
               paste(collapse = ", "))
      }
    }
  ]

  # Check position
  if (!is.integer(summary_stats$pos))
    stop("Position must be an integer")
  if (!summary_stats[!is.na(pos), all(0 < pos)])
    stop("Position must be a positive integer")

  if (summary_stats[
    !is.na(pos),
    any(pos <= 0)
  ])
    stop("Non-positive pos in summary statistics")

  # Check alleles
  summary_stats[
    !stringr::str_detect(ref, "^[ACGT]+$"),
    if (.N != 0)
      stop("Invalid ref, e.g., ",
           head(ref) |>
             paste(collapse = ", "))
  ]
  summary_stats[
    !stringr::str_detect(alt, "^[ACGT]+$"),
    if (.N != 0)
      stop("Invalid alt, e.g., ",
           head(alt) |>
             paste(collapse = ", "))
  ]
  summary_stats[
    ref == alt,
    if (.N != 0)
      stop("Some alleles are identical for ref and alt")
  ]

  # Check p-value
  if (!is.numeric(summary_stats$pval))
    stop("Non-numeric pval variable")
  if (any(summary_stats[, pval < 0]))
    stop("Values of pval cannot be below 0")
  if (any(summary_stats[, 1 < pval]))
    stop("Values of pval cannot exceed 1")

  # Check effect size
  if (!is.numeric(summary_stats$effect))
    stop("Non-numeric effect variable")

  # Check effect SE
  if (!is.numeric(summary_stats$effect_se))
    stop("Non-numeric effect_se variable")
  if (any(summary_stats$effect_se == 0))
    stop("Some effect_se are 0")

  # if (data.table:::selfrefok(dt) != 1) {
  #   warning("Internal self reference doesn't match")
  # }

  invisible(summary_stats)
}

set_as_summary_stats <- function(
  dt
) {
  stopifnot(data.table::is.data.table(dt))
  stopifnot(data.table:::selfrefok(dt)==1)
  data.table::setattr(dt, "class", unique(c("summary_stats", class(dt))))
}


#' @export
is_summary_stats <- function(object) inherits(object, "summary_stats")

#' @export
empty_summary_stats <- function(
  build = NULL
) {
  summary_stats <- data.table::data.table(
    variant_id = character(0),
    chr = character(0),
    pos = integer(0),
    ref = character(0),
    alt = character(0),
    effect = numeric(0),
    effect_se = numeric(0),
    pval = numeric(0)
  ) |>
    set_as_summary_stats()
  if (!is.null(build)) {
    assign_current_as_build(summary_stats, build)
  }
  summary_stats
}

#' @export
new_summary_stats <- function(
  dt,
  build = get_build(dt),
  all_colnames = summary_stats_column_names,
  complete_missing_stats = TRUE,
  ignore_warning = FALSE
) {
  if (!data.table::is.data.table(dt)) {
    summary_stats <- data.table::as.data.table(dt)
    set_as_summary_stats(summary_stats)
  } else {
    summary_stats <- data.table::copy(dt)
    set_as_summary_stats(summary_stats)
  }

  assign_standardized_names(
    summary_stats,
    all_colnames = all_colnames
  )

  if (nrow(summary_stats) == 0) {
    if (!ignore_warning)
      warning("No variants in data.")

    return(
      empty_summary_stats(build = build)
    )
  }

  if (complete_missing_stats) {
    complete_missing_statistics(summary_stats)
  }

  fix_missing_chr_pos_ref_alt(summary_stats)

  capitalize_alleles(summary_stats)

  if (!("variant_id" %in% names(summary_stats))) {
    assign_variant_id_from_chr_pos_ref_alt(summary_stats, build = build)
  }

  assign_current_as_build(summary_stats, build)

  data.table::setcolorder(summary_stats, c("variant_id", "chr", "pos", "ref", "alt"))

  validate_summary_stats(summary_stats)

  summary_stats
}

complete_missing_statistics <- function(
    summary_stats
) {
  if (!"effect" %in% colnames(summary_stats)) {
    if (!"odds_ratio" %in% colnames(summary_stats)
        | !"pval" %in% colnames(summary_stats)) {
      stop("No effect size column found.")
    }
    warning("No effect size column found, approximating standard effects from ",
            "odds ratio (PMID: 11113947), standard errors from the p-value.")
    summary_stats[
      ,
      effect := log(odds_ratio) * sqrt(3)/pi
    ][
      ,
      effect_se := -effect / qnorm(pval/2)
    ][]
  } else if (!"effect_se" %in% colnames(summary_stats)) {
    if (!"pval" %in% colnames(summary_stats)) {
      stop("No p-value column found, unable to infer standard error.")
    }
    summary_stats[
      ,
      effect_se := abs(effect / qnorm(pval/2))
    ][]
  } else if (!"pval" %in% colnames(summary_stats)) {
    summary_stats[
      ,
      pval := 2 * pnorm(-abs(effect/effect_se))
    ][]
  }
  invisible(summary_stats)
}

fix_missing_chr_pos_ref_alt <- function(
    summary_stats
) {
  if (all(c("chr", "pos", "ref", "alt") %in% colnames(summary_stats))) {
    assign_correct_chr(summary_stats)
    summary_stats[, pos := as.integer(pos)]
    return(invisible(summary_stats))
  }
  if (all(is_correct_variant_id_format(summary_stats$variant_id))) {
    return(invisible(assign_chr_pos_ref_alt_from_variant_id(summary_stats)))
  }
  assign_correct_chr(summary_stats)
  infer_and_assign_position(summary_stats)
  # Some MVP summary statistics do not directly list the effect allele...
  infer_and_assign_ref_alt_alleles(summary_stats)
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
assign_variant_id_from_chr_pos_ref_alt <- function(summary_stats, ...) {
  summary_stats[
    ,
    variant_id := get_variant_id_from_chr_pos_ref_alt(summary_stats, ...)
  ][]
  invisible(summary_stats)
}
infer_and_assign_position <- function(summary_stats) {
  if ("pos" %in% colnames(summary_stats)) {
    if (!is.numeric(summary_stats$pos))
      summary_stats[, pos := as.integer(pos)][]
    return(invisible(summary_stats))
  }

  # In some CHARGE files, there is no position column but it is encoded in other
  # columns
  summary_stats[
    ,
    pos := infer_position(summary_stats)
  ][] |>
    invisible()
}


capitalize_alleles <- function(summary_stats) {
  summary_stats[
    ,
    ref := stringr::str_to_upper(ref)
  ][
    ,
    alt := stringr::str_to_upper(alt)
  ][]
}


assign_chr_pos_ref_alt_from_variant_id <- function(summary_stats) {
  parsed_variant_ids <- parse_variant_id(summary_stats$variant_id)

  new_columns <- setdiff(names(parsed_variant_ids),
                         names(summary_stats))

  if (length(new_columns) == 0)
    return(invisible(summary_stats))

  summary_stats[
    ,
    (new_columns) := parsed_variant_ids[, mget(new_columns)]
  ][]
  invisible(summary_stats)
}


assign_current_as_build <- function(
    summary_stats,
    current_build = get_build(summary_stats),
    columns = c("variant_id", "id", "chr", "pos", "start", "end") |>
      intersect(names(summary_stats)),
    force = FALSE
) {
  if (has_build_version(summary_stats, current_build, columns = columns) & !force) {
    data.table::setattr(summary_stats, "build", current_build)
    return(invisible(summary_stats))
  }

  if (!"variant_id" %in% names(summary_stats)) {
    assign_variant_id_from_chr_pos_ref_alt(summary_stats, build = current_build)
  }

  summary_stats[
    ,
    get_build_colnames(current_build, columns) := mget(columns)
  ][]
  data.table::setattr(summary_stats, "build", current_build)
  invisible(summary_stats)
}


infer_and_assign_ref_alt_alleles <- function(
    summary_stats
) {
  if (all(c("ref", "alt") %in% colnames(summary_stats))) {
    return(summary_stats)
  }

  alleles <- infer_ref_alt_alleles(summary_stats)

  summary_stats[
    ,
    c("ref", "alt") := alleles[, .(ref, alt)]
  ][]

  if ("assumed_ea_col" %in% colnames(alleles)) {
    summary_stats[, assumed_ea_col := alleles$assumed_ea_col][]
  }

  invisible(summary_stats)
}
