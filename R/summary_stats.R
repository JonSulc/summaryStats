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
    c("variant_id", "chr", "pos", "ref", "alt", "effect", "pval", "effect_se"),
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
  if (!is.numeric(summary_stats$pos))
    stop("Position must be numeric")
  if (!summary_stats[, all(pos %% 1 == 0)])
    stop("Non-integer pos variable")
  if (any(summary_stats$pos <= 0))
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

new_summary_stats <- function(
  dt,
  build = get_build(dt),
  all_colnames = summary_stats_column_names,
  complete_missing_stats = TRUE,
  ignore_warning = FALSE,
  required_columns = c("chr", "pos", "ref", "alt", "effect", "pval")
) {
  if (!data.table::is.data.table(dt)) {
    summary_statistics <- data.table::as.data.table(dt)
    set_as_summary_stats(summary_statistics)
  } else {
    summary_statistics <- data.table::copy(dt)
    set_as_summary_stats(summary_statistics)
  }

  assign_standardized_names(
    summary_statistics,
    all_colnames = all_colnames
  )

  if (nrow(summary_statistics) == 0) {
    if (!ignore_warning)
      warning("No variants in data.")
    missing_columns <- setdiff(
      required_columns,
      colnames(summary_statistics)
    )
    if (length(missing_columns) != 0)
      summary_statistics[
        ,
        (missing_columns) :=
          data.table::data.table(
            chr = character(0),
            pos = numeric(0),
            ref = character(0),
            alt = character(0),
            effect = numeric(0),
            pval = numeric(0)
          )[, .SD, .SDcols = missing_columns]
      ]
  }

  if (complete_missing_stats) {
    complete_missing_statistics(summary_statistics)
  }

  fix_missing_chr_pos_ref_alt(summary_statistics)

  # Some MVP summary statistics do not directly list the effect allele...
  infer_ref_alt_alleles(summary_statistics)
  capitalize_alleles(summary_statistics)

  if (!any(c("chr", "pos", "ref", "alt") %in% names(summary_statistics))) {
    assign_chr_pos_ref_alt_from_variant_id(summary_statistics)
  }

  assign_correct_pos(summary_statistics)
  assign_correct_chr(summary_statistics)

  if (!("variant_id" %in% names(summary_statistics))) {
    assign_variant_id_from_chr_pos_ref_alt(summary_statistics, build = build)
  }

  assign_current_as_build(summary_statistics, build)

  data.table::setcolorder(summary_statistics, c("variant_id", "chr", "pos", "ref", "alt"))

  summary_statistics
}

assign_standardized_names <- function(
    summary_stats,
    all_colnames = summary_stats_column_names,
    ...
) {
  old_names <- get_standard_colnames(
    names(summary_stats),
    all_colnames = all_colnames,
    ...
  )

  if (length(unlist(old_names)) == 0) return(invisible(summary_stats))

  data.table::setnames(
    summary_stats,
    old = unlist(old_names),
    new = names(unlist(old_names)),
    skip_absent = TRUE
  )
}

complete_missing_statistics <- function(
    summary_stats
) {
  if (!"effect" %in% colnames(summary_stats)) {
    if (!"odds_ratio" %in% colnames(summary_stats)
        | !"pval" %in% colnames(summary_stats)) {
      stop("No effect size column found, unable to use as summary statsistics.")
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
      effect_se := -abs(effect) / qnorm(pval/2)
    ][]
  } else if (!"pval" %in% colnames(summary_stats)) {
    summary_stats[
      ,
      pval := 2 * pnorm(-abs(effect/effect_se))
    ][]
  }
}

fix_missing_chr_pos_ref_alt <- function(
    summary_stats
) {
  if (all(c("chr", "pos", "ref", "alt") %in% colnames(summary_stats))) {
    return(invisible())
  }
  if (variant_id_is_correct_format(summary_stats)) {
    return(invisible(assign_chr_pos_ref_alt_from_variant_id(summary_stats)))
  }
  infer_position(summary_stats)
  infer_ref_alt_alleles(summary_stats)
}

variant_id_is_correct_format <- function(
    summary_stats
) {
  if (!"variant_id" %in% colnames(summary_stats))
    return(FALSE)

  grepl("^(chr[0-9xX]+)_([0-9]+)_([^_]+)_([^_]+)_.*$", summary_stats$variant_id[1])
}


# In some CHARGE files, there is no position column but it is encoded in other
# columns
infer_position <- function(summary_stats) {
  if ("pos" %in% colnames(summary_stats)) {
    return(invisible(summary_stats))
  }
  position_encoding_columns <- intersect(
    summary_stats_extra_column_names$position_encoding$column_name,
    names(summary_stats)
  )
  if (length(position_encoding_columns) == 0) {
    stop("Unable to infer position from existing columns")
  }

  summary_stats[
    ,
    pos := stringr::str_match(
      get(position_encoding_columns[1]),
      summary_stats_extra_column_names$position_encoding[
        position_encoding_columns[1],
        regex,
        on = "column_name"
      ]
    )[, 2] |>
      as.integer()
  ][] |>
    invisible()
}


# In some MVP summary statistics, alleles are listed as 1 and 2 with a separate
# column indicating which is the effect allele (already renamed alt in
# loading). The ref allele is set as the other one here.
infer_ref_alt_alleles <- function(
  summary_stats
) {
  if ("ref" %in% colnames(summary_stats)) {
    return(invisible(summary_stats))
  }
  # ref is inferred from variant_id at a later stage
  if ("variant_id" %in% colnames(summary_stats)) {
    return()
  }

  cat("Reference and alternate allele columns not found, attempting to infer... ")
  allele_1 <- intersect(
    summary_stats_extra_column_names$allele_1,
    colnames(summary_stats)
  )
  allele_2 <- intersect(
    summary_stats_extra_column_names$allele_2,
    colnames(summary_stats)
  )
  if (length(allele_1) == 0 |
      length(allele_2) == 0) {
    stop("Failed to determine reference and alternate alleles.")
  }

  if (length(allele_1) != 1) {
    warning("More than one column labeled as allele 1:\n",
            paste(allele_1, collapse = ", "), "\n",
            "Using ", allele_1[1], " as allele 1, assuming ", allele_1[2],
            " is the effect (coded) allele")
    setnames(summary_stats, old = allele_1[2], new = "alt")
    allele_1 <- allele_1[1]
  }
  if (length(allele_2) != 1) {
    warning("More than one column labeled as allele 2:\n",
            paste(allele_2, collapse = ", "), "\n",
            "Using ", allele_2[1], " as allele 2, assuming ", allele_2[2],
            " is the effect (coded) allele")
    setnames(summary_stats, old = allele_2[2], new = "alt")
    allele_2 <- allele_2[1]
  }

  if (!"alt" %in% names(summary_stats)) {
    message("No effect allele specified")
    warning("No effect allele specified in summary statistics file,",
            " assuming effect allele is ", allele_2)
    summary_stats[
      ,
      assumed_ea_col := allele_2
    ][]
    data.table::setnames(
      summary_stats,
      old = c(allele_1, allele_2),
      new = c("ref", "alt")
    )
    return(summary_stats)
  }
  summary_stats[
    ,
    ref := fifelse(
      get(allele_1) == alt,
      get(allele_2),
      get(allele_1)
    )
  ][]
  cat("Success!\n")
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


get_build <- function(summary_stats) {
  if (nrow(summary_stats) == 0) {
    if ("variant_id_b38" %in% names(summary_stats)) return("b38")
    if ("variant_id_b37" %in% names(summary_stats)) return("b37")
    if ("variant_id_b36" %in% names(summary_stats)) return("b36")
    stop("Build could not be inferred from empty summary_stats.")
  }

  if (!("variant_id" %in% names(summary_stats)))
    stop("No variant_id column, unable to infer genome assembly.")

  build <- summary_stats[
    !is.na(variant_id)
  ][
    1,
    stringr::str_match(variant_id, "_(b[0-9]+)$")
  ][, 2]

  if (!is.na(build)) return(build)

  for (build in c("b38", "b37", "b36")) {
    if (paste0("variant_id_", build) %in% names(summary_stats)
        & all(summary_stats[, .(chr, pos, ref, alt)] ==
              summary_stats[, mget(paste0(c("chr_", "pos_", "ref_", "alt_"), build))],
              na.rm = TRUE))
      return(build)
  }

  stop("Unable to infer genome assembly from variant ids.")
}

assign_current_as_build <- function(
    summary_stats,
    current_build = get_build(dt),
    columns = c("variant_id", "id", "chr", "pos", "start", "end") |>
      intersect(names(dt)),
    force = FALSE
) {
  if (has_build_version(summary_stats, current_build, columns = columns) & !force) {
    return(invisible(summary_stats))
  }

  if (!"variant_id" %in% names(summary_stats)) {
    assign_variant_id_from_chr_pos_ref_alt(summary_stats, build = current_build)
  }

  summary_stats[
    ,
    get_build_colnames(current_build, columns) := mget(columns)
  ][]
  invisible(summary_stats)
}

has_build_version <- function(summary_stats, build, columns) {
  all(
    get_build_colnames(build, columns = columns) %in% names(summary_stats)
  )
}

get_build_colnames <- function(
    build,
    columns = c("chr", "pos", "variant_id", "ref", "alt")
) {
  sprintf("%s_%s",
          columns,
          build)
}

assign_correct_pos <- function(summary_stats) {
  if (!is.numeric(summary_stats$pos))
    summary_stats[, pos := as.numeric(pos)][]
  invisible(summary_stats)
}
assign_correct_chr <- function(summary_stats) {
  if (is.numeric(summary_stats$chr)) {
    summary_stats[
      ,
      chr := as.character(chr) # Force conversion in case DT is empty
    ][
      ,
      chr := paste0("chr", chr)
    ][]
  } else {
    summary_stats[
      !substr(chr, 1, 3) == "chr",
      chr := paste0("chr", chr)
    ][]
  }
  invisible(summary_stats)
}


get_standard_colnames <- function(
    data_colnames,
    all_colnames = summary_stats_column_names,
    required_columns = NULL
) {
  summary_stats_colnames <- lapply(
    all_colnames,
    intersect,
    data_colnames
  )

  if (any(sapply(summary_stats_colnames[required_columns], length) == 0)) {
    warning("Some required columns are missing: ",
            summary_stats_colnames[required_columns][
              sapply(summary_stats_colnames[required_columns], length) == 0
            ] |> names() |> paste(collapse = ", "))
  }

  if (any(sapply(summary_stats_colnames, length) > 1)) {
    warning(sprintf(
      "Argument 'all_colnames' matched multiple columns for %s, using first match.",
      do.call(paste,
              list(names(all_colnames)[sapply(summary_stats_colnames, length) > 1],
                   sep = ", "))))
    summary_stats_colnames <- lapply(summary_stats_colnames, \(x) x[1])
  }

  summary_stats_colnames[!is.na(summary_stats_colnames)]
}
