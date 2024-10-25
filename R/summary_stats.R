#' @import data.table
#' @importFrom stringr str_match str_to_upper

as_summary_stats <- function(
    summary_stats,
    build = get_build(summary_stats),
    all_colnames = summary_stats_column_names,
    complete_missing_stats = TRUE,
    ignore_warning = FALSE
) {
  if (!data.table::is.data.table(summary_stats)) {
    summary_stats <- data.table::as.data.table(summary_stats)
  } else {
    summary_stats <- data.table::copy(summary_stats)
  }

  assign_standardized_names(summary_stats)

  if (nrow(summary_stats) == 0) {
    if (!ignore_warning)
      warning("No variants in data.")
    missing_columns <- setdiff(
      c("chr", "pos", "ref", "alt", "effect", "effect_se"),
      colnames(summary_stats)
    )
    if (length(missing_columns) != 0)
      summary_stats[
        ,
        (missing_columns) :=
          data.table::data.table(
            chr = character(0),
            pos = numeric(0),
            ref = character(0),
            alt = character(0),
            effect = numeric(0),
            effect_se = numeric(0)
          )[, .SD, .SDcols = missing_columns]
      ]
  }

  if (complete_missing_stats) {
    complete_missing_statistics(summary_stats)
  }

  fix_missing_chr_pos_ref_alt(summary_stats)

  # Some MVP summary statistics do not directly list the effect allele...
  infer_ref_alt_alleles(summary_stats)
  capitalize_alleles(summary_stats)

  if (!any(c("chr", "pos", "ref", "alt") %in% names(summary_stats))) {
    assign_chr_pos_ref_alt_from_variant_id(summary_stats)
  }

  assign_correct_pos(summary_stats)
  assign_correct_chr(summary_stats)

  if (!("variant_id" %in% names(summary_stats))) {
    assign_variant_id_from_chr_pos_ref_alt(summary_stats, build = build)
  }

  assign_current_as_build(summary_stats, build)

  summary_stats
}

assign_standardized_names <- function(
    summary_stats,
    all_colnames = summary_stats_column_names
) {
  old_names <- get_colnames(summary_stats,
                            all_colnames = all_colnames)

  if (length(unlist(old_names)) == 0) return()

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
      stop("No effect size column found, unable to use summary statistics for MR")
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
    summary_stats_extra_column_names()$position_encoding$column_name,
    names(summary_stats)
  )
  if (length(position_encoding_columns) == 0) {
    stop("Unable to infer position from existing columns")
  }

  summary_stats[
    ,
    pos := stringr::str_match(
      get(position_encoding_columns[1]),
      summary_stats_extra_column_names()$position_encoding[
        position_encoding_columns[1],
        regex,
        on = "column_name"
      ]
    )[, 2] |>
      as.integer()
  ][] |>
    invisible()
}


# In MVP summary statistics, alleles are listed as 1 and 2 with a separate
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
    summary_stats_extra_column_names()$allele_1,
    colnames(summary_stats)
  )
  allele_2 <- intersect(
    summary_stats_extra_column_names()$allele_2,
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
    stop("Build could not be inferred from empty variant_table.")
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
    current_build = get_build(summary_stats),
    columns = c("variant_id", "chr", "pos", "ref", "alt"),
    force = FALSE
) {
  if (has_build_version(summary_stats, current_build) & !force) {
    return(invisible(summary_stats))
  }
  summary_stats[
    ,
    get_build_colnames(current_build, columns) := mget(columns)
  ][]
  invisible(summary_stats)
}

has_build_version <- function(summary_stats, build) {
  all(
    get_build_colnames(build) %in% names(summary_stats)
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


get_colnames <- function(summary_stats,
                         all_colnames = summary_stats_column_names) {
  summary_stats_colnames <- lapply(
    all_colnames,
    intersect,
    names(summary_stats)
  )

  if (length(summary_stats_colnames$effect) == 0) {
    warning("Missing 'effect' column")
  }
  if (length(summary_stats_colnames$variant_id) == 0
      & (length(summary_stats_colnames$chr) == 0
         | length(summary_stats_colnames$pos) == 0)) {
    warning("Columns required for MR are missing:\n",
            paste(
              c("chr", "pos")[sapply(summary_stats_colnames[c("chr", "pos")], length) == 0],
              collapse = ", "
            ))
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
