#' @import data.table
#' @importFrom stringr str_match

is_correct_variant_id_format <- function(
    variant_id
) {
  if (is.null(variant_id))
    return(FALSE)

  grepl("^(chr[0-9xX]+)_([0-9]+)_([ACGT]+)_([ACGT]+)_.*$", variant_id)
}
get_variant_id_from_chr_pos_ref_alt <- function(summary_stats,
                                                build = "b38") {
  if (!all(c("chr", "pos", "alt", "ref") %in% colnames(summary_stats))) {
    stop("Missing columns to create variant IDs: ",
         paste(c("chr", "pos", "alt", "ref")[
           !c("chr", "pos", "alt", "ref") %in% colnames(summary_stats)
         ],
         collapse = ", "))
  }
  summary_stats[
    ,
    fifelse(
      is.na(chr) | is.na(pos),
      NA_character_,
      do.call(sprintf,
              c(.("%s_%s_%s_%s_%s"), .SD, .(build)))
    ),
    .SDcols = c("chr", "pos", "ref", "alt")
  ]
}
parse_variant_id <- function(variant_id) {
  if (!all(is_correct_variant_id_format(variant_id))) {
    stop("Unable to parse missing information from variant_id")
  }
  parsed <- stringr::str_match(
    variant_id,
    "^(chr[^_]+)_([0-9]+)_([acgtACGT]+)_([acgtACGT]+)_[^_]+$"
  )[, -1, drop = FALSE] |>
    data.table::as.data.table() |>
    setNames(c("chr", "pos", "ref", "alt"))
  parsed[
    ,
    pos := as.integer(pos)
  ][]
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


get_build <- function(summary_stats) {
  if (!is.null(attr(summary_stats, "build"))) return(attr(summary_stats, "build"))
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

  # for (build in c("b38", "b37", "b36")) {
  #   if (paste0("variant_id_", build) %in% names(summary_stats)
  #       & all(summary_stats[, .(chr, pos, ref, alt)] ==
  #             summary_stats[, mget(paste0(c("chr_", "pos_", "ref_", "alt_"), build))],
  #             na.rm = TRUE))
  #     return(build)
  # }

  stop("Unable to infer genome assembly from variant ids.")
}


has_build_version <- function(
  summary_stats,
  build,
  columns = c("variant_id", "chr", "pos", "start", "end") |>
    intersect(names(summary_stats))
) {
  all(
    get_build_colnames(build, columns = columns) %in% names(summary_stats)
  )
}
get_build_colnames <- function(
    build,
    columns = c("variant_id", "chr", "pos")
) {
  sprintf("%s_%s",
          columns,
          build)
}


# In some MVP summary statistics, alleles are listed as 1 and 2 with a separate
# column indicating which is the effect allele (already renamed alt in
# loading). The ref allele is set as the other one here.
infer_ref_alt_alleles <- function(
  summary_stats
) {
  if (all(c("ref", "alt") %in% colnames(summary_stats))) {
    return(summary_stats[, .(ref, alt)])
  }

  if ("alt" %in% colnames(summary_stats)) {
    return(sort_a1_a2_into_ref_alt(summary_stats))
  }

  allele1 <- intersect(
    summary_stats_extra_column_names$allele1,
    colnames(summary_stats)
  )
  allele2 <- intersect(
    summary_stats_extra_column_names$allele2,
    colnames(summary_stats)
  )

  if (1 <= length(allele1) & 1 <= length(allele2)) {
    allele1 <- allele1[1]
    allele2 <- allele2[1]
  } else if (2 <= length(allele1)) {
    allele2 <- allele1[2]
    allele1 <- allele1[1]
  } else if (2 <= length(allele2)) {
    allele1 <- allele2[1]
    allele2 <- allele2[2]
  } else {
    stop("Failed to identify reference and alternate alleles.")
  }
  warning("No effect allele identified, assuming reference allele is '",
          allele1, "' and alternate (effect) allele is '", allele2, "'.")

  summary_stats[
    ,
    .(
      ref = get(allele1),
      alt = get(allele2),
      assumed_ea_col = allele2
    )
  ]
}
sort_a1_a2_into_ref_alt <- function(
  summary_stats
) {
  allele1 <- intersect(
    summary_stats_extra_column_names$allele1,
    colnames(summary_stats)
  )
  allele2 <- intersect(
    summary_stats_extra_column_names$allele2,
    colnames(summary_stats)
  )
  if (1 <= length(allele1) & 1 <= length(allele2)) {
    allele1 <- allele1[1]
    allele2 <- allele2[1]
  } else if (2 <= length(allele1)) {
    allele2 <- allele1[2]
    allele1 <- allele1[1]
  } else if (2 <= length(allele2)) {
    allele1 <- allele2[1]
    allele2 <- allele2[2]
  } else {
    stop("Failed to identify reference and alternate alleles.")
  }
  summary_stats[
    ,
    .(ref = data.table::fifelse(
      get(allele1) == alt,
      get(allele2),
      get(allele1)),
      alt = alt)
  ]
}


infer_position <- function(summary_stats) {
  position_encoding_columns <- intersect(
    summary_stats_extra_column_names$position_encoding$column_name,
    names(summary_stats)
  )
  if (length(position_encoding_columns) == 0) {
    stop("Unable to infer position from existing columns")
  }

  summary_stats[
    ,
    stringr::str_match(
      get(position_encoding_columns[1]),
      summary_stats_extra_column_names$position_encoding[
        position_encoding_columns[1],
        regex,
        on = "column_name"
      ]
    )[, 2] |>
      as.integer()
  ]
}
