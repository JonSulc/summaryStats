#' @import data.table
#' @importFrom stringr str_match

get_variant_id_from_chr_pos_ref_alt <- function(summary_stats,
                                                build = "b38") {
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
assign_variant_id_from_chr_pos_ref_alt <- function(summary_stats, ...) {
  summary_stats[
    ,
    variant_id := get_variant_id_from_chr_pos_ref_alt(summary_stats, ...)
  ][]
  invisible(summary_stats)
}

parse_variant_id <- function(variant_id) {
  stringr::str_match(
    variant_id,
    "^(chr[0-9xX]+)_([0-9]+)_([^_]+)_([^_]+)_.*$"
  ) |>
    data.table::as.data.table() |>
    setNames(c("variant_id", "chr", "pos", "ref", "alt")) |>
    assign_correct_pos()
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
