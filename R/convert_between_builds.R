#' @import data.table
#' @importFrom rtracklayer liftOver import.chain

# LiftOver chain files for the conversion between builds
chain_b37_b38 <- "~/rcp_storage/common/Users/sulc/data/ucsc/hg19ToHg38.over.chain" |>
  rtracklayer::import.chain()
chain_b38_b37 <- "~/rcp_storage/common/Users/sulc/data/ucsc/hg38ToHg19.over.chain" |>
  rtracklayer::import.chain()
chain_b36_b38 <- "~/rcp_storage/common/Users/sulc/data/ucsc/hg18ToHg38.over.chain" |>
  rtracklayer::import.chain()
# No file found for b38 -> b36, have to go through b37...
chain_b37_b36 <- "~/rcp_storage/common/Users/sulc/data/ucsc/hg19ToHg18.over.chain" |>
  rtracklayer::import.chain()

#' @export
convert_to_build <- function(
    summary_stats,
    build = "b38",
    current_build = get_build(summary_stats),
    chain = NULL,
    force = FALSE
) {
  assign_current_as_build(summary_stats, current_build = current_build)

  if (current_build == build) {
    return(invisible(summary_stats))
  }

  if (has_build_version(summary_stats, build) & !force) {
    summary_stats[
      ,
      c("variant_id", "chr", "pos") :=
        mget(get_build_colnames(build, c("variant_id", "chr", "pos")))
    ][]
    return(invisible(summary_stats))
  }

  # There's no chain file for hg38 -> 36
  if (is.null(chain)) {
    if (current_build == "b38" & build == "b36") {
      convert_to_build(summary_stats, current_build = "b38", build = "b37")
      convert_to_build(summary_stats, current_build = "b37", build = "b36")
      return(
        invisible(summary_stats)
      )
    }
    chain <- get_chain(from = current_build, to = build)
  }

  converted_positions <- calculate_converted_positions(summary_stats,
                                                       chain,
                                                       build)
  summary_stats[
    ,
    names(converted_positions) := converted_positions
  ][]
  assign_current_as_build(summary_stats, build)

  invisible(summary_stats)
}


get_chain <- function(from, to) {
  if (!paste(from, to) %in% c("b38 b37",
                              "b37 b38", "b37 b36",
                              "b36 b38", "b36 b37")) {
    stop("Attempting non-supported conversion from ", from, " to ", to)
  }
  get(sprintf("chain_%s_%s", from, to))
}


calculate_converted_positions <- function(
    summary_stats,
    chain,
    build
) {
  if (nrow(summary_stats) == 0) {
    return(
      data.table::data.table(
        variant_id = character(0),
        chr = character(0),
        pos = numeric(0),
        ref = character(0),
        alt = character(0)
      )
    )
  }
  if (all(summary_stats[, is.na(chr)])) {
    return(summary_stats)
  }

  positions <- summary_stats[, .(chr, pos, ref, alt)]
  positions[
    is.na(chr) | is.na(pos),
    c("chr", "pos") := NA
  ][
    ,
    # Positions are changed when applying by chr, need to be able to resort them
    index := seq_len(nrow(positions))
  ][]
  positions <- positions[
    ,
    c(.SD,
      .(start = pos,
        end   = pos))
  ][
    !is.na(chr),
    data.table::foverlaps(
      .SD,
      get_chain_dt(chain, chr),
      mult = "first"
    )[
      ,
      .(
        pos = data.table::fifelse(
          rev,
          end - offset - (i.start - start),
          i.start - offset
        ),
        new_chr = chr,
        ref = ref,
        alt = alt,
        index = index
      )
    ],
    by = chr
  ][
    ,
    .(chr = new_chr,
      pos = pos,
      ref = ref,
      alt = alt,
      index = index)
  ][
    ,
    c(.(variant_id = get_variant_id_from_chr_pos_ref_alt(.SD, build = build)),
      .SD)
  ][
    data.table::data.table(index = seq_len(nrow(positions))) |>
      cbind(positions[, .(ref, alt)]),
    on = c("index", "ref", "alt")
  ][
    ,
    -"index"
  ]

  if (any(is.na(positions$variant_id))) {
    warning(positions[is.na(variant_id), .N],
            " SNP(s) have no mapped position in build ",
            build)
  }

  positions
}

get_chain_dt <- function(chain, chr) {
  data.table::as.data.table(chain[[chr]]@ranges) |>
    cbind(data.table::data.table(
      offset  = chain[[chr]]@offset,
      chr = rep(chain[[chr]]@space,
                chain[[chr]]@length),
      rev = rep(chain[[chr]]@reversed,
                chain[[chr]]@length)
    )) |>
    setkey(start, end)
}
