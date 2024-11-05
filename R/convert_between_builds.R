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


get_chain <- function(from, to) {
  if (!paste(from, to) %in% c("b38 b37",
                              "b37 b38", "b37 b36",
                              "b36 b38", "b36 b37")) {
    stop("Attempting non-supported conversion from ", from, " to ", to)
  }
  get(sprintf("chain_%s_%s", from, to))
}


convert_to_genomic_range <- function(
  summary_stats
) {
  summary_stats[
    ,
    .(start = pos,
      end   = pos) |>
      c(.SD),
    .SDcols = -"pos"
  ]
}


convert_genomic_range_to_summary_stats <- function(
  genomic_range,
  current_build
) {
  if (any(!is.na(genomic_range$chr))) {
    if (genomic_range[, !all(start == end)]) {
      stop("'start' and 'end' values do not match.")
    }
  }

  genomic_range[
    ,
    .(pos = start) |>
      c(.SD),
    .SDcols = -c("start", "end")
  ] |>
    assign_variant_id_from_chr_pos_ref_alt(build = current_build) |>
    assign_current_as_build(current_build = current_build)
}


calculate_converted_positions <- function(
    query,
    target_build,
    current_build = get_build(query),
    chain = NULL,
    ...
) {
  if (nrow(query) == 0) {
    return(
      query
    )
  }
  if (all(query[, is.na(chr)])) {
    return(query)
  }
  if (is.numeric(query$chr)) {
    query <- copy(query)
    query[, chr := paste0("chr", chr)]
  }

  if (!all(c("start", "end") %chin% names(query))) {
    return(
      convert_to_genomic_range(query) |>
        calculate_converted_positions(target_build = target_build,
                                      current_build = current_build,
                                      chain = chain) |>
        convert_genomic_range_to_summary_stats(target_build)
    )
  }

  if (all(get_build_colnames(
    build = target_build,
    columns = c("chr", "start", "end")
  ) %chin% colnames(query))) {
    if (!all(get_build_colnames(
      build = current_build,
      columns = c("chr", "start", "end")
    ) %chin% colnames(query))) {
      query[
        ,
        get_build_colnames(
          build = current_build,
          columns = c("chr", "start", "end")
        ) := c("chr", "start", "end")
      ]
    }
    query[
      ,
      c("chr", "start", "end") := get_build_colnames(
        build = target_build,
        columns = c("chr", "start", "end")
      )
    ]
    return(query)
  }

  if (is.null(chain)) {
    # There's no chain for b38 -> b36...
    if (current_build == "b38" & target_build == "b36") {
      return(
        calculate_converted_positions(
          query,
          target_build = "b37",
          current_build = "b38",
          ...
        ) |>
          calculate_converted_positions(
            target_build = "b36",
            current_build = "b37"
          )
      )
    } else {
      chain <- get(
        sprintf("chain_%s_%s",
                current_build,
                target_build)
      )
    }
  }

  column_order <- data.table::copy(colnames(query))

  query <- data.table::copy(query)[
    is.na(chr) | is.na(start) | is.na(end),
    c("chr", "start", "end") := NA
  ][
    ,
    # Positions are changed when applying by chr, need to be able to resort them
    index := seq_len(nrow(query))
  ][]

  if (any(is.na(query$chr))) {
    warning("Dropping ", query[is.na(chr), .N], " SNPs with invalid chr")
  }

  query <- query[
    !is.na(chr),
    data.table::foverlaps(
      .SD,
      get_chain_dt(chain, chr),
      ...
    )[
      ,
      .(
        start = data.table::fifelse(
          rev,
          start - offset + end - data.table::fifelse(end < i.end, end, i.end),
          data.table::fifelse(start < i.start, i.start, start) - offset
        ),
        end = data.table::fifelse(
          rev,
          start - offset + end - data.table::fifelse(start < i.start, i.start, start),
          data.table::fifelse(end < i.end, end, i.end) - offset
        ),
        new_chr = chr
      ) |>
        c(.SD),
      .SDcols = -c("chr", "start", "end", "width", "offset", "rev", "i.start", "i.end")
    ],
    by = chr
  ][
    ,
    chr := new_chr
  ][
    ,
    new_chr := NULL
  ][]

  setorder(query, index, chr, start)

  query[
    ,
    index := NULL
  ][]

  setcolorder(query, column_order)

  if (any(is.na(query$chr))) {
    warning(query[is.na(chr), .N],
            " SNP(s) have no mapped position in build ",
            target_build)
  }

  query
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
