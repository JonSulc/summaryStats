#' @import data.table

# If `x` argument of sample is length 1 integer, will sample from 1:x instead
rsample <- function(
  x,
  size,
  replace = TRUE,
  ...
) {
  if (length(x) == 1) return(rep(x, size))
  sample(x, size, replace = replace, ...)
}

fake_dt <- function(
  nrows = 500,
  chr = 1,
  start = 1,
  end = 12345,
  ref = c("A", "C", "G", "T"),
  alt = c("A", "C", "G", "T"),
  pval = NULL,
  effect = NULL
) {
  if (is.null(pval)) pval <- runif(nrows)
  if (is.null(effect)) effect <- rnorm(nrows)
  if (is.numeric(chr)) chr <- paste0("chr", chr)
  dt <- data.table::data.table(
    chr = rsample(chr, nrows),
    pos = rsample(start:end, nrows, replace = FALSE),
    ref = sample(ref, nrows, replace = TRUE),
    pval = pval,
    effect = effect
  )
  dt[
    ,
    alt := sample(setdiff(alt, ref), .N, replace = TRUE),
    by = ref
  ][]
  dt[order(chr, pos), .(chr, pos, ref, alt, effect, pval)]
}

fake_stats <- function(
  ...,
  build = "b38"
) {
  summary_stats <- fake_dt(...)[
    ,
    c(.(variant_id = get_variant_id_from_chr_pos_ref_alt(.SD, build = build)),
      .SD,
      .(effect_se = abs(effect / qnorm(pval/2))))
  ]
  set_as_summary_stats(summary_stats)
  summary_stats
}
