library(quickcheck)

test_that("Variant ID construction/parsing works", {
  quickcheck::for_all(
    summary_stats = quickcheck::data.table_(
      chr = integer_bounded(1, 22),
      pos = quickcheck::integer_positive(),
      ref = character_letters(),
      alt = character_letters(),
      effect = quickcheck::double_(),
      pval = double_bounded(1e-10, 1)
    ),
    property = function(summary_stats) {
      summary_stats <- new_summary_stats(summary_stats, build = "b38")
      summary_stats[
        ,
        variant_id := get_variant_id_from_chr_pos_ref_alt(summary_stats)
      ]
      expect_equal(parse_variant_id(summary_stats[, variant_id]),
                   summary_stats[, .(variant_id, chr, pos, ref, alt)])
    }
  )
})
