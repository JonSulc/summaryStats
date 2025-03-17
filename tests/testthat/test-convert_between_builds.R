query1 <- data.table::data.table(
  chr        = "chr3",
  start      = 10042,
  end        = 10515
) |>
  new_genomic_ranges() |>
  suppressWarnings()

query2 <- data.table::data.table(
  chr        = "chr3",
  start      = 10519,
  end        = 12378
) |>
  new_genomic_ranges() |>
  suppressWarnings()

test_that("Converting to the same build doesn't change the object", {
  gr <- new_genomic_ranges(query1) |>
    suppressWarnings()
  expect_equal(
    calculate_converted_positions(gr, "b36", "b36"),
    gr
  )
})

test_that("Converting does not change the original object", {
  query1_bkp <- data.table::copy(query1)
  calculate_converted_positions(
    query1,
    target_build = "b37",
    current_build = "b38",
    chain = chain_b38_b37
  )
  expect_equal(query1, query1_bkp)
})

test_that("Conversion matches rtracklayer::liftOver", {
  for (query in list(query1, query2)) {
    expect_equal(
      calculate_converted_positions(
        query,
        target_build = "b37",
        current_build = "b38",
        chain = chain_b38_b37
      )[
        ,
        .(chr, start, end)
      ],
      data.table::as.data.table(
        rtracklayer::liftOver(
          query[
            ,
            GenomicRanges::GRanges(
              chr,
              IRanges::IRanges(start, end)
            )
          ],
          chain = chain_b38_b37
        )
      )[
        order(start),
        .(chr   = as.character(seqnames),
          start = as.numeric(start),
          end   = as.numeric(end))
      ][
        order(chr)
      ]
    )
  }
})

test_that("Converting multiple ranges preserves order", {
  expect_equal(
    calculate_converted_positions(
      rbind(query1, query2),
      target_build = "b37",
      current_build = "b38",
      chain = chain_b38_b37
    ),
    rbind(
      calculate_converted_positions(
        query1,
        target_build = "b37",
        current_build = "b38",
        chain = chain_b38_b37
      ),
      calculate_converted_positions(
        query2,
        target_build = "b37",
        current_build = "b38",
        chain = chain_b38_b37
      )
    )
  )
})

summary_stats <- data.table::data.table(
  chr = 3,
  pos = 10042,
  ref = "A",
  alt = "C",
  effect = 1,
  pval = .5
) |>
  new_summary_stats(build = "b38") |>
  suppressWarnings()

test_that("Converting to genomic ranges works", {
  expect_true(
    all(c("start", "end") %chin% names(convert_to_genomic_range(summary_stats)))
  )
})

test_that("Converting summary stats works", {
  expect_equal(
    calculate_converted_positions(
      summary_stats,
      target_build = "b37",
      current_build = "b38"
    )[
      ,
      c(chr[[1]], pos[[1]])
    ],
    c("chr1", 249240559)
  )
})

test_that("Converting with NAs works", {
  expect_no_error(
    data.table::data.table(
      chr = NA,
      pos = NA,
      ref = "A",
      alt = "C",
      effect = 1,
      pval = 1
    ) |>
      as_summary_stats(build = "b38") |>
      calculate_converted_positions(chain_b38_b37,
                                    target_build = "b37",
                                    current_build = "b38") |>
      suppressWarnings()
  )
  expect_no_error(
    data.table::data.table(
      chr = c(NA, "chr1"),
      pos = c(NA, 1),
      ref = "A",
      alt = "C",
      effect = 1,
      pval = 1
    ) |>
      as_summary_stats(build = "b38") |>
      calculate_converted_positions(target_build = "b37",
                                    current_build = "b38") |>
      suppressWarnings()
  )
})

test_that("Missing SNPs in build conversion are handled properly", {
  # This SNP does not have an equivalent in b37, should raise warning
  summary_stats <- data.table(
    variant_id = "chr5_49861645_A_G_b38",
    effect = 1,
    pval = .5
  ) |>
    as_summary_stats() |>
    suppressWarnings()
  expect_warning(capture.output(calculate_converted_positions(summary_stats,
                                                              target_build = "b37",
                                                              current_build = "b38")),
                 "no mapped pos")
})
