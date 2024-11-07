#' @import data.table

#' @title Retrieve summary statistics for a genome region from dbGaP on LISP
#'   server
#' @description Easy method to retrieve summary statistics from MVP and CHARGE
#'   on LISP server.
#' @param chromosome Specifies the chromosome for the region of interest,
#'   specified either as the number alone (i.e., 1 or "1") or preceded by
#'   "chr" (i.e., "chr1"). Can also be a data.frame or data.table, in which
#'   case it must have (at least) three columns: "chr", "start", "end", each
#'   row of which specifies a different region of interest.
#' @param start Starting position of the region of interest.
#' @param end End position of the region of interest.
#' @param traits Vector of traits of interest, either keywords (e.g.,
#'   "ejection fraction") or trait IDs (e.g., "HFpEF_EUR")
#' @param source Cohort or consortium from which the data is to be loaded. Must
#'   be one of "mvp" or "charge".
#' @param build Specifies the genome build for the genomic regions of interest.
#'   Can be one of "b38" (GRCh38/hg38), "b37" (GRCh37/hg19) or "b36"
#'   (NCBI36/hg18).
#'
#' @returns A data.table containing the summary statistics of any of the
#'   variants in the region(s) of interest loaded from the files corresponding
#'   to the requested traits.
#' @examples
#' get_dbgap_associations(chromosome = 3,
#'                        start = 123456,
#'                        end = 133456,
#'                        traits = "ejection fraction",
#'                        source = "mvp")
#'
#' @export
get_dbgap_associations <- function(
  chromosome,
  start,
  end,
  traits,
  source = c("mvp", "charge"),
  build = "b38"
) {
  source <- match.arg(source, tolower(source))
  new_genomic_ranges(
    chromosome = chromosome,
    start      = start,
    end        = end,
    build      = build
  ) |>
    get_dbgap_associations_from_genomic_ranges(
      traits        = traits,
      source        = source,
      current_build = build
    )
}


get_dbgap_associations_from_genomic_ranges <- function(
    genomic_ranges,
    traits,
    source,
    current_build = "b38",
    ...
) {
  if (nrow(genomic_ranges) == 0) {
    stop("No data to be extracted from ", toupper(source), " data.")
  }
  traits <- find_dbgap_matching_traits(traits, source = source)
  traits[
    ,
    filename := file.path(get(sprintf("data_path_%s_prefix", source)), filename)
  ]

  if (any(traits$build == "")) {
    traits[
      (build == ""),
      build := infer_build_from_file(
        filename,
        source = source
      ),
      by = filename
    ]
    update_charge_trait_builds(traits)
  }

  genomic_ranges <- setNames(nm = unique(traits$build[traits$build != ""])) |>
    lapply(
      \(target_build) {
        calculate_converted_positions(genomic_ranges,
                                      target_build = target_build,
                                      current_build = current_build)
      }
    )

  traits[
    build != "",
    traits_data := .(.(
      load_dbgap_data(
        filename,
        genomic_ranges = genomic_ranges[[build]],
        source = source,
        build = build
      ) |>
        calculate_converted_positions(target_build = current_build,
                                      current_build = build)
    )),
    by = .(build, filename)
  ][
    build == "",
    traits_data := .(.(
      data.table::data.table()
    ))
  ][]

  traits[
    ,
    # Build needs to be dropped, otherwise is taken as argument in later
    # function calls from within the DT
    -"build"
  ][
    ,
    .(
      trait_data = .(c(
        .SD,
        traits_data[[1]]
      ))
    ),
    .SDcols = -"traits_data",
    by = "filename"
  ][
    ,
    rbindlist(trait_data, fill = TRUE)
  ][
    # Where there are no SNPs in common, the data is filled with NA
    !is.na(variant_id)
  ]
}


load_dbgap_data <- function(
    filename,
    genomic_ranges,
    source = c("mvp", "charge"),
    build = "",
    ...
) {
  source <- match.arg(source)

  if (build == "") {
    build <- infer_build_from_file(filename, source = source)
  }

  cat("Loading", filename, "\n")

  summary_stats <- data.table::fread(
    cmd = cl_awk_from_genomic_ranges(
      filename,
      source,
      genomic_ranges,
      ...
    )
  ) |>
    as_summary_stats(build = build) |>
    # Some MVP files have duplicate SNP entries (e.g., mvp-padxt2d)
    unique()

  summary_stats
}


find_dbgap_matching_traits <- function(
    keyword,
    source
) {
  traits <- data.table::data.table(keyword = keyword)

  found_traits <- data.table::rbindlist(
    list(
      get(sprintf("all_%s_traits", source))[
        traits,
        c(.(keyword = filename),
          .SD),
        on = c(filename = "keyword"), nomatch = 0
      ],
      get(sprintf("all_%s_traits", source))[
        traits,
        c(.(keyword = id),
          .SD),
        on = c(id = "keyword"), nomatch = 0
      ]
    )
  )

  if (all(keyword %in% found_traits$keyword)){
    return(found_traits)
  }

  found_traits <- found_traits[
    traits[
      !found_traits$keyword,
      lapply(keyword, \(key) get(sprintf("all_%s_traits", source))[
        grep(key, trait, ignore.case = TRUE),
        {
          if (.N == 0) {
            warning(sprintf("No %s traits found for keyword %s",
                            toupper(source),
                            key))
          } else {
            cat(sprintf("Found %s trait(s) for keyword '%s':",
                        toupper(source),
                        key),
                paste(trait, collapse = ", "),
                "\n")
          }
          .SD |> c(.(keyword = key))
        }
      ]) |> data.table::rbindlist(),
      on = "keyword"
    ],
    on = colnames(found_traits)
  ]

  if (nrow(found_traits) == 0 | all(is.na(found_traits$filename))) {
    warning("No ", toupper(source), " outcomes found matching any keyword.")
    return()
  }

  if (any(is.na(found_traits$filename)))
    warning("Some traits were not found: ",
            paste(found_traits[is.na(filename), keyword], collapse = ", "))

  if (any(found_traits$build == "")) {
    # Only currently an issue with CHARGE data
    stopifnot(source == "charge")
    cat("Some builds are unknown, inferring...\n")
    found_traits[
      build == "" & !is.na(filename),
      c("build", "build_match") := infer_build_from_file(
        file.path(path_prefix, filename),
        source = source
      ),
      by = filename
    ][]
    cat("Updating lookup table\n")
    get(sprintf("update_%s_trait_builds", source))(found_traits)
  }

  found_traits[
    !is.na(filename),
    .SD,
    .SDcols = setdiff(colnames(found_traits), c("keyword", "build_match"))
  ]
}
