#' @import data.table
#' @importFrom jsonlite fromJSON
#' @importFrom stringr str_match

# RefSNP data is only loaded if necessary and saved in an internal environment
internal_data <- new.env(parent = emptyenv())
# RefSNP data for the inference of genome build
data_path_refsnp_db22 <- "~/rcp_storage/common/Users/sulc/data/dbsnp/chr22.csv"

read_json_file <- function(filename, nlines = -1L, ...) {
  readLines(filename, n = nlines, ...) |>
    lapply(parse_single_line_to_snp) |>
    data.table::rbindlist()
}


parse_single_line_to_snp <- function(single_line, ...) {
  jsonlite::fromJSON(single_line) |>
    parse_single_snp_json(...)
}


parse_single_snp_json <- function(single_snp_json, output_db_file = NULL) {
  is_genomic <- are_assemblies_from_seq_id_traits_genomic(
    single_snp_json$primary_snapshot_data$placements_with_allele$placement_annot$seq_id_traits_by_assembly
  )
  db <- data.table::data.table(
    build = single_snp_json$primary_snapshot_data$placements_with_allele$placement_annot$seq_id_traits_by_assembly[is_genomic] |>
      sapply(\(seq_id_traits) {
        stringr::str_match(seq_id_traits$assembly_name, "^[^.]+") |>
          stringr::str_replace("^GRCh", "b")
      }),
    chr = single_snp_json$primary_snapshot_data$placements_with_allele$alleles[is_genomic] |>
      lapply(\(alleles) alleles$allele$spdi$seq_id |>
               get_chr_from_seq_id()),
    pos = single_snp_json$primary_snapshot_data$placements_with_allele$alleles[is_genomic] |>
      lapply(\(alleles) alleles$allele$spdi$position),
    ref = single_snp_json$primary_snapshot_data$placements_with_allele$alleles[is_genomic] |>
      lapply(\(alleles) alleles$allele$spdi$deleted_sequence),
    alt = single_snp_json$primary_snapshot_data$placements_with_allele$alleles[is_genomic] |>
      lapply(\(alleles) alleles$allele$spdi$inserted_sequence)
  )
  if (nrow(db) != 0) {
    db <- db[
      ,
      c(build = build, chr = chr, pos = pos, ref = ref, alt = alt),
      by = seq_len(nrow(db))
    ][
      ,
      .(build, chr, pos, ref, alt)
    ][
      ref != alt
    ] |>
      unique()
  }

  if (is.null(output_db_file)) return(db)
  data.table::fwrite(db, output_db_file, append = TRUE)
  NULL
}


are_assemblies_from_seq_id_traits_genomic <- function(seq_id_traits) {
  sapply(seq_id_traits, \(traits){
    length(traits) != 0 &
      all(nrow(traits) != 0) &
      all(traits$is_chromosome)
  })
}


get_chr_from_seq_id <- function(seq_id) {
  chr <- stringr::str_match(seq_id, "^NC_([0-9]+)[.]")[, 2] |>
    as.numeric()
  paste0("chr", chr)
}


load_genomic_ranges_from_file <- function(filename,
                                          chr = 22,
                                          nsnps = 10000,
                                          source = c("mvp", "charge")) {
  variants_dt <- filtered_fread_on_chr(filename,
                                       chr = chr,
                                       source = source,
                                       nsnps = nsnps) |>
    as_variant_table(build = "b38")
  variants_dt[
    ,
    .(chr, pos, ref, alt)
  ]
}


infer_build_from_file <- function(filename,
                                  nsnps = 10000,
                                  source = c("mvp", "charge"),
                                  refsnp_db = refsnp_db22()) {
  cat("Inferring build from", filename, "\n")
  genomic_ranges <- load_genomic_ranges_from_file(filename,
                                                  nsnps = nsnps,
                                                  source = source)

  if (nrow(genomic_ranges) == 0) {
    warning("Unable to infer build from file.")
    return()
  }

  build <- infer_build(genomic_ranges,
                       refsnp_db,
                       return_build_match = TRUE)
  build <- build[which.max(build_match)]
  build[
    ,
    sprintf("\t-> Build %s, %.2f%% match\n", build, 100*build_match) |>
      cat()
  ]
  build$build
}


infer_build <- function(genomic_ranges,
                        refsnp_db,
                        return_build_match = FALSE) {
  genomic_ranges_builds <- refsnp_db[
    genomic_ranges,
    on = c("chr", "pos", "ref", "alt")
  ][
    is.na(build),
    build := "b36"
  ][
    ,
    .(n = .N,
      build_match = .N/nrow(genomic_ranges)),
    by = build
  ]

  if (!return_build_match) {
    return(
      genomic_ranges_builds[
        which.max(n),
        build
      ]
    )
  }

  genomic_ranges_builds[
    ,
    .(build, build_match)
  ]
}


extract_db <- function(filename,
                       batch_size = 100,
                       output_db_file = stringr::str_replace(filename,
                                                             "json[.]bz2$",
                                                             "csv"),
                       skip_lines = -1) {
  if (file.exists(output_db_file) & skip_lines == -1) {
    stop("File already exists")
  }

  con <- file(filename, "r")
  on.exit(close(con))

  if (skip_lines != -1) {
    for (i in seq_len(skip_lines))
      readLines(con, n = 1)
  }

  while (TRUE) {
    lines <- readLines(con, n = 1)
    if ( length(lines) == 0 ) {
      break
    }
    purrr::walk(lines,
                parse_single_line_to_snp,
                output_db_file = output_db_file)
  }
}

# Used for testing and debugging
extract_line_number <- function(filename, skip, n = 3) {
  con <- file(filename, "r")
  on.exit(close(con))
  for (i in seq_len(skip)) readLines(con, n = 1)
  readLines(con, n = n)
}


refsnp_db22 <- function() {
  if (!"refsnp_db22" %in% names(internal_data)) {
    cat("Loading SNP reference for genome inference...\n")
    internal_data[["refsnp_db22"]] <- data.table::fread(data_path_refsnp_db22())
  }
  internal_data[["refsnp_db22"]]
}
