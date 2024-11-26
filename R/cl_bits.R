#' @import data.table
#' @importFrom stringr str_match str_detect

file_is_csv <- function(filename) {
  dt_output <- data.table::fread(filename, nrows = 1, verbose = TRUE) |>
    capture.output() |>
    stringr::str_match("sep='([^']+)'")
  dt_output[!is.na(dt_output[, 1]), 2][1] == ","
}

are_values_quoted <- function(
  filename
) {
  data.table::fread(filename, nrows = 1, quote = "") |>
    names() |>
    stringr::str_detect("\"") |>
    all()
}

check_quotes <- function(
  value,
  value_needs_to_be_quoted
) {
  if (!value_needs_to_be_quoted) return(value)
  paste0("\\\"", value, "\\\"")
}

cl_read_dbgap_file_template <- function(filename, source) {
  if (source == "mvp") {
    return(
      sprintf("zcat %s | awk%s '%%s'",
              filename,
              ifelse(file_is_csv(filename),
                     " -F','",
                     ""))
    )
  }
  sprintf("awk%s '%%s' %s",
          ifelse(file_is_csv(filename),
                 " -F','",
                 ""),
          filename)
}

cl_find_column_indices <- function(filename, source) {
  if (source == "mvp") {
    cat("The following 'gzip: stdout: Broken pipe' is normal, don't @ me:")
  }
  summary_stats <- data.table::fread(
    cmd = sprintf("%s %s%s",
                  ifelse(source == "mvp", "zcat ", "head -n 2"),
                  filename,
                  ifelse(source == "mvp", " | head -n 2", ""))
  )

  column_names <- get_colnames(summary_stats)

  setNames(nm = names(column_names)) |>
    lapply(\(column_name) {
    list(name  = column_names[[column_name]],
         index = which(colnames(summary_stats) == column_names[[column_name]]),
         prefix = {
           if (column_name == "chr" && grepl("^chr", summary_stats[[column_names[[column_name]]]])) {
             return("chr")
           }
           ""
         })
  })
}

cl_variable_is <- function(
  column_names,
  variable_name,
  variable_values,
  values_are_quoted
) {
  sprintf(
    "$%s == \"%s\" || $%s == \"%s\"",
    column_names[[variable_name]]$index,
    check_quotes(column_names[[variable_name]]$name, values_are_quoted),
    column_names[[variable_name]]$index,
    check_quotes(variable_values, values_are_quoted)
  )
}

cl_variable_is_within <- function(
  column_names,
  variable_name,
  minimum_values,
  maximum_values
) {
  c(
    sprintf(
      "$%s == \"%s\"",
      column_names[[variable_name]]$index,
      check_quotes(column_names[[variable_name]]$name, values_are_quoted)
    ),
    sprintf(
      "%s <= $%s && $%s <= %s",
      minimum_values, column_names[[variable_name]]$index,
      column_names[[variable_name]]$index, maximum_values
    )
  ) |>
    paste(collapse = " || ")
}

cl_strip_chr <- function(
  chr
) {
  if (grepl("^chr", chr)) {
    return(substring(chr, 4))
  }
  chr
}

cl_variant_is_within_genomic_ranges <- function(
    column_names,
    genomic_ranges,
    values_are_quoted
) {
  genomic_ranges[
    ,
    .(cl = cl_variant_is_within_single_genomic_range(
      column_names,
      cl_strip_chr(chr),
      start, end,
      values_are_quoted = values_are_quoted
    )),
    by = chr
  ]$cl
}

cl_variant_is_within_single_genomic_range <- function(
    column_names,
    chr,
    starts,
    ends,
    values_are_quoted
) {
  sprintf(
    "(%s) && (%s)",
    cl_variable_is(column_names,
                   "chr", chr,
                   values_are_quoted = values_are_quoted),
    cl_variable_is_within(column_names,
                          "pos", starts, ends)
  )
}

cl_awk_from_genomic_ranges <- function(
  filename,
  source,
  genomic_ranges
) {
  column_names <- cl_find_column_indices(filename, source)
  sprintf(
    cl_read_dbgap_file_template(filename, source),
    cl_variant_is_within_genomic_ranges(
      column_names,
      genomic_ranges,
      values_are_quoted = are_values_quoted(filename)
    )
  )
}
