#' @import data.table

#' @export
all_mvp_traits <- data.table::fread(
  "~/rcp_storage/common/Users/sulc/data/mvp/mvp_traits.tsv"
)
data_path_mvp_prefix <- "~/Databases/MVP/release"

all_charge_traits_path <- "~/rcp_storage/common/Users/sulc/data/charge/CHARGE_lookup_table.csv"
#' @export
all_charge_traits <- data.table::fread(
  all_charge_traits_path
)
data_path_charge_prefix <- "~/Databases/CHARGE/authorized_data/deflated_organized/submission"
