# 01_import_data.R
# Read the five yearly raw CSV exports (2001-2005) and bind into one table.
# NOTE: 2002-2005 files contain duplicate column names (Date, PH1..PH10 appear
# twice: summer & autumn). We use default check.names=TRUE so duplicates become
# "Date.1", "PH1.1", ... and then select a UNIFORM set of columns. Missing
# columns in a given year are filled with NA (safe union).
# Source data is NEVER modified (principle 3) -- we only read and combine.
source("00_config.R")

raw_files <- list.files(RAW, pattern = "data\\.csv$", full.names = TRUE)
message("Found raw files:")
print(basename(raw_files))

read_one <- function(f) {
  # default check.names = TRUE standardises duplicate headers
  d <- read.csv(f, fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  yr <- as.integer(str_extract(basename(f), "\\d{4}"))
  d$Year <- yr
  # Uniform columns kept across all years (absent -> NA after bind)
  keep <- c("Year", "Line", "Rep", "Mean_Summer", "FW", "DW", "Mean_Autumn",
            "PH_Autumn")
  miss <- setdiff(keep, names(d))
  for (m in miss) d[[m]] <- NA
  # Coerce all kept columns to character BEFORE binding so that any
  # year-to-year type mismatch (e.g. Rep stored as int vs char) cannot
  # trigger a vctrs incompatibility error. Numeric coercion happens in 02.
  out <- d[keep]
  out[] <- lapply(out, as.character)
  out
}

raw_list <- lapply(raw_files, read_one)
raw_all  <- bind_rows(raw_list)

saveRDS(raw_all, file.path(INT, "raw_all.rds"))

cat("\nTotal rows imported:", nrow(raw_all), "\n")
cat("Rows per year:\n")
print(table(raw_all$Year))
cat("\nUniform columns retained:\n")
print(names(raw_all))
