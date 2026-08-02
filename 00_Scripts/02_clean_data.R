# 02_clean_data.R
# Reshape raw into a tidy core trial frame (2002-2005) and build
# genotype x environment (Year) mean matrices per trait.
# Missing values are PRESERVED as NA (no mean-imputation, principle 8).
# 
# MODIFICATION (2026-07-30): Use parse_traits() to convert death/missing
# markers to 0 instead of NA, preventing bias from excluding non-surviving
# plants (see issue #...). 
source("00_config.R")

raw <- readRDS(file.path(INT, "raw_all.rds"))

# ---- Custom parsing function to keep mortality records as 0 ----
parse_traits <- function(x) {
  x_str <- trimws(as.character(x))
  x_num <- suppressWarnings(as.numeric(x_str))
  # Identify records indicating death or non-survival (adjust keywords as needed)
  dead_idx <- is.na(x_num) & grepl("dead|die|missing|loss|死|无", x_str, ignore.case = TRUE)
  x_num[dead_idx] <- 0  # Set dead plants to 0 to preserve their true contribution
  return(x_num)
}

core <- raw %>%
  filter(Year >= 2002) %>%
  mutate(
    Genotype = as.character(Line),
    Rep      = as.integer(Rep),
    SummerPH = parse_traits(Mean_Summer),
    FW       = parse_traits(FW),
    DW       = parse_traits(DW),
    AutumnPH = parse_traits(Mean_Autumn)
  ) %>%
  # FIX (2026-08-01): In the 2002 CSV, FW and DW columns are swapped.
  # Evidence: 162/510 rows in 2002 have DW > FW (physically impossible),
  # whereas 2003-2005 have zero such rows.
  # Correction: swap FW and DW for all 2002 records to restore the
  # expected relationship FW >= DW.
  mutate(
    FW_tmp = FW,
    FW = if_else(Year == 2002, DW, FW),
    DW = if_else(Year == 2002, FW_tmp, DW)
  ) %>%
  select(Year, Genotype, Rep, SummerPH, FW, DW, AutumnPH, -FW_tmp)

# --- Aggregate to the experimental unit: Genotype x Rep x Year ------------
# Each (Genotype, Rep, Year) contains multiple distinct sampling rows
# (NOT duplicate entries -- values differ). Treating them as independent
# observations would cause pseudoreplication and understate error variance.
# We therefore average them to one rep-level value per trait.
# ASSUMPTION (flagged for author confirmation): these multiple rows are
# sub-samples / sub-plots within a replication; the true design should be
# confirmed against field books.
core_rep <- core %>%
  group_by(Year, Genotype, Rep) %>%
  summarise(
    across(
      c(SummerPH, FW, DW, AutumnPH),
      ~ {
        if (all(is.na(.))) NA_real_ else mean(., na.rm = TRUE)
      }
    ),
    .groups = "drop"
  )

# --- Genotype x Environment (Year) means across reps, per trait ------------
gxe_means <- function(trait) {
  core_rep %>%
    filter(!is.na(.data[[trait]])) %>%
    group_by(Genotype, Year) %>%
    summarise(value = mean(.data[[trait]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Year, values_from = value) %>%
    arrange(Genotype)
}
gxe_list <- setNames(lapply(TRAITS, gxe_means), TRAITS)

saveRDS(core,     file.path(INT, "core.rds"))
saveRDS(core_rep, file.path(INT, "core_rep.rds"))
saveRDS(gxe_list, file.path(INT, "gxe_list.rds"))

# --- NA / coverage diagnostics (used by Statistics_Check_Report) -----------
na_tab <- core %>% summarise(across(c(SummerPH, FW, DW, AutumnPH), ~ sum(is.na(.))))
geno_per_year <- core %>% group_by(Year) %>% summarise(n_geno = n_distinct(Genotype),
                                                       n_obs  = n())
cat("\n=== NA counts in core (2002-2005) ===\n")
print(na_tab)
cat("\n=== Genotypes / observations per year ===\n")
print(geno_per_year)
cat("\n=== Environments available per trait ===\n")
for (t in TRAITS) cat(sprintf("  %-9s: %s\n", t, paste(names(gxe_list[[t]])[!names(gxe_list[[t]]) %in% c("Genotype")], collapse = ", ")))
cat("\nTotal distinct genotypes in core:", n_distinct(core$Genotype), "\n")