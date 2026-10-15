##### 07. Assemble Supplementary Data 1 as a multi-tab workbook #####

##### Setup #####

rm(list = ls())
gc()

setwd("C:/Users/simon/Dropbox/2026-01 Toxoplasma XL-MS/R-Project-Toxoplasma-XL-MS/")
source("./script/PACKAGES.R")
lapply(list.files("./script/functions/", pattern = "\\.R$", full.names = TRUE), source)

##### Read the pipeline outputs #####

# Each tab is read from the stable `_latest` output of the script that produced it, so the
# workbook always reflects the most recent complete run of the pipeline.
csms        <- read_excel("./outputs/Combined_Data_CSMs_latest.xlsx")
crosslinks  <- read_excel("./outputs/Combined_Data_Crosslinks_latest.xlsx")
indexed     <- read_excel("./outputs/Combined_Data_Crosslinks_Disambiguated_Indexed_latest.xlsx")
residues    <- read_excel("./outputs/Combined_Data_Residue_Pairs_latest.xlsx")
ppis        <- read_excel("./outputs/Combined_Data_PPIs_latest.xlsx")

##### Score cutoffs and target / decoy counts #####

# Counts and FDR are computed from the tables themselves, so tab F cannot drift from the
# data it describes. FDR = decoys / targets, matching report_fdr().
count_level <- function(data, level, cutoffs) {
  tibble(
    `Aggregation Level` = level,
    `Score Cutoffs`     = cutoffs,
    `Total Count`       = nrow(data),
    `Target Count`      = sum(!data$`Is Decoy`),
    `Decoy Count`       = sum(data$`Is Decoy`)) |>
    mutate(`Estimated FDR` = `Decoy Count` / `Target Count`)
}

# The cutoff descriptions are entered manually from the messages printed by 01, 02, 05 and
# 06. The PPI gate is the conjunction of the crosslink-level cutoffs, which the input
# crosslinks have already passed, and the OR gate applied at PPI level.
score_cutoff_summary <- bind_rows(
  count_level(csms, "Crosslink-spectrum match",
              "Score \u226520 AND DeltaScore \u226510"),
  count_level(crosslinks, "Crosslink",
              "Score \u226539.6 AND DeltaScore \u226510.4"),
  count_level(residues, "Residue pair",
              "Score \u226540.4 AND DeltaScore \u226510.8"),
  count_level(ppis, "Protein-protein interaction",
              paste("(Score \u226539.6 AND DeltaScore \u226510.4) AND",
                    "(Score \u226564.1 OR DeltaScore \u226524.7)")))

print(score_cutoff_summary)

##### Assemble the workbook #####

# Decoys are retained in tabs A to D, since they are what the reported FDR is estimated
# from. Tab E is the final interaction network, so decoys are dropped there.
supplementary_data_1 <- list(
  Supplementary_Data_1A = csms,
  Supplementary_Data_1B = crosslinks,
  Supplementary_Data_1C = indexed,
  Supplementary_Data_1D = residues |> filter(!`Is Decoy`) |> select(-`Is Decoy`),
  Supplementary_Data_1E = ppis  |> filter(!`Is Decoy`) |> select(-`Is Decoy`),
  Supplementary_Data_1F = score_cutoff_summary)

##### Save #####

openxlsx::write.xlsx(
  supplementary_data_1,
  file = paste0("./outputs/Supplementary_Data_1_", Sys.Date(), ".xlsx"))
