##### Aggregate crosslinks to residue pairs and re-filter to 1% FDR #####

##### Setup #####

rm(list=ls())
gc()
setwd("C:/Users/simon/Dropbox/R-Project-Crosslinking/")
source("./script/PACKAGES.R")
lapply(list.files("./script/functions/", pattern = "\\.R$", full.names = TRUE), source)

crosslinks <- read_excel("./outputs/Combined_Data_Crosslinks_Disambiguated_Indexed_2026-02-24.xlsx")

##### Aggregate crosslinks to unique residue pairs #####

# Note: A/B columns are already in alphanumeric protein order from the previous
# indexing script, so no pmin/pmax is needed — grouping directly on Leading Protein
# A/B and Crosslink Position A/B gives the canonical residue pair identity.
residues <- crosslinks %>%
  # Consider only unambiguous residue pairs
  filter(
    `# Accessions A` == 1 & `# Accessions B` == 1,
    `# Crosslink Positions A` == 1 & `# Crosslink Positions B` == 1) %>%
  # Calculate total # CSMs and unique crosslinks supporting each residue pair
  group_by(`Leading Protein A`, `Crosslink Position A`, `Leading Protein B`, `Crosslink Position B`) %>%
  mutate(
    `# CSMs`             = sum(`# CSMs`),
    `# Unique Crosslinks` = n()) %>%
  ungroup() %>%
  # Keep the highest scoring crosslink for each residue pair
  arrange(desc(`Identification Score`)) %>%
  distinct(`Leading Protein A`, `Crosslink Position A`, `Leading Protein B`, `Crosslink Position B`, .keep_all = TRUE)

##### Identify optimal score cutoffs to obtain 1% residue-pair-level FDR #####

residues_score_cutoffs <- Identify_Score_Cutoffs(
  residues, gate = "AND", fdr.level = 0.01,
  min.score = 20, max.score = 50,
  min.delta.score = 10, max.delta.score = 50)

# Row 1 maximises target residue pairs at <= 1% FDR
if (nrow(residues_score_cutoffs) == 0) stop("No valid score cutoffs found — check input data or relax threshold ranges.")
MIN_SCORE <- residues_score_cutoffs[1, ]$score
MIN_DELTA <- residues_score_cutoffs[1, ]$delta.score

# Filter to 1% FDR
residues <- residues %>%
  filter(`Identification Score` >= MIN_SCORE, `Identification DeltaScore` >= MIN_DELTA)

# Check FDR
residues %>%
  group_by(`Is Decoy`) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(proportion = n / sum(n)) %>%
  print()

##### Calculate local FDR #####

residues <- calculate_local_fdr(residues)

ggplot(residues, aes(x = `Identification Score`, y = `Local FDR`)) +
  geom_point()

##### Save data #####

# Save target residue pairs only (decoys excluded from final output)
openxlsx::write.xlsx(
  filter(residues, !`Is Decoy`),
  file = paste0("./outputs/Combined_Data_Residue_Pairs_", Sys.Date(), ".xlsx"))
