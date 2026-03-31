##### Aggregate inter-protein crosslinks to protein-protein interactions and re-filter to 1% FDR #####

##### Setup #####

rm(list=ls())
gc()
setwd("C:/Users/simon/Dropbox/R-Project-Crosslinking/")
source("./script/PACKAGES.R")
lapply(list.files("./script/functions/", pattern = "\\.R$", full.names = TRUE), source)

crosslinks <- read_excel("./outputs/Combined_Data_Crosslinks_Disambiguated_Indexed_2026-02-24.xlsx")

##### Aggregate crosslinks to PPIs #####

# Note: pmin/pmax on Leading Protein A/B is intentional here — unlike the residue
# pair script (where A/B order was already canonical), PPI aggregation needs to
# collapse e.g. A=Prot1/B=Prot2 and A=Prot2/B=Prot1 to the same interaction.
ppis <- crosslinks %>%
  # Keep only unambiguous inter-protein crosslinks
  filter(`Leading Protein A` != `Leading Protein B`) %>%
  filter(`# Accessions A` == 1 & `# Accessions B` == 1) %>%
  # Canonicalise protein order and build interaction identifier
  mutate(
    `Interaction Identifier` = paste(
      pmin(`Leading Protein A`, `Leading Protein B`),
      pmax(`Leading Protein A`, `Leading Protein B`),
      sep = "-"),
    `Protein A` = pmin(`Leading Protein A`, `Leading Protein B`),
    `Protein B` = pmax(`Leading Protein A`, `Leading Protein B`)) %>%
  group_by(`Is Decoy`, `Interaction Identifier`, `Protein A`, `Protein B`) %>%
  summarise(
    # PPI score = sum of contributing crosslink scores
    `Identification Score`      = sum(`Identification Score`),
    `Identification DeltaScore` = sum(`Identification DeltaScore`),
    `# CSMs`                    = sum(`# CSMs`),
    `# Crosslinks`              = n(),
    `Experiments Detected`      = paste0(
      sort(unique(na.omit(unlist(str_split(`Experiments Detected`, pattern=";"))))),
      collapse = ";"),
    `Disambiguated Accessions A` = paste0(
      unique(na.omit(unlist(str_split(`Disambiguated Accessions A`, pattern=";")))),
      collapse = ";"),
    `Disambiguated Accessions B` = paste0(
      unique(na.omit(unlist(str_split(`Disambiguated Accessions B`, pattern=";")))),
      collapse = ";"),
    .groups = "drop")

# Check FDR before filtering
ppis %>%
  group_by(`Is Decoy`) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(freq = n / sum(n)) %>%
  print()

##### Identify and apply optimal score cutoffs for 1% PPI-level FDR #####

# OR gate used here on the basis of target vs decoy distribution
ppis_score_cutoffs <- Identify_Score_Cutoffs(
  ppis, gate = "OR", fdr.level = 0.01,
  min.score = 20, max.score = 100,
  min.delta.score = 10, max.delta.score = 50)

if (nrow(ppis_score_cutoffs) == 0) stop("No valid score cutoffs found — check input data or relax threshold ranges.")
MIN_SCORE <- ppis_score_cutoffs[1, ]$score
MIN_DELTA <- ppis_score_cutoffs[1, ]$delta.score

# OR filter matches the OR gate used in cutoff search
ppis <- ppis %>%
  filter(`Identification Score` >= MIN_SCORE | `Identification DeltaScore` >= MIN_DELTA)

# Check FDR after filtering
ppis %>%
  group_by(`Is Decoy`) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(proportion = n / sum(n)) %>%
  print()

##### Calculate local FDR #####

ppis <- calculate_local_fdr(ppis, score_col="Identification Score", decoy_col="Is Decoy", bw=10, n=1024)

ggplot(ppis, aes(x=`Identification Score`, y=`Local FDR`)) +
  geom_point() +
  scale_x_log10()

ggplot(ppis, aes(x=`Identification Score`)) +
  geom_density() +
  scale_x_log10()

##### Discard decoys and save #####

ppis <- ppis %>%
  filter(!`Is Decoy`)

# How many proteins are in the final PPI network?
message("Proteins in PPI network: ", length(unique(c(ppis$`Protein A`, ppis$`Protein B`))))

# Save PPIs
openxlsx::write.xlsx(ppis, file=paste0("./outputs/Combined_Data_PPIs_", Sys.Date(), ".xlsx"))
