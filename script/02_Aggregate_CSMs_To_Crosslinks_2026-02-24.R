##### Aggregate CSMs to unique crosslinks and re-filter to 1% #####
#####
##### Setup #####

rm(list=ls())
gc()
setwd("C:/Users/simon/Dropbox/R-Project-Crosslinking/")
source("./script/PACKAGES.R")
lapply(list.files("./script/functions/", pattern = "\\.R$", full.names = TRUE), source)
conflicts_prefer(base::union)

csms <- read_excel("./outputs/Combined_Data_CSMs_2026-02-24.xlsx")

#####
##### Aggregate combined data to crosslinks #####

# Aggregate to crosslinks and check FDR without filtering
crosslinks_unfiltered <- csms %>%
  mutate(Crosslink_ID=paste(pmin(`Sequence A`, `Sequence B`), pmax(`Sequence A`, `Sequence B`), sep="_")) %>%
  arrange(desc(`Identification Score`)) %>%
  distinct(Crosslink_ID, .keep_all = T)
crosslinks_unfiltered %>%
  group_by(`Is Decoy`) %>%
  summarise(n=n()) %>%
  mutate(proportion=n/sum(n))

# Identify optimal score cutoffs to obtain 1% crosslink-level FDR
crosslinks_score_cutoffs <- Identify_Score_Cutoffs(crosslinks_unfiltered, gate = "AND", fdr.level=0.01, min.score = 20, max.score = 50, min.delta.score = 10, max.delta.score = 50)
# Cutoffs derived from combined data and applied consistently to all subsets.
# Row 1 is the combination maximising target crosslinks at <= 1% FDR.
if (nrow(crosslinks_score_cutoffs) == 0) stop("No valid score cutoffs found — check input data or relax threshold ranges.")
MIN_SCORE <- crosslinks_score_cutoffs[1, ]$score
MIN_DELTA <- crosslinks_score_cutoffs[1, ]$delta.score

# Create a function to aggregate CSMs to crosslinks
aggregate_csms_to_crosslinks <- function(csms_data, min_score, min_delta) {
  csms_data %>%
    #Identify unique peptide sequence pairs
    mutate(Crosslink_ID = paste(
      pmin(`Sequence A`, `Sequence B`),
      pmax(`Sequence A`, `Sequence B`),
      sep="_")) %>%
    #Count how many CSMs support the sequence pair, which experiments the CSMs come from 
    group_by(Crosslink_ID) %>%
    mutate(
      `# CSMs` = n(),
      `Experiments Detected`=paste0(sort(unique(Experiment)), collapse=";")) %>%
    ungroup() %>%
    # Keep the CSM with the best score for each unique peptide pair
    arrange(desc(`Identification Score`)) %>%
    distinct(Crosslink_ID, .keep_all = T) %>%
    select(-Crosslink_ID) %>%
    #Apply score filters for 1% crosslink-level FDR
    filter(
      `Identification Score`>=min_score,
      `Identification DeltaScore`>=min_delta) %>%
    #Keep only useful columns
    select(
      `Identification Score`, `Identification DeltaScore`, `Is Decoy`, `# CSMs`, `Experiments Detected`,
      `Sequence A`, `Sequence Length A`, `Leading Protein A`, `Accessions A`=`Accessions Protein A`, `# Accessions A`,`Crosslink Positions A`, `Peptide Position A`, `# Crosslink Positions A`,
      `Sequence B`, `Sequence Length B`, `Leading Protein B`, `Accessions B`=`Accessions Protein B`, `# Accessions B`,`Crosslink Positions B`, `Peptide Position B`, `# Crosslink Positions B`)
}

#Aggregate CSMs to crosslinks based on unique peptide sequence pairs, filter to 1% FDR
crosslinks <- aggregate_csms_to_crosslinks(csms, min_score = MIN_SCORE, min_delta = MIN_DELTA)
openxlsx::write.xlsx(crosslinks, file=paste0("./outputs/Combined_Data_Crosslinks_", Sys.Date(), ".xlsx"))

#####
##### Aggregate data subsets #####

crosslinks_n2_cavitation <- aggregate_csms_to_crosslinks(csms %>% filter(Experiment %in% c("N2 cavitation 1", "N2 cavitation 2")), min_score = MIN_SCORE, min_delta = MIN_DELTA)
openxlsx::write.xlsx(crosslinks_n2_cavitation, file=paste0("./outputs/N2_Lysate_Crosslinks_", Sys.Date(), ".xlsx"))

crosslinks_fixed <- aggregate_csms_to_crosslinks(csms %>% filter(Experiment %in% c("Fixed cells 1", "Fixed cells 2")), min_score = MIN_SCORE, min_delta = MIN_DELTA)
openxlsx::write.xlsx(crosslinks_fixed, file=paste0("./outputs/Fixed_Cells_Crosslinks_", Sys.Date(), ".xlsx"))

crosslinks_n2_cavitation_1 <- aggregate_csms_to_crosslinks(csms %>% filter(Experiment=="N2 cavitation 1"), min_score = MIN_SCORE, min_delta = MIN_DELTA)
openxlsx::write.xlsx(crosslinks_n2_cavitation_1, file=paste0("./outputs/N2_Lysate_1_Crosslinks_", Sys.Date(), ".xlsx"))

crosslinks_n2_cavitation_2 <- aggregate_csms_to_crosslinks(csms %>% filter(Experiment=="N2 cavitation 2"), min_score = MIN_SCORE, min_delta = MIN_DELTA)
openxlsx::write.xlsx(crosslinks_n2_cavitation_2, file=paste0("./outputs/N2_Lysate_2_Crosslinks_", Sys.Date(), ".xlsx"))

crosslinks_fixed_1 <- aggregate_csms_to_crosslinks(csms %>% filter(Experiment=="Fixed cells 1"), min_score = MIN_SCORE, min_delta = MIN_DELTA)
openxlsx::write.xlsx(crosslinks_fixed_1, file=paste0("./outputs/Fixed_Cells_1_Crosslinks_", Sys.Date(), ".xlsx"))

crosslinks_fixed_2 <- aggregate_csms_to_crosslinks(csms %>% filter(Experiment=="Fixed cells 2"), min_score = MIN_SCORE, min_delta = MIN_DELTA)
openxlsx::write.xlsx(crosslinks_fixed_2, file=paste0("./outputs/Fixed_Cells_2_Crosslinks_", Sys.Date(), ".xlsx"))
