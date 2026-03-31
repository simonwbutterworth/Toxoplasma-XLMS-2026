##### Get crosslinked residue position (since PD assignments are inconsistent) #####

##### Setup #####

rm(list=ls())
gc()
setwd("C:/Users/simon/Dropbox/R-Project-Crosslinking/")
source("./script/PACKAGES.R")
lapply(list.files("./script/functions/", pattern = "\\.R$", full.names = TRUE), source)

crosslinks <- read_excel("./outputs/Combined_Data_Crosslinks_Disambiguated_2026-02-24.xlsx")
sequences  <- read.csv("inputs/ToxoDB_TGRH88_Protein_Sequences.csv")

##### For decoy crosslinks, assign position based only on identified peptide sequence (1-indexing) #####

# Peptide Position A/B are renamed to Crosslink Position A/B here (no protein sequence lookup needed).
# Note: Peptide Position A/B are therefore not present in select() below — they were consumed by rename().

crosslinks_decoys <- crosslinks %>%
  filter(`Is Decoy` == TRUE) %>%
  # Rename peptide positions as crosslink positions (decoys have no protein sequence to look up)
  rename(
    `Crosslink Position A` = `Peptide Position A`,
    `Crosslink Position B` = `Peptide Position B`) %>%
  # Identify crosslinked residue from peptide sequence
  mutate(
    `Crosslinked Residue A` = substr(x=`Sequence A`, start=`Crosslink Position A`, stop=`Crosslink Position A`),
    `Crosslinked Residue B` = substr(x=`Sequence B`, start=`Crosslink Position B`, stop=`Crosslink Position B`)) %>%
  # Fix N-term crosslinks (position 0 indicates N-terminal crosslink)
  # true=1L to match integer type of Crosslink Position columns
  mutate(
    `Crosslinked Residue A` = if_else(`Crosslink Position A` == 0L, true="N-term", false=`Crosslinked Residue A`),
    `Crosslinked Residue B` = if_else(`Crosslink Position B` == 0L, true="N-term", false=`Crosslinked Residue B`),
    `Crosslink Position A`  = if_else(`Crosslink Position A` == 0L, true=1L,       false=`Crosslink Position A`),
    `Crosslink Position B`  = if_else(`Crosslink Position B` == 0L, true=1L,       false=`Crosslink Position B`)) %>%
  select(
    -`Crosslink Positions A`,
    -`Crosslink Positions B`) %>%
  relocate(c(`Crosslink Position A`, `Crosslinked Residue A`), .after = "# Crosslink Positions A")


##### For target crosslinks, determine crosslinked residue position from protein sequence #####

crosslinks_target <- crosslinks %>%
  filter(`Is Decoy` == FALSE)

# For peptide A:
crosslinks_target <- crosslinks_target %>%
  left_join(sequences, by = c("Leading Protein A" = "Gene_ID")) %>%
  rename(`Protein Sequence A` = Protein_Sequence) %>%
  mutate(
    # Locate peptide within full protein sequence (string = full protein, pattern = peptide)
    `Sequence Position A`  = str_locate(string=`Protein Sequence A`, pattern=`Sequence A`)[, 1],
    `Crosslink Position A`  = as.integer(`Sequence Position A`) + as.integer(`Peptide Position A`) - 1L,
    `Crosslinked Residue A` = substr(`Protein Sequence A`, `Crosslink Position A`, `Crosslink Position A`),
    `Is Lysine`             = `Crosslinked Residue A` == "K") %>%
  # Fix: N-terminal crosslink on start methionine (peptide starts at position 1, crosslink at position 0)
  mutate(
    `Crosslinked Residue A` = if_else(`Sequence Position A` == 1L & `Crosslink Position A` == 0L, true="N-term", false=`Crosslinked Residue A`),
    `Crosslink Position A`  = if_else(`Sequence Position A` == 1L & `Crosslink Position A` == 0L, true=1L,       false=`Crosslink Position A`)) %>%
  # Fix: N-terminal crosslink on residue after start methionine (peptide starts at 2, crosslink at 1)
  mutate(
    `Crosslinked Residue A` = if_else(`Sequence Position A` == 2L & `Crosslink Position A` == 1L, true="N-term", false=`Crosslinked Residue A`),
    `Crosslink Position A`  = if_else(`Sequence Position A` == 2L & `Crosslink Position A` == 1L, true=2L,       false=`Crosslink Position A`)) %>%
  # If peptide occurs multiple times in protein or accession is ambiguous, position is unresolvable
  mutate(
    `Crosslinked Residue A` = if_else(`# Crosslink Positions A` > 1L, true=NA_character_, false=`Crosslinked Residue A`),
    `Crosslink Position A`  = if_else(`# Crosslink Positions A` > 1L, true=NA_integer_,   false=`Crosslink Position A`)) %>%
  select(
    -`Crosslink Positions A`,
    -`Peptide Position A`,
    -`Protein Sequence A`,
    -`Sequence Position A`,
    -`Is Lysine`) %>%
  relocate(c(`Crosslink Position A`, `Crosslinked Residue A`), .after = "# Crosslink Positions A")

# For peptide B:
crosslinks_target <- crosslinks_target %>%
  left_join(sequences, by = c("Leading Protein B" = "Gene_ID")) %>%
  rename(`Protein Sequence B` = Protein_Sequence) %>%
  mutate(
    `Sequence Position B`   = str_locate(string=`Protein Sequence B`, pattern=`Sequence B`)[, 1],
    `Crosslink Position B`  = as.integer(`Sequence Position B`) + as.integer(`Peptide Position B`) - 1L,
    `Crosslinked Residue B` = substr(`Protein Sequence B`, `Crosslink Position B`, `Crosslink Position B`),
    `Is Lysine`             = `Crosslinked Residue B` == "K") %>%
  # Fix: N-terminal crosslink on start methionine
  mutate(
    `Crosslinked Residue B` = if_else(`Sequence Position B` == 1L & `Crosslink Position B` == 0L, true="N-term", false=`Crosslinked Residue B`),
    `Crosslink Position B`  = if_else(`Sequence Position B` == 1L & `Crosslink Position B` == 0L, true=1L,       false=`Crosslink Position B`)) %>%
  # Fix: N-terminal crosslink on residue after start methionine
  mutate(
    `Crosslinked Residue B` = if_else(`Sequence Position B` == 2L & `Crosslink Position B` == 1L, true="N-term", false=`Crosslinked Residue B`),
    `Crosslink Position B`  = if_else(`Sequence Position B` == 2L & `Crosslink Position B` == 1L, true=2L,       false=`Crosslink Position B`)) %>%
  # If peptide occurs multiple times in protein or accession is ambiguous, position is unresolvable
  mutate(
    `Crosslinked Residue B` = if_else(`# Crosslink Positions B` > 1L, true=NA_character_, false=`Crosslinked Residue B`),
    `Crosslink Position B`  = if_else(`# Crosslink Positions B` > 1L, true=NA_integer_,   false=`Crosslink Position B`)) %>%
  select(
    -`Crosslink Positions B`,
    -`Peptide Position B`,
    -`Protein Sequence B`,
    -`Sequence Position B`,
    -`Is Lysine`) %>%
  relocate(c(`Crosslink Position B`, `Crosslinked Residue B`), .after = "# Crosslink Positions B")

# Recombine target and decoy crosslinks
crosslinks <- bind_rows(crosslinks_decoys, crosslinks_target)


#####
##### Organise A and B in protein alphanumeric order #####

crosslinks <- crosslinks %>%
  mutate(swap = `Leading Protein A` > `Leading Protein B`)

# Identify paired A/B columns
a_cols <- names(crosslinks)[str_ends(names(crosslinks), " A")]
b_cols <- str_replace(a_cols, " A$", " B")

# Sanity check: every A column must have a matching B column
stopifnot(all(b_cols %in% names(crosslinks)))

# Swap all A/B pairs where swap == TRUE, using temporary columns to avoid overwrite conflicts
for (i in seq_along(a_cols)) {
  a   <- a_cols[i]
  b   <- b_cols[i]
  tmp <- paste0("__tmp__", i)
  
  if (!identical(class(crosslinks[[a]]), class(crosslinks[[b]]))) {
    warning("Column type mismatch before swap: ", a, " vs ", b)
  }
  
  crosslinks[[tmp]] <- crosslinks[[a]]
  crosslinks[[a]]   <- if_else(crosslinks$swap, crosslinks[[b]],   crosslinks[[a]])
  crosslinks[[b]]   <- if_else(crosslinks$swap, crosslinks[[tmp]], crosslinks[[b]])
  crosslinks[[tmp]] <- NULL
}

crosslinks <- crosslinks %>% select(-swap)

# Create interaction identifier (A and B are now in alphanumeric order)
crosslinks <- crosslinks %>%
  mutate(Protein_Interaction_Identifier = paste(`Leading Protein A`, `Leading Protein B`, sep="-"))


#####
##### Save results #####

openxlsx::write.xlsx(
  crosslinks,
  file = paste0("./outputs/Combined_Data_Crosslinks_Disambiguated_Indexed_", Sys.Date(), ".xlsx"))
