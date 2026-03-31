##### Combine all CSM data, filter to 1% FDR, update protein accessions to gene IDs #####

##### Setup #####

rm(list=ls())
gc()
setwd("C:/Users/simon/Dropbox/R-Project-Crosslinking/")
source("./script/PACKAGES.R")
lapply(list.files("./script/functions/", pattern = "\\.R$", full.names = TRUE), source)

MIN_SCORE <- 20
MIN_DELTA <- 10

##### Import and merge CSMs data for all experiments #####

# Import data
csms_nitro_rep1 <- read_excel("./data/20250314_N2_Lysate_DSSO_SEC_HPLC_Crosslink_MS2_Scans.xlsx")   %>% mutate(Experiment="N2 cavitation 1")
csms_nitro_rep2 <- read_excel("./data/20250501_N2_Lysate_DSSO_SEC_HPLC_Crosslink_MS2_Scans.xlsx")   %>% mutate(Experiment="N2 cavitation 2")
csms_fixed_rep1 <- read_excel("./data/20250523_Fixed_Cells_DSSO_SEC_HPLC_Crosslink_MS2_Scans.xlsx") %>% mutate(Experiment="Fixed cells 1")
csms_fixed_rep2 <- read_excel("./data/20250606_Fixed_Cells_DSSO_SEC_HPLC_Crosslink_MS2_Scans.xlsx") %>% mutate(Experiment="Fixed cells 2")

# Merge based on shared columns
shared.colnames <- Reduce(intersect, list(
  colnames(csms_nitro_rep1),
  colnames(csms_nitro_rep2),
  colnames(csms_fixed_rep1),
  colnames(csms_fixed_rep2)))
csms <- 
  bind_rows(
    csms_nitro_rep1[shared.colnames],
    csms_nitro_rep2[shared.colnames],
    csms_fixed_rep1[shared.colnames],
    csms_fixed_rep2[shared.colnames]) %>%
  # Fill in accessions for decoys with a dummy ID
  mutate(
    `Leading Protein A`    = if_else(`Is Decoy`, true=paste0("DECOY_", `Sequence A`), false=`Leading Protein A`),
    `Leading Protein B`    = if_else(`Is Decoy`, true=paste0("DECOY_", `Sequence B`), false=`Leading Protein B`),
    `Accessions Protein A` = if_else(`Is Decoy`, true=paste0("DECOY_", `Sequence A`), false=`Accessions Protein A`),
    `Accessions Protein B` = if_else(`Is Decoy`, true=paste0("DECOY_", `Sequence B`), false=`Accessions Protein B`)) %>%
  # Remove whitespace from accessions
  mutate(across(.cols=c(`Accessions Protein A`, `Accessions Protein B`), ~gsub(" ", "", .))) %>%
  # Count the number of unique accession and crosslink site positions for each crosslink
  mutate(
    `# Crosslink Positions A` = lengths(strsplit(as.character(`Crosslink Positions A`), ";")),
    `# Crosslink Positions B` = lengths(strsplit(as.character(`Crosslink Positions B`), ";")),
    `# Accessions A`= lengths(lapply(strsplit(as.character(`Accessions Protein A`), ";"), unique)),
    `# Accessions B`= lengths(lapply(strsplit(as.character(`Accessions Protein B`), ";"), unique)),
    `Sequence Length A`=nchar(`Sequence A`),
    `Sequence Length B`=nchar(`Sequence B`),
    `Is Ambiguous`=`# Accessions A`>1 | `# Accessions B`>1 | `# Crosslink Positions A`>1 | `# Crosslink Positions B`>1) 

##### Filter CSMs for minimum score cutoffs, count number and FDR ####

# Apply CSM score filters and check FDR
csms <- csms %>%
  filter(`Identification Score`>=MIN_SCORE & `Identification DeltaScore`>=MIN_DELTA)

# Count number of CSMs and check FDR of combined data
csms %>%
  group_by(`Is Decoy`) %>%
  summarise(n=n()) %>%
  mutate(proportion=n/sum(n))

# Count number of CSMs and check FDR for each experiment
csms %>%
  group_by(Experiment, `Is Decoy`) %>%
  summarise(n=n()) %>%
  mutate(proportion=n/sum(n))

##### Update protein accessions to gene IDs

# Update accession IDs
id_mapping <- read.csv("./inputs/ToxoDB_TGRH88_Aliases.csv") %>%
  select(Gene_ID, Alias) %>%
  filter(grepl("^KAF", Alias)) %>%
  mutate(Alias=as.character(Alias))
csms <- csms %>%
  # Update Leading Protein A
  left_join(id_mapping, by=c("Leading Protein A"="Alias")) %>% 
  mutate(`Leading Protein A`=coalesce(Gene_ID, `Leading Protein A`)) %>%
  select(-Gene_ID) %>%
  # Update Leading Protein B
  left_join(id_mapping, by=c("Leading Protein B"="Alias")) %>% 
  mutate(`Leading Protein B`=coalesce(Gene_ID, `Leading Protein B`)) %>%
  select(-Gene_ID) %>%
  # Update Accessions Protein A
  separate_rows(`Accessions Protein A`, sep = ";") %>%
  left_join(id_mapping, by = c("Accessions Protein A" = "Alias")) %>%
  mutate(`Accessions Protein A` = coalesce(Gene_ID, `Accessions Protein A`)) %>%
  select(-Gene_ID) %>%
  group_by(across(-`Accessions Protein A`)) %>%
  summarise(`Accessions Protein A` = paste(unique(`Accessions Protein A`), collapse = ";"), .groups = "drop") %>%
  # Update Accessions Protein B
  separate_rows(`Accessions Protein B`, sep = ";") %>%
  left_join(id_mapping, by = c("Accessions Protein B" = "Alias")) %>%
  mutate(`Accessions Protein B` = coalesce(Gene_ID, `Accessions Protein B`)) %>%
  select(-Gene_ID) %>%
  group_by(across(-`Accessions Protein B`)) %>%
  summarise(`Accessions Protein B` = paste(unique(`Accessions Protein B`), collapse = ";"), .groups = "drop")

##### Save CSMs #####

openxlsx::write.xlsx(csms, file=paste0("./outputs/Combined_Data_CSMs_", Sys.Date(), ".xlsx"))

