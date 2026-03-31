##### Disambiguate crosslinks #####

##### Setup #####

rm(list=ls())
gc()
setwd("C:/Users/simon/Dropbox/R-Project-Crosslinking/")
source("./script/PACKAGES.R")
lapply(list.files("./script/functions/", pattern = "\\.R$", full.names = TRUE), source)

crosslinks <- read_excel("./outputs/Combined_Data_Crosslinks_2026-02-24.xlsx")
output_file <- "Combined_Data"

# Helper: summarise proportion of ambiguous non-decoy crosslinks
# Note: explicit print() required when running as a sourced script
summarise_ambiguity <- function(crosslinks) {
  crosslinks %>%
    filter(!`Is Decoy`) %>%
    mutate(`Is Ambiguous` = `# Accessions A` > 1 | `# Accessions B` > 1) %>%
    group_by(`Is Ambiguous`) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(freq = n / sum(n)) %>%
    print()
}

##### Disambiguate crosslinks #####

# Check proportion of ambiguous crosslinks before disambiguation
summarise_ambiguity(crosslinks)

# Step 1: Disambiguate redundant accessions with no uniquely identifying crosslinks
crosslinks_disambig <- Disambiguate_Redundant_Accessions(crosslinks = crosslinks, force_retain = "TGRH88_057810")
summarise_ambiguity(crosslinks_disambig)

# Step 2: Collapse homolog groups to a representative member (first accession in each vector)
# List names are labels only; the first accession in each vector is used as the
# representative. Order is intentional — place the canonical/most abundant member first.
crosslinks_disambig <- Disambiguate_Homolog_Groups(
  crosslinks_disambig,
  homolog.list = list(
    # Group 1: Redundant accessions with more complicated patterns
    TGRH88_000170 = c("TGRH88_000170", "TGRH88_000150", "TGRH88_000130"),
    TGRH88_030830 = c("TGRH88_030830", "TGRH88_030810", "TGRH88_030820", "TGRH88_030800", "TGRH88_030770"),
    # Group 2: Identical and near-identical homologs that cannot be distinguished
    TGRH88_061120 = c("TGRH88_061120", "TGRH88_061080"),
    TGRH88_063130 = c("TGRH88_063130", "TGRH88_025950"),
    ROP29         = c("TGRH88_087780", "TGRH88_084170", "TGRH88_084110"),
    Beta_Tubulin  = c("TGRH88_011050", "TGRH88_043000", "TGRH88_000700"),
    TGRH88_058620 = c("TGRH88_058620", "TGRH88_057560"),
    TGRH88_024020 = c("TGRH88_024020", "TGRH88_024010"),
    SRS59KB       = c("TGRH88_064840", "TGRH88_064830"),
    GAPM1A        = c("TGRH88_036510", "TGRH88_036500"),
    TGRH88_056760 = c("TGRH88_056760", "TGRH88_056750", "TGRH88_056720"),
    TGRH88_016110 = c("TGRH88_016110", "TGRH88_016100", "TGRH88_016090"),
    # Group 3: Homolog groups assigned to the most abundant member (most uniquely identifying crosslinks)
    EF1gamma      = c("TGRH88_057220", "TGRH88_075760", "TGRH88_077590"),
    Histone_H3    = c("TGRH88_069360", "TGRH88_059090"),
    ROP7_ROP4     = c("TGRH88_020200", "TGRH88_020190", "TGRH88_020180"),
    NTPase_II     = c("TGRH88_066180", "TGRH88_066160", "TGRH88_066170"),
    UbC           = c("TGRH88_057810", "TGRH88_014510"), # polyUbC retained as leading protein; TGRH88_014510 is RPL40-UbC fusion
    TGRH88_016370 = c("TGRH88_016370", "TGRH88_057800"),
    Histone_H2Babv= c("TGRH88_016930", "TGRH88_021270", "TGRH88_064730"),
    HSP70_BIP     = c("TGRH88_078290", "TGRH88_064630", "TGRH88_051090"),
    EF1_GEF       = c("TGRH88_045970", "TGRH88_058460"),
    # Group 4: Manually curated accessions to collapse
    RPL7          = c("TGRH88_053800", "TGRH88_012560"),
    RPL27         = c("TGRH88_068330", "TGRH88_008040"),
    RPL27A        = c("TGRH88_069140", "TGRH88_052510"),
    MyoA_MyoC     = c("TGRH88_041990", "TGRH88_073760"),
    PSMC5_PSMC6   = c("TGRH88_047560", "TGRH88_069520"),
    RPL22         = c("TGRH88_031200", "TGRH88_020510"),
    RPL36         = c("TGRH88_032200", "TGRH88_021330")))
summarise_ambiguity(crosslinks_disambig)

# Step 3: Disambiguate likely intra-protein crosslinks where one peptide is unambiguous
crosslinks_disambig <- Disambiguate_Intra_Crosslinks(crosslinks_disambig)
summarise_ambiguity(crosslinks_disambig)

# Inspect remaining ambiguous crosslinks
crosslinks_ambig_remaining <- crosslinks_disambig %>%
  filter(!`Is Decoy`) %>%
  filter(`# Accessions A` > 1 | `# Accessions B` > 1)

# Of the unambiguous non-decoy crosslinks, how many are inter- vs. intra-protein?
crosslinks_disambig %>%
  filter(!`Is Decoy`) %>%
  filter(`# Accessions A` == 1 & `# Accessions B` == 1) %>%
  mutate(`Is Inter` = `Leading Protein A` != `Leading Protein B`) %>%
  group_by(`Is Inter`) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(freq = n / sum(n)) %>%
  print()

##### Save outputs #####

# Save unambiguous disambiguated crosslinks
openxlsx::write.xlsx(
  filter(crosslinks_disambig, `# Accessions A` == 1 & `# Accessions B` == 1),
  file = paste0("./outputs/", output_file, "_Crosslinks_Disambiguated_", Sys.Date(), ".xlsx"),
  na.string = "NA")

# Save all proteins with crosslinks (leading proteins + any accessions removed during disambiguation)
crosslinked_proteins <- bind_rows(
  crosslinks_disambig %>% filter(!`Is Decoy`) %>% select(Gene_ID = `Leading Protein A`),
  crosslinks_disambig %>% filter(!`Is Decoy`) %>% select(Gene_ID = `Leading Protein B`),
  data.frame(Gene_ID = crosslinks_disambig %>% filter(!`Is Decoy`) %>%
               pull(`Disambiguated Accessions A`) %>% na.omit() %>% str_split(";") %>% unlist()),
  data.frame(Gene_ID = crosslinks_disambig %>% filter(!`Is Decoy`) %>%
               pull(`Disambiguated Accessions B`) %>% na.omit() %>% str_split(";") %>% unlist())) %>%
  distinct() %>%
  filter(!is.na(Gene_ID))
openxlsx::write.xlsx(
  crosslinked_proteins,
  file = paste0("./outputs/", output_file, "_Crosslinked_Proteins_", Sys.Date(), ".xlsx"))
