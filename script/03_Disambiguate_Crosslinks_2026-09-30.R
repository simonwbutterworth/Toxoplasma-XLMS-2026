##### 03. Disambiguate crosslink protein accessions #####

##### Setup #####

rm(list = ls())
gc()

setwd("C:/Users/simon/Dropbox/2026-01 Toxoplasma XL-MS/R-Project-Toxoplasma-XL-MS/")
source("./script/PACKAGES.R")
lapply(list.files("./script/functions/", pattern = "\\.R$", full.names = TRUE), source)
ME49_ORTHOLOGS <- "./inputs/TGRH88_TO_TGGT1_TGME49_Syntenic_OrthoMCLDB7_Groups.csv"

crosslinks <- read_excel("./outputs/Combined_Data_Crosslinks_latest.xlsx")

# Disambiguation rewrites the leading protein and accession columns, so the ME49 IDs
# carried over from 01 no longer describe the retained accessions. They are dropped here
# and re-mapped from the disambiguated IDs at the end of this script.
crosslinks <- crosslinks |>
  select(-matches("ME49"))

# Proportion of non-decoy crosslinks whose protein of origin is still ambiguous. This is
# protein-level ambiguity only: a crosslink can be unambiguous in protein but still
# ambiguous in which residue was linked (`# Crosslink Positions` > 1), which is resolved
# separately in 04. An explicit print() is needed when running as a sourced script.
summarise_ambiguity <- function(crosslinks, stage) {
  crosslinks |>
    filter(!`Is Decoy`) |>
    summarise(
      stage       = stage,
      unambiguous = sum(`# Accessions A` == 1 & `# Accessions B` == 1),
      ambiguous   = sum(`# Accessions A` > 1 | `# Accessions B` > 1)) |>
    mutate(freq_ambiguous = ambiguous / (ambiguous + unambiguous)) |>
    print()
  
  invisible(crosslinks)
}

##### Disambiguate #####

summarise_ambiguity(crosslinks, "before disambiguation")

# Step 1: discard accessions with no uniquely identifying crosslinks.
# TGRH88_057810 (poly-ubiquitin) is force-retained: it is genuinely expressed but its
# peptides are shared with the RPL40-UbC fusion, so it has no unambiguous crosslinks.
crosslinks_disambig <- Disambiguate_Redundant_Accessions(
  crosslinks, force_retain = "TGRH88_057810")

summarise_ambiguity(crosslinks_disambig, "after redundant accessions")

# Step 2: collapse homolog groups to a representative member. The first accession of each
# vector is the representative and list names are labels only, so the order within each
# vector is deliberate: the canonical or most abundant member is placed first.
crosslinks_disambig <- Disambiguate_Homolog_Groups(
  crosslinks_disambig,
  homolog.list = list(
    # Group 1: redundant accessions with more complicated patterns
    TGRH88_000170  = c("TGRH88_000170", "TGRH88_000150", "TGRH88_000130"),
    TGRH88_030830  = c("TGRH88_030830", "TGRH88_030810", "TGRH88_030820",
                       "TGRH88_030800", "TGRH88_030770"),
    # Group 2: identical and near-identical homologs that cannot be distinguished
    TGRH88_061120  = c("TGRH88_061120", "TGRH88_061080"),
    TGRH88_063130  = c("TGRH88_063130", "TGRH88_025950"),
    ROP29          = c("TGRH88_087780", "TGRH88_084170", "TGRH88_084110"),
    Beta_Tubulin   = c("TGRH88_011050", "TGRH88_043000", "TGRH88_000700"),
    TGRH88_058620  = c("TGRH88_058620", "TGRH88_057560"),
    TGRH88_024020  = c("TGRH88_024020", "TGRH88_024010"),
    SRS59KB        = c("TGRH88_064840", "TGRH88_064830"),
    GAPM1A         = c("TGRH88_036510", "TGRH88_036500"),
    TGRH88_056760  = c("TGRH88_056760", "TGRH88_056750", "TGRH88_056720"),
    TGRH88_016110  = c("TGRH88_016110", "TGRH88_016100", "TGRH88_016090"),
    # Group 3: homolog groups assigned to the most abundant member, defined as the member
    # with the most uniquely identifying crosslinks
    EF1gamma       = c("TGRH88_057220", "TGRH88_075760", "TGRH88_077590"),
    Histone_H3     = c("TGRH88_069360", "TGRH88_059090"),
    ROP7_ROP4      = c("TGRH88_020200", "TGRH88_020190", "TGRH88_020180"),
    NTPase_II      = c("TGRH88_066180", "TGRH88_066160", "TGRH88_066170"),
    # Poly-ubiquitin retained as the leading protein; TGRH88_014510 is the RPL40-UbC fusion
    UbC            = c("TGRH88_057810", "TGRH88_014510"),
    TGRH88_016370  = c("TGRH88_016370", "TGRH88_057800"),
    Histone_H2Babv = c("TGRH88_016930", "TGRH88_021270", "TGRH88_064730"),
    HSP70_BIP      = c("TGRH88_078290", "TGRH88_064630", "TGRH88_051090"),
    EF1_GEF        = c("TGRH88_045970", "TGRH88_058460"),
    # Group 4: manually curated accessions to collapse
    RPL7           = c("TGRH88_053800", "TGRH88_012560"),
    RPL27          = c("TGRH88_068330", "TGRH88_008040"),
    RPL27A         = c("TGRH88_069140", "TGRH88_052510"),
    MyoA_MyoC      = c("TGRH88_041990", "TGRH88_073760"),
    PSMC5_PSMC6    = c("TGRH88_047560", "TGRH88_069520"),
    RPL22          = c("TGRH88_031200", "TGRH88_020510"),
    RPL36          = c("TGRH88_032200", "TGRH88_021330")))

summarise_ambiguity(crosslinks_disambig, "after homolog groups")

# Step 3: assign likely intra-protein crosslinks where one peptide is unambiguous
crosslinks_disambig <- Disambiguate_Intra_Crosslinks(crosslinks_disambig)

summarise_ambiguity(crosslinks_disambig, "after intra-protein")

##### Inspect #####

# Crosslinks still ambiguous after all three strategies, for manual curation. Candidates
# here are what should be added to homolog.list above.
crosslinks_ambig_remaining <- crosslinks_disambig |>
  filter(!`Is Decoy`, `# Accessions A` > 1 | `# Accessions B` > 1)

message(nrow(crosslinks_ambig_remaining), " crosslinks remain ambiguous.")

# Of the unambiguous non-decoy crosslinks, how many are inter- versus intra-protein?
crosslinks_disambig |>
  filter(!`Is Decoy`, `# Accessions A` == 1, `# Accessions B` == 1) |>
  count(`Is Inter` = `Leading Protein A` != `Leading Protein B`) |>
  mutate(freq = n / sum(n)) |>
  print()

##### Re-map ME49 gene IDs #####

# Unambiguous crosslinks only, so `# Accessions A/B` are all 1 and are dropped along with
# the accession lists they counted. Decoys are carried through the rest of the pipeline so
# that residue-pair-level (05) and PPI-level (06) FDR can be estimated at those levels;
# crosslink-level filtering in 02 does not constrain either.
crosslinks_out <- crosslinks_disambig |>
  filter(`# Accessions A` == 1, `# Accessions B` == 1) |>
  select(-`Accessions A`, -`Accessions B`, -`# Accessions A`, -`# Accessions B`) |>
  mutate(across(c(`Leading Protein A`, `Leading Protein B`,
                  `Disambiguated Accessions A`, `Disambiguated Accessions B`),
                \(x) to_me49(x, ME49_ORTHOLOGS),
                .names = "{str_replace(.col, ' (?=[AB]$)', ' ME49 ')}"))

##### Save outputs #####

openxlsx::write.xlsx(
  crosslinks_out,
  file = paste0("./outputs/Combined_Data_Crosslinks_Disambiguated_", Sys.Date(), ".xlsx"))
openxlsx::write.xlsx(
  crosslinks_out,
  file = paste0("./outputs/Combined_Data_Crosslinks_Disambiguated_latest.xlsx"))

# Every protein implicated by a crosslink: the retained leading proteins, plus the
# accessions that disambiguation removed, which remain candidates rather than the assigned
# protein. Decoys are excluded, since this is a list of real proteins.
crosslinked_gene_ids <- crosslinks_disambig |>
  filter(!`Is Decoy`) |>
  select(`Leading Protein A`, `Leading Protein B`,
         `Disambiguated Accessions A`, `Disambiguated Accessions B`) |>
  pivot_longer(everything(), values_to = "Gene_ID") |>
  pull(Gene_ID) |>
  str_split(";") |>
  unlist() |>
  na.omit() |>
  unique()

crosslinked_proteins <- tibble(Gene_ID = crosslinked_gene_ids) |>
  mutate(ME49_Gene_ID = to_me49(Gene_ID, ME49_ORTHOLOGS))

openxlsx::write.xlsx(
  crosslinked_proteins,
  file = paste0("./outputs/Combined_Data_Crosslinked_Proteins_", Sys.Date(), ".xlsx"))