##### 01. Import and merge CSMs, map accessions to gene IDs, filter to minimum scores #####

##### Setup #####

rm(list = ls())
gc()

setwd("C:/Users/simon/Dropbox/2026-01 Toxoplasma XL-MS/R-Project-Toxoplasma-XL-MS/")
source("./script/PACKAGES.R")
lapply(list.files("./script/functions/", pattern = "\\.R$", full.names = TRUE), source)

MIN_SCORE <- 20
MIN_DELTA <- 10

##### Import and merge CSMs for all experiments #####

# One entry per experiment: label -> file. New experiments are added here only.
csm_files <- c(
  "N2 cavitation 1" = "20250314_N2_Lysate_DSSO_SEC_HPLC_Crosslink_MS2_Scans.xlsx",
  "N2 cavitation 2" = "20250501_N2_Lysate_DSSO_SEC_HPLC_Crosslink_MS2_Scans.xlsx",
  "Fixed cells 1"   = "20250523_Fixed_Cells_DSSO_SEC_HPLC_Crosslink_MS2_Scans.xlsx",
  "Fixed cells 2"   = "20250606_Fixed_Cells_DSSO_SEC_HPLC_Crosslink_MS2_Scans.xlsx")

# Read each XlinkX export, tag it with its experiment, and stack. bind_rows takes the union
# of columns, so a column missing from one export arrives as NA for that experiment rather
# than being dropped from all of them.
csms <- imap_dfr(csm_files, ~ read_excel(file.path("./data", .x)) |> mutate(Experiment = .y))

# Stop if an export is missing a column the pipeline depends on. `Peptide Position` is
# required by 02, which defines unique crosslinks by crosslinker site rather than by
# peptide sequence alone. `Crosslink localisation probability` is required by the residue
# indexing diagnostics in 04.
required_cols <- c(
  "Is Decoy", "Identification Score", "Identification DeltaScore",
  "Sequence A", "Sequence B",
  "Peptide Position A", "Peptide Position B",
  "Crosslink localisation probability A", "Crosslink localisation probability B",
  "Modifications A", "Modifications B",
  "Leading Protein A", "Leading Protein B",
  "Accessions Protein A", "Accessions Protein B",
  "Crosslink Positions A", "Crosslink Positions B")

missing_cols <- setdiff(required_cols, colnames(csms))
if (length(missing_cols)) {
  stop("Missing column(s) in CSM exports: ", paste(missing_cols, collapse = ", "))
}

# Drop uninformative columns. `Checked` is FALSE and `Identification Usage` is
# "Unconsidered" for every CSM. `Unmodified mass Peptide A MH+ [Da]` has no Peptide B
# counterpart in the export, so it cannot be swapped with the rest of the A/B pair in 02.
csms <- csms |>
  select(-any_of(c(
    "Checked",
    "Identification Usage",
    "Unmodified mass Peptide A MH+ [Da]")))

##### Clean accessions #####

csms <- csms |>
  mutate(
    # Decoys carry no accession. A dummy ID derived from the peptide sequence lets them
    # survive the grouping, joining and counting steps, so they remain available for FDR
    # estimation at every level of the pipeline.
    across(c(`Leading Protein A`, `Accessions Protein A`),
           ~ if_else(`Is Decoy`, paste0("DECOY_", `Sequence A`), .x)),
    across(c(`Leading Protein B`, `Accessions Protein B`),
           ~ if_else(`Is Decoy`, paste0("DECOY_", `Sequence B`), .x)),
    # Strip whitespace from the semicolon-delimited accession lists
    across(c(`Accessions Protein A`, `Accessions Protein B`), ~ gsub(" ", "", .x)))

##### Map protein accessions from KAF aliases to ToxoDB gene IDs #####

# Alias -> Gene_ID lookup, restricted to KAF (GenBank protein) aliases
alias_tbl <- read.csv("./inputs/ToxoDB_TGRH88_Aliases.csv") |>
  select(Gene_ID, Alias) |>
  filter(grepl("^KAF", Alias)) |>
  mutate(Alias = as.character(Alias))

# The mapping must be 1:1, otherwise a lookup would silently return only the first hit
stopifnot(!anyDuplicated(alias_tbl$Alias))

alias_map <- setNames(alias_tbl$Gene_ID, alias_tbl$Alias)

# Translate a semicolon-delimited accession string element by element, preserving order,
# length and duplicates.
#
# Duplicates are retained because `Accessions Protein X` is positionally paired with
# `Crosslink Positions X`: the nth accession corresponds to the nth crosslink position. The
# same protein appears twice when a peptide maps to it at two crosslinkable sites.
# De-duplicating the accessions would break that pairing, and the disambiguation functions
# in 03 — which take the index of a matched accession and use it to pull the corresponding
# position — would return the wrong residue. Unmatched accessions (decoys, contaminants)
# are passed through unchanged.
translate_accessions <- function(x) {
  vapply(strsplit(x, ";"),
         function(a) paste(coalesce(unname(alias_map[a]), a), collapse = ";"),
         character(1))
}

# translate_accessions would convert an NA into the literal string "NA"
stopifnot(!anyNA(csms$`Accessions Protein A`), !anyNA(csms$`Accessions Protein B`),
          !anyNA(csms$`Leading Protein A`),    !anyNA(csms$`Leading Protein B`))

csms <- csms |>
  mutate(across(c(`Leading Protein A`, `Leading Protein B`,
                  `Accessions Protein A`, `Accessions Protein B`),
                translate_accessions))

##### Count accessions and crosslink positions per CSM #####

# Counted after the gene-ID mapping so the counts describe the accessions actually stored
# in the table.
csms <- csms |>
  mutate(
    # Positions are counted as-is, not de-duplicated: two entries for the same accession
    # represent genuine ambiguity about which residue within that protein was crosslinked.
    `# Crosslink Positions A` = lengths(strsplit(as.character(`Crosslink Positions A`), ";")),
    `# Crosslink Positions B` = lengths(strsplit(as.character(`Crosslink Positions B`), ";")),
    # Accessions are de-duplicated for counting, since a protein listed twice is still one
    # protein. The duplicate is retained in the string to preserve pairing with positions.
    `# Accessions A`    = lengths(lapply(strsplit(`Accessions Protein A`, ";"), unique)),
    `# Accessions B`    = lengths(lapply(strsplit(`Accessions Protein B`, ";"), unique)),
    `Sequence Length A` = nchar(`Sequence A`),
    `Sequence Length B` = nchar(`Sequence B`),
    # A CSM is ambiguous if either the protein of origin or the crosslinked residue within
    # that protein is uncertain
    `Is Ambiguous` = `# Accessions A` > 1 | `# Accessions B` > 1 |
      `# Crosslink Positions A` > 1 | `# Crosslink Positions B` > 1)

# Assert the invariant relied on by 03 and 04: the accession list and the crosslink
# position list are the same length on both sides of every CSM.
check_pairing <- function(accessions, positions, side) {
  bad <- lengths(strsplit(accessions, ";")) != lengths(strsplit(as.character(positions), ";"))
  if (any(bad)) {
    stop(sum(bad), " CSM(s) have mismatched accession / crosslink-position lists on side ", side)
  }
}

check_pairing(csms$`Accessions Protein A`, csms$`Crosslink Positions A`, "A")
check_pairing(csms$`Accessions Protein B`, csms$`Crosslink Positions B`, "B")

##### Filter CSMs to minimum score cutoffs and report FDR #####

csms <- csms |>
  filter(`Identification Score` >= MIN_SCORE, `Identification DeltaScore` >= MIN_DELTA)

# FDR = decoys / targets
report_fdr(csms)               # combined data
report_fdr(csms, Experiment)   # per experiment

##### Annotate CSMs with ME49 gene IDs #####

# ME49 identifiers are added alongside the RH88 IDs, which remain the key used throughout
# the pipeline. The new columns keep the A/B suffix last so that they stay aligned with
# their own side when the A/B pairs are swapped into canonical order in 02.
csms_annotated <- csms |>
  mutate(across(c(`Leading Protein A`, `Leading Protein B`,
                  `Accessions Protein A`, `Accessions Protein B`),
                \(x) to_me49(x, "./inputs/TGRH88_TO_TGGT1_TGME49_Syntenic_OrthoMCLDB7_Groups.csv"),
                .names = "{str_replace(.col, ' (?=[AB]$)', ' ME49 ')}"))

##### Save CSMs #####

# Dated copy for provenance, plus a stable filename so downstream scripts do not need
# their input path edited on each run.
openxlsx::write.xlsx(csms_annotated,
                     file = paste0("./outputs/Combined_Data_CSMs_", Sys.Date(), ".xlsx"))
openxlsx::write.xlsx(csms_annotated,
                     file = "./outputs/Combined_Data_CSMs_latest.xlsx")