##### 02. Aggregate CSMs to unique crosslinks and re-filter to 1% crosslink-level FDR #####

##### Setup #####

rm(list = ls())
gc()

# Force C collation. The canonical A/B ordering below compares peptide sequences with `>`,
# which in R follows the locale's collation rules, so the same CSM could be ordered
# differently on another machine and give a different Crosslink_ID. C collation is stable.
Sys.setlocale("LC_COLLATE", "C")

setwd("C:/Users/simon/Dropbox/2026-01 Toxoplasma XL-MS/R-Project-Toxoplasma-XL-MS/")
source("./script/PACKAGES.R")
lapply(list.files("./script/functions/", pattern = "\\.R$", full.names = TRUE), source)

csms <- read_excel("./outputs/Combined_Data_CSMs_latest.xlsx")

##### Canonical A/B ordering #####

# XlinkX assigns A and B arbitrarily, so the same physical crosslink can be reported with
# its peptides in either order. Each CSM is sorted so that A is always the smaller
# sequence, with ties broken on crosslink site. Without this, one crosslink splits into
# two Crosslink_IDs.
csms <- swap_ab_columns(
  csms,
  csms$`Sequence A` > csms$`Sequence B` |
    (csms$`Sequence A` == csms$`Sequence B` &
       csms$`Peptide Position A` > csms$`Peptide Position B`))

##### Crosslink definition #####

# Peptide pair plus the crosslink site within each peptide. The site is essential, since
# the same pair of peptides can be linked between different residues. Position 0 denotes a
# peptide N-terminal link.
csms <- csms |>
  mutate(Crosslink_ID = paste(`Sequence A`, `Peptide Position A`,
                              `Sequence B`, `Peptide Position B`, sep = "|"))

##### Derive crosslink-level score cutoffs #####

# One row per crosslink carrying its best-scoring CSM, before any score filter, so the
# cutoffs are derived at crosslink level on the same data they will be applied to.
crosslinks_unfiltered <- csms |>
  arrange(desc(`Identification Score`)) |>
  distinct(Crosslink_ID, .keep_all = TRUE)

# FDR = decoys / targets
report_fdr(crosslinks_unfiltered)

# Cutoffs are derived once on the combined data and applied unchanged to every subset, so
# all subsets remain directly comparable. Row 1 maximises target crosslinks at <= 1% FDR.
crosslinks_score_cutoffs <- Identify_Score_Cutoffs(
  crosslinks_unfiltered, gate = "AND", fdr.level = 0.01,
  min.score = 20, max.score = 50, min.delta.score = 10, max.delta.score = 50)

if (nrow(crosslinks_score_cutoffs) == 0) {
  stop("No valid score cutoffs found — check input data or relax threshold ranges.")
}

MIN_SCORE <- crosslinks_score_cutoffs[1, ]$score
MIN_DELTA <- crosslinks_score_cutoffs[1, ]$delta.score
message("Crosslink-level cutoffs: score >= ", MIN_SCORE, ", delta score >= ", MIN_DELTA)

##### Aggregate CSMs to crosslinks #####

# Takes CSMs, or a subset of them, rather than an already-collapsed frame: aggregating from
# CSM level is what gives `# CSMs` and `Experiments Detected` their values, and it keeps
# Crosslink_ID identical across the combined data and every subset.
aggregate_csms_to_crosslinks <- function(csms_data, min_score, min_delta) {
  stopifnot("Crosslink_ID" %in% names(csms_data))
  
  csms_data |>
    # Count supporting CSMs and the experiments they came from. Counted before the score
    # filter, so `# CSMs` is total spectral support for the crosslink rather than support
    # above threshold.
    group_by(Crosslink_ID) |>
    mutate(
      `# CSMs`               = n(),
      `Experiments Detected` = paste(sort(unique(Experiment)), collapse = ";")) |>
    ungroup() |>
    # Keep the best-scoring CSM for each unique crosslink
    arrange(desc(`Identification Score`)) |>
    distinct(Crosslink_ID, .keep_all = TRUE) |>
    # Apply the crosslink-level 1% FDR cutoffs
    filter(
      `Identification Score`      >= min_score,
      `Identification DeltaScore` >= min_delta) |>
    # Keep only the columns used downstream. `Accessions X` is renamed here because the
    # disambiguation functions in 03 expect that name.
    select(
      Crosslink_ID, `Identification Score`, `Identification DeltaScore`, `Is Decoy`,
      `# CSMs`, `Experiments Detected`,
      `Sequence A`, `Sequence Length A`, `Peptide Position A`, `Modifications A`,
      `Crosslink localisation probability A`,
      `Crosslink Positions A`, `# Crosslink Positions A`,
      `Leading Protein A`, `Leading Protein ME49 A`,
      `Accessions A`      = `Accessions Protein A`,
      `Accessions ME49 A` = `Accessions Protein ME49 A`, `# Accessions A`,
      `Sequence B`, `Sequence Length B`, `Peptide Position B`, `Modifications B`,
      `Crosslink localisation probability B`,
      `Crosslink Positions B`, `# Crosslink Positions B`,
      `Leading Protein B`, `Leading Protein ME49 B`,
      `Accessions B`      = `Accessions Protein B`,
      `Accessions ME49 B` = `Accessions Protein ME49 B`, `# Accessions B`)
}

crosslinks <- aggregate_csms_to_crosslinks(csms, MIN_SCORE, MIN_DELTA)

# Spectral support must account for every CSM that passed into the aggregation
stopifnot(sum(crosslinks$`# CSMs`) <= nrow(csms))

report_fdr(crosslinks)

openxlsx::write.xlsx(crosslinks,
                     file = paste0("./outputs/Combined_Data_Crosslinks_", Sys.Date(), ".xlsx"))
openxlsx::write.xlsx(crosslinks,
                     file = "./outputs/Combined_Data_Crosslinks_latest.xlsx")

##### Aggregate data subsets #####

subsets <- list(
  N2_Lysate     = c("N2 cavitation 1", "N2 cavitation 2"),
  Fixed_Cells   = c("Fixed cells 1", "Fixed cells 2"),
  N2_Lysate_1   = "N2 cavitation 1",
  N2_Lysate_2   = "N2 cavitation 2",
  Fixed_Cells_1 = "Fixed cells 1",
  Fixed_Cells_2 = "Fixed cells 2")

stopifnot(all(unlist(subsets) %in% csms$Experiment))

crosslink_subsets <- imap(subsets, function(experiments, name) {
  out <- aggregate_csms_to_crosslinks(
    csms |> filter(Experiment %in% experiments),
    min_score = MIN_SCORE, min_delta = MIN_DELTA)
  openxlsx::write.xlsx(out,
                       file = paste0("./outputs/", name, "_Crosslinks_", Sys.Date(), ".xlsx"))
  out
})

# Crosslink counts and FDR per subset
bind_rows(
  report_fdr(crosslinks) |> mutate(subset = "Combined"),
  imap_dfr(crosslink_subsets, ~ report_fdr(.x) |> mutate(subset = .y))) |>
  select(subset, targets, decoys, FDR) |>
  print()