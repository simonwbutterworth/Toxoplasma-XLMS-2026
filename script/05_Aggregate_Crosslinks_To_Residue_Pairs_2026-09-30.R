##### 05. Aggregate crosslinks to unique residue pairs and re-filter to 1% FDR #####

##### Setup #####

rm(list = ls())
gc()

setwd("C:/Users/simon/Dropbox/2026-01 Toxoplasma XL-MS/R-Project-Toxoplasma-XL-MS/")
source("./script/PACKAGES.R")
lapply(list.files("./script/functions/", pattern = "\\.R$", full.names = TRUE), source)

crosslinks <- read_excel("./outputs/Combined_Data_Crosslinks_Disambiguated_Indexed_latest.xlsx")

##### Aggregate crosslinks to unique residue pairs #####

# A/B are already in canonical order from 04, by leading protein and then by residue
# position, so grouping directly on the four columns gives the residue pair identity and no
# pmin/pmax is needed. The position tie-break in 04 is what makes this true for
# intra-protein crosslinks; without it, 100-200 and 200-100 would be separate pairs.
#
# Decoy positions are in peptide coordinates, since decoys have no protein sequence, so
# decoy residue pairs are keyed on DECOY_<sequence> plus peptide position. They are
# comparable only with each other, which is all the FDR estimate requires.
residues <- crosslinks |>
  # Keep only crosslinks whose residue is resolved. Tested on the resolved position itself
  # rather than on `# Crosslink Positions`: since 04 resolves against the protein sequence,
  # XlinkX's position count is no longer a reliable proxy, and grouping on an NA position
  # would pool every unresolved crosslink in a protein into one spurious pair.
  filter(!is.na(`Crosslink Position A`), !is.na(`Crosslink Position B`)) |>
  group_by(`Leading Protein A`, `Crosslink Position A`,
           `Leading Protein B`, `Crosslink Position B`) |>
  mutate(
    # Spectral and crosslink support pooled across every crosslink for this residue pair.
    # Counted before the score filter below, so these are total support rather than
    # support above threshold.
    `# CSMs`              = sum(`# CSMs`),
    `# Unique Crosslinks` = n(),
    # Union of experiments across all contributing crosslinks, not only the representative
    `Experiments Detected` = paste(
      sort(unique(na.omit(unlist(str_split(`Experiments Detected`, ";"))))),
      collapse = ";")) |>
  ungroup() |>
  # Keep the highest-scoring crosslink as the representative for each residue pair
  arrange(desc(`Identification Score`)) |>
  distinct(`Leading Protein A`, `Crosslink Position A`,
           `Leading Protein B`, `Crosslink Position B`, .keep_all = TRUE) |>
  mutate(Residue_Pair_ID = paste0(`Leading Protein A`, ":", `Crosslink Position A`, "-",
                                  `Leading Protein B`, ":", `Crosslink Position B`))

report_fdr(residues)

##### Local FDR #####

# Computed before the score filter: the density estimate needs the full target and decoy
# score distributions, and after filtering to 1% FDR too few decoys remain to estimate one.
#
# This is a different quantity from the tail FDR used to set the cutoffs. Local FDR is the
# posterior probability that a single residue pair at a given score is incorrect, so its
# correct form is decoy_density / (decoy_density + target_density).
residues <- calculate_local_fdr(residues)

local_fdr_plot <- ggplot(residues, aes(x = `Identification Score`, y = `Local FDR`,
                                       colour = `Is Decoy`)) +
  geom_point(size = 0.5, alpha = 0.5)

print(local_fdr_plot)

##### Identify score cutoffs for 1% residue-pair-level FDR #####

residues_score_cutoffs <- Identify_Score_Cutoffs(
  residues, gate = "AND", fdr.level = 0.01,
  min.score = 20, max.score = 50,
  min.delta.score = 10, max.delta.score = 50)

if (nrow(residues_score_cutoffs) == 0) {
  stop("No valid score cutoffs found — check input data or relax threshold ranges.")
}

# Row 1 maximises target residue pairs at <= 1% FDR
MIN_SCORE <- residues_score_cutoffs[1, ]$score
MIN_DELTA <- residues_score_cutoffs[1, ]$delta.score
message("Residue-pair cutoffs: score >= ", MIN_SCORE, ", delta score >= ", MIN_DELTA)

residues <- residues |>
  filter(`Identification Score` >= MIN_SCORE, `Identification DeltaScore` >= MIN_DELTA)

report_fdr(residues)

##### Inspect #####

# Residue pairs linking the same residue in the same protein: two copies of a homodimer
# crosslinked at equivalent positions. These are real, but worth knowing the count.
n_self_pairs <- residues |>
  filter(!`Is Decoy`,
         `Leading Protein A` == `Leading Protein B`,
         `Crosslink Position A` == `Crosslink Position B`) |>
  nrow()

message("Self-residue pairs (possible homodimers): ", n_self_pairs)

# Intra- versus inter-protein
residues |>
  filter(!`Is Decoy`) |>
  count(`Is Inter` = `Leading Protein A` != `Leading Protein B`) |>
  mutate(freq = n / sum(n)) |>
  print()

##### Save data #####

openxlsx::write.xlsx(
  residues,
  file = paste0("./outputs/Combined_Data_Residue_Pairs_", Sys.Date(), ".xlsx"))
openxlsx::write.xlsx(
  residues,
  file = "./outputs/Combined_Data_Residue_Pairs_latest.xlsx")