##### 06. Aggregate inter-protein crosslinks to PPIs and re-filter to 1% FDR #####

##### Setup #####

rm(list = ls())
gc()

# Force C collation. Protein accessions are compared with `<` and `>` below, which in R
# follows the locale's collation rules. C collation is stable across machines.
Sys.setlocale("LC_COLLATE", "C")

setwd("C:/Users/simon/Dropbox/2026-01 Toxoplasma XL-MS/R-Project-Toxoplasma-XL-MS/")
source("./script/PACKAGES.R")
lapply(list.files("./script/functions/", pattern = "\\.R$", full.names = TRUE), source)

crosslinks <- read_excel("./outputs/Combined_Data_Crosslinks_Disambiguated_Indexed_latest.xlsx")

# Collapse a set of semicolon-delimited lists into one sorted, de-duplicated list.
# Returns NA rather than "" when there is nothing to report.
merge_lists <- function(x) {
  values <- sort(unique(na.omit(unlist(str_split(x, ";")))))
  values <- values[values != ""]
  if (!length(values)) NA_character_ else paste(values, collapse = ";")
}

##### Aggregate crosslinks to PPIs #####

# A/B are already in canonical protein order from 04, so no pmin/pmax is needed. Applying
# one here would misalign the accession columns: `Disambiguated Accessions A` refers to the
# original A side, so reordering the protein columns alone would separate them.
stopifnot(with(crosslinks, all(`Leading Protein A` <= `Leading Protein B`)))

ppis <- crosslinks |>
  # Inter-protein crosslinks only; protein-level unambiguity is guaranteed by 03
  filter(`Leading Protein A` != `Leading Protein B`) |>
  mutate(`Interaction Identifier` = paste(`Leading Protein A`, `Leading Protein B`,
                                          sep = "-")) |>
  # ME49 IDs are a function of the leading protein, so grouping on them cannot split a group
  group_by(`Is Decoy`, `Interaction Identifier`,
           `Protein A`      = `Leading Protein A`,
           `Protein B`      = `Leading Protein B`,
           `Protein A ME49` = `Leading Protein ME49 A`,
           `Protein B ME49` = `Leading Protein ME49 B`) |>
  summarise(
    # PPI score = sum of contributing crosslink scores. This makes the score scale with
    # crosslink multiplicity: target PPIs accumulate crosslinks, whereas decoy PPIs (keyed
    # on DECOY_<peptide sequence>) almost never do, so the separation between the two
    # distributions is partly structural rather than a difference in confidence. The
    # PPI-level FDR below is therefore optimistic.
    `Identification Score`      = sum(`Identification Score`),
    `Identification DeltaScore` = sum(`Identification DeltaScore`),
    `# CSMs`                    = sum(`# CSMs`),
    `# Crosslinks`              = n(),
    # Distinct residue pairs supporting the interaction. Several crosslinks hitting one
    # residue pair is weaker evidence than the same number hitting several.
    `# Residue Pairs` = n_distinct(
      paste(`Crosslink Position A`, `Crosslink Position B`)[
        !is.na(`Crosslink Position A`) & !is.na(`Crosslink Position B`)]),
    `Experiments Detected`            = merge_lists(`Experiments Detected`),
    `Disambiguated Accessions A`      = merge_lists(`Disambiguated Accessions A`),
    `Disambiguated Accessions B`      = merge_lists(`Disambiguated Accessions B`),
    `Disambiguated Accessions ME49 A` = merge_lists(`Disambiguated Accessions ME49 A`),
    `Disambiguated Accessions ME49 B` = merge_lists(`Disambiguated Accessions ME49 B`),
    .groups = "drop")

# PPIs supported only by crosslinks whose residues could not be resolved
message(sum(ppis$`# Residue Pairs` == 0 & !ppis$`Is Decoy`),
        " target PPIs have no residue-resolved supporting crosslink.")

report_fdr(ppis)

##### Local FDR #####

# Computed before the score filter: the density estimate needs the full target and decoy
# distributions, and too few decoys survive filtering to support one.
ppis <- calculate_local_fdr(ppis, score_col = "Identification Score",
                            decoy_col = "Is Decoy", bw = 10, n = 1024)

ppi_local_fdr_plot <- ggplot(ppis, aes(x = `Identification Score`, y = `Local FDR`,
                                       colour = `Is Decoy`)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_x_log10() +
  simon_theme

print(ppi_local_fdr_plot)

ppi_score_plot <- ggplot(ppis, aes(x = `Identification Score`, colour = `Is Decoy`)) +
  geom_density() +
  scale_x_log10() +
  simon_theme

print(ppi_score_plot)

##### Identify and apply score cutoffs for 1% PPI-level FDR #####

# An OR gate is used here on the basis of the target and decoy score distributions, rather
# than the AND gate applied at crosslink and residue-pair level.
ppis_score_cutoffs <- Identify_Score_Cutoffs(
  ppis, gate = "OR", fdr.level = 0.01,
  min.score = 20, max.score = 100,
  min.delta.score = 10, max.delta.score = 50)

if (nrow(ppis_score_cutoffs) == 0) {
  stop("No valid score cutoffs found — check input data or relax threshold ranges.")
}

MIN_SCORE <- ppis_score_cutoffs[1, ]$score
MIN_DELTA <- ppis_score_cutoffs[1, ]$delta.score
message("PPI cutoffs (OR gate): score >= ", MIN_SCORE, " OR delta score >= ", MIN_DELTA)

# OR filter, matching the OR gate used in the cutoff search
ppis <- ppis |>
  filter(`Identification Score` >= MIN_SCORE | `Identification DeltaScore` >= MIN_DELTA) %>%
  arrange(`Is Decoy`)

report_fdr(ppis)

##### Save #####

ppis_targets <- ppis |> filter(!`Is Decoy`)

message("PPIs: ", nrow(ppis_targets),
        " | proteins in network: ",
        length(unique(c(ppis_targets$`Protein A`, ppis_targets$`Protein B`))))

openxlsx::write.xlsx(ppis,
                     file = paste0("./outputs/Combined_Data_PPIs_", Sys.Date(), ".xlsx"))
openxlsx::write.xlsx(ppis,
                     file = "./outputs/Combined_Data_PPIs_latest.xlsx")

