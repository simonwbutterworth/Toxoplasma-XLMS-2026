##### 04. Convert peptide crosslink sites to protein residue positions #####

##### Setup #####

rm(list = ls())
gc()

# Force C collation. The canonical A/B ordering below compares accessions with `>`, which
# in R follows the locale's collation rules, so the same crosslink could be ordered
# differently on another machine. C collation is stable.
Sys.setlocale("LC_COLLATE", "C")

setwd("C:/Users/simon/Dropbox/2026-01 Toxoplasma XL-MS/R-Project-Toxoplasma-XL-MS/")
source("./script/PACKAGES.R")
lapply(list.files("./script/functions/", pattern = "\\.R$", full.names = TRUE), source)

crosslinks <- read_excel("./outputs/Combined_Data_Crosslinks_Disambiguated_latest.xlsx")
sequences  <- read.csv("./inputs/ToxoDB_TGRH88_Protein_Sequences.csv")

# A repeated Gene_ID would duplicate crosslinks on lookup
stopifnot(!anyDuplicated(sequences$Gene_ID))

protein_seq <- setNames(sequences$Protein_Sequence, sequences$Gene_ID)

##### Locate the crosslinked residue in the protein #####

# The protein position is derived from the protein sequence rather than taken from XlinkX:
#
#     crosslink position = peptide start in protein + peptide position - 1
#
# For an N-terminal link (Peptide Position 0) it is simply the peptide start, since the
# link is on the peptide's first residue.
#
# XlinkX's own `Crosslink Positions` column is used only to break ties, because:
#   * it is 0-indexed — the sequence-derived position is exactly XlinkX + 1 in all 108,632
#     comparable CSMs, and reading it as 1-based lands on an essentially random residue;
#   * it is the sentinel -1 for every N-terminal link.
# Where a peptide occurs more than once in its protein, the XlinkX position is what picks
# the correct occurrence.
locate_crosslink <- function(peptide, peptide_position, protein, xlinkx_positions) {
  if (is.na(protein) || is.na(peptide) || is.na(peptide_position)) return(NA_integer_)
  
  starts <- str_locate_all(protein, fixed(peptide))[[1]][, "start"]
  if (!length(starts)) return(NA_integer_)
  
  # One candidate protein position per occurrence of the peptide
  candidates <- if (peptide_position == 0L) starts else starts + peptide_position - 1L
  if (length(candidates) == 1L) return(as.integer(candidates))
  
  # Peptide repeats: keep the candidate XlinkX also reports, converting 0- to 1-based.
  # N-terminal links carry the -1 sentinel and cannot be resolved this way, so they fall
  # through to NA, which is correct since the occurrence is genuinely unknown.
  reported <- suppressWarnings(as.integer(unlist(str_split(xlinkx_positions, ";")))) + 1L
  hits     <- unique(intersect(candidates, reported))
  if (length(hits) == 1L) as.integer(hits) else NA_integer_
}

# Index one side, A or B
index_side <- function(crosslinks, side) {
  peptide  <- crosslinks[[paste0("Sequence ", side)]]
  pep_pos  <- as.integer(crosslinks[[paste0("Peptide Position ", side)]])
  protein  <- unname(protein_seq[crosslinks[[paste0("Leading Protein ", side)]]])
  xlinkx   <- as.character(crosslinks[[paste0("Crosslink Positions ", side)]])
  is_decoy <- crosslinks$`Is Decoy`
  
  position <- vapply(
    seq_along(peptide),
    function(i) locate_crosslink(peptide[i], pep_pos[i], protein[i], xlinkx[i]),
    integer(1))
  
  # Decoys have no protein sequence, so they stay in peptide coordinates. They are retained
  # so that residue-pair (05) and PPI-level (06) FDR can be estimated; their positions are
  # not comparable with target positions, only with each other. An N-terminal link (peptide
  # position 0) is placed on residue 1.
  position <- if_else(is_decoy, pmax(pep_pos, 1L), position)
  
  residue <- case_when(
    pep_pos == 0L ~ "N-term",
    is_decoy      ~ str_sub(peptide, pmax(pep_pos, 1L), pmax(pep_pos, 1L)),
    TRUE          ~ str_sub(protein, position, position))
  residue[is.na(position)] <- NA_character_
  
  crosslinks[[paste0("Crosslink Position ",  side)]] <- position
  crosslinks[[paste0("Crosslinked Residue ", side)]] <- residue
  crosslinks
}

crosslinks <- index_side(crosslinks, "A")
crosslinks <- index_side(crosslinks, "B")

##### Check the indexing #####

# Every DSSO crosslink must be on a lysine or a protein N-terminus. This is the strongest
# validation available for the position formula and the 1-based indexing.
for (side in c("A", "B")) {
  residue <- crosslinks[[paste0("Crosslinked Residue ", side)]]
  wrong   <- !crosslinks$`Is Decoy` & !is.na(residue) & !residue %in% c("K", "N-term")
  if (any(wrong)) {
    warning(sum(wrong), " target crosslink(s) on side ", side,
            " are on neither a lysine nor an N-terminus: ",
            paste(unique(residue[wrong]), collapse = ", "))
  }
}

# Sides whose residue could not be resolved: peptide absent from the protein, protein
# absent from the sequence file, or a repeated peptide that XlinkX does not pin down.
# These keep NA positions and will not form residue pairs in 05.
crosslinks |>
  filter(!`Is Decoy`) |>
  summarise(
    unresolved_A     = sum(is.na(`Crosslink Position A`)),
    unresolved_B     = sum(is.na(`Crosslink Position B`)),
    missing_sequence = sum(!`Leading Protein A` %in% names(protein_seq) |
                             !`Leading Protein B` %in% names(protein_seq))) |>
  print()

# Disagreements with XlinkX, after correcting its 0-based indexing. Reported only: the
# protein sequence is the authority for where a peptide sits. Two known causes, both
# XlinkX quirks rather than errors in the indexing here:
#   * peptide starts at protein position 2 — XlinkX numbers a peptide beginning just after
#     the initiator methionine as though it began at position 1 (0/515 agree in the raw
#     data at start 2, versus ~100% at starts 1, 3 and beyond);
#   * repetitive peptides with two candidate lysines (e.g. VGGKVGGK), where XlinkX's
#     `Crosslink Positions` names a different lysine than its own `Peptide Position`.
#     These carry a low localisation probability.
diagnose_xlinkx_disagreement <- function(crosslinks, side) {
  # Extract before building the tibble: inside tibble(), a column named `side` would shadow
  # the argument and be recycled to nrow(), breaking the paste0() lookups below.
  peptide  <- crosslinks[[paste0("Sequence ", side)]]
  protein  <- unname(protein_seq[crosslinks[[paste0("Leading Protein ", side)]]])
  position <- crosslinks[[paste0("Crosslink Position ", side)]]
  residue  <- crosslinks[[paste0("Crosslinked Residue ", side)]]
  loc_prob <- crosslinks[[paste0("Crosslink localisation probability ", side)]]
  reported <- str_split(as.character(crosslinks[[paste0("Crosslink Positions ", side)]]), ";")
  
  tibble(
    peptide_side  = side,
    is_decoy      = crosslinks$`Is Decoy`,
    residue       = residue,
    position      = position,
    peptide_start = str_locate(protein, fixed(peptide))[, "start"],
    localised     = loc_prob >= 0.95,
    agrees        = mapply(function(p, r) p %in% (suppressWarnings(as.integer(r)) + 1L),
                           position, reported)) |>
    filter(!is_decoy, !is.na(position), residue != "N-term", !agrees) |>
    count(peptide_side,
          cause = case_when(peptide_start == 2 ~ "peptide starts after initiator Met",
                            !localised         ~ "site not confidently localised",
                            TRUE               ~ "unexplained"))
}

bind_rows(diagnose_xlinkx_disagreement(crosslinks, "A"),
          diagnose_xlinkx_disagreement(crosslinks, "B")) |>
  print()

##### Order A and B canonically #####

# Ordered by leading protein, then by residue position. The position tie-break matters:
# without it, intra-protein crosslinks keep whatever order XlinkX gave them, so the same
# link appears as both 100-200 and 200-100 and splits into two residue pairs in 05.
crosslinks <- swap_ab_columns(
  crosslinks,
  crosslinks$`Leading Protein A` > crosslinks$`Leading Protein B` |
    (crosslinks$`Leading Protein A` == crosslinks$`Leading Protein B` &
       coalesce(crosslinks$`Crosslink Position A` > crosslinks$`Crosslink Position B`, FALSE)))

crosslinks <- crosslinks |>
  mutate(Protein_Interaction_Identifier = paste(`Leading Protein A`, `Leading Protein B`,
                                                sep = "-")) |>
  # Drop the XlinkX position lists now they have been used: they are 0-indexed and would be
  # misread alongside the 1-indexed `Crosslink Position` columns.
  select(-`Crosslink Positions A`, -`Crosslink Positions B`,
         -`Peptide Position A`, -`Peptide Position B`) |>
  relocate(`Crosslink Position A`, `Crosslinked Residue A`,
           .after = `# Crosslink Positions A`) |>
  relocate(`Crosslink Position B`, `Crosslinked Residue B`,
           .after = `# Crosslink Positions B`)

##### Save results #####

openxlsx::write.xlsx(
  crosslinks,
  file = paste0("./outputs/Combined_Data_Crosslinks_Disambiguated_Indexed_",
                Sys.Date(), ".xlsx"))
openxlsx::write.xlsx(
  crosslinks,
  file = "./outputs/Combined_Data_Crosslinks_Disambiguated_Indexed_latest.xlsx")

##### Count number of crosslinks per protein #####

# Inter-protein crosslinks are counted for both partners; intra-protein crosslinks are
# counted once, by dropping the duplicated B side.
link_count <- crosslinks |>
  filter(!`Is Decoy`) |>
  mutate(Type = if_else(`Leading Protein A` == `Leading Protein B`,
                        "Intra_Links", "Inter_Links")) |>
  pivot_longer(c(`Leading Protein A`, `Leading Protein B`), values_to = "Protein") |>
  filter(!(Type == "Intra_Links" & name == "Leading Protein B")) |>
  count(Protein, Type) |>
  pivot_wider(names_from = Type, values_from = n, values_fill = 0) |>
  mutate(Total_Links = Intra_Links + Inter_Links)

openxlsx::write.xlsx(
  link_count,
  file = "./outputs/Combined_Data_Crosslinks_Per_Protein_latest.xlsx")