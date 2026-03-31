##### Disambiguate accessions with no uniquely identifying crosslinks #####

# For this step, we assume that any protein accessions that do not have any
# non-ambiguous crosslinks are not expressed (or not expressed highly enough
# to detect any crosslinks) and therefore can be discarded from consideration,
# with the ambiguous crosslink assigned to the nonredundant accession. This is
# only implemented if the ambiguous side of the crosslink has a single
# nonredundant accession, as more complex cases require manual inspection.
#
# Arguments:
#   crosslinks   : data frame of crosslinks (output of aggregate_csms_to_crosslinks)
#   force_retain : optional character vector of accessions to always treat as
#                  nonredundant regardless of whether they have unambiguous
#                  crosslinks (e.g. force_retain = "TGRH88_057810" for poly-ubiquitin)
#
# TODO: add optional complete.crosslink.set parameter to allow disambiguation
#       using a larger reference set when processing a subset of crosslinks.

Disambiguate_Redundant_Accessions <- function(crosslinks, force_retain = character(0)) {
  
  # Validate inputs
  required_cols <- c(
    "Accessions A", "Accessions B",
    "# Accessions A", "# Accessions B",
    "Crosslink Positions A", "Crosslink Positions B",
    "# Crosslink Positions A", "# Crosslink Positions B",
    "Leading Protein A", "Leading Protein B")
  missing_cols <- setdiff(required_cols, colnames(crosslinks))
  if (length(missing_cols)) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  # Identify accessions with uniquely identifying (unambiguous) crosslinks,
  # plus any accessions explicitly requested to retain
  nonredundant_accessions <- unique(c(
    crosslinks$`Accessions A`[crosslinks$`# Accessions A` == 1],
    crosslinks$`Accessions B`[crosslinks$`# Accessions B` == 1],
    force_retain))
  
  # Initialise logging columns if not already present (supports iterative calls)
  if (!"Disambiguated A"            %in% colnames(crosslinks)) crosslinks$`Disambiguated A`            <- FALSE
  if (!"Disambiguated B"            %in% colnames(crosslinks)) crosslinks$`Disambiguated B`            <- FALSE
  if (!"Disambiguated Accessions A" %in% colnames(crosslinks)) crosslinks$`Disambiguated Accessions A` <- NA_character_
  if (!"Disambiguated Accessions B" %in% colnames(crosslinks)) crosslinks$`Disambiguated Accessions B` <- NA_character_
  
  # Pre-split all accession and position columns once outside rowwise()
  # to avoid repeated str_split calls per row
  accnsA_list <- strsplit(crosslinks$`Accessions A`,        ";")
  accnsB_list <- strsplit(crosslinks$`Accessions B`,        ";")
  posnsA_list <- strsplit(crosslinks$`Crosslink Positions A`, ";")
  posnsB_list <- strsplit(crosslinks$`Crosslink Positions B`, ";")
  
  crosslinks_disambig <- crosslinks %>%
    mutate(.row = row_number()) %>%
    rowwise() %>%
    mutate(
      # Retrieve pre-split vectors for this row
      accnsA = list(accnsA_list[[.row]]),
      accnsB = list(accnsB_list[[.row]]),
      posnsA = list(posnsA_list[[.row]]),
      posnsB = list(posnsB_list[[.row]]),
      
      # Which accessions on each side are in the nonredundant set
      nonambig_A = list(accnsA[accnsA %in% nonredundant_accessions]),
      nonambig_B = list(accnsB[accnsB %in% nonredundant_accessions]),
      
      # Flag rows to collapse: currently ambiguous but exactly one nonredundant accession exists
      Collapse_A = (`# Accessions A` > 1) & (length(nonambig_A) == 1),
      Collapse_B = (`# Accessions B` > 1) & (length(nonambig_B) == 1),
      
      # Index of the nonredundant accession within the full accession list
      # (used to pull the matching crosslink position)
      Match_Index_A = if (Collapse_A) match(nonambig_A[[1]], accnsA) else NA_integer_,
      Match_Index_B = if (Collapse_B) match(nonambig_B[[1]], accnsB) else NA_integer_,
      
      # Update accession and leading protein columns
      `Accessions A`      = if (Collapse_A) nonambig_A[[1]] else `Accessions A`,
      `Accessions B`      = if (Collapse_B) nonambig_B[[1]] else `Accessions B`,
      `Leading Protein A` = if (Collapse_A) nonambig_A[[1]] else `Leading Protein A`,
      `Leading Protein B` = if (Collapse_B) nonambig_B[[1]] else `Leading Protein B`,
      `# Accessions A`    = if (Collapse_A) 1L              else `# Accessions A`,
      `# Accessions B`    = if (Collapse_B) 1L              else `# Accessions B`,
      
      # Update crosslink position counts:
      # if positions == accessions (one position per accession, no within-protein ambiguity),
      # collapse to the single matching position
      `# Crosslink Positions A` = if (Collapse_A & length(posnsA) == length(accnsA)) 1L else `# Crosslink Positions A`,
      `# Crosslink Positions B` = if (Collapse_B & length(posnsB) == length(accnsB)) 1L else `# Crosslink Positions B`,
      
      # Pull the matched crosslink position as a plain string
      `Crosslink Positions A` = if (Collapse_A & `# Crosslink Positions A` == 1L) {
        posnsA[[Match_Index_A]]
      } else {
        paste(`Crosslink Positions A`, collapse = ";")
      },
      `Crosslink Positions B` = if (Collapse_B & `# Crosslink Positions B` == 1L) {
        posnsB[[Match_Index_B]]
      } else {
        paste(`Crosslink Positions B`, collapse = ";")
      },
      
      # Log whether disambiguation occurred and which accessions were removed
      `Disambiguated A` = `Disambiguated A` | Collapse_A,
      `Disambiguated B` = `Disambiguated B` | Collapse_B,
      `Disambiguated Accessions A` = as.character(if (Collapse_A) {
        removed <- setdiff(accnsA, nonambig_A[[1]])
        paste(c(na.omit(`Disambiguated Accessions A`), removed), collapse = ";")
      } else {
        `Disambiguated Accessions A`
      }),
      `Disambiguated Accessions B` = as.character(if (Collapse_B) {
        removed <- setdiff(accnsB, nonambig_B[[1]])
        paste(c(na.omit(`Disambiguated Accessions B`), removed), collapse = ";")
      } else {
        `Disambiguated Accessions B`
      })) %>%
    
    # Drop helper columns
    select(-accnsA, -accnsB, -posnsA, -posnsB,
           -nonambig_A, -nonambig_B,
           -Collapse_A, -Collapse_B,
           -Match_Index_A, -Match_Index_B,
           -.row) %>%
    ungroup()
  
  # Summary
  n_disambig_A <- sum(crosslinks_disambig$`Disambiguated A`, na.rm = TRUE)
  n_disambig_B <- sum(crosslinks_disambig$`Disambiguated B`, na.rm = TRUE)
  message("Disambiguated ", n_disambig_A, " Protein A and ", n_disambig_B, " Protein B accessions.")
  
  return(crosslinks_disambig)
}
