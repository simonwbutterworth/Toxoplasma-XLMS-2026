##### Disambiguate crosslinks by favouring intra-protein links #####

# Intra-protein crosslinks can be assumed to be more probable than inter-protein
# crosslinks. In the case that one peptide of a crosslink has a single unique
# protein accession, but the other peptide has multiple possible accessions that
# include the unique accession of the first peptide, we assume this crosslink is
# most likely an intra-link on the uniquely identified protein and assign the
# ambiguous peptide to that protein.

Disambiguate_Intra_Crosslinks <- function(crosslinks) {
  
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
  
  # Initialise logging columns if not already present (supports iterative calls)
  if (!"Disambiguated A"            %in% colnames(crosslinks)) crosslinks$`Disambiguated A`            <- FALSE
  if (!"Disambiguated B"            %in% colnames(crosslinks)) crosslinks$`Disambiguated B`            <- FALSE
  if (!"Disambiguated Accessions A" %in% colnames(crosslinks)) crosslinks$`Disambiguated Accessions A` <- NA_character_
  if (!"Disambiguated Accessions B" %in% colnames(crosslinks)) crosslinks$`Disambiguated Accessions B` <- NA_character_
  
  # Pre-split accession and position columns once outside rowwise()
  accnsA_list <- strsplit(crosslinks$`Accessions A`,          ";")
  accnsB_list <- strsplit(crosslinks$`Accessions B`,          ";")
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
      
      # Collapse A: B side is unambiguous and its accession appears among A's accessions
      # Collapse B: A side is unambiguous and its accession appears among B's accessions
      # Note: use %in% rather than str_detect to avoid partial-string false matches
      # (e.g. "TGRH88_001" would match "TGRH88_0011" with str_detect)
      Collapse_A = (`# Accessions B` == 1) & (`# Accessions A` > 1) & (`Leading Protein B` %in% accnsA),
      Collapse_B = (`# Accessions A` == 1) & (`# Accessions B` > 1) & (`Leading Protein A` %in% accnsB),
      
      # Index of the unambiguous accession within the ambiguous list
      Match_Index_A = if (Collapse_A) {
        idx <- which(accnsA == `Leading Protein B`)
        if (length(idx) > 0) idx[[1]] else NA_integer_
      } else { NA_integer_ },
      Match_Index_B = if (Collapse_B) {
        idx <- which(accnsB == `Leading Protein A`)
        if (length(idx) > 0) idx[[1]] else NA_integer_
      } else { NA_integer_ },
      
      # Capture the original Leading Protein values before any mutation
      # to avoid the A update corrupting the B collapse condition (and vice versa)
      .orig_LP_A = `Leading Protein A`,
      .orig_LP_B = `Leading Protein B`,
      
      # Update accessions and leading proteins to the unambiguous accession
      `Leading Protein A` = if (Collapse_A) .orig_LP_B else `Leading Protein A`,
      `Leading Protein B` = if (Collapse_B) .orig_LP_A else `Leading Protein B`,
      `Accessions A`      = if (Collapse_A) .orig_LP_B else `Accessions A`,
      `Accessions B`      = if (Collapse_B) .orig_LP_A else `Accessions B`,
      `# Accessions A`    = if (Collapse_A) 1L         else `# Accessions A`,
      `# Accessions B`    = if (Collapse_B) 1L         else `# Accessions B`,
      
      # If positions == accessions (one position per accession), collapse count to 1
      `# Crosslink Positions A` = if (Collapse_A & length(posnsA) == length(accnsA)) 1L else `# Crosslink Positions A`,
      `# Crosslink Positions B` = if (Collapse_B & length(posnsB) == length(accnsB)) 1L else `# Crosslink Positions B`,
      
      # Pull the matched crosslink position as a plain string
      `Crosslink Positions A` = if (Collapse_A & `# Crosslink Positions A` == 1L &
                                    !is.na(Match_Index_A) & Match_Index_A <= length(posnsA)) {
        posnsA[[Match_Index_A]]
      } else {
        paste(posnsA, collapse = ";")
      },
      `Crosslink Positions B` = if (Collapse_B & `# Crosslink Positions B` == 1L &
                                    !is.na(Match_Index_B) & Match_Index_B <= length(posnsB)) {
        posnsB[[Match_Index_B]]
      } else {
        paste(posnsB, collapse = ";")
      },
      
      # Log whether disambiguation occurred and which accessions were removed
      `Disambiguated A` = `Disambiguated A` | Collapse_A,
      `Disambiguated B` = `Disambiguated B` | Collapse_B,
      `Disambiguated Accessions A` = as.character(if (Collapse_A) {
        removed <- setdiff(accnsA, .orig_LP_B)
        paste(c(na.omit(`Disambiguated Accessions A`), removed), collapse = ";")
      } else {
        `Disambiguated Accessions A`
      }),
      `Disambiguated Accessions B` = as.character(if (Collapse_B) {
        removed <- setdiff(accnsB, .orig_LP_A)
        paste(c(na.omit(`Disambiguated Accessions B`), removed), collapse = ";")
      } else {
        `Disambiguated Accessions B`
      })) %>%
    
    # Drop helper columns
    select(-accnsA, -accnsB, -posnsA, -posnsB,
           -Collapse_A, -Collapse_B,
           -Match_Index_A, -Match_Index_B,
           -.row, -.orig_LP_A, -.orig_LP_B) %>%
    ungroup()
  
  # Summary
  n_disambig_A <- sum(crosslinks_disambig$`Disambiguated A`, na.rm = TRUE)
  n_disambig_B <- sum(crosslinks_disambig$`Disambiguated B`, na.rm = TRUE)
  message("Disambiguated ", n_disambig_A, " Protein A and ", n_disambig_B, " Protein B accessions.")
  
  return(crosslinks_disambig)
}
