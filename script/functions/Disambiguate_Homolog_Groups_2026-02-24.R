##### Collapse manually defined homolog groups in crosslink data #####

# Some crosslinks are ambiguous because the accessions are identical or near-
# identical homologs. In this case, the accessions cannot and need not be
# distinguished and can be collapsed to a representative member.
#
# Arguments:
#   crosslinks    : data frame of crosslinks (output of aggregate_csms_to_crosslinks
#                   or Disambiguate_Redundant_Accessions)
#   homolog.list  : named or unnamed list of character vectors, each defining a
#                   homolog group. The first element of each vector is used as
#                   the representative accession.
#                   Example: list(c("TGRH88_001234", "TGRH88_001235"),
#                                 c("TGRH88_005678", "TGRH88_005679"))

Disambiguate_Homolog_Groups <- function(crosslinks, homolog.list) {
  
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
  if (!is.list(homolog.list) || length(homolog.list) == 0) {
    stop("homolog.list must be a non-empty list of character vectors.")
  }
  
  # Initialise logging columns if not already present (supports iterative calls)
  if (!"Disambiguated A"            %in% colnames(crosslinks)) crosslinks$`Disambiguated A`            <- FALSE
  if (!"Disambiguated B"            %in% colnames(crosslinks)) crosslinks$`Disambiguated B`            <- FALSE
  if (!"Disambiguated Accessions A" %in% colnames(crosslinks)) crosslinks$`Disambiguated Accessions A` <- NA_character_
  if (!"Disambiguated Accessions B" %in% colnames(crosslinks)) crosslinks$`Disambiguated Accessions B` <- NA_character_
  
  # Iterate through each homolog group
  for (protein.group in homolog.list) {
    
    replacement <- protein.group[1]
    
    # Pre-split accession and position columns once per iteration
    accnsA_list <- strsplit(crosslinks$`Accessions A`,          ";")
    accnsB_list <- strsplit(crosslinks$`Accessions B`,          ";")
    posnsA_list <- strsplit(crosslinks$`Crosslink Positions A`, ";")
    posnsB_list <- strsplit(crosslinks$`Crosslink Positions B`, ";")
    
    crosslinks <- crosslinks %>%
      mutate(.row = row_number()) %>%
      rowwise() %>%
      mutate(
        # Retrieve pre-split vectors for this row
        accnsA = list(accnsA_list[[.row]]),
        accnsB = list(accnsB_list[[.row]]),
        posnsA = list(posnsA_list[[.row]]),
        posnsB = list(posnsB_list[[.row]]),
        
        # Collapse if: currently ambiguous, all accessions are within the homolog
        # group, and the replacement accession is among them
        Collapse_A = (`# Accessions A` > 1) & all(accnsA %in% protein.group) & (replacement %in% accnsA),
        Collapse_B = (`# Accessions B` > 1) & all(accnsB %in% protein.group) & (replacement %in% accnsB),
        
        # Index of the replacement value within the accession list
        Match_Index_A = if (Collapse_A) match(replacement, accnsA) else NA_integer_,
        Match_Index_B = if (Collapse_B) match(replacement, accnsB) else NA_integer_,
        
        # Update accession and leading protein to the representative member
        `Accessions A`      = if (Collapse_A) replacement else `Accessions A`,
        `Accessions B`      = if (Collapse_B) replacement else `Accessions B`,
        `Leading Protein A` = if (Collapse_A) replacement else `Leading Protein A`,
        `Leading Protein B` = if (Collapse_B) replacement else `Leading Protein B`,
        `# Accessions A`    = if (Collapse_A) 1L          else `# Accessions A`,
        `# Accessions B`    = if (Collapse_B) 1L          else `# Accessions B`,
        
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
          removed <- setdiff(accnsA, replacement)
          paste(c(na.omit(`Disambiguated Accessions A`), removed), collapse = ";")
        } else {
          `Disambiguated Accessions A`
        }),
        `Disambiguated Accessions B` = as.character(if (Collapse_B) {
          removed <- setdiff(accnsB, replacement)
          paste(c(na.omit(`Disambiguated Accessions B`), removed), collapse = ";")
        } else {
          `Disambiguated Accessions B`
        })) %>%
      
      ungroup() %>%
      select(-accnsA, -accnsB, -posnsA, -posnsB,
             -Collapse_A, -Collapse_B,
             -Match_Index_A, -Match_Index_B,
             -.row)
    
    message("Processed homolog group: ", paste(protein.group, collapse = " / "))
  }
  
  # Summary
  n_disambig_A <- sum(crosslinks$`Disambiguated A`, na.rm = TRUE)
  n_disambig_B <- sum(crosslinks$`Disambiguated B`, na.rm = TRUE)
  message("Total disambiguated: ", n_disambig_A, " Protein A and ", n_disambig_B, " Protein B accessions.")
  
  return(crosslinks)
}
