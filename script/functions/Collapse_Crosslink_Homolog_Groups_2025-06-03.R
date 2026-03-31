# Collapse manually defined homolog groups in crosslink data #


Collapse_Crosslink_Homolog_Groups <- function(crosslinks, homolog.list) {

# Create a column to record whether a homolog group has been collapsed
crosslinks <- crosslinks %>%
  mutate(
    `Homolog Disambiguated A`=F,
    `Homolog Disambiguated B`=F)
  
# Iterate through each protein homolog group
for (protein.group in homolog.list) {
  
  # Assign the first member of the group as the replacement value
  replacement <- protein.group[1]
  
  crosslinks <- crosslinks %>%
    rowwise() %>%
    # If the # accessions is >1 but all accessions are within the homolog group, we can disambiguate the crosslink by assigning to the leading group member:
    mutate(
      # Separate the accessions and crosslink positions 
      accnsA = list(strsplit(as.character(`Accessions A`), ";")[[1]]),
      accnsB = list(strsplit(as.character(`Accessions B`), ";")[[1]]),
      posnsA = list(strsplit(as.character(`Crosslink Positions A`), ";")[[1]]),
      posnsB = list(strsplit(as.character(`Crosslink Positions B`), ";")[[1]]),
      # Conditions to collapse: accession is ambiguous, all possible accessions are within the group, replacement value is among possible accessions
      Collapse_A = (`# Accessions A` > 1) && all(accnsA %in% protein.group) && (replacement %in% accnsA),
      Collapse_B = (`# Accessions B` > 1) && all(accnsB %in% protein.group) && (replacement %in% accnsB),
      # If conditions met, determine the position of the replacement value in the accessions and positions 
      Match_Index_A = if (Collapse_A) which(accnsA == replacement) else NA_integer_,
      Match_Index_B = if (Collapse_B) which(accnsB == replacement) else NA_integer_,
      # If conditions met, set leading protein and accession to the replacement value and make a record of this
      `Accessions A` = if (Collapse_A) replacement else `Accessions A`,
      `Accessions B` = if (Collapse_B) replacement else `Accessions B`,
      `Leading Protein A` = if (Collapse_A) replacement else `Leading Protein A`,
      `Leading Protein B` = if (Collapse_B) replacement else `Leading Protein B`,
      `# Accessions A` = if (Collapse_A) 1 else `# Accessions A`,
      `# Accessions B` = if (Collapse_B) 1 else `# Accessions B`,
      `Homolog Disambiguated A` = if (Collapse_A) TRUE else `Homolog Disambiguated A`,
      `Homolog Disambiguated B` = if (Collapse_B) TRUE else `Homolog Disambiguated B`,
      # If conditions met AND crosslink positions = # accessions (i.e. no multiple matches within same protein), set both to 1
      `# Crosslink Positions A` = if (Collapse_A && length(posnsA) == length(accnsA)) 1 else `# Crosslink Positions A`,
      `# Crosslink Positions B` = if (Collapse_B && length(posnsB) == length(accnsB)) 1 else `# Crosslink Positions B`,
      # If conditions met AND updated # crosslink positions = 1, pull the match index position as the crosslink position
      `Crosslink Positions A` = if (Collapse_A && `# Crosslink Positions A` == 1) posnsA[Match_Index_A] else `Crosslink Positions A`,
      `Crosslink Positions B` = if (Collapse_B && `# Crosslink Positions B` == 1) posnsB[Match_Index_B] else `Crosslink Positions B`
    ) %>%
    ungroup() %>%
    # Drop helper cols
    select(-accnsA, -accnsB, -posnsA, -posnsB, -Collapse_A, -Collapse_B, -Match_Index_A, -Match_Index_B)
  
  # Check progress
  print(protein.group)
  }
  return(crosslinks)
    
}
  
      
     
