library(tidyverse)

# Identifies optimal score and delta score thresholds to achieve a target FDR.
# Identify_AND_Gate_Score_Cutoffs: CSM must meet BOTH score and delta score thresholds.
# Identify_OR_Gate_Score_Cutoffs:  CSM must meet EITHER score or delta score threshold.
# Both are wrappers around the shared inner function Identify_Score_Cutoffs().

Identify_Score_Cutoffs <- function(
    df,
    gate = c("AND", "OR"),
    fdr.level = 0.01,
    min.score = 20,
    max.score = 100,
    min.delta.score = 10,
    max.delta.score = 100,
    interval = 0.1) {
  
  gate <- match.arg(gate)
  
  # Check for required columns
  required_cols <- c("Is Decoy", "Identification Score", "Identification DeltaScore")
  missing_cols  <- setdiff(required_cols, colnames(df))
  if (length(missing_cols)) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  # Extract vectors once for performance
  score_vals <- df$`Identification Score`
  delta_vals <- df$`Identification DeltaScore`
  is_decoy   <- df$`Is Decoy`
  
  # Build threshold grid
  scores <- seq(min.score,       max.score,       by = interval)
  deltas <- seq(min.delta.score, max.delta.score, by = interval)
  grid   <- expand.grid(score = scores, delta.score = deltas)
  
  # Vectorised FDR calculation over grid
  message("Computing FDR over ", nrow(grid), " threshold combinations (gate = ", gate, ")...")
  
  if (gate == "AND") {
    counts <- mapply(function(s, d) {
      sel    <- score_vals >= s & delta_vals >= d
      total  <- sum(sel)
      decoy  <- sum(sel & is_decoy)
      c(total = total, decoy = decoy)
    }, grid$score, grid$delta.score)
  } else {
    counts <- mapply(function(s, d) {
      sel    <- score_vals >= s | delta_vals >= d
      total  <- sum(sel)
      decoy  <- sum(sel & is_decoy)
      c(total = total, decoy = decoy)
    }, grid$score, grid$delta.score)
  }
  
  grid$total  <- counts["total", ]
  grid$decoy  <- counts["decoy", ]
  grid$target <- grid$total - grid$decoy
  grid$fdr    <- ifelse(grid$total > 0, grid$decoy / grid$total, 0)
  
  message("Done.")
  
  # Filter to target FDR and sort by most targets retained, then lowest thresholds
  result <- grid %>%
    filter(fdr <= fdr.level) %>%
    arrange(desc(target), score, delta.score)
  
  if (nrow(result) == 0) {
    warning(
      "No threshold combination achieved FDR <= ", fdr.level, " (gate = ", gate, "). ",
      "Consider relaxing min.score or min.delta.score."
    )
  }
  
  return(result)
}
