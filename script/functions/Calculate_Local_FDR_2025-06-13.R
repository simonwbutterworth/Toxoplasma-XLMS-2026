library(tidyverse)

# Calculate local FDR based on identification scores
calculate_local_fdr <- function(df,
                                score_col = "Identification Score",
                                decoy_col = "Is Decoy",
                                bw = 5,
                                n = 1024) {
  # Extract scores
  scores <- df[[score_col]]
  is_decoy <- df[[decoy_col]]
  
  scores.target <- scores[!is_decoy]
  scores.decoy  <- scores[is_decoy]
  
  # Define the density range
  min_score <- min(scores, na.rm = TRUE)
  max_score <- max(scores, na.rm = TRUE)
  
  # Estimate densities
  density.target <- density(scores.target, from = min_score, to = max_score, n = n, bw = bw)
  density.decoy  <- density(scores.decoy,  from = min_score, to = max_score, n = n, bw = bw)
  
  # Mixing proportions
  pi.target <- length(scores.target) / nrow(df)
  pi.decoy  <- length(scores.decoy) / nrow(df)
  
  # Interpolation function
  interpolate.density <- function(x, density) {
    approx(density$x, density$y, xout = x, rule = 2)$y
  }
  
  # Compute interpolated densities
  density_target <- interpolate.density(scores, density.target)
  density_decoy  <- interpolate.density(scores, density.decoy)
  
  # Compute local FDR
  df$`Local FDR` <- (pi.decoy * density_decoy) /
    (pi.decoy * density_decoy + pi.target * density_target)
  
  return(df)
}
