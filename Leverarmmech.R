rm(list=ls()) 
setwd("/Users/raquelcardenas/Desktop/prog/allspecimens")

library(ggplot2)
library(dplyr)
library(tidyr)


files <- list.files("/Users/raquelcardenas/Desktop/prog/allspecimens", pattern = "\\.csv$", full.names = TRUE)

calculate_ma_grid <- function(file) {
  landmarks <- read.csv(file, skip = 1, header = TRUE)
  
  # Species code from filename (characters 5 and 6)
  species_code <- substr(basename(file), 5, 6)
  species_levels <- c("Eg", "Hs", "Ns", "Om", "Pt", "Dk")  
  species_labels <- c("Echinosorex gymnura", "Hylomys suillus", "Neotetracus sinensis", "Otohylomys megalotis", "Podogymnura truei", "Deinogalerix koenigswaldi")
  species <- factor(species_code, levels = species_levels, labels = species_labels)
  
  # TMJ midpoint
  tmj <- (landmarks[4, c("X", "Y", "Z")] + landmarks[5, c("X", "Y", "Z")]) / 2
  
  # Out-levers (bite points)
  out_levers <- sapply(list(
    Incisor = landmarks[13, c("X", "Y", "Z")],
    Molar = landmarks[17, c("X", "Y", "Z")]
  ), function(bite) sqrt(sum((tmj - bite)^2)))
  
  # In-levers (muscle attachment points)
  in_levers <- sapply(list(
    `Superficial masseter - Dorsal` = landmarks[8, c("X", "Y", "Z")],
    `Superficial masseter - Ventral` = landmarks[9, c("X", "Y", "Z")],
    `Deep masseter` = landmarks[18, c("X", "Y", "Z")],
    Temporalis = landmarks[2, c("X", "Y", "Z")]
  ), function(muscle) sqrt(sum((tmj - muscle)^2)))
  
  #  MAs calc
  ma_values <- c()
  for (muscle_name in names(in_levers)) {
    for (bite_name in names(out_levers)) {
      key <- paste("MA", muscle_name, bite_name, sep = "_")
      ma_values[key] <- in_levers[muscle_name] / out_levers[bite_name]
    }
  }
  
  # Returning as named vector with species
  return(c(ma_values, species = as.character(species)))
}



results_matrix <- t(sapply(files, function(f) {
  tryCatch(calculate_ma_grid(f), error = function(e) rep(NA, 9))  # 8 MA + 1 species
}))

results_df <- data.frame(results_matrix, stringsAsFactors = FALSE)
results_df$specimen <- basename(files)
results_df$species <- as.factor(results_df$species)

# Color palette
species_colors <- c(
  "Echinosorex gymnura" = "coral1", 
  "Hylomys suillus" = "goldenrod1",
  "Neotetracus sinensis" = "chartreuse1",
  "Otohylomys megalotis" = "blueviolet",
  "Podogymnura truei" = "deepskyblue",
  "Deinogalerix koenigswaldi" = "blue"
)

results_df$color <- species_colors[as.character(results_df$species)]

results_long <- results_df %>%
  pivot_longer(cols = starts_with("MA_"),
               names_to = "muscle_bite",
               values_to = "MA")

# Adding this to separate columsn
results_long <- results_long %>%
  separate(muscle_bite, into = c("MA_tag", "muscle", "bite"), sep = "_") %>%
  select(-MA_tag)  #


results_long$MA <- as.numeric(results_long$MA)

# New plot
library(ggplot2)
unique_muscles <- unique(results_long$muscle)
  for (m in unique_muscles) {
  p <- ggplot(filter(results_long, muscle == m), aes(x = bite, y = MA, fill = species)) +
    geom_boxplot(outlier.shape = 21,      
               outlier.fill = "white",   #
               outlier.colour = "black", 
               outlier.size = 2) +
    scale_fill_manual(values = species_colors) +
    labs(
      title = paste("Mechanical Advantage (MA) -", m),
      x = "Biting scenario",
      y = "Mechanical Advantage"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)  #
}





