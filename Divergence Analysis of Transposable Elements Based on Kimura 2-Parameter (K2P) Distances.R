# Data Organization
# Species name vector
species_names <- c("Aed","Arg", "Cap", "Dap", "Dro", "Eri", "Hom", "Hya", "Lep", "Pen", "Pen.m", "Pen.c", "Por", "Pro", "Tac", "Tig")

# Process data for each species using a loop
for (species in species_names) {
  # Dynamically construct variable name for raw data
  data_var_name <- paste0(species, ".TE_divsum")
  
  # Dynamically construct variable name for processed data
  summarized_var_name <- paste0(species, "_summarized")
  
  # Retrieve corresponding dataframe from global environment
  species_data <- get(data_var_name, envir = .GlobalEnv)
  
  # Modify classification: unify categories containing "SINE" to "SINE"
  species_data$Classification <- ifelse(grepl("SINE", species_data$Classification), "SINE", species_data$Classification)
  
  # Group by Div and Classification, then sum percent values
  summarized_data <- species_data %>% 
    group_by(Div, Classification) %>%
    summarise(percent = sum(percent), .groups = 'drop')  %>%
    filter(Classification != "Unknown")
  
  # Store processed dataframe back to global environment
  assign(summarized_var_name, summarized_data, envir = .GlobalEnv)
}

library(ggsci)
library(scales)

# Create color mapping
categories <- c("DNA.DTA", "DNA.DTC", "DNA.DTH", "DNA.DTT","DNA.DTM",
                "DNA.Helitron", "LINE.unknown", "LTR.Copia", "LTR.Gypsy",
                "LTR.unknown", "LTR.unknown_BEL", "MITE.DTA", "MITE.DTC",
                "MITE.DTH", "MITE.DTM", "MITE.DTT", "Penelope", "polinton",
                "SINE", "TIR.Kolobok", "TIR.PiggyBac", "TIR.Sola2", "TIR.Tc1_Mariner",
                "Unknown")

igv_colors_24 <- pal_igv("default")(24)
color_map <- setNames(igv_colors_24, categories)

# Visualization template for each species (example with Aed)
Aed_TE.DIV <- ggplot(data = Aed_summarized, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence (%)", 
       y = "TE content (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right") +
  guides(fill = guide_legend(nrow = 12,  # Set number of rows
                             ncol = 2))  # Set number of columns

# Example modification for unknown elements plot
ggplot(data = Aed_unknown_data, aes(x = Div, y = percent)) +
  geom_col(aes(fill = Classification), width = 0.8, color = "white") +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_fill_lancet() +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura 2‐Parameter divergence (%)", 
       y = "TE content (%)", fill = NULL) + 
  theme_classic() +
  theme(legend.position = "right") +
  guides(fill = guide_legend(nrow = 12,  # Set number of rows
                             ncol = 2))  # Set number of columns
