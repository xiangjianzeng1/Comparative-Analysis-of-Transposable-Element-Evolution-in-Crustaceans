# Load required packages
library(tidyverse)
library(rtracklayer)
library(dplyr)

# Create empty dataframe to store genome sizes for each species
genome_sizes = data.frame(Species = character(), Genome_Size = numeric())

# Define species list
species_list1 <- c("Aedes_aegypti", "Drosophila_melanogaster", "Daphnia_pulex", "Tigriopus_copepod", 
                   "Lepeophtheirus_salmonis", "Hyalella_azteca", "Penaeus_chinensis", "Penaeus_monodon", 
                   "Penaeus_vannamei", "Portunus_trituberculatus", "Procambarus_clarkii", "Eriocheir_sinensis", 
                   "Homarus_americanus", "Argiope_bruennichi")

# Calculate genome size for each species
for (species in species_list1) {
  # Get genome size from corresponding fai file
  genome_size = sum(read.table(file = paste0("12.ref/", species, ".genome.fa.fai"))$V2)
  
  # Add current species and genome size to dataframe
  genome_sizes = rbind(genome_sizes, data.frame(Species = species, Genome_Size = genome_size))
}

# Replace underscores with spaces in species names
genome_sizes <- genome_sizes %>%
  mutate(Species = gsub("_", " ", Species))

# Create empty dataframe for TE annotation data
allTE = data.frame()

# Define species codes for path construction
species_list23 <- c("Aed", "Dro", "Dap", "Tig", 
                    "Lep", "Hya", "Pen.c", "Pen.m", 
                    "Pen", "Por", "Pro", "Eri", 
                    "Hom", "Arg")

# Process each species' TE annotations
for (species in species_list23) {
  # Build file path
  te_file_path = paste0("23.EDTA/", species, "/", species, ".genome.fa.mod.EDTA.TEanno.gff3")
  
  # Check file existence
  if (file.exists(te_file_path)) {
    # Read TE annotation file
    te_data = import(te_file_path) %>%
      as_tibble() %>%
      mutate(Species = species)
    
    # Add current species' TE data to main dataframe
    allTE = rbind(allTE, te_data)
  } else {
    message(paste("Warning: File not found for species", species))
  }
}

# Process merged data and add genome sizes
allTE_processed = allTE %>%
  mutate(Classification0 = str_split(Classification, "/", simplify = TRUE)[, 1]) %>%
  left_join(genome_sizes, by = "Species")

# Remove unnecessary columns
allTE_processed_filtered = allTE_processed %>%
  select(-phase, -TSD, -TIR, -motif, -tsd, -Parent)

# Plotting
library(tidyverse)
library(ggplot2)
library(dplyr)

# Calculate total TE size per species
te_size <- allTE_summ %>%
  group_by(Species) %>%
  summarise(TE_size = sum(size))

# Merge genome size and TE size data
genome_te_data <- left_join(te_size, genome_sizes, by = "Species")

genome_te_data_sorted_size <- genome_te_data[order(genome_te_data$Genome_Size), ]

# Create genome size vs TE content plot
library(ggrepel)
p1 <- ggplot(genome_te_data_sorted_size, aes(x = log2(Genome_Size), y = log2(TE_size))) +
  geom_point(color = "black", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "#8074AC", size = 1, fill = "#0055A5") +
  geom_text_repel(aes(label = Species), size = 3, box.padding = 0.35, point.padding = 0.5,
                  segment.color = 'grey50') +
  theme_classic() +
  labs(x = "Genome size (log2(bp))", y = "TE content (log2(bp))") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

p1

# Calculate correlation coefficient
log_genome_size <- log2(genome_te_data_sorted_size$Genome_Size)
log_te_size <- log2(genome_te_data_sorted_size$TE_size)

# Compute Pearson correlation
correlation_coefficient <- cor(log_genome_size, log_te_size, method = "pearson")

# Print results
print(correlation_coefficient)


# Plot superfamily classification
# Data processing
superfamily_data <- allTE_summ_Class %>%
  mutate(Classification = ifelse(str_detect(Classification, "/"), 
                               str_split(Classification, "/", simplify = TRUE)[, 2], 
                               Classification)) %>% # Keep content after "/", retain all if no "/"
  group_by(Species, Classification, Classification0) %>%
  summarise(size = sum(size, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(proportion = size / sum(size)) %>%
  ungroup() 

superfamily_data$Species <- factor(superfamily_data$Species, levels = custom_species_order)

library(RColorBrewer)
# Define color palette count
colourCount <- length(unique(superfamily_data$Classification))

# Create plot
p3 <- ggplot(superfamily_data, aes(x = Species, y = proportion, 
                                 fill = as.factor(Classification), 
                                 colour = Classification0)) +
  geom_bar(stat = "identity", position = "stack", size = 0.7, width = 0.9) +
  theme_bw() +
  labs(title = "Proportion of Transposable Elements (Superfamily) in Sixteen Species", 
       x = NULL, y = NULL) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colourCount)) +  # Custom fill colors for Classification
  scale_colour_manual(values = class0_colors, name = "Classification") +  # Border colors for Classification0
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),  # Adjust legend key size
        legend.text = element_text(size = 8)) +
  guides(fill = guide_legend(order = 1, title = NULL),  # Remove Classification legend title
         colour = guide_legend(order = 2, title = "Classification"))  # Adjust Classification0 legend

p3

# Correlation between TE family counts and genome size
# Add genome size information

# Remove specific species
genome_sizes <- genome_sizes %>%
  filter(Species != "Capitulum mitella")

allTE_summ_Class <- allTE_summ_Class %>%
  filter(Species != "Tigriopus californicus") %>%
  left_join(genome_sizes, by = "Species")

# Prepare count data
te_count_data <- allTE_summ_Class %>%
  select(Species, Classification, count, Genome_Size) %>%
  mutate(log_count = log10(count + 1), 
         log_genome_size = log10(Genome_Size))

# Calculate correlations for counts
correlation_results_count <- te_count_data %>%
  filter(!is.na(log_count) & !is.na(log_genome_size)) %>% # Remove missing values
  group_by(Classification) %>%
  summarize(
    correlation = if(n() >= 3) cor(log_count, log_genome_size, method = "pearson") else NA,
    p_value = if(n() >= 3) cor.test(log_count, log_genome_size, method = "pearson")$p.value else NA
  ) %>%
  ungroup() %>%
  filter(Classification != "TIR/Sola2",
         Classification != "pararetrovirus")

# Prepare size data
te_size_data <- allTE_summ_Class %>%
  select(Species, Classification, size, Genome_Size) %>%
  mutate(log_size = log10(size + 1), 
         log_genome_size = log10(Genome_Size))

# Calculate correlations for sizes
correlation_results_size <- te_size_data %>%
  filter(!is.na(log_size) & !is.na(log_genome_size)) %>%
  group_by(Classification) %>%
  summarize(
    correlation = if(n() >= 3) cor(log_size, log_genome_size, method = "pearson") else NA,
    p_value = if(n() >= 3) cor.test(log_size, log_genome_size, method = "pearson")$p.value else NA
  ) %>%
  ungroup()

# Create bubble plot
library(ggplot2)
library(ggrepel)

# Filter significant results
correlation_results_count <- correlation_results_count[correlation_results_count$correlation > 0.5, ]

# Plot bubble chart
ggplot(correlation_results_count, 
       aes(x = reorder(Classification, correlation), 
           y = correlation, 
           size = -log10(p_value), 
           color = correlation)) +
  geom_text_repel(aes(label = sprintf("r = %.2f\np = %.2e", correlation, p_value)),
                  size = 2.5, 
                  box.padding = 0.35, 
                  point.padding = 0.5,
                  fontface = "bold") +
  geom_point(alpha = 0.5) +
  scale_size_continuous(range = c(3, 15)) +
  scale_color_gradient2(midpoint = 0.5, 
                       low = "green", 
                       mid = "blue", 
                       high = "red") +
  scale_y_continuous(limits = c(0.5, 1)) +
  labs(x = NULL,
       size = "-log10(p-value)",
       color = "Correlation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 6),
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_line(color = "gray", linetype = "dashed"))






# Plot correlation between TE count and genome size

TE_count <- te_count_data %>%
  group_by(Species) %>%
  summarise(total_count = sum(count)) %>%
  left_join(genome_sizes, by = "Species")

ggplot(TE_count, aes(x = log10(Genome_Size), y = log10(total_count))) +
  geom_point(color = "black", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "#8074AC", size = 1, fill = "#0055A5") +
  geom_text_repel(aes(label = Species), size = 3, box.padding = 0.2, point.padding = 0.1,
                  segment.color = 'grey50') +
  theme_classic() +
  labs(x = "genome size (log10(bp))", y = "TE content (log10(count))") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
p4

# Data preparation
log_genome_size1 <- log10(TE_count$Genome_Size)
log_te_size1 <- log10(TE_count$total_count)

correlation_test1 <- cor.test(log_genome_size1, log_te_size1, method = "pearson")

# Output correlation coefficient and P-value
correlation_coefficient1 <- correlation_test1$estimate
p_value1 <- correlation_test1$p.value
print(correlation_coefficient1) # 0.9052379
print(p_value1) # 2.083162e-05

######## TE size analysis
log_genome_size1 <- log10(TE_size$Genome_Size)
log_te_size1 <- log10(TE_size$total_size)

correlation_test1 <- cor.test(log_genome_size1, log_te_size1, method = "pearson")

# Output correlation results
correlation_coefficient1 <- correlation_test1$estimate
p_value1 <- correlation_test1$p.value
print(correlation_coefficient1) # 0.9052379
print(p_value1) # 2.083162e-05

###########################################################
# Plot TE family size correlation with genome size
te_size_data <- allTE_summ_Class %>%
  select(Species, Classification, size, Genome_Size) %>%
  mutate(log_size = log10(size + 1), log_genome_size = log10(Genome_Size))

TE_size <- te_size_data %>%
  group_by(Species) %>%
  summarise(total_size = sum(size)) %>%
  left_join(genome_sizes, by = "Species")

ggplot(TE_size, aes(x = log10(Genome_Size), y = log10(total_size))) +
  geom_point(color = "black", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "#8074AC", size = 1, fill = "#0055A5") +
  geom_text_repel(aes(label = Species), size = 3, box.padding = 0.2, point.padding = 0.1,
                  segment.color = 'grey50') +
  theme_classic() +
  labs(x = "genome size (log10(bp))", y = "TE size (log10(bp))") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

te_size_data <- allTE_summ_Class %>%
  select(Species, Classification, size, Genome_Size) %>%
  mutate(log_size = log10(size + 1), log_genome_size = log10(Genome_Size)) %>%
  select(Species, Classification, log_size, log_genome_size)

# Calculate correlations for size data
correlation_results_size <- te_size_data %>%
  filter(!is.na(log_size) & !is.na(log_genome_size)) %>% # Remove missing values
  group_by(Classification) %>%
  summarize(
    correlation = if(n() >= 3) cor(log_size, log_genome_size, method = "pearson") else NA,
    p_value = if(n() >= 3) cor.test(log_size, log_genome_size, method = "pearson")$p.value else NA
  ) %>%
  ungroup() %>%
  filter(Classification != "TIR/Sola2") %>%
  filter(Classification != "pararetrovirus")

# Print results
print(correlation_results_size)

# Create bubble plot
library(ggplot2)
library(ggrepel)
correlation_results_size <- correlation_results_size[correlation_results_size$correlation > 0.5, ]

# Create bubble plot with labels
ggplot(correlation_results_size, aes(x = reorder(Classification, correlation), y = correlation, size = -log10(p_value), color = correlation)) +
  geom_text_repel(aes(label = sprintf("r = %.2f\np = %.2e", correlation, p_value)), size = 2.5, box.padding = 0.1, point.padding = 0.1,
                  fontface = "bold") +
  geom_point(alpha = 0.5) +
  scale_size_continuous(range = c(3, 15), guide = guide_legend(override.aes = list(alpha = 0.5))) +
  scale_color_gradient2(midpoint = 0.5, low = "green", mid = "blue", high = "red", guide = guide_colorbar(barwidth = 1, barheight = 2.5)) +
  scale_y_continuous(limits = c(0.5, 1)) +
  labs(x = NULL,
       size = "-log10(p-value)",
       color = "Correlation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(0.5, "cm"), # Adjust legend key size
        legend.text = element_text(size = 6), # Adjust legend text size
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_line(color = "gray", linetype = "dashed"))

########## TE count analysis
# Create bubble plot for count data
library(ggplot2)
library(ggrepel)

correlation_results_count <- correlation_results_count[correlation_results_count$correlation > 0.5, ] 

# Create bubble plot with labels
ggplot(correlation_results_count, aes(x = reorder(Classification, correlation), y = correlation, size = -log10(p_value), color = correlation)) +
  geom_text_repel(aes(label = sprintf("r = %.2f\np = %.2e", correlation, p_value)), size = 2.5, box.padding = 0.1, point.padding = 0.1,
                  fontface = "bold") +
  geom_point(alpha = 0.5) +
  scale_size_continuous(range = c(3, 15), guide = guide_legend(override.aes = list(alpha = 0.5))) +
  scale_color_gradient2(midpoint = 0.5, low = "green", mid = "blue", high = "red", guide = guide_colorbar(barwidth = 1, barheight = 2.5)) +
  scale_y_continuous(limits = c(0.5, 1)) +
  labs(x = NULL,
       size = "-log10(p-value)",
       color = "Correlation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(0.5, "cm"), # Adjust legend key size
        legend.text = element_text(size = 6), # Adjust legend text size
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_line(color = "gray", linetype = "dashed"))
