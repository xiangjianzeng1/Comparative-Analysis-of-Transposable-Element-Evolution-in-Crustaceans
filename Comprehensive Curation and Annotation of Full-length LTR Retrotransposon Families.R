# Complete LTR analysis
# Data processing

# Load necessary R packages
library(dplyr)
library(readr)
library(rtracklayer)  # For importing GFF files
library(stringr)

# Define species list to process
species_list <- c("Aed", "Arg", "Dap", "Eri","Dro", "Hom", "Hya", "Lep", "Pen", "Pen.m", "Pen.c", "Por", "Pro")

# Initialize list to store dataframes for each species
species_data_frames <- list()

# Batch process data for each species and store results in list
for (species in species_list) {
  # Import and process GFF files for each species
  intactLTR_gff_file <- paste0("31.intactLTR/", species, "/intact.LTR.gff3")
  intactLTR_gff <- import(intactLTR_gff_file)
  intactLTR_gff <- as_tibble(intactLTR_gff) %>%
    select(seqnames, start, end, width, strand, ID, Name, ltr_identity, motif, tsd)
  
  # Read and process LTR_age data for each species
  intactLTR_age_file <- paste0("23.EDTA/", species, "/", species, ".genome.fa.mod.EDTA.raw/LTR/", species, ".genome.fa.mod.pass.list")
  intactLTR_age <- read_table(intactLTR_age_file) %>%
    select(LTR_loc = `#LTR_loc`, Insertion_Time) %>%
    distinct()
  
  # Process intactLTR and join with intactLTR_age
  intactLTR <- mutate(intactLTR_gff, 
                      LTR_loc = str_c(seqnames, ":", start, "..", end)) %>%
    inner_join(intactLTR_age, by = "LTR_loc")
  
  # Read classification information for each species
  intactLTR_cls_file <- paste0("31.intactLTR/", species, "/intact.LTR.fa.rexdb.cls.tsv")
  intactLTR_cls <- read_delim(intactLTR_cls_file, 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
  
  # Perform left join and select relevant columns
  intactLTR_final <- left_join(intactLTR, 
                               intactLTR_cls, 
                               by = c("ID" = "#TE")) %>%
    select(ID, Name, seqnames, start, end, width, ltr_identity, motif, tsd, Order, 
           Superfamily, Clade, Complete, Domains, Insertion_Time)
  
  # Store dataframe in list with species name
  species_data_name <- paste0(species, ".species_data")
  species_data_frames[[species_data_name]] <- intactLTR_final
}

# Combine all species data
all.species_data <- bind_rows(species_data_frames, .id = "Species")

# Filter data
all.species_data.filter <- filter(all.species_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown" ,"NA", "Tc1_Marine" ,"PiggyBac")))

# Select relevant columns
selected_data <- select(all.species_data.filter, Species, Superfamily, Insertion_Time)

# Modify species names
selected_data <- selected_data %>%
  mutate(Species = case_when(
    grepl("Aed", Species) ~ "(A.aeg)",
    grepl("Arg", Species) ~ "(A.bru)",
    grepl("Dap", Species) ~ "(D.pul)",
    grepl("Eri", Species) ~ "(E.sin)",
    grepl("Hom", Species) ~ "(H.ame)",
    grepl("Hya", Species) ~ "(H.azt)",
    grepl("Lep", Species) ~ "(L.sal)",
    grepl("Pen\\.c", Species) ~ "(P.chi)",
    grepl("Pen\\.m", Species) ~ "(P.mon)",
    grepl("Pen", Species) ~ "(P.van)",
    grepl("Por", Species) ~ "(P.tri)",
    grepl("Pro", Species) ~ "(P.cla)",
    grepl("Dro", Species) ~ "(D.mel)",
    TRUE ~ Species
  )) 

############################
# Visualization
library(ggplot2)
library(ggsci)
library(dplyr)

# Convert insertion time to million years
selected_data <- selected_data %>%
  mutate(Insertion_Time_mya = Insertion_Time / 1e6) %>%
  filter(!is.na(Superfamily) & !(Superfamily %in% c("NA", "mixture")))

# Set custom ordering
custom_species_order <- c("E.sin", "P.tri", "P.cla", "H.ame", "P.van", "P.chi", "P.mon", "H.azt", "L.sal", "D.pul", "D.mel", "A.aeg", "A.bru")
custom_supmerfamily_order <- c("mixture","Bel-Pao", "Copia",  "Gypsy")
selected_data$Species <- factor(selected_data$Species, levels = custom_species_order)
selected_data$Superfamily <- factor(selected_data$Superfamily, levels = custom_supmerfamily_order)

# Define color palette
violin_color <- c( "Gypsy" = "#D55740" , "Copia" = "#6CB9D2", "Bel-Pao" = "#479E88", "mixture"= "#f4a01a")

# Create violin plot
selected_data<- selected_data %>%
  filter(Insertion_Time_mya <= 3)

ggplot(selected_data, aes(x = Species, y = Insertion_Time_mya, fill = Superfamily)) +
  geom_violin(trim = T , size = 0.1) +
  scale_fill_manual(values = violin_color) +
  facet_grid(Superfamily ~ ., scales = "fixed") +
  labs(x = "Species", y = "Insertion Time (mya)") +
  scale_y_reverse(limits = c(3, 0), breaks = c(0, 1,2,3)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.1),
        legend.position = "top",
        panel.spacing = unit(0.3, "lines")) +
  guides(fill = guide_legend(title = NULL ,ncol = 6))

##############################
# Alternative visualization
# Calculate median insertion times
ordered_species <- selected_data %>%
  group_by(Species) %>%
  summarise(median_insertion_time = median(Insertion_Time_mya, na.rm = TRUE)) %>%
  arrange(median_insertion_time) %>%
  pull(Species)

# Order factors
selected_data$Species <- factor(selected_data$Species, levels = ordered_species)

# Create alternative violin plot
ggplot(selected_data, aes(x = Species, y = Insertion_Time_mya, fill = Superfamily)) +
  geom_violin(trim = TRUE) +
  scale_fill_npg() +
  facet_wrap(~ Superfamily, scales = "free_y") +
  labs(title = "Distribution of Insertion Times by Species and TE Superfamily",
       x = "Species",
       y = "Insertion Time (mya)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###############################
# Access data by species
species_data_frames[["Arg.species_data"]]  # Access data for Arg species

# Visualization parameters
superfamily_colors <- c(
  "Bel-Pao" = "#1D4787",  # Blue
  "Copia" = "#FF0000",    # Red
  "Gypsy" = "#69B657"    # Green
)

# Create density plots for each species
for (species_name in names(species_data_frames)) {
  species_data <- species_data_frames[[species_name]]
  species_filtered <- filter(species_data, !is.na(Superfamily) & !(Superfamily %in% c("unknown", "NA")))
  
  plot <- ggplot(species_filtered, aes(x = Insertion_Time, color = Superfamily)) +
    geom_density(aes(y = ..density..), adjust=1) +
    scale_x_continuous(name = 'Insertion Time (Mys)', limits = c(min_time, max_time)) +
    scale_color_manual(values = superfamily_colors) +
    theme_classic() +
    theme(legend.position = "right") +
    ggtitle(paste("Insertion Time Density for", species_name))
  
  print(plot)
}
