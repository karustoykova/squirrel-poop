# ==========================================
# 1. LOAD LIBRARIES
# ==========================================

library("phyloseq")   # For microbiome data analysis
library("cowplot")    # For combining plots
library("vegan")      # For rarefaction
library("ggplot2")
library("reshape2")


# ==========================================
# 2. LOAD FILES AND EXTRACT DATA
# ==========================================

# Load phyloseq object
phyloseq_obj <- readRDS("data/phyloseq_object.rds")

# Remove control sample
phyloseq_obj <- subset_samples(phyloseq_obj, sample_names(phyloseq_obj) != "squirrelsample_Rural_110")

# Extract key components
otu_data <- otu_table(phyloseq_obj)
head(otu_data)

taxonomy_data <- tax_table(phyloseq_obj)
head(taxonomy_data)

metadata <- sample_data(phyloseq_obj)
head(metadata)

# Summarize abundance of each taxon across all samples
taxa_sums <- taxa_sums(phyloseq_obj)
print(taxa_sums)

# Change ASV names to a sequential number (e.g., ASV_1, ASV_2, etc.)
asv_names <- paste0("ASV_", 1:ntaxa(phyloseq_obj))
tail(asv_names)

# Replace ASV names in the phyloseq object with the new sequential names
taxa_names(phyloseq_obj) <- asv_names
tail(taxa_names(phyloseq_obj))


# Sample names, based on ID number
meta <- as.data.frame(sample_data(phyloseq_obj))
new_sample_names <- paste0("squirrel_sampleID", seq_len(nrow(meta)))
sample_names(phyloseq_obj) <- new_sample_names


# ==========================================
# 3. COUNT NORMALIZATION
# ==========================================

# Normalize read depth across all samples using rarefaction
phyloseq_obj_rare <- rarefy_even_depth(phyloseq_obj, rngseed = 4865)

# Check number of reads after rarefaction
cat("Reads per sample after rarefaction:\n")
head(sample_sums(phyloseq_obj_rare))


# Transform counts into proportions (relative abundance)
phyloseq_obj_prop <- transform_sample_counts(phyloseq_obj_rare, function(x){ x / sum(x) })

# Verify relative abundance transformation
cat("Total reads after relative abundance transformation:\n")
head(sample_sums(phyloseq_obj_prop))


# ==========================================
# 4. SUMMARIZE RESULTS TO HIGHTER TAX LEVEL
# ==========================================

# Summarize to genus level
phyloseq_obj_genus <- tax_glom(phyloseq_obj, taxrank = "Genus" )
phyloseq_obj_rare_genus <- tax_glom(phyloseq_obj_rare, taxrank = "Genus")
phyloseq_obj_prop_genus <- tax_glom(phyloseq_obj_prop, taxrank = "Genus")

# Summarize to family level
phyloseq_obj_rare_family <- tax_glom(phyloseq_obj_rare, taxrank = "Family")
phyloseq_obj_prop_family <- tax_glom(phyloseq_obj_prop, taxrank = "Family")
phyloseq_obj_rare_family <- tax_glom(phyloseq_obj_rare, taxrank = "Family")

# Summarize to phylum level
phyloseq_obj_rare_phylum <- tax_glom(phyloseq_obj_rare, taxrank = "Phylum")
phyloseq_obj_prop_phylum <- tax_glom(phyloseq_obj_prop, taxrank = "Phylum")
phyloseq_obj_rare_phylum <- tax_glom(phyloseq_obj_rare, taxrank = "Phylum")


# ==========================================
# 5. SUBSET DATA
# ==========================================

# Grey squirrels: Campus, Peri-Urban, Rural, Suburb
grey_campus <- subset_samples(phyloseq_obj, host_color == "Grey" & habitat_class == "Campus")
grey_peri_urban <- subset_samples(phyloseq_obj, host_color == "Grey" & habitat_class == "Peri-Urban")
grey_rural <- subset_samples(phyloseq_obj, host_color == "Grey" & habitat_class == "Rural")
grey_suburb <- subset_samples(phyloseq_obj, host_color == "Grey" & habitat_class == "Suburb")

# Intermediate squirrels: Campus, Peri-Urban, Rural, Suburb
intermediate_campus <- subset_samples(phyloseq_obj, host_color == "Intermediate" & habitat_class == "Campus")
intermediate_peri_urban <- subset_samples(phyloseq_obj, host_color == "Intermediate" & habitat_class == "Peri-Urban")
intermediate_rural <- subset_samples(phyloseq_obj, host_color == "Intermediate" & habitat_class == "Rural")
intermediate_suburb <- subset_samples(phyloseq_obj, host_color == "Intermediate" & habitat_class == "Suburb")

# Black squirrels: Campus, Peri-Urban, Rural, Suburb
black_campus <- subset_samples(phyloseq_obj, host_color == "Black" & habitat_class == "Campus")
black_peri_urban <- subset_samples(phyloseq_obj, host_color == "Black" & habitat_class == "Peri-Urban")
black_rural <- subset_samples(phyloseq_obj, host_color == "Black" & habitat_class == "Rural")
black_suburb <- subset_samples(phyloseq_obj, host_color == "Black" & habitat_class == "Suburb")


# ==========================================
# 6. HEATMAP OF TOP 10 GENERA ABUNDANCE
# ==========================================

# Abundance matrix, calculate total abundance for each genus
abundance_matrix <- as.data.frame(otu_table(phyloseq_obj_genus))

# Assign genus names to columns
genus_names <- as.character(tax_table(phyloseq_obj_genus)[, "Genus"])
colnames(abundance_matrix) <- genus_names

# Calculate total abundance for each genus
total_abundance <- colSums(abundance_matrix)

# Identify top 10 genera
top_10_genera <- names(sort(total_abundance, decreasing = TRUE))[1:10]
print(top_10_genera)

# Filter to keep only the top 10 genera in abundance
abundance_top10 <- abundance_matrix[, colnames(abundance_matrix) %in% top_10_genera]

# Add sample names to the abundance data
abundance_top10$Sample <- rownames(abundance_top10)

# Melt abundance data for ggplot
abundance_long <- psmelt(phyloseq_obj_prop_genus)
abundance_long <- abundance_long[abundance_long$Genus %in% top_10_genera, ]

# Plot heatmap for top 10 genera
heatmap_top10 <- ggplot(abundance_long, aes(x = Sample, y = Genus, fill = Abundance)) + 
  geom_tile() + 
  scale_fill_gradientn(colors = c("#36004d", "#27517b", "#20727b", "#209473", "#fce51e"),
                       name = "Relative Abundance") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) + 
  ggtitle("Relative Abundance of Top 10 Genera by Host Color") 

# Display heatmap
heatmap_top10

# Save heatmap
ggsave(
  filename = "heatmap_top10_genera.png", 
  plot = heatmap_top10, 
  width = 10, 
  height = 8, 
  dpi = 300
)
