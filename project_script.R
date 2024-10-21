# ==========================================
# 1. LOAD LIBRARIES
# ==========================================

library("phyloseq")         # For microbiome data analysis
library("cowplot")          # For combining plots
library("vegan")            # For rarefaction
library("ggplot2")          # For visualization
library("reshape2")         # For reshaping data between wide and long formats
library("dplyr")            # For data manipulation
library("tidyr")            # For tidying data
library("ggrepel")          # For adding labels to ggplot2 plots
library("MicrobiomeStat")   # For statistical analysis of microbiome data
library("writexl")          # For writing data frames to Excel files


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


# ==========================================
# 7. ALPHA DIVERSITY
# ==========================================

#Variable with the alpha diversity indexes to generate
alpha_ind <- c("Observed", "Shannon", "Fisher", "Simpson", "Chao1")

#Generate Alpha diversity
alpha_tab <- estimate_richness(phyloseq_obj, 
                               measures = alpha_ind)

# Variable with the column names to add from the metadata
cols_to_add_alpha <- c("habitat_class", "Run", "host_color")

# Turn metadata into a data frame
metadata_df <- data.frame(metadata)

# Combine alpha diversity results with metadata
alpha_tab_comb <- bind_cols(alpha_tab, metadata_df[cols_to_add_alpha])
head(alpha_tab_comb)


# Loop over each alpha diversity index BY COLOR
for (alpha in alpha_ind) {
  
  # Perform the Kruskal-Wallis test for the current alpha diversity index
  kruskal_test_results <- kruskal.test(reformulate("host_color", response = alpha), data = alpha_tab_comb)
  
  # Print the results
  cat("\nKruskal-Wallis test for", alpha, "grouped by host_color:\n")
  print(kruskal_test_results)
}


# Loop over each alpha diversity index BY HABITAT
for (alpha in alpha_ind) {
  
  # Perform the Kruskal-Wallis test for the current alpha diversity index
  kruskal_test_results <- kruskal.test(reformulate("habitat_class", response = alpha), data = alpha_tab_comb)
  
  # Print the results
  cat("\nKruskal-Wallis test for", alpha, "grouped by habitat class:\n")
  print(kruskal_test_results)
}


#=========================================
# 8. VISUALZATION OF ALPHA DIVERSITY
#=========================================

# Pivot data to use in plots
alpha_tab_comb_long <- pivot_longer(data = alpha_tab_comb, 
                                    cols = alpha_ind, 
                                    names_to = "Index", 
                                    values_to = "Value")

# Violin plot for all alpha diversities grouped by color
alpha_diversity_color_plot <- ggplot(data = alpha_tab_comb_long, 
                                     aes(x = host_color, 
                                         y = Value,
                                         fill = host_color)) + 
  geom_jitter(height = 0, 
              width = 0.25, 
              alpha = 0.75, 
              aes(color = host_color)) +
  geom_violin(trim = FALSE, 
              alpha = 0.5, 
              aes(color = host_color)) +
  facet_wrap("Index", scales = "free") + 
  theme_light() +
  labs(x = "Host Color", 
       title = "Alpha Diversity by Host Color") +
  theme(strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 9),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, 
                                  face = "bold"),
        axis.title = element_text(face = "bold"))

# Display plot
alpha_diversity_color_plot

# Save plot
ggsave("alpha_diversity_color_plot.png",
       plot = alpha_diversity_color_plot,
       width = 12,
       height = 8,
       units = "in",
       dpi = 300)

# Violin plot for all alpha diversities grouped by habitat
alpha_diversity_habitat_plot <- ggplot(data = alpha_tab_comb_long, 
                                       aes(x = habitat_class, 
                                           y = Value, 
                                           fill = habitat_class)) + 
  geom_jitter(height = 0, 
              width = 0.25, 
              alpha = 0.75, 
              aes(color = habitat_class)) +
  geom_violin(trim = FALSE, 
              alpha = 0.5, 
              aes(color = habitat_class)) +
  facet_wrap("Index", scales = "free") + 
  theme_light() +
  labs(x = "Habitat Class", 
       title = "Alpha Diversity by Habitat Class") +
  theme(strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 9),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, 
                                  face = "bold"),
        axis.title = element_text(face = "bold"))

# Display plot
alpha_diversity_habitat_plot

# Save plot
ggsave("alpha_diversity_habitat_plot.png",
       plot = alpha_diversity_habitat_plot,
       width = 12,
       height = 8,
       units = "in",
       dpi = 300)


#=========================================
# 9.BETA DIVERSITY
#=========================================

# Subset the phyloseq object to include only Black squirrels
black_squirrels <- subset_samples(phyloseq_obj, host_color == "Black")

# Subset the phyloseq object to include only Grey squirrels
grey_squirrels <- subset_samples(phyloseq_obj, host_color == "Grey")

# Subset the phyloseq object to include only Intermediate squirrels
intermediate_squirrels <- subset_samples(phyloseq_obj, host_color == "Intermediate")


# Distance matrix using Jaccard method
jacc_dist_black <- distance(black_squirrels, method = "jaccard", binary = TRUE)
jacc_dist_grey <- distance(grey_squirrels, method = "jaccard", binary = TRUE)
jacc_dist_intermediate <- distance(intermediate_squirrels, method = "jaccard", binary = TRUE)

# Distance matrix using Bray-Curtis method
bray_dist_black <- distance(black_squirrels, method = "bray")
bray_dist_grey <- distance(grey_squirrels, method = "bray")
bray_dist_intermediate <- distance(intermediate_squirrels, method = "bray")


# Ordination using Jaccard method
jacc_ord_black <- ordinate(black_squirrels, method = "MDS", distance = jacc_dist_black)
jacc_ord_grey <- ordinate(grey_squirrels, method = "MDS", distance = jacc_dist_grey)
jacc_ord_intermediate <- ordinate(intermediate_squirrels, method = "MDS", distance = jacc_dist_intermediate)

# Ordination using Bray-Curtis method
bray_ord_black <- ordinate(black_squirrels, method = "MDS", distance = bray_dist_black)
bray_ord_grey <- ordinate(grey_squirrels, method = "MDS", distance = bray_dist_grey)
bray_ord_intermediate <- ordinate(intermediate_squirrels, method = "MDS", distance = bray_dist_intermediate)


# Jaccard plot for Black squirrels
jacc_plot_black <- plot_ordination(black_squirrels, jacc_ord_black, color = "habitat_class") +
  theme_classic() +
  stat_ellipse(alpha = 0.25) +
  ggtitle("Jaccard - Black Squirrels") +
  labs(color = "Habitat Class") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Jaccard plot for Grey squirrels
jacc_plot_grey <- plot_ordination(grey_squirrels, jacc_ord_grey, color = "habitat_class") +
  theme_classic() +
  stat_ellipse(alpha = 0.25) +
  ggtitle("Jaccard - Grey Squirrels") +
  labs(color = "Habitat Class") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Jaccard plot for Intermediate squirrels
jacc_plot_intermediate <- plot_ordination(intermediate_squirrels, jacc_ord_intermediate, color = "habitat_class") +
  theme_classic() +
  stat_ellipse(alpha = 0.25) +
  ggtitle("Jaccard - Intermediate Squirrels") +
  labs(color = "Habitat Class") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Bray-Curtis plot for Black squirrels
bray_plot_black <- plot_ordination(black_squirrels, bray_ord_black, color = "habitat_class") +
  theme_classic() +
  stat_ellipse(alpha = 0.25) +
  ggtitle("Bray-Curtis - Black Squirrels") +
  labs(color = "Habitat Class") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Bray-Curtis plot for Grey squirrels
bray_plot_grey <- plot_ordination(grey_squirrels, bray_ord_grey, color = "habitat_class") +
  theme_classic() +
  stat_ellipse(alpha = 0.25) +
  ggtitle("Bray-Curtis - Grey Squirrels") +
  labs(color = "Habitat Class") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Bray-Curtis plot for Intermediate squirrels
bray_plot_intermediate <- plot_ordination(intermediate_squirrels, bray_ord_intermediate, color = "habitat_class") +
  theme_classic() +
  stat_ellipse(alpha = 0.25) +
  ggtitle("Bray-Curtis - Intermediate Squirrels") +
  labs(color = "Habitat Class") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )


# Combine Jaccard plots
combined_jacc_plot <- plot_grid(jacc_plot_black, jacc_plot_grey, jacc_plot_intermediate, 
                                ncol = 1, align = "v", labels = c("A", "B", "C"))

# Combine Bray-Curtis plots
combined_bray_plot <- plot_grid(bray_plot_black, bray_plot_grey, bray_plot_intermediate, 
                                ncol = 1, align = "v", labels = c("A", "B", "C"))


# Save combined Jaccard plot
ggsave(filename = "combined_jacc_plot.png", plot = combined_jacc_plot, width = 10, height = 12, dpi = 300)

# Save combined Bray-Curtis plot
ggsave(filename = "combined_bray_plot.png", plot = combined_bray_plot, width = 10, height = 12, dpi = 300)


#=========================================
# 10.DIFFERENTIAL ABUNDANCE
#=========================================

# Transpose the OTU table if taxa are not rows
t_phyloseq_obj <- phyloseq::t(phyloseq_obj)

# Convert phyloseq object to mstat objects
mstat_obj <- mStat_convert_phyloseq_to_data_obj(t_phyloseq_obj)

# Check if host_color and habitat_class exist
print(colnames(mstat_obj$meta.dat))

# Create interaction variable in meta.dat
mstat_obj$meta.dat$host_habitat <- interaction(mstat_obj$meta.dat$host_color, 
                                               mstat_obj$meta.dat$habitat_class)

# Perform differential abundance analysis
test_results <- generate_taxa_test_single(
  data.obj = mstat_obj,
  group.var = "host_habitat",
  adj.vars = "habitat_class",
  feature.dat.type = "count",
  feature.level = c("Genus"),
  prev.filter = 0.1,
  abund.filter = 0.0001
)

# Clean and format the test results for visualization
formatted_results <- lapply(test_results$Genus, function(x) {
  x[c("Coefficient", "SE", "P.Value", "Adjusted.P.Value", "Mean.Abundance", "Prevalence")] <- 
    round(x[c("Coefficient", "SE", "P.Value", "Adjusted.P.Value", "Mean.Abundance", "Prevalence")], 3)
  
  colnames(x) <- c("Variable", "Coefficient", "Standard Error", "P-Value", "Adjusted P-Value", "Mean Abundance", "Prevalence")
  return(x)
})

# Combine into a data frame
final_results <- do.call(rbind, formatted_results)

# Ensure P.Value is numeric
final_results$P.Value <- as.numeric(as.character(final_results$`P-Value`))

# Calculate log2 fold change and -log10(p-value)
final_results$Log2FC <- log2(final_results$Coefficient)
final_results$NegLog10P <- -log10(final_results$P.Value)

# Determine significance
final_results$Significant <- final_results$`Adjusted P-Value` < 0.05


# Save the final results to a Excel file (for better readability)
# write_xlsx(final_results, "differential_abundance_results.xlsx")


# Define the y-axis limit
y_limit <- 5

# Volcano plot with labels only for significant genera
volcano_plot <- ggplot(final_results, aes(x = Log2FC, y = NegLog10P, color = `Mean Abundance`)) +
  geom_point(aes(shape = Significant), alpha = 0.7, size = 4) +
  geom_label_repel(data = final_results_filtered %>% filter(Significant == TRUE), 
                   aes(label = Variable), 
                   size = 3, 
                   show.legend = FALSE,
                   segment.color = 'grey50',
                   nudge_x = 0.1,
                   direction = 'both',
                   box.padding = 0.5,
                   max.overlaps = Inf) +  # Allow unlimited overlaps
  scale_color_gradientn(colors = c("#36004d", "#27517b", "#20727b", "#209473", "#fce51e"), 
                        name = "Mean Abundance") +  # Gradient scale for abundance
  labs(title = "Volcano Plot of Differential Abundance",
       x = "Log2 Fold Change",
       y = "-Log10 P-Value") +
  coord_cartesian(ylim = c(0, y_limit)) +
  theme_minimal() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),  # Make the title bold
        legend.position = "top")

# Print the volcano plot
print(volcano_plot)

# Save the volcano plot
ggsave("volcano_plot_differential_abundance.png", 
       plot = volcano_plot, 
       width = 10, 
       height = 6, 
       dpi = 300)
