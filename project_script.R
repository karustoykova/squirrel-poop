# ==========================================
# 1. LOAD LIBRARIES
# ==========================================

library("phyloseq")   # For microbiome data analysis
library("cowplot")    # For combining plots


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

# WIP
# - count normalization, summarization at diff taxonomic levels, etc.

