#####################################################################################
########################__ORDINATION_&_RELATIVE_ABUNDANCE__##########################
#####################################################################################

library(phyloseq)
library(dplyr)
library(ggplot2)
library(vegan)
library(ape)
library(ggsci)
library(microViz)
library(tidyr)
library(speedyseq)

# Load data from nf-core ampliseq
physeq <- readRDS("dada2_phyloseq.rds")

# Now we read in the phylo tree computed by qiime2
tree <- read.tree("tree.nwk")

# And add it to our phyloseq object
dada2_phyloseq1 <- merge_phyloseq(physeq, phy_tree(tree))

#####################################################################################
############################____ALPHA__DIVERSITY____#################################

# Define the character strings and columns for filtering
filter_column1 <- "Timepoint"
filter_string1 <- "Baseline"

filter_column2 <- "Treatment"
filter_string2 <- "PAN"

# Filter the phyloseq object
filtered_physeq <- subset_samples(dada2_phyloseq1, 
                                    get(filter_column1) == filter_string1 & 
                                      get(filter_column2) == filter_string2)

# Print the filtered sample data
print(sample_data(filtered_physeq))

unfilt_data <- prune_species(speciesSums(filtered_physeq) > 0, filtered_physeq)

# Calculate richness measures
richness_data <- estimate_richness(unfilt_data)

# Now make the plots using this data in GraphPad Prism
write.csv(richness_data, file = "PAN_baseline_alpha_diversity.csv", row.names = FALSE)

#####################################################################################
#########################_____BETA__DIVERSITY______##################################

# Remove OTUs that appear less than 5 times in all samples 
da = dada2_phyloseq1
wh0 = genefilter_sample(da, filterfun_sample(function(x) x > 5), A=0.5*nsamples(da))
da1 = prune_taxa(wh0, da)

# Transform to even sampling depth
da1 = transform_sample_counts(da1, function(x) 1E6 * x/sum(x))

# Keep 5 most abundant taxa in all samples
phylum.sum = tapply(taxa_sums(da1), tax_table(da1)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
da1 = prune_taxa((tax_table(da1)[, "Phylum"] %in% top5phyla), da1)

# See what your metadata options are
metadata_vars <- colnames(sample_data(da1))
print(metadata_vars)

##############################################################################################
##############___________SHAM__and__CCI__at__baseline__VEHICLE__VS__PAN__________#############

# This section was repeated for each timepoint withe the independent variable being the probiotics

# Adjust "Timepoint" and "baseline" to match your actual column and value names
da1_baseline <- subset_samples(da1, Timepoint == "Baseline")

#Compute weighted unifrac distances
wunifrac_dist <- distance(da1_baseline, method = "wunifrac")

# Plot PCoA
pcoa_res <- ordinate(da1_baseline, method = "PCoA", distance = "wunifrac")
p1 <- plot_ordination(da1_baseline, pcoa_res, color = "Treatment") +
  geom_point(size = 3) +
  ggtitle("Effect of Pan-probiotics at baseline") +
  theme_minimal() + 
  stat_ellipse(type = "t", linetype = 2)

p1

# Perform dist matrix
distance_matrix <- phyloseq::distance(da1_baseline, method = "wunifrac")

# Perform PERMANOVA
grouping_variable <- "Treatment"
permanova_result1 <- adonis2(distance_matrix ~ get(grouping_variable), as(sample_data(da1_baseline), "data.frame"))

# SAve PERMANOVA results
write.csv(permanova_result1, file = "all_baseline_permanova.csv", row.names = FALSE)

######################################################################################
#############################___CUTSOM__COLOR__PALETTES__#############################

genus_palette <- tax_palette(
  data = clean_data, rank = "Genus", n = 25, pal = "greenArmytage",
  add = c(Other = "purple")
)

phylum_palette <- tax_palette(
  data = clean_data, rank = "Phylum", n = 10, pal = "greenArmytage",
  add = c(Other = "purple")
)

###########################______PHYLUM______#####################################

good_data <- subset_taxa(dada2_phyloseq1, Kingdom == "Bacteria")

clean_data <- phyloseq_validate(good_data) %>%
  tax_fix() %>%
  tax_agg(rank = "Genus")

sham_veh_phylum <- clean_data %>%
  ps_filter(Timepoint == "Baseline", Treatment == "PAN", .keep_all_taxa = TRUE)

sham_veh_phylum %>%
  comp_barplot("Phylum", n_taxa = 5, merge_other = TRUE, label = NULL, palette = phylum_palette) +
  facet_wrap(vars(Treatment), scales = "free") + 
  ggtitle(
    "All groups at 2 weeks",
  ) +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

######################################################################################
##################################_____GENUS_____#####################################

good_data <- subset_taxa(dada2_phyloseq1, Kingdom == "Bacteria")

clean_data <- phyloseq_validate(good_data) %>%
  tax_fix() %>%
  tax_agg(rank = "Genus")

# ba male veh
ba_m_veh_genus <- clean_data %>%
  ps_filter(Timepoint == "35d", Treatment == "Vehicle", Sex == "M", Injury == "Sham", .keep_all_taxa = TRUE)

ba_m_veh_genus %>%
  comp_barplot("Genus", n_taxa = 10, merge_other = TRUE, label = NULL, palette = genus_palette) +
  facet_wrap(vars(Treatment), scales = "free") +
  ggtitle("6 weeks Vehicle Sham Males") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold")) +
  scale_x_discrete(labels = ba_m_veh_genus@sam_data$MouseID)

# ba female veh
ba_f_veh_genus <- clean_data %>%
  ps_filter(Timepoint == "35d", Treatment == "Vehicle", Sex == "F", Injury == "Sham", .keep_all_taxa = TRUE)

ba_f_veh_genus %>%
  comp_barplot("Genus", n_taxa = 10, merge_other = TRUE, label = NULL, palette = genus_palette) +
  facet_wrap(vars(Treatment), scales = "free") +
  ggtitle("6 weeks Vehicle Sham Females") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold")) +
  scale_x_discrete(labels = ba_f_veh_genus@sam_data$MouseID)

# ba male pan
ba_m_pan_genus <- clean_data %>%
  ps_filter(Timepoint == "35d", Treatment == "PAN", Sex == "M", Injury == "Sham", .keep_all_taxa = TRUE)

ba_m_pan_genus %>%
  comp_barplot("Genus", n_taxa = 10, merge_other = TRUE, label = NULL, palette = genus_palette) +
  facet_wrap(vars(Treatment), scales = "free") +
  ggtitle("6 weeks Pan Sham Males") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold")) +
  scale_x_discrete(labels = ba_m_pan_genus@sam_data$MouseID)

# ba female pan
ba_f_pan_genus <- clean_data %>%
  ps_filter(Timepoint == "35d", Treatment == "PAN", Sex == "F", Injury == "Sham", .keep_all_taxa = TRUE)

ba_f_pan_genus %>%
  comp_barplot("Genus", n_taxa = 10, merge_other = TRUE, label = NULL, palette = genus_palette) +
  facet_wrap(vars(Treatment), scales = "free") +
  ggtitle("6 weeks Pan Sham Females") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold")) +
  scale_x_discrete(labels = ba_f_pan_genus@sam_data$MouseID)

## ba cci mice

# ba male veh
ba_m_veh_genus <- clean_data %>%
  ps_filter(Timepoint == "35d", Treatment == "Vehicle", Sex == "M", Injury == "CCI", .keep_all_taxa = TRUE)

ba_m_veh_genus %>%
  comp_barplot("Genus", n_taxa = 10, merge_other = TRUE, label = NULL, palette = genus_palette) +
  facet_wrap(vars(Treatment), scales = "free") +
  ggtitle("6 weeks Vehicle CCI Males") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold")) +
  scale_x_discrete(labels = ba_m_veh_genus@sam_data$MouseID)

# ba female veh
ba_f_veh_genus <- clean_data %>%
  ps_filter(Timepoint == "35d", Treatment == "Vehicle", Sex == "F", Injury == "CCI", .keep_all_taxa = TRUE)

ba_f_veh_genus %>%
  comp_barplot("Genus", n_taxa = 10, merge_other = TRUE, label = NULL, palette = genus_palette) +
  facet_wrap(vars(Treatment), scales = "free") +
  ggtitle("6 weeks Vehicle CCI Females") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold")) +
  scale_x_discrete(labels = ba_f_veh_genus@sam_data$MouseID)

# ba male pan
ba_m_pan_genus <- clean_data %>%
  ps_filter(Timepoint == "35d", Treatment == "PAN", Sex == "M", Injury == "CCI", .keep_all_taxa = TRUE)

ba_m_pan_genus %>%
  comp_barplot("Genus", n_taxa = 10, merge_other = TRUE, label = NULL, palette = genus_palette) +
  facet_wrap(vars(Treatment), scales = "free") +
  ggtitle("6 weeks Pan CCI Males") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold")) +
  scale_x_discrete(labels = ba_m_pan_genus@sam_data$MouseID)

# ba female pan
ba_f_pan_genus <- clean_data %>%
  ps_filter(Timepoint == "35d", Treatment == "PAN", Sex == "F", Injury == "CCI", .keep_all_taxa = TRUE)

ba_f_pan_genus %>%
  comp_barplot("Genus", n_taxa = 10, merge_other = TRUE, label = NULL, palette = genus_palette) +
  facet_wrap(vars(Treatment), scales = "free") +
  ggtitle("6 weeks Pan CCI Females") +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold")) +
  scale_x_discrete(labels = ba_f_pan_genus@sam_data$MouseID)

###############################

sham_veh_genus %>%
  comp_barplot(
    tax_level = "Genus",
    label = "Sex", # name an alternative variable to label axis
    n_taxa = 15, # give more taxa unique colours
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = TRUE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5", # is the default (use NA to remove outlines)
  ) +
  coord_flip()

htmp <- sham_veh_genus %>%
  ps_mutate(Cage = as.character(Cage)) %>%
  tax_transform("log2", add = 1, chain = TRUE) %>%
  comp_heatmap(
    taxa = tax_top(sham_veh_genus, n = 30), grid_col = NA, name = "Log2p",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    colors = heat_palette(palette = viridis::turbo(11)),
    row_names_side = "left", row_dend_side = "right", sample_side = "bottom",
    sample_anno = sampleAnnotation(
      Cage = anno_sample_cat(
        var = "Cage", col = c(C1 = "red", C2 = "blue", C3 = "green", C4 = "orange",
                                 C5 = "purple", C6 = "grey", C7 = "yellow", C8 = "red4",
                                 C9 = "blue4", C10 = "green4", C11 = "orange4", C12 = "purple4"),
        box_col = NA, legend_title = "Cage", size = grid::unit(4, "mm")
      )
    )
  )

ComplexHeatmap::draw(
  object = htmp, annotation_legend_list = attr(htmp, "AnnoLegends"),
  merge_legends = TRUE
)

#####################################################################################
##################################____F__/__B___#####################################

df <- dada2_phyloseq1 %>% tax_glom(taxrank = "Phylum") %>% 
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>% select(OTU, Phylum, Sample, Abundance) %>% 
  spread(Sample, Abundance) 

write.table(df, file = "PAN_relabund_phylum.csv", sep = ",", row.names = F, col.names = T)

# Add cages to metadata
######################################################
cages = read.csv("cages.csv", header = FALSE)

sample_data(dada2_phyloseq1)$Cage <- cages

# This is your desired column data, assuming 'cages' is already defined
cages_data <- cages 

# Correct way to rename and assign the column data
sample_data(dada2_phyloseq1)[, "Cage"] <- cages_data

# Verify the change
colnames(sample_data(dada2_phyloseq1))

View(sample_data(dada2_phyloseq1))
