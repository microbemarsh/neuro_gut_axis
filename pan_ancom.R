#####################################################################
#################__ANCOMBC2_PAN_PROBIOTICS_&_TBI__###################
#####################################################################

library(phyloseq)
library(ANCOMBC)
library(tidyverse)
library(ape)
library(DT)

# Load in my data
raw_physeq <- readRDS("dada2_phyloseq.rds")

# Now we read in the phylo tree computed by qiime2
tree <- read.tree("tree.nwk")

# And add it to our phyloseq object
physeq <- merge_phyloseq(raw_physeq, phy_tree(tree))

# Remove ASVs that break ancombc
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

# ASVs that are no bueno, definitely a better way to do but this is what ANCOMBC2 recommends
badTaxa = c("3ada00c2495bf99d9e0b37663a7ebb41", "93ca769ae201bb1896db00a095254f31", "18bd5640bde80098cbe9a16f97e3b5f7", 
            "452afb6e36f623626dbf48eb303f997b", "9a50c310c570fa39a4245715c680fbfa", "1d7834c6774e38d56ed36b94d6a1b8e6", 
            "70b75647201517cd92b0c9a47eef2531", "ef77162ed5af70274d0a9fe4a4b95ef7", "0dc730c264aa9465ad446195387bea6e", 
            "ed943af1e01a86e2032a7568e8a109d3", "07161dd5f7582d23aac3cc0e165fa987", "a663ce952e845f471b7befe3da01fb1a", 
            "6f1cdfb7e25de2cf3a6cdd75731e7558", "67589da16f2325aca684ac0cf9db2805", "a91fba08c52f5225fd2da894957ab8e3", 
            "8c63b6c34307dddd8477ba30b45e30aa", "ae9b4a948a038f5bedbeebdf2177a910", "6a5228691ba2110e180528e52c03843c", 
            "bee9b4cd24d11a477b134915a0288fea", "28b46ba937b35199648c5e8e5ea411fe", "b7e351cda153dd1ef71bbffa7beadb92", 
            "8a7180f642c6cfe389bf0054632bb492", "748589cabd8a372b95c297827614582c", "e3995b4c960b199a3e144b98b3a47e88", 
            "de8fbb97462d3d1c2776d50fea981ab2", "334087f1453114b9fc76b2ea1ea114cc", "5eeba974b1b592dbd8181c46f372042a",
            "7bc32f5149c74ede9d8a4d5f682b86c6", "c083bea0eee804e67bd856ce40ebbb2a", "20fea129e8d10bb146f6ae389365d306", 
            "99394bbb620a1e1173b2af8c1db030b8", "e69ea0a521268081d7e2ea27a5069024", "6c1205d2d9d213533bd025d14fdbe335", 
            "674d126a03d95767489a10fef7ee078f", "e6e7fcb4277f3886d18fdde923a3db0a", "c1c07a84fd1f4a1df1b02fbdfcb54f93", 
            "c88dcd003681ddd02cf6c6ea4fdbb868", "fda7994e8f4e23aeac39c12f5cb77c1a", "1993b6a31dcf285bfd508a079f2b6455", 
            "db1a86b10760c3f8e5338855178ee352", "2ce3694ccb7606936fb6a1c2952fd045", "7c982b154373fc90fa1a5fea49d730f2", 
            "238e1dd03080f33c38106d760c691f27", "5bcb0d939e56a24e647f439aedc605e9", "3243402b15098ef6e6fb7c24e81757ea", 
            "c5b3a9d0da1b353b2e2842070fb8f010", "61ebe6bc68581c0b726665f67081674b", "17e8f557c793f2fb778b3acb73bcc5aa", 
            "3e0562875a3e7e8064b04f8a61fe7fea", "0deb81c1c562dbac476b2cf1471ebf48", "d5e9ebc22ca9b62b1fffe7782e365747", 
            "a6f1dda94db6fd6e6761a502c190e04f", "0c06326fb7a30a29d1065665bc51f744", "74679f3f81f03cd57d914304ce5fdbd1", 
            "61c37c1e5cf6b2f180632911a3ef2775", "495a59a0aa0ced822fdfae190b48d603", "63a206a073822a020585f225ce0baa21", 
            "a4d598c4b0ca87e01e75e6ef01c34d96", "7e3e77021c0de79c4a484a6b934962ff", "4e6fc05a624a1bb3dc9f2aa03cca6fdf", 
            "48271e1ce614331799b1c5a1edcc84a3", "f7ad339b9d39d000a2d90e2bd743b47a", "8322218563b33f3b8214730d9f8b06fd", 
            "dddd0458189606240dc7fe2473a38fd3", "ae1bf7d98157e5a108ee01a3d8a18dfd", "b79953dbfda019f098fc9c54093407d2", 
            "d1eba4439565931cb8626cfa2baeff5d", "29c16acd3ee197be8d0e5c4af4c5767c", "239a693a3e61c35d13d9c57cb1f8ae46", 
            "baff21923b7226f744174447c6421d74", "7b4e1a2e5ea9ed735669c8bb6d8d0e6d", "5661cfcf5a71630963f6b22e472d9dd6",
            "2b3b64ed13bbe0e244f1555509db42d4")

# Cleaned up physeq for ancombc
clean_physeq = pop_taxa(physeq, badTaxa)

# Phyloseq to ancombc
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(clean_physeq)

# subset to groups to compare
tse2 = tse[, tse$Treatment == "PAN" & tse$Injury == "Sham"]

# Set seed for ancombc
set.seed(123)

# Run ancombc2
output3 = ancombc2(data = tse2, assay_name = "counts", tax_level = "Species",
                   fix_formula = "Timepoint", rand_formula = NULL,
                   p_adj_method = "holm", pseudo_sens = TRUE,
                   prv_cut = 0.2, lib_cut = 1000, s0_perc = 0.05,
                   group = "Timepoint", struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                   iter_control = list(tol = 1e-2, max_iter = 20, 
                                       verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl(),
                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                   trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                               nrow = 2, 
                                                               byrow = TRUE),
                                                        matrix(c(-1, 0, 1, -1),
                                                               nrow = 2, 
                                                               byrow = TRUE),
                                                        matrix(c(1, 0, 1, -1),
                                                               nrow = 2, 
                                                               byrow = TRUE)),
                                        node = list(2, 2, 1),
                                        solver = "ECOS",
                                        B = 100))
# How many structural zeroes?
tab_zero = output3$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")

# Print output to dataframe
res_prim3 = output3$res

# SAve the ancombc stats for each group comparison
write.csv(res_prim3, file = "ancom_timepoint_PAN_sham.csv", row.names = FALSE)
