#################################################################
library(readr)
library(ggpicrust2)
library(tidyverse)
library(patchwork)
library(ggprism)
library(KEGGREST)
library(ggh4x)

meta = read.delim("PAN_Metadata.tsv")

kegg_abundance <- ko2kegg_abundance(file = "~/Documents/HMRI/PAN_paper/picrust/pred_metagenome_unstrat.tsv") # Or use data(kegg_abundance)
# Subset data for comparison groups
veh_sham_meta = meta %>%
  filter(Treatment == "Vehicle" & Injury == "Sham")

print(veh_sham_meta$sampleID)

veh_sham = kegg_abundance %>%
  select("Villapol_284_B", "Villapol_285_B", "Villapol_286_B", "Villapol_287_B",
         "Villapol_288_B", "Villapol_289_B", "Villapol_290_B", "Villapol_291_B",  
         "Villapol_292_B", "Villapol_293_B", "Villapol_35d_284", "Villapol_35d_285",
         "Villapol_35d_286", "Villapol_35d_287", "Villapol_35d_288", "Villapol_35d_289",
         "Villapol_35d_290", "Villapol_35d_291", "Villapol_35d_292", "Villapol_35d_293")

daa_results_df <- pathway_daa(veh_sham, metadata = veh_sham_meta, group = "Timepoint", daa_method = "LinDA")
daa_annotated_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df, ko_to_kegg = TRUE)
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = veh_sham,
                      daa_results_df = daa_annotated_results_df,
                      Group = veh_sham_meta$Timepoint,
                      ko_to_kegg = TRUE,
                      p_values_threshold = 0.05,
                      order = "pathway_class",
                      select = NULL,
                      p_value_bar = TRUE,
                      colors = NULL,
                      x_lab = "pathway_name")


p

# If that dont work
daa_results_df1 <- daa_results_df %>% filter(p_adjust < 0.05) %>% slice(1:29)

# Annotate pathway results using KO to KEGG conversion
daa_annotated_results_df1 <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df1, ko_to_kegg = TRUE)

# Generate pathway error bar plot
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = veh_sham, daa_results_df = daa_annotated_results_df1, 
                      Group = veh_sham_meta$Timepoint, p_values_threshold = 0.05, 
                      order = "pathway_class", select = NULL, ko_to_kegg = TRUE, 
                      p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")

p

#################################################################################################
#         _________VEH____SHAM___ABOVE____          ________VEH______CCI_____BELOW_______
################################################################################################

meta = read.delim("PAN_Metadata.tsv")

kegg_abundance <- ko2kegg_abundance(file = "~/Documents/HMRI/PAN_paper/picrust/pred_metagenome_unstrat.tsv")

# Subset data for comparison groups
veh_cci_meta = meta %>%
  filter(Treatment == "Vehicle" & Injury == "CCI")

print(veh_cci_meta$sampleID)

veh_cci = kegg_abundance %>%
  select("Villapol_254_B", "Villapol_255_B", "Villapol_256_B", "Villapol_257_B", "Villapol_258_B",  
 "Villapol_259_B", "Villapol_260_B", "Villapol_261_B", "Villapol_262_B", "Villapol_263_B",  
"Villapol_264_B", "Villapol_265_B", "Villapol_266_B", "Villapol_267_B", "Villapol_268_B",  
"Villapol_269_B", "Villapol_270_B", "Villapol_271_B", "Villapol_272_B", "Villapol_273_B",  
"Villapol_35d_254", "Villapol_35d_255", "Villapol_35d_256", "Villapol_35d_257", "Villapol_35d_258",
"Villapol_35d_259", "Villapol_35d_261", "Villapol_35d_262", "Villapol_35d_263", "Villapol_35d_264",
"Villapol_35d_265", "Villapol_35d_266", "Villapol_35d_267", "Villapol_35d_268", "Villapol_35d_269",
"Villapol_35d_270", "Villapol_35d_271", "Villapol_35d_272", "Villapol_35d_273")

daa_results_df <- pathway_daa(veh_cci, metadata = veh_cci_meta, group = "Timepoint", daa_method = "LinDA")
daa_annotated_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df, ko_to_kegg = TRUE)
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = veh_cci,
                      daa_results_df = daa_annotated_results_df,
                      Group = veh_cci_meta$Timepoint,
                      ko_to_kegg = TRUE,
                      p_values_threshold = 0.05,
                      order = "pathway_class",
                      select = NULL,
                      p_value_bar = TRUE,
                      colors = NULL,
                      x_lab = "pathway_name")


p

# If that dont work
daa_results_df1 <- daa_results_df %>% filter(p_adjust < 0.05) %>% slice(1:29)

# Annotate pathway results using KO to KEGG conversion
daa_annotated_results_df1 <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df1, ko_to_kegg = TRUE)

# Generate pathway error bar plot
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = veh_cci, daa_results_df = daa_annotated_results_df1, 
                      Group = veh_cci_meta$Timepoint, p_values_threshold = 0.05, 
                      order = "pathway_class", select = NULL, ko_to_kegg = TRUE, 
                      p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")

p

#################################################################################################
#         _________VEH____CCI___ABOVE____          ________PAN______SHAM_____BELOW_______
################################################################################################

meta = read.delim("PAN_Metadata.tsv")

kegg_abundance <- ko2kegg_abundance(file = "~/Documents/HMRI/PAN_paper/picrust/pred_metagenome_unstrat.tsv")

# Subset data for comparison groups
pan_sham_meta = meta %>%
  filter(Treatment == "PAN" & Injury == "Sham")

print(pan_sham_meta$sampleID)

pan_sham = kegg_abundance %>%
  select("Villapol_274_B", "Villapol_275_B", "Villapol_276_B", "Villapol_277_B",  
 "Villapol_278_B", "Villapol_279_B", "Villapol_280_B", "Villapol_281_B",
 "Villapol_282_B", "Villapol_283_B", "Villapol_35d_274", "Villapol_35d_275",
"Villapol_35d_276", "Villapol_35d_277", "Villapol_35d_278", "Villapol_35d_279",
"Villapol_35d_280", "Villapol_35d_281", "Villapol_35d_282", "Villapol_35d_283")

daa_results_df <- pathway_daa(pan_sham, metadata = pan_sham_meta, group = "Timepoint", daa_method = "LinDA")
daa_annotated_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df, ko_to_kegg = TRUE)
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = pan_sham,
                      daa_results_df = daa_annotated_results_df,
                      Group = pan_sham_meta$Timepoint,
                      ko_to_kegg = TRUE,
                      p_values_threshold = 0.05,
                      order = "pathway_class",
                      select = NULL,
                      p_value_bar = TRUE,
                      colors = NULL,
                      x_lab = "pathway_name")


p

# If that dont work
daa_results_df1 <- daa_results_df %>% filter(p_adjust < 0.05) %>% slice(1:29)

# Annotate pathway results using KO to KEGG conversion
daa_annotated_results_df1 <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df1, ko_to_kegg = TRUE)

# Generate pathway error bar plot
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = pan_sham, daa_results_df = daa_annotated_results_df1, 
                      Group = pan_sham_meta$Timepoint, p_values_threshold = 0.05, 
                      order = "pathway_class", select = NULL, ko_to_kegg = TRUE, 
                      p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")

p

#################################################################################################
#         _________PAN____SHAM___ABOVE____          ________PAN______CCI_____BELOW_______
################################################################################################

meta = read.delim("PAN_Metadata.tsv")

kegg_abundance <- ko2kegg_abundance(file = "~/Documents/HMRI/PAN_paper/picrust/pred_metagenome_unstrat.tsv")

# Subset data for comparison groups
pan_cci_meta = meta %>%
  filter(Treatment == "PAN" & Injury == "CCI")

print(pan_cci_meta$sampleID)

pan_cci = kegg_abundance %>%
  select("Villapol_234_B", "Villapol_236_B", "Villapol_237_B", "Villapol_238_B", "Villapol_239_B",  
         "Villapol_240_B", "Villapol_241_B", "Villapol_242_B", "Villapol_243_B", "Villapol_244_B",  
         "Villapol_245_B", "Villapol_246_B", "Villapol_247_B", "Villapol_248_B", "Villapol_249_B", 
         "Villapol_250_B", "Villapol_251_B", "Villapol_252_B", "Villapol_253_B", "Villapol_35d_234",
         "Villapol_35d_236", "Villapol_35d_237", "Villapol_35d_238", "Villapol_35d_239", "Villapol_35d_240",
         "Villapol_35d_241", "Villapol_35d_242", "Villapol_35d_243", "Villapol_35d_244", "Villapol_35d_245",
         "Villapol_35d_246", "Villapol_35d_247", "Villapol_35d_248", "Villapol_35d_249", "Villapol_35d_250",
         "Villapol_35d_251", "Villapol_35d_252", "Villapol_35d_253")

daa_results_df <- pathway_daa(pan_cci, metadata = pan_cci_meta, group = "Timepoint", daa_method = "LinDA")
daa_annotated_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df, ko_to_kegg = TRUE)
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = pan_cci,
                      daa_results_df = daa_annotated_results_df,
                      Group = pan_cci_meta$Timepoint,
                      ko_to_kegg = TRUE,
                      p_values_threshold = 0.05,
                      order = "pathway_class",
                      select = NULL,
                      p_value_bar = TRUE,
                      colors = NULL,
                      x_lab = "pathway_name")


p

##############################################################################
################_____Trying__with__MetaCyc__pathways_____#####################

metacyc = read.delim("METACYC_path_abun_unstrat_descrip.tsv")

pan_cci_meta = meta %>%
  filter(Treatment == "PAN" & Injury == "CCI")

df = metacyc[,-2]

pan_cci = df %>%
  select("pathway","Villapol_234_B", "Villapol_236_B", "Villapol_237_B", "Villapol_238_B", "Villapol_239_B",  
         "Villapol_240_B", "Villapol_241_B", "Villapol_242_B", "Villapol_243_B", "Villapol_244_B",  
         "Villapol_245_B", "Villapol_246_B", "Villapol_247_B", "Villapol_248_B", "Villapol_249_B", 
         "Villapol_250_B", "Villapol_251_B", "Villapol_252_B", "Villapol_253_B", "Villapol_35d_234",
         "Villapol_35d_236", "Villapol_35d_237", "Villapol_35d_238", "Villapol_35d_239", "Villapol_35d_240",
         "Villapol_35d_241", "Villapol_35d_242", "Villapol_35d_243", "Villapol_35d_244", "Villapol_35d_245",
         "Villapol_35d_246", "Villapol_35d_247", "Villapol_35d_248", "Villapol_35d_249", "Villapol_35d_250",
         "Villapol_35d_251", "Villapol_35d_252", "Villapol_35d_253") %>%
  column_to_rownames(var = "pathway")
         
daa_results_df <- pathway_daa(pan_cci, 
                              metadata = pan_cci_meta, 
                              group = "Timepoint", 
                              daa_method = "LinDA")

daa_results_df1 = daa_results_df %>% 
  filter(p_values < 0.02) 

daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc",
                                               daa_results_df = daa_results_df1,
                                               ko_to_kegg = FALSE)

# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = pan_cci,
                      daa_results_df = daa_annotated_results_df ,
                      Group = pan_cci_meta$Timepoint,
                      ko_to_kegg = FALSE,
                      p_values_threshold = 0.05,
                      order = "group",
                      select = NULL,
                      p_value_bar = TRUE,
                      colors = NULL,
                      x_lab = "description")

p





# If that dont work
daa_results_df1 <- daa_results_df %>% filter(p_adjust < 0.05) %>% slice(1:29)

# Annotate pathway results using KO to KEGG conversion
daa_annotated_results_df1 <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df1, ko_to_kegg = TRUE)

# Generate pathway error bar plot
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = pan_cci, daa_results_df = daa_annotated_results_df1, 
                      Group = pan_cci_meta$Timepoint, p_values_threshold = 0.05, 
                      order = "pathway_class", select = NULL, ko_to_kegg = TRUE, 
                      p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")

p



p2 <- pathway_heatmap(kegg_abundance %>% rownames_to_column("feature") %>% filter(feature %in% daa_annotated_sub_method_results_df1$feature) %>% column_to_rownames("feature"), metadata, "traitment")