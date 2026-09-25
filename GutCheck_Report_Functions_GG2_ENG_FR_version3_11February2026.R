library(tidyverse)
library(lessR)
library(stringi)
library(magick)
library(vegan)
library(reshape2)
library(ecodist)

## Script containing the main functions currently used for report generation 


## Version February 11th, 2025

# !!! Major change : this version uses reference quartile table, rather than a reference BIOM table as input for all functions

# Other changes: 
## Added interleaving for the phylum figures to potentially help with overlaps & also adjusted height, width parameters 
### minor change Feb 18th, 2026 - probiotic bacteria recommendation added and short_report figures added
### minor change March 18th, 2026 - GutPotentialDF now inner_joined with combined_DF so that individual Vitamin scores are included in Excel
### minor change April 9th, 2026 - added Short Report Excel file 
###                             - in cleanGenus function Genus == " g__" added a white space
### minor change May 12th, 2026 - update how column names are read for alpha-diversity
### minor change June 10th, 2026 - updated short figure colour and removed labels 

### Required input files 
#   reference_valuesAGP_GG2_allsamples_genus 
#   benefBac_AGP_quantiles_by_genus_GG2_summedLactobacilli
#   biomtable 



### Main functions list



### Functions that are used-within the main functions 

# cleanGenusNames(inputFilePrefix, InputTable)
# createGenusTable(inputFilePrefix, InputTable) #similar to above, in the future can combine them  
# calculateSymptomScoreGenus(inputFilePrefix, UpBacteriaList, DownBacteriaList, InputTable, ReferenceQuartileTable)
# calculateProcessScoreGenus(inputFilePrefix, UpBacteriaList, DownBacteriaList, InputTable, ReferenceQuartileTable)
# 
#






#reference_valuesAGP_GG2_allsamples_genus <-  read_tsv("C:/Users/farid/OneDrive/Desktop/NucliqBiologics/Reference_datasets/AGP/GG2/AGP_quantiles_by_genus_GG2.tsv.txt")

reference_valuesAGP_GG2_allsamples_genus <- read_tsv(file.path(getwd(), "AGP_quantiles_by_genus_GG2.tsv.txt"))

reference_valuesAGP_GG2_allsamples_genus <- reference_valuesAGP_GG2_allsamples_genus %>% 
  add_row(Genus = "Escherichia", q25_percent = 0, q50_percent = 0, q75_percent = 0.1) %>%  #E.coli normally represents around 0.1% of the gut flora in humans (Al-Zyoud et al. 2019)  
  add_row(Genus = "Klebsiella", q25_percent = 0, q50_percent = 0.0156, q75_percent = 0.031) %>%   # Klebsiella in adults (median: 0.0156, range: 0–0.031) (Liao et al. 2024) 
  add_row(Genus = "Lacticaseibacillus", q25_percent = 0, q50_percent = 0, q75_percent = 0)  %>%  # if it is present at all it is already more than in the ref. pop. as it is typically not detected 
  add_row(Genus = "Weizmannia", q25_percent = 0, q50_percent = 0, q75_percent = 0)  %>%  # if it is present at all it is already more than in the ref. pop. as it is typically not detected 
  add_row(Genus = "Lactiplantibacillus", q25_percent = 0, q50_percent = 0, q75_percent = 0)  # if it is present at all it is already more than in the ref. pop. as it is typically not detected 


# FYI no Enterococcus in AGP women only cohort, but it is present in the general cohort, that's why it is not manually added 


reference_valuesAGP_GG2_allsamples_genus




#benefBac_reference <-  read_tsv("C:/Users/farid/OneDrive/Desktop/NucliqBiologics/Reference_datasets/AGP/GG2/benefBac_AGP_quantiles_by_genus_GG2_summedLactobacilli.txt")
benefBac_reference <- read_tsv(file.path(getwd(), "benefBac_AGP_quantiles_by_genus_GG2_summedLactobacilli.txt"))
benefBac_reference





### Functions for data pre-processing

# Helper: resolve alpha-diversity.tsv path regardless of whether it is a file (new fix)
# or a directory (legacy buggy Java behaviour where a folder was created with that name)
get_alpha_diversity_path <- function() {
  path_as_dir <- "./alpha-diversity.tsv"
  if (dir.exists(path_as_dir)) {
    # Old Java bug created a directory – read the file inside it
    inner <- file.path(path_as_dir, "alpha-diversity.tsv")
    if (file.exists(inner)) return(inner)
  }
  # Normal case: it is a plain file
  return(path_as_dir)
} 


cleanGenusNames <- function(inputFilePrefix, InputTable) {
  #how to use: instead of InputTable, specify biomtable 
  #this will output a NON-Aggregated Table 
  
  sampleID <- paste(inputFilePrefix)  
  
  sampleID_table_all <- InputTable %>% 
    select(sampleID, taxonomy) %>% 
    separate(taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")  
  #filter_at(sampleID, ~ . != 0)
  
  
  
  sampleID_table_all %>% 
    select(sampleID, Genus, Family) %>% 
    mutate(Genus = gsub("_[A-Z]$", '', Genus)) %>%      ## note the change
    mutate(Genus = gsub("_[A-Z]_.*$", '', Genus)) %>%  
    mutate(Genus = gsub("_[0-9].*$", '', Genus)) %>% 
    mutate(Family = gsub("_[A-Z]$", '', Family)) %>%      ## note the change
    mutate(Family = gsub("_[A-Z]_.*$", '', Family)) %>%  
    mutate(Family = gsub("_[0-9].*$", '', Family)) %>% 
    mutate(Genus = if_else(Genus == "g" | Genus == " g__" | is.na(Genus), paste("unclassified", Family), Genus)) %>%  #note that unknown genera with known family shown as unclassified family name
    select(-Family) %>% 
    group_by(Genus) %>% 
    summarise(RA = sum(.data[[sampleID]])) %>% 
    arrange(desc(RA)) %>% 
    mutate_at("Genus", str_replace, " ", "") %>% 
    mutate_at("Genus", str_replace, "\\[", "") %>% 
    mutate_at("Genus", str_replace, "\\]", "") %>%
    mutate_at("Genus", str_replace, "g__$", "Other") %>%
    mutate_at("Genus", str_replace, ".*human", "Other") %>%
    mutate_at("Genus", str_replace, "[f,g]__", "") %>% 
    mutate_at("Genus", str_replace, "g_$", "Other") %>% 
    replace_na(list(Genus = "Other")) %>%  #to replace NA have to replace_na not mutate
    mutate_at("Genus", str_replace, "unclassifiedNA", "Other")
  
}






createGenusTable <- function(inputFilePrefix, InputTable)  {
  
  
  #how to use: instead of InputTable, specify biomtable 
  #this will output Aggregated Table 
  
  sampleID <- paste(inputFilePrefix)  
  
  sampleID_table_all <- InputTable %>% 
    select(sampleID, taxonomy) %>% 
    separate(taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")  
  #filter_at(sampleID, ~ . != 0)
  
  
  
  genus <- sampleID_table_all %>% 
    select(sampleID, Genus, Family) %>% 
    mutate(Genus = gsub("_[A-Z]$", '', Genus)) %>%     
    mutate(Genus = gsub("_[A-Z][A-Z]$", '', Genus)) %>% #remove 2letter codes March2025
    mutate(Genus = gsub("_[A-Z]_.*$", '', Genus)) %>%   
    mutate(Family = gsub("_[A-Z]$", '', Family)) %>%      
    mutate(Family = gsub("_[A-Z]_.*$", '', Family)) %>%  
    mutate(Family = gsub("_[0-9].*$", '', Family)) %>% 
    mutate(Genus = if_else(Genus == "g" | Genus == "g__" | is.na(Genus), paste("unclassified", Family), Genus)) %>%  #unknown genera with known family shown as unclassified family name
    select(-Family) %>% 
    group_by(Genus) %>% 
    summarise(RA = sum(.data[[sampleID]])) %>% 
    arrange(desc(RA)) %>% 
    mutate_at("Genus", str_replace, " ", "") %>% 
    mutate_at("Genus", str_replace, "\\[", "") %>% 
    mutate_at("Genus", str_replace, "\\]", "") %>%
    mutate_at("Genus", str_replace, "g__$", "Other") %>%
    mutate_at("Genus", str_replace, ".*human", "Other") %>%
    mutate_at("Genus", str_replace, "[f,g]__", "") %>% 
    mutate_at("Genus", str_replace, "g_$", "Other") %>% 
    mutate(Genus = gsub("_[0-9].*$", '', Genus)) %>% #this is placed towards the end to make sure that genus names starting with a number do not get replaced UNTIL g__ has been removed 
    replace_na(list(Genus = "Other")) %>%  #to replace NA have to replace_na not mutate
    mutate_at("Genus", str_replace, "unclassifiedNA", "Other")
  
  all_genus <- aggregate(RA ~ Genus, genus, sum)
  return(all_genus)
}





#make sure sample val and q1 and q3 comparison are correct

#If this is a symptom, then 

#For bacteria pos.ass with Symptom: sample > q3_val: counts as point i.e. bacteria that increases symptom too high  
#For bacteria neg.ass with symptom: sample < q1_val: counts as point i.e. bacteria that decreases symptom too low 
#more points: means more likelihood of symptom 

#same for inflammation related scores etc

# !!! But if using to calculate things such Muscle Strength or Aerobic Endurance, logic is reverse - to get higher score there needs to be MORE bacteria pos. associated with these things and less bacteria neg. associated with these things  

calculateSymptomScoreGenus <- function(inputFilePrefix, UpBacteriaList, DownBacteriaList, InputTable, ReferenceQuartileTable) {
  
  
  all_genus <- cleanGenusNames(inputFilePrefix, InputTable) #this is not aggregated 
  genus_ag <- aggregate(RA ~ Genus, all_genus, sum) ### this step sums up the values for all genera, incl. Other, shouldn't contain dupl.rows 
  
  UpVec <- c() #Take the vector outside of the loop, so that whatever is output from a loop gets appended to it and the end result is one vector with multiple elements 
  for(i in UpBacteriaList) { 
    q3_val <- ReferenceQuartileTable %>% filter(Genus == i) %>%  select(q75_percent)
    q3_val <- q3_val/100
    sample_val <- subset(genus_ag, Genus == i)$RA
    sample_val <- ifelse(length(sample_val) == 0, 0, sample_val)
    val_comparison <- sample_val >= q3_val
    UpVec <- c(UpVec, val_comparison)
  }
  
  DownVec <- c()
  for (i in DownBacteriaList) { 
    q1_val <- ReferenceQuartileTable %>%  filter(Genus == i) %>%  select(q25_percent)
    q1_val <- q1_val/100 
    sample_val <- subset(genus_ag, Genus == i)$RA
    sample_val <- ifelse(length(sample_val) == 0, 0, sample_val)
    val_comparison <- sample_val < q1_val
    DownVec <- c(DownVec, val_comparison)
  }
  
  SumVec <- c(UpVec, DownVec)
  total_genus_associated_with_symptom <- sum(SumVec , na.rm = TRUE)/length(SumVec)
  
  return(total_genus_associated_with_symptom)
  
}














### Logic explanation

##### Example: Lactobacillus and Roseburia pos.assoc. with SCORE
#####          Bacteroides neg.asssoc. with SCORE
#####   if sample has Lactobacillus and Roseburia > q1 of ref.value ==> this adds up to 2 points 
#####   if sample has Bacteroides < q3 of ref.value (within normal value) ==> this adds up to 1 point 
#####   total points 3/3 

####    but if there is Bacteroides > q3 (more than normal ref.value) ==> 0 points 
####    total points 2/3 

##### NOTE: this will calculate instances of TRUE in a list 
##### testlist <- c(TRUE, FALSE, FALSE)
##### sum(testlist , na.rm = TRUE)/length(testlist) #should give 0.33 


calculateProcessScoreGenus <- function(inputFilePrefix, UpBacteriaList, DownBacteriaList, InputTable, ReferenceQuartileTable) {
  
  
  all_genus <- cleanGenusNames(inputFilePrefix, InputTable) #this is not aggregated 
  genus_ag <- aggregate(RA ~ Genus, all_genus, sum) ### this step sums up the values for all genera, incl. Other, shouldn't contain dupl.rows 
  
  UpVec <- c() #Take the vector outside of the loop, so that whatever is output from a loop gets appended to it and the end result is one vector with multiple elements 
  for(i in UpBacteriaList) { 
    q1_val <- ReferenceQuartileTable %>% filter(Genus == i) %>%  select(q25_percent)
    q1_val <- q1_val/100
    sample_val <- subset(genus_ag, Genus == i)$RA
    sample_val <- ifelse(length(sample_val) == 0, 0, sample_val)
    val_comparison <- sample_val >= q1_val ## not how this is different from Symptom 
    UpVec <- c(UpVec, val_comparison)
  }
  
  DownVec <- c()
  for (i in DownBacteriaList) { 
    q3_val <- ReferenceQuartileTable %>%  filter(Genus == i) %>%  select(q75_percent)
    q3_val <- q3_val/100 
    sample_val <- subset(genus_ag, Genus == i)$RA
    sample_val <- ifelse(length(sample_val) == 0, 0, sample_val)
    val_comparison <- sample_val < q3_val
    DownVec <- c(DownVec, val_comparison)
  }
  
  SumVec <- c(UpVec, DownVec)
  total_genus_associated_with_Score <- sum(SumVec , na.rm = TRUE)/length(SumVec)
  
  return(total_genus_associated_with_Score)
  
}














#######################################
######## MAIN FUNCTIONS ###############
#######################################


#_________________ getSampleTaxonomicInfo_GG2 __________________________

getSampleTaxonomicInfo_GG2_v3 <- function(inputFilePrefix, InputTable, language = "ENG") {
  
  palette_colours <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17","#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02","#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#E5D8BD", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#E6F5C9", "#FFF2AE", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3","#FF7F00","#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62","#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7","#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#FFED6F")
  
  #colours previously used: "#D9D9D9", "#F2F2F2", "#F4CAE4",  "#F1E2CC",  "#FFFFCC","#FDDAEC"
  
  output_dir <- paste("./",inputFilePrefix,sep="")
  
  dir.create(output_dir)
  
  
  
  sampleID <- paste(inputFilePrefix)  
  
  sampleID_table_all <- biomtable %>% 
    select(sampleID, taxonomy) %>% 
    separate(taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% 
    filter_at(sampleID, ~ . != 0)
  
  sampleID_table_all
  
  
  
  ### Phylum 
  
  
  
  
  
  all_phylum <-  
    sampleID_table_all %>% 
    select(sampleID, Phylum) %>% 
    group_by(Phylum)  %>% 
    mutate(Phylum = gsub("_[A-Z]$", '', Phylum)) %>%      ## note the change
    mutate(Phylum = gsub("_[A-Z]_.*$", '', Phylum)) %>%   ## note the change
    summarise(RA = sum(.data[[sampleID]])) %>% 
    mutate_at("Phylum", str_replace, " ", "") %>% 
    mutate_at("Phylum", str_replace, "\\[", "") %>% 
    mutate_at("Phylum", str_replace, "\\]", "") %>% 
    mutate_at("Phylum", str_replace, "p__$", "Other") %>% ## note the change 
    replace_na(list(Phylum = "Other")) %>%                ## note the change 
    mutate_at("Phylum", str_replace, "p__", "")
  
  
  all_phylum <- if(language == "ENG")       # This DF is exported to Taxonomy file 
  { 
    all_phylum <- all_phylum 
  } else {
    
    all_phylum <- all_phylum %>% 
      mutate(Phylum = gsub("Other", 'Les autres phyla', Phylum))   
  }
  
  
  ### Now making a dataframe used for the graph 
  
  phylum_ref <- all_phylum
  i1 <- phylum_ref$RA <= 0.01
  phylum_ref$Phylum[i1] <- "Other"
  
  phylum_ref <- if(language == "ENG")
  { 
    phylum_ref <- phylum_ref 
  } else {
    phylum_ref <- phylum_ref %>% 
      mutate(Phylum = gsub("Other", 'Les autres phyla', Phylum))  
    
  }
  
  phylum_ag <- aggregate(RA ~ Phylum, phylum_ref, sum)
  phylum_ag$Percentage <- paste0(phylum_ag$Phylum,' [',100*round(phylum_ag$RA,10),'%',']')
  
  
  
  phylum_ag_sorted <- phylum_ag %>% arrange(RA)
  
  n <- nrow(phylum_ag_sorted)
  half_n <- ceiling(n / 2)
  
  df_small <- head(phylum_ag_sorted, half_n)
  df_large <- tail(phylum_ag_sorted, n - half_n)
  
  
  
  interleave_rows <- function(df1, df2) {
    idx <- order(c(1:nrow(df1), 1:nrow(df2)))
    combined_df <- rbind(df1, df2)[idx, ]
    return(combined_df)
  }
  
  phylum_ag_interleaved <- interleave_rows(df_small, df_large) #%>%  filter(RA > 0.01) # If Phylum % is below 1% --> do not display? 
  print(phylum_ag_interleaved)
  
  
  
  
  ### note the change in phylum names and colour assignments!
  color_assgn <- c("Actinobacteriota"="#FC8D62",
                   "Bacteroidota"="#A6D854",
                   "Campylobacterota"="#BEBADA",
                   "Cyanobacteria"="#E78AC3",
                   "Deferribacterota"="#FFFFB3",
                   "Deinococcota"="#8DD3C7",
                   "Desulfobacterota"="#E5C494",
                   "Elusimicrobiota"="#80B1D3",
                   "Eremiobacterota"="#CCEBC5",
                   "Firmicutes"="#BC80BD",
                   "Fusobacteriota"="#B3B3B3",
                   "Methanobacteriota"="#8DA0CB",
                   "Les autres phyla" = "#D9D9D9",
                   "Other" = "#D9D9D9",
                   "Patescibacteria"="#FFED6F",
                   "Proteobacteria"="#FFD92F",
                   "Spirochaetota"="#FCCDE5",
                   "Synergistota"="#FB8072",
                   "Thermoplasmatota"="#66C2A5",
                   "Verrucomicrobiota"="#FDB462")
  
  
  color_assgn_df <- as.data.frame(color_assgn) %>%
    rownames_to_column(var = "Phylum")
  
  
  color_for_sampleDF <- inner_join(color_assgn_df, phylum_ag_interleaved, by = "Phylum")
  print(color_for_sampleDF)
  
  color_for_sample <- inner_join(phylum_ag_interleaved, color_assgn_df, by = "Phylum") %>%
    #arrange(Phylum) %>%
    select(color_assgn)
  
  
  assign("overwriteable_phylum_ag_in_globalenv", phylum_ag_interleaved, envir = .GlobalEnv)
  assign("overwriteable_colours_in_globalenv", color_for_sample, envir = .GlobalEnv)
  
  
  donutplotname <- paste(output_dir,"/",inputFilePrefix,"_phylumDonutChart", sep = "")
  
  assign("overwriteable_pathToPlot_in_globalenv", donutplotname, envir = .GlobalEnv)
  
  
  
  #Then create the PieChart
  
  ### Note the updated values_size and width
  
  PieChart(x = Phylum, y = RA, data = overwriteable_phylum_ag_in_globalenv, 
           fill = as.matrix(overwriteable_colours_in_globalenv$color_assgn), 
           values = "%",
           values_digits = 1, #round to 2 decimal place
           values_position = "in",        #values_size = 0.78,
           values_size = 0.90,
           values_color = "black",
           clockwise = TRUE,
           init_angle = 0,
           main = NULL, width = 10, height = 8, pdf_file = overwriteable_pathToPlot_in_globalenv) #height = 8
  
  donutPDF_path <- paste(output_dir,"/",inputFilePrefix,"_phylumDonutChart.pdf", sep = "")
  donutPNG_path <- paste(output_dir,"/",inputFilePrefix,"_phylumDonutChart.png", sep = "")
  pdf_convert(donutPDF_path, format = "png", pages = NULL, filenames = donutPNG_path, dpi = 300, opw = "", upw = "", verbose = TRUE)
  
  
  
  ###Genus
  
  all_genus <- createGenusTable(inputFilePrefix, biomtable)
  
  all_genus <- if (language == "ENG") { 
    all_genus <- all_genus
    
  } else {
    
    all_genus <- all_genus %>% 
      mutate(Genus = gsub("unclassified", "Les", Genus)) %>% 
      mutate(Genus = gsub("eae", "eae non classifiés", Genus)) %>% 
      mutate(Genus = gsub("Other", "Les autres bactéries", Genus))
  }
  
  
  
  
  genus_ref <- all_genus #make a copy, then replace values below threshold in genus_ref
  i6 <- genus_ref$RA <= 0.003 # keep only 0.3% and higher  
  genus_ref$Genus[i6] <- if (language == "ENG") { 
    genus_ref$Genus[i6] <- "Other"
    
  } else {
    genus_ref$Genus[i6] <- "Les autres bactéries"
  }
  
  genus_ag <- aggregate(RA ~ Genus, genus_ref, sum)
  genus_ag$Percentage <- paste0(genus_ag$Genus,' [',100*round(genus_ag$RA,10),'%',']')
  
  
  
  
  
  genus_plot <- ggplot(genus_ag, aes(x=RA, y= reorder(Genus,RA), fill=Percentage)) + 
    geom_bar(width = 1, stat = "identity") + 
    theme_classic() + 
    #coord_polar("y", start=0) +
    theme_classic() + 
    scale_fill_manual(values = palette_colours) + 
    guides(fill=guide_legend(title="")) +
    theme(legend.key.size = unit(1,"line")) +
    labs(caption = "Most abundant bacterial genera present in your sample") +
    ylab("") + xlab("") +
    geom_text(aes(label=paste0(round(RA,4)*100,"%")), hjust = 0.5 - sign(genus_ag$RA)/2, size = 3.5) +
    theme(legend.position="none") +
    scale_x_continuous(limits = c(0, max(genus_ag$RA + 0.05))) +  theme(axis.title.x=element_blank(),
                                                                        axis.text.x=element_blank(),
                                                                        axis.ticks.x=element_blank())
  
  
  genus_plot <- if (language == "ENG") { 
    genus_plot <- genus_plot 
  } else {
    genus_plot <- genus_plot + labs(caption = "Les bactéries les plus abondantes présentes dans votre échantillon") 
    
  }
  
  
  plotname_genus <- paste(output_dir,"/",inputFilePrefix,"_genus.png", sep = "")
  
  ggsave(plotname_genus, plot = last_plot())
  
  ### Output Table 
  
  all_genus <- all_genus %>%   ## dataframe to export 
    arrange(desc(RA)) %>% 
    filter(RA != 0) %>% 
    mutate(Percentage = round(RA*100,6)) %>% 
    select(Genus, Percentage)
  
  table_name_all_genus <- paste(output_dir,"/",inputFilePrefix,"_all_genera.txt", sep = "")
  write.table(all_genus, table_name_all_genus, sep = "\t", quote = FALSE, row.names = FALSE)
  
  
  #getting taxonomic data into a dataframe for Levels 7, 5, 4, 3 
  
  ## SPECIES 
  
  species <- 
    sampleID_table_all %>% 
    select(sampleID, Species) %>% 
    mutate(Species = gsub("_[0-9]*$", '', Species)) %>% #remove _num from species name
    mutate(Species = gsub("_[A-Z]$", '', Species)) %>% #remove _alpha from species name 
    mutate(Species = gsub("_[A-Z]_", ' ', Species)) %>%      ## remove _alpha_ codes from genus
    mutate(Species = gsub("_[A-Z][A-Z] ", ' ', Species)) %>%      ## remove _double alpha codes from genus #March2025 
    mutate(Species = gsub("_[A-Z] ", ' ', Species)) %>%      ## remove _alpha codes from genus 
    mutate(Species = gsub("_[A-Z]_[0-9]* ", ' ', Species)) %>%  #in genus remove alpha fol by num
    mutate(Species = gsub(" [0-9]* ", ' ', Species)) %>% # in genus remove num codes of any length
    mutate(Species = gsub("\\[", '', Species)) %>% 
    mutate(Species = gsub("\\]", '', Species)) %>% 
    mutate(Species = gsub("s__$", NA, Species)) %>%
    mutate(Species = gsub("s_$", NA, Species)) %>%
    mutate(Species = gsub("s__", '', Species)) %>% 
    mutate(Species = gsub("^ ", '', Species)) %>% #remove left hand space
    replace_na(list(Species = "Other (species classification not available)")) %>% #to replace NA use replace_na not mutate
    group_by(Species) %>% 
    summarise(RA = sum(.data[[sampleID]])) %>% 
    arrange(desc(RA))
  
  species_ref <- aggregate(RA ~ Species, species , sum)
  
  
  i7 <- species_ref$RA <= 0.001 #this gets only species with >0.1% abundance (in case a plot is needed and there is limited space)
  species_ref$Species[i7] <- "Other"
  species_ag <- aggregate(RA ~ Species, species_ref, sum)
  
  species_ag$Percentage <- paste0(species_ag$Species,' [',100*round(species_ag$RA,10),'%',']')
  
  #the code below shows all species present in a given sample: which is what needs to be exported 
  all_species <- aggregate(RA ~ Species, species, sum) %>%   ## dataframe to export 
    arrange(desc(RA)) %>% 
    mutate(Percentage = round(RA*100,6)) %>% 
    select(Species, Percentage)
  
  
  all_species <- if (language == "ENG") {
    all_species <- all_species 
  } else {
    all_species <- all_species %>% 
      mutate(Species = gsub("Other", "Les autres bactéries", Species)) %>% 
      mutate(Species = gsub("species classification not available", "la classification des espèces n'est pas disponible", Species)) 
    
  }  
  
  
  
  
  
  table_name_all_species <- paste(output_dir,"/",inputFilePrefix,"_all_species.txt", sep = "")
  
  write.table(all_species, table_name_all_species, sep = "\t", quote = FALSE, row.names = FALSE)
  
  
  
  ## FAMILY 
  
  family <- 
    sampleID_table_all %>% 
    select(sampleID, Family) %>% 
    mutate(Family = gsub("_[A-Z]$", '', Family)) %>%      
    mutate(Family = gsub("_[A-Z]_.*$", '', Family)) %>%  
    mutate(Family = gsub("_[0-9].*$", '', Family)) %>% 
    group_by(Family)  %>%
    summarise(RA = sum(.data[[sampleID]])) %>% 
    mutate_at("Family", str_replace, " ", "") %>% 
    mutate_at("Family", str_replace, "\\[", "") %>% 
    mutate_at("Family", str_replace, "\\]", "") %>% 
    mutate_at("Family", str_replace, "f__$", "Other") %>%
    mutate_at("Family", str_replace, ".*NA", "Other") %>% 
    mutate_at("Family", str_replace, "f__", "")
  
  
  
  all_family <- aggregate(RA ~ Family, family, sum) %>%    ## dataframe to export
    arrange(desc(RA)) %>% 
    mutate(Percentage = round(RA*100,6)) %>% 
    select(Family, Percentage)
  
  ## CLASS
  
  
  class_level <- 
    sampleID_table_all %>% 
    select(sampleID, Class) %>% 
    group_by(Class)  %>%
    summarise(RA = sum(.data[[sampleID]])) %>% 
    mutate_at("Class", str_replace, " ", "") %>% 
    mutate_at("Class", str_replace, "\\[", "") %>% 
    mutate_at("Class", str_replace, "\\]", "") %>% 
    mutate_at("Class", str_replace, "c__$", "Other") %>%
    mutate_at("Class", str_replace, ".*NA", "Other") %>% 
    mutate_at("Class", str_replace, "c__", "")
  
  all_class <- aggregate(RA ~ Class, class_level, sum) %>%     ## dataframe to export
    arrange(desc(RA)) %>% 
    mutate(Percentage = round(RA*100,6)) %>% 
    select(Class, Percentage)
  
  
  
  ## ORDER
  
  order_level <- 
    sampleID_table_all %>% 
    select(sampleID, Order) %>% 
    group_by(Order)  %>%
    summarise(RA = sum(.data[[sampleID]])) %>% 
    mutate_at("Order", str_replace, " ", "") %>% 
    mutate_at("Order", str_replace, "\\[", "") %>% 
    mutate_at("Order", str_replace, "\\]", "") %>% 
    mutate_at("Order", str_replace, "o__$", "Other") %>%
    mutate_at("Order", str_replace, ".*NA", "Other") %>% 
    mutate_at("Order", str_replace, "o__", "")
  
  all_order <- aggregate(RA ~ Order, order_level, sum) %>%     ## dataframe to export
    arrange(desc(RA)) %>% 
    mutate(Percentage = round(RA*100,6)) %>% 
    select(Order, Percentage)
  
  
  table_name_all_taxonomy_xls <- paste(output_dir,"/",inputFilePrefix,"_taxonomy.xlsx", sep = "")
  
  sheets <- list("Phylum" = all_phylum, "Class" = all_class, "Order" = all_order, "Family" = all_family, "Genus" = all_genus, "Species" = all_species) #assume sheet1 and sheet2 are data frames
  write_xlsx(sheets, table_name_all_taxonomy_xls)
}













getGutPotentialScores_v3_GG2 <- function(inputFilePrefix, ReferenceGenusTable = AGP_biomtable_Genus_clean_agg_GG2, Template = "Gutcheck", ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus) {
  
  
  output_dir <- paste("./",inputFilePrefix,"/GutPotential", sep="")
  
  dir.create(output_dir)
  
  Template <- paste(Template)
  
  sampleID <- paste(inputFilePrefix)  
  all_genus <- createGenusTable(inputFilePrefix, biomtable)
  
  genus_ref <- aggregate(RA ~ Genus, all_genus, sum)
  genus_ag <- aggregate(RA ~ Genus, genus_ref, sum)
  
  
  genus_ag$Percentage <- paste0(genus_ag$Genus,' [',100*round(genus_ag$RA,10),'%',']')
  
  
  
  ##___ FUNCTIONS used within this function _____
  
  make_rainbow_gauge_plot <- function(InputDataFrame, InputScore, InputLabel, InputColorGroup) {
    
    #xmax=1.3 changed from xmax=2 on March 1st, 2024
    #color hex codes also changed on March 1st, 2024
    #geom_text(aes(x = 0, y = 0, label = label, colour=group), size=12) size changed from 7 to 12 
    gauge_plot <- ggplot(InputDataFrame, aes(fill = InputColorGroup, ymax = InputScore*2, ymin = 0, xmax = 1.3, xmin = 1)) +
      geom_rect(aes(ymax=2, ymin=0, xmax=1.3, xmin=1), fill = ifelse(InputScore <0.33, "#FF7F7F",
                                                                     ifelse(InputScore>=0.33 & InputScore<0.65, "#fad8ac","#c0edbb"))) +
      geom_rect() + 
      coord_polar(theta = "y",  start=0) + xlim(c(0, 2.5)) + ylim(c(0,2)) +
      geom_text(aes(x = 0, y = 0, label = InputLabel, colour=InputColorGroup), size=12)  + 
      theme_void() +
      scale_fill_manual(values = c("red"="#ff1616", "orange"="#e7c536", "green"="#7dd956")) +
      scale_colour_manual(values = c("red"="#ff1616", "orange"="#e7c536", "green"="#7dd956")) +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank()) +
      guides(fill=FALSE) +
      guides(colour=FALSE) 
    
    return(gauge_plot)
    
  }
  
  
  
  
  make_rainbow_gauge_plot_short <- function(InputDataFrame, InputScore, InputLabel, InputColorGroup, FigureTitleLabel) {
    
    
    gauge_plot_short <- ggplot(InputDataFrame, aes(fill = InputColorGroup, ymax = InputScore*2, ymin = 0, xmax = 1.3, xmin = 1)) +
      geom_rect(aes(ymax=2, ymin=0, xmax=1.3, xmin=1), fill =  "white") +
      geom_rect() + 
      coord_polar(theta = "y",  start=0) + xlim(c(0, 2.5)) + ylim(c(0,2)) +
      geom_text(aes(x = 0, y = 0, label = InputLabel, colour=group), size=17, colour = "white", fontface = "bold.italic")  + 
      geom_text(aes(x=1.7, y=1, label= FigureTitleLabel), size=9, colour = "white", fontface = "bold") +
      #facet_wrap(~title, ncol = 5) +
      theme_void() +
      scale_fill_manual(values = c("red"="#ff1616", "orange"="#e7c536", "green"="#7dd956")) +
      scale_colour_manual(values = c("red"="#ff1616", "orange"="#e7c536", "green"="#7dd956")) +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank()) +
      guides(fill=FALSE) +
      guides(colour=FALSE) 
    # theme(plot.background = element_rect(fill = "#130083"))
    
    return(gauge_plot_short)
    
  }
  
  
  make_slider_score_plot <- function(InputScore, InputScoreLabel) {
    
    slider <- data.frame("values"=seq(0, 1, by = 0.01))
    
    slider_plot <- ggplot(slider, aes(x=values, y = 0))  +
      geom_line(aes(color = values),  size=12, lineend = "round") + 
      geom_point(mapping = aes(x=InputScore),size=11, shape = 16, color = "white") +
      scale_color_gradient(low = "light blue", high = "dark blue") + 
      theme_void() +
      xlab("") +
      ylab("") +
      scale_y_continuous(limits = c(-1, 1), breaks = 1) + 
      annotate('text', x = InputScore, y = 0, label = InputScoreLabel, color =  "black", fontface = 'italic', size = 3.2) +
      theme(legend.position="none")
    
    if (Template == "Gutcheck") {
      slider_plot <- slider_plot 
      
    } else {
      slider_plot <- slider_plot +   scale_color_gradient(low = "#e0c9f0", high = "#662C8F")   
      
    }
    
    return(slider_plot)
    
  }
  
  
  
  
  
  
  ### Gut Potential Sections indicated below
  
  #______________SCFA_Production_________________________________
  
  
  Genus_Butyrate_RAup <- c("Blautia", "Roseburia",
                           "Faecalibacterium",
                           "Anaerostipes",
                           "Clostridium",
                           "Eubacterium",
                           "Coprococcus",
                           "Ruminococcus",
                           "Agathobaculum", #added to DB April 2024
                           "Gemmiger",         #new
                           "Anaerobutyricum", #new
                           "Agathobacter")    #new
  
  #removed: Subdoligranulum and Actinomyces
  
  Genus_Butyrate_RAdown <- c("")
  
  #Calculate score at Genus level  
  
  
  total_genus_associated_with_Butyrate <- calculateProcessScoreGenus(inputFilePrefix, Genus_Butyrate_RAup, Genus_Butyrate_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  Butyrate_score_column <- total_genus_associated_with_Butyrate*100
  Butyrate_df <- data.frame(Butyrate_score_column) %>% 
    mutate(group=ifelse(Butyrate_score_column <33, "red", ifelse(Butyrate_score_column>=33 & Butyrate_score_column<65, "orange","green")),
           label=paste0(round(Butyrate_score_column,0), "%"), 
           Butyrate_col=paste0(round(Butyrate_score_column,2)))
  
  #print(Butyrate_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_ButyrateProduction.txt", sep = "")
  write.table(Butyrate_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  GP_slider_plot <- make_slider_score_plot(Butyrate_df$Butyrate_score_column/100, Butyrate_df$label)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_ButyrateProduction.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  
  Genus_Propionate_RAup <- c("Bacteroides",
                             "Blautia",
                             "Roseburia",
                             "Clostridium",
                             "Eubacterium",
                             "Coprococcus",
                             "Prevotella",
                             "Alistipes",
                             "Phascolarctobacterium",
                             "Parabacteroides",
                             "Akkermansia",
                             "Lactobacillus",
                             "Veillonella",
                             "Dialister",
                             "Ruminococcus")
  
  Genus_Propionate_RAdown <- c("") 
  
  
  total_genus_associated_with_Propionate <- calculateProcessScoreGenus(inputFilePrefix, Genus_Propionate_RAup, Genus_Propionate_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  Propionate_score_column <- total_genus_associated_with_Propionate*100
  Propionate_df <- data.frame(Propionate_score_column) %>% 
    mutate(group=ifelse(Propionate_score_column <33, "red", ifelse(Propionate_score_column>=33 & Propionate_score_column<65, "orange","green")),
           label=paste0(round(Propionate_score_column,0), "%"), 
           Propionate_col=paste0(round(Propionate_score_column,2)))
  
  #print(Propionate_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_PropionateProduction.txt", sep = "")
  write.table(Propionate_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  GP_slider_plot <- make_slider_score_plot(Propionate_df$Propionate_score_column/100, Propionate_df$label)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_PropionateProduction.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  Genus_Acetate_RAup <- c("Bacteroides",
                          "Blautia",
                          "Roseburia",
                          "Clostridium",
                          "Eubacterium",
                          "Coprococcus",
                          "Bifidobacterium",
                          "Prevotella",
                          "Akkermansia",
                          "Parabacteroides",
                          "Alistipes",
                          "Dorea",
                          "Lactobacillus", 
                          "Veillonella",
                          "Streptococcus",
                          "Ruminococcus", #new
                          "Agathobacter")  #new
  
  
  
  Genus_Acetate_RAdown <- c("")
  
  
  total_genus_associated_with_Acetate <- calculateProcessScoreGenus(inputFilePrefix, Genus_Acetate_RAup, Genus_Acetate_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  Acetate_score_column <- total_genus_associated_with_Acetate*100
  Acetate_df <- data.frame(Acetate_score_column) %>% 
    mutate(group=ifelse(Acetate_score_column <33, "red", ifelse(Acetate_score_column>=33 & Acetate_score_column<65, "orange","green")),
           label=paste0(round(Acetate_score_column,0), "%"), 
           Acetate_col=paste0(round(Acetate_score_column,2)))
  
  #print(Acetate_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_AcetateProduction.txt", sep = "")
  write.table(Acetate_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  GP_slider_plot <-make_slider_score_plot(Acetate_df$Acetate_score_column/100, Acetate_df$label)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_AcetateProduction.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  
  Genus_SCFA_RAup <- c(Genus_Acetate_RAup, Genus_Butyrate_RAup, Genus_Propionate_RAup) #to keep only unique elements from each vector can use purrr:reduce(list(a,b,c)), union) but that will change how bacteria are weighted
  
  Genus_SCFA_RAdown <- c("")
  
  total_genus_associated_with_SCFA <- calculateProcessScoreGenus(inputFilePrefix, Genus_SCFA_RAup, Genus_SCFA_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  SCFA_score_column <- total_genus_associated_with_SCFA*100
  SCFA_df <- data.frame(SCFA_score_column) %>% 
    mutate(group=ifelse(SCFA_score_column <33, "red", ifelse(SCFA_score_column>=33 & SCFA_score_column<65, "orange","green")),
           label=paste0(round(SCFA_score_column,0), "%"), 
           SCFA_col=paste0(round(SCFA_score_column,2)))
  
  SCFA_label <- SCFA_df$label
  SCFA_group <- SCFA_df$group
  
  #print(SCFA_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_SCFAProduction.txt", sep = "")
  write.table(SCFA_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  make_rainbow_gauge_plot(SCFA_df, SCFA_score_column/100, SCFA_label, SCFA_group)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_SCFAProduction.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  
  #______________Metabolic_Health________________________
  
  
  Genus_Lactose_RAup <- c("Bifidobacterium", "Lactobacillus","Lactococcus", "Faecalibacterium", "Roseburia", "Agathobacter")
  
  
  
  Genus_Lactose_RAdown <- c("")
  
  
  total_genus_associated_with_Lactose <- calculateProcessScoreGenus(inputFilePrefix, Genus_Lactose_RAup, Genus_Lactose_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  Lactose_score_column <- total_genus_associated_with_Lactose*100
  Lactose_df <- data.frame(Lactose_score_column) %>% 
    mutate(group=ifelse(Lactose_score_column <33, "red", ifelse(Lactose_score_column>=33 & Lactose_score_column<65, "orange","green")),
           label=paste0(round(Lactose_score_column,0), "%"), 
           Lactose_col=paste0(round(Lactose_score_column,2)))
  
  #print(Lactose_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_LactoseProduction.txt", sep = "")
  write.table(Lactose_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  GP_slider_plot <-make_slider_score_plot(Lactose_df$Lactose_score_column/100, Lactose_df$label)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_LactoseMetabolism.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  Genus_Carbohydrate_RAup <- c("Akkermansia",
                               "Alistipes",
                               "Anaerostipes",
                               "Bacteroides",
                               "Bifidobacterium",
                               "Blautia",
                               "Clostridium",
                               "Coprococcus",
                               "Dialister",
                               "Dorea",
                               "Eubacterium",
                               "Faecalibacterium",
                               "Lactobacillus",
                               "Parabacteroides",
                               "Phascolarctobacterium",
                               "Prevotella",
                               "Roseburia",
                               "Ruminococcus",
                               "Streptococcus",
                               "Anaerobutyricum")
  
  Genus_Carbohydrate_RAdown <- c("")
  
  
  total_genus_associated_with_Carbohydrate <- calculateProcessScoreGenus(inputFilePrefix, Genus_Carbohydrate_RAup, Genus_Carbohydrate_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  Carbohydrate_score_column <- total_genus_associated_with_Carbohydrate*100
  Carbohydrate_df <- data.frame(Carbohydrate_score_column) %>% 
    mutate(group=ifelse(Carbohydrate_score_column <33, "red", ifelse(Carbohydrate_score_column>=33 & Carbohydrate_score_column<65, "orange","green")),
           label=paste0(round(Carbohydrate_score_column,0), "%"), 
           Carbohydrate_col=paste0(round(Carbohydrate_score_column,2)))
  
  #print(Carbohydrate_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_CarbohydrateProduction.txt", sep = "")
  write.table(Carbohydrate_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  GP_slider_plot <-make_slider_score_plot(Carbohydrate_df$Carbohydrate_score_column/100, Carbohydrate_df$label)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_CarbohydrateMetabolism.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  Genus_Protein_RAup <- c("Bacillus",
                          "Clostridium",
                          "Staphylococcus",
                          "Streptococcus",
                          "Bacteroides") #new addition
  
  
  
  Genus_Protein_RAdown <- c("")
  
  
  total_genus_associated_with_Protein <- calculateProcessScoreGenus(inputFilePrefix, Genus_Protein_RAup, Genus_Protein_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  Protein_score_column <- total_genus_associated_with_Protein*100
  Protein_df <- data.frame(Protein_score_column) %>% 
    mutate(group=ifelse(Protein_score_column <33, "red", ifelse(Protein_score_column>=33 & Protein_score_column<65, "orange","green")),
           label=paste0(round(Protein_score_column,0), "%"), 
           Protein_col=paste0(round(Protein_score_column,2)))
  
  #print(Protein_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_ProteinProduction.txt", sep = "")
  write.table(Protein_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  GP_slider_plot <-make_slider_score_plot(Protein_df$Protein_score_column/100, Protein_df$label)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_ProteinMetabolism.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  Genus_Gluten_RAup <- c(
    "Bacillus",
    "Bifidobacterium",
    "Clostridium",
    "Lactobacillus",
    "Prevotella",
    "Staphylococcus",
    "Stenotrophomonas",
    "Streptococcus"
  )
  
  Genus_Gluten_RAdown <- c("")
  
  
  total_genus_associated_with_Gluten <- calculateProcessScoreGenus(inputFilePrefix, Genus_Gluten_RAup, Genus_Gluten_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  Gluten_score_column <- total_genus_associated_with_Gluten*100
  Gluten_df <- data.frame(Gluten_score_column) %>% 
    mutate(group=ifelse(Gluten_score_column <33, "red", ifelse(Gluten_score_column>=33 & Gluten_score_column<65, "orange","green")),
           label=paste0(round(Gluten_score_column,0), "%"), 
           Gluten_col=paste0(round(Gluten_score_column,2)))
  
  #print(Gluten_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_GlutenProduction.txt", sep = "")
  write.table(Gluten_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  GP_slider_plot <-make_slider_score_plot(Gluten_df$Gluten_score_column/100, Gluten_df$label)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_GlutenMetabolism.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  Genus_Fat_RAup <- c("Anaerostipes",
                      "Akkermansia",
                      "Alistipes",
                      "Bacteroides",
                      "Bifidobacterium",
                      "Blautia",
                      "Clostridium",
                      "Coprococcus",
                      "Dialister",
                      "Dorea",
                      "Eubacterium",
                      "Faecalibacterium",
                      "Lactobacillus",
                      "Parabacteroides",
                      "Phascolarctobacterium",
                      "Prevotella",
                      "Roseburia",
                      "Ruminococcus",
                      "Streptococcus",
                      "Agathobaculum", #new
                      "Gemmiger",       #new
                      "Anaerobutyricum", #new
                      "Agathobacter")   #new
  
  #Escherichia" removed, as it is filtered from AGP
  #Subdoligranulum removed not in updated taxonomy
  #"Actinomyces" removed: low frequency
  #"Veillonella" removed, no references 
  
  
  Genus_Fat_RAdown <- c("")
  
  
  total_genus_associated_with_Fat <- calculateProcessScoreGenus(inputFilePrefix, Genus_Fat_RAup, Genus_Fat_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  Fat_score_column <- total_genus_associated_with_Fat*100
  Fat_df <- data.frame(Fat_score_column) %>% 
    mutate(group=ifelse(Fat_score_column <33, "red", ifelse(Fat_score_column>=33 & Fat_score_column<65, "orange","green")),
           label=paste0(round(Fat_score_column,0), "%"), 
           Fat_col=paste0(round(Fat_score_column,2)))
  
  #print(Fat_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_FatProduction.txt", sep = "")
  write.table(Fat_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  GP_slider_plot <-make_slider_score_plot(Fat_df$Fat_score_column/100, Fat_df$label)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_FatMetabolism.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  Genus_OverallMetabolism_RAup <- c(Genus_Lactose_RAup, Genus_Fat_RAup, Genus_Carbohydrate_RAup, Genus_Protein_RAup, Genus_Gluten_RAup)
  Genus_OverallMetabolism_RAdown <- c("")
  
  total_genus_associated_with_OverallMetabolism <- calculateProcessScoreGenus(inputFilePrefix, Genus_OverallMetabolism_RAup, Genus_OverallMetabolism_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  OverallMetabolism_score_column <- total_genus_associated_with_OverallMetabolism*100
  OverallMetabolism_df <- data.frame(OverallMetabolism_score_column) %>% 
    mutate(group=ifelse(OverallMetabolism_score_column <33, "red", ifelse(OverallMetabolism_score_column>=33 & OverallMetabolism_score_column<65, "orange","green")),
           label=paste0(round(OverallMetabolism_score_column,0), "%"), 
           OverallMetabolism_col=paste0(round(OverallMetabolism_score_column,2)))
  
  OverallMetabolism_label <- OverallMetabolism_df$label
  OverallMetabolism_group <- OverallMetabolism_df$group
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_OverallMetabolismProduction.txt", sep = "")
  write.table(OverallMetabolism_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  
  make_rainbow_gauge_plot(OverallMetabolism_df, OverallMetabolism_score_column/100, OverallMetabolism_label, OverallMetabolism_group)
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_OverallMetabolism.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  make_rainbow_gauge_plot_short(OverallMetabolism_df, OverallMetabolism_score_column/100, OverallMetabolism_label, OverallMetabolism_group, "")
  plotname_GutProcess_short <- paste(output_dir,"/",inputFilePrefix,"_OverallMetabolism_shortReport.png", sep = "")
  ggsave(plotname_GutProcess_short, plot = last_plot())
  
  
  Genus_VitaminB1_RAup <-  c(
    "Bacteroides",
    "Prevotella",
    "Clostridium",
    "Lactobacillus",
    "Ruminococcus",
    "Bifidobacterium",
    "Fusobacterium"
  )
  Genus_VitaminB1_RAdown <- c("")
  
  
  total_genus_associated_with_VitaminB1 <- calculateProcessScoreGenus(inputFilePrefix, Genus_VitaminB1_RAup, Genus_VitaminB1_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  VitaminB1_score_column <- total_genus_associated_with_VitaminB1*100
  VitaminB1_df <- data.frame(VitaminB1_score_column) %>%
    mutate(group=ifelse(VitaminB1_score_column <33, "red", ifelse(VitaminB1_score_column>=33 & VitaminB1_score_column<65, "orange","green")),
           label=paste0(round(VitaminB1_score_column,0), "%"),
           VitaminB1_col=paste0(round(VitaminB1_score_column,2)))
  #print(VitaminB1_df)
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_VitaminB1Production.txt", sep = "")
  write.table(VitaminB1_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  GP_slider_plot <-make_slider_score_plot(VitaminB1_df$VitaminB1_score_column/100, VitaminB1_df$label)
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_VitaminB1Synthesis.png", sep = "")
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  Genus_VitaminB2_RAup <- c(
    "Bacteroides",
    "Prevotella",
    "Clostridium",
    "Lactobacillus",
    "Ruminococcus"
  )
  Genus_VitaminB2_RAdown <- c("")
  total_genus_associated_with_VitaminB2 <- calculateProcessScoreGenus(inputFilePrefix, Genus_VitaminB2_RAup, Genus_VitaminB2_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  VitaminB2_score_column <- total_genus_associated_with_VitaminB2*100
  VitaminB2_df <- data.frame(VitaminB2_score_column) %>%
    mutate(group=ifelse(VitaminB2_score_column <33, "red", ifelse(VitaminB2_score_column>=33 & VitaminB2_score_column<65, "orange","green")),
           label=paste0(round(VitaminB2_score_column,0), "%"),
           VitaminB2_col=paste0(round(VitaminB2_score_column,2)))
  #print(VitaminB2_df)
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_VitaminB2Production.txt", sep = "")
  write.table(VitaminB2_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  GP_slider_plot <-make_slider_score_plot(VitaminB2_df$VitaminB2_score_column/100, VitaminB2_df$label)
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_VitaminB2Synthesis.png", sep = "")
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  Genus_VitaminB3_RAup <- c(
    "Bacteroides",
    "Prevotella",
    "Clostridium",
    "Ruminococcus",
    "Bifidobacterium",
    "Fusobacterium"
  )
  #"Helicobacter" removed not present in AGP
  Genus_VitaminB3_RAdown <- c("")
  total_genus_associated_with_VitaminB3 <- calculateProcessScoreGenus(inputFilePrefix, Genus_VitaminB3_RAup, Genus_VitaminB3_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  VitaminB3_score_column <- total_genus_associated_with_VitaminB3*100
  VitaminB3_df <- data.frame(VitaminB3_score_column) %>%
    mutate(group=ifelse(VitaminB3_score_column <33, "red", ifelse(VitaminB3_score_column>=33 & VitaminB3_score_column<65, "orange","green")),
           label=paste0(round(VitaminB3_score_column,0), "%"),
           VitaminB3_col=paste0(round(VitaminB3_score_column,2)))
  
  #print(VitaminB3_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_VitaminB3Production.txt", sep = "")
  
  write.table(VitaminB3_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  GP_slider_plot <-make_slider_score_plot(VitaminB3_df$VitaminB3_score_column/100, VitaminB3_df$label)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_VitaminB3Synthesis.png", sep = "")
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  Genus_VitaminB5_RAup <- c("Bacteroides", "Prevotella", "Ruminococcus")
  Genus_VitaminB5_RAdown <- c("")
  total_genus_associated_with_VitaminB5 <- calculateProcessScoreGenus(inputFilePrefix, Genus_VitaminB5_RAup, Genus_VitaminB5_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  VitaminB5_score_column <- total_genus_associated_with_VitaminB5*100
  VitaminB5_df <- data.frame(VitaminB5_score_column) %>%
    mutate(group=ifelse(VitaminB5_score_column <33, "red", ifelse(VitaminB5_score_column>=33 & VitaminB5_score_column<65, "orange","green")),
           label=paste0(round(VitaminB5_score_column,0), "%"),
           VitaminB5_col=paste0(round(VitaminB5_score_column,2)))
  #print(VitaminB5_df)
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_VitaminB5Production.txt", sep = "")
  write.table(VitaminB5_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  GP_slider_plot <-make_slider_score_plot(VitaminB5_df$VitaminB5_score_column/100, VitaminB5_df$label)
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_VitaminB5Synthesis.png", sep = "")
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  Genus_VitaminB6_RAup <- c("Bacteroides",
                            "Prevotella",
                            "Bifidobacterium",
                            "Collinsella")
  Genus_VitaminB6_RAdown <- c("")
  total_genus_associated_with_VitaminB6 <- calculateProcessScoreGenus(inputFilePrefix, Genus_VitaminB6_RAup, Genus_VitaminB6_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  VitaminB6_score_column <- total_genus_associated_with_VitaminB6*100
  VitaminB6_df <- data.frame(VitaminB6_score_column) %>%
    mutate(group=ifelse(VitaminB6_score_column <33, "red", ifelse(VitaminB6_score_column>=33 & VitaminB6_score_column<65, "orange","green")),
           label=paste0(round(VitaminB6_score_column,0), "%"),
           VitaminB6_col=paste0(round(VitaminB6_score_column,2)))
  #print(VitaminB6_df)
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_VitaminB6Production.txt", sep = "")
  write.table(VitaminB6_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  GP_slider_plot <-make_slider_score_plot(VitaminB6_df$VitaminB6_score_column/100, VitaminB6_df$label)
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_VitaminB6Synthesis.png", sep = "")
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  Genus_VitaminB7_RAup <- c("Bacteroides",
                            "Lactobacillus",
                            "Fusobacterium",
                            "Campylobacter")
  Genus_VitaminB7_RAdown <- c("")
  total_genus_associated_with_VitaminB7 <- calculateProcessScoreGenus(inputFilePrefix, Genus_VitaminB7_RAup, Genus_VitaminB7_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  VitaminB7_score_column <- total_genus_associated_with_VitaminB7*100
  VitaminB7_df <- data.frame(VitaminB7_score_column) %>%
    mutate(group=ifelse(VitaminB7_score_column <33, "red", ifelse(VitaminB7_score_column>=33 & VitaminB7_score_column<65, "orange","green")),
           label=paste0(round(VitaminB7_score_column,0), "%"),
           VitaminB7_col=paste0(round(VitaminB7_score_column,2)))
  
  
  #print(VitaminB7_df)
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_VitaminB7Production.txt", sep = "")
  write.table(VitaminB7_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  GP_slider_plot <-make_slider_score_plot(VitaminB7_df$VitaminB7_score_column/100, VitaminB7_df$label)
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_VitaminB7Synthesis.png", sep = "")
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  
  Genus_VitaminB9_RAup <- c("Bacteroides",
                            "Prevotella",
                            "Clostridium",
                            "Lactobacillus",
                            "Lactococcus",
                            "Streptococcus",
                            "Bifidobacterium",
                            "Fusobacterium"
  )
  #Salmonella removed
  Genus_VitaminB9_RAdown <- c("")
  total_genus_associated_with_VitaminB9 <- calculateProcessScoreGenus(inputFilePrefix, Genus_VitaminB9_RAup, Genus_VitaminB9_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  VitaminB9_score_column <- total_genus_associated_with_VitaminB9*100
  VitaminB9_df <- data.frame(VitaminB9_score_column) %>%
    mutate(group=ifelse(VitaminB9_score_column <33, "red", ifelse(VitaminB9_score_column>=33 & VitaminB9_score_column<65, "orange","green")),
           label=paste0(round(VitaminB9_score_column,0), "%"),
           VitaminB9_col=paste0(round(VitaminB9_score_column,2)))
  #print(VitaminB9_df)
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_VitaminB9Production.txt", sep = "")
  write.table(VitaminB9_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  GP_slider_plot <-make_slider_score_plot(VitaminB9_df$VitaminB9_score_column/100, VitaminB9_df$label)
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_VitaminB9Synthesis.png", sep = "")
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  
  Genus_VitaminB12_RAup <- c("Bacteroides",
                             "Bacillus",
                             "Prevotella",
                             "Clostridium",
                             "Faecalibacterium",
                             "Ruminococcus",
                             "Lactobacillus",
                             "Bifidobacterium",
                             "Fusobacterium")
  Genus_VitaminB12_RAdown <- c("")
  total_genus_associated_with_VitaminB12 <- calculateProcessScoreGenus(inputFilePrefix, Genus_VitaminB12_RAup, Genus_VitaminB12_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  VitaminB12_score_column <- total_genus_associated_with_VitaminB12*100
  VitaminB12_df <- data.frame(VitaminB12_score_column) %>%
    mutate(group=ifelse(VitaminB12_score_column <33, "red", ifelse(VitaminB12_score_column>=33 & VitaminB12_score_column<65, "orange","green")),
           label=paste0(round(VitaminB12_score_column,0), "%"),
           VitaminB12_col=paste0(round(VitaminB12_score_column,2)))
  #print(VitaminB12_df)
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_VitaminB12Production.txt", sep = "")
  write.table(VitaminB12_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  GP_slider_plot <-make_slider_score_plot(VitaminB12_df$VitaminB12_score_column/100, VitaminB12_df$label)
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_VitaminB12Synthesis.png", sep = "")
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  
  Genus_VitaminK2_RAup <- c("Bacillus",
                            "Bacteroides",
                            "Eubacterium",
                            "Veillonella",
                            "Enterococcus",
                            "Lactococcus")
  #Propionibacterium (renamed to Cutibacterium) removed due to low freq in AGP
  #Mycobacterium removed due to low freq in AGP
  #Escherichia and Enterobacteria not present in AGP dataset after bloom filtering
  Genus_VitaminK2_RAdown <- c("")
  total_genus_associated_with_VitaminK2 <- calculateProcessScoreGenus(inputFilePrefix, Genus_VitaminK2_RAup, Genus_VitaminK2_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  
  VitaminK2_score_column <- total_genus_associated_with_VitaminK2*100
  VitaminK2_df <- data.frame(VitaminK2_score_column) %>%
    mutate(group=ifelse(VitaminK2_score_column <33, "red", ifelse(VitaminK2_score_column>=33 & VitaminK2_score_column<65, "orange","green")),
           label=paste0(round(VitaminK2_score_column,0), "%"),
           VitaminK2_col=paste0(round(VitaminK2_score_column,2)))
  
  #print(VitaminK2_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_VitaminK2Production.txt", sep = "")
  write.table(VitaminK2_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  GP_slider_plot <-make_slider_score_plot(VitaminK2_df$VitaminK2_score_column/100, VitaminK2_df$label)
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_VitaminK2Synthesis.png", sep = "")
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  Genus_OverallVitaminSynthesis_RAup <- c(Genus_VitaminB1_RAup, Genus_VitaminB2_RAup, Genus_VitaminB3_RAup, Genus_VitaminB5_RAup, Genus_VitaminB6_RAup, Genus_VitaminB7_RAup, Genus_VitaminB9_RAup, Genus_VitaminB12_RAup, Genus_VitaminK2_RAup)
  Genus_OverallVitaminSynthesis_RAdown <- c("")
  
  
  total_genus_associated_with_OverallVitaminSynthesis <- calculateProcessScoreGenus(inputFilePrefix, Genus_OverallVitaminSynthesis_RAup, Genus_OverallVitaminSynthesis_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  OverallVitaminSynthesis_score_column <- total_genus_associated_with_OverallVitaminSynthesis*100
  OverallVitaminSynthesis_df <- data.frame(OverallVitaminSynthesis_score_column) %>%
    mutate(group=ifelse(OverallVitaminSynthesis_score_column <33, "red", ifelse(OverallVitaminSynthesis_score_column>=33 & OverallVitaminSynthesis_score_column<65, "orange","green")),
           label=paste0(round(OverallVitaminSynthesis_score_column,0), "%"),
           OverallVitaminSynthesis_col=paste0(round(OverallVitaminSynthesis_score_column,2)))
  
  OverallVitaminSynthesis_label <- OverallVitaminSynthesis_df$label
  OverallVitaminSynthesis_group <- OverallVitaminSynthesis_df$group
  
  #print(OverallVitaminSynthesis_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_OverallVitaminSynthesisProduction.txt", sep = "")
  write.table(OverallVitaminSynthesis_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  make_rainbow_gauge_plot(OverallVitaminSynthesis_df, OverallVitaminSynthesis_score_column/100, OverallVitaminSynthesis_label, OverallVitaminSynthesis_group)
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_OverallVitaminSynthesis.png", sep = "")
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  Genus_OverallImmuneScore_RAup <- c(Genus_VitaminB1_RAup, Genus_VitaminB2_RAup, Genus_VitaminB3_RAup, Genus_VitaminB5_RAup, Genus_VitaminB6_RAup, Genus_VitaminB7_RAup, Genus_VitaminB9_RAup, Genus_VitaminB12_RAup, Genus_VitaminK2_RAup, Genus_Acetate_RAup, Genus_Butyrate_RAup, Genus_Propionate_RAup)
  Genus_OverallImmuneScore_RAdown <- c("")
  
  
  
  total_genus_associated_with_OverallImmuneScore <- calculateProcessScoreGenus(inputFilePrefix, Genus_OverallImmuneScore_RAup, Genus_OverallImmuneScore_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  OverallImmuneScore_score_column <- total_genus_associated_with_OverallImmuneScore*100
  OverallImmuneScore_df <- data.frame(OverallImmuneScore_score_column) %>% 
    mutate(group=ifelse(OverallImmuneScore_score_column <33, "red", ifelse(OverallImmuneScore_score_column>=33 & OverallImmuneScore_score_column<65, "orange","green")),
           label=paste0(round(OverallImmuneScore_score_column,0), "%"), 
           OverallImmuneScore_col=paste0(round(OverallImmuneScore_score_column,2)))
  
  OverallImmuneScore_label <- OverallImmuneScore_df$label
  OverallImmuneScore_group <- OverallImmuneScore_df$group
  
  #print(OverallImmuneScore_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_OverallImmuneScoreProduction.txt", sep = "")
  write.table(OverallImmuneScore_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  make_rainbow_gauge_plot(OverallImmuneScore_df, OverallImmuneScore_score_column/100, OverallImmuneScore_label, OverallImmuneScore_group)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_OverallImmuneScore.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  make_rainbow_gauge_plot_short(OverallImmuneScore_df, OverallImmuneScore_score_column/100, OverallImmuneScore_label, OverallImmuneScore_group, "")
  plotname_GutProcess_short <- paste(output_dir,"/",inputFilePrefix,"_OverallImmuneScore_shortReport.png", sep = "")
  ggsave(plotname_GutProcess_short, plot = last_plot())
  
  
  #_________________________Gut-Brain-Axis______________________
  
  Genus_Serotonin_RAup <- c("Streptococcus", 
                            "Enterococcus", 
                            "Roseburia",
                            "Corynebacterium")
  
  #"Corynebacterium" previously not used due to low freq 
  
  
  Genus_Serotonin_RAdown <- c("")
  
  
  
  total_genus_associated_with_Serotonin <- calculateProcessScoreGenus(inputFilePrefix, Genus_Serotonin_RAup, Genus_Serotonin_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  Serotonin_score_column <- total_genus_associated_with_Serotonin*100
  Serotonin_df <- data.frame(Serotonin_score_column) %>% 
    mutate(group=ifelse(Serotonin_score_column <33, "red", ifelse(Serotonin_score_column>=33 & Serotonin_score_column<65, "orange","green")),
           label=paste0(round(Serotonin_score_column,0), "%"), 
           Serotonin_col=paste0(round(Serotonin_score_column,2)))
  
  #print(Serotonin_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_SerotoninProduction.txt", sep = "")
  write.table(Serotonin_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  GP_slider_plot <-make_slider_score_plot(Serotonin_df$Serotonin_score_column/100, Serotonin_df$label)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_SerotoninSynthesis.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  
  
  Genus_GABA_RAup <- c("Lactobacillus", 
                       "Lactococcus", 
                       "Enterococcus", 
                       "Pediococcus", 
                       "Bifidobacterium", 
                       "Parabacteroides", 
                       "Streptococcus")
  
  Genus_GABA_RAdown <-c("")
  
  Genus_GABA_RAdown <- c("")
  
  
  total_genus_associated_with_GABA <- calculateProcessScoreGenus(inputFilePrefix, Genus_GABA_RAup, Genus_GABA_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  GABA_score_column <- total_genus_associated_with_GABA*100
  GABA_df <- data.frame(GABA_score_column) %>% 
    mutate(group=ifelse(GABA_score_column <33, "red", ifelse(GABA_score_column>=33 & GABA_score_column<65, "orange","green")),
           label=paste0(round(GABA_score_column,0), "%"), 
           GABA_col=paste0(round(GABA_score_column,2)))
  
  #print(GABA_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_GABAProduction.txt", sep = "")
  write.table(GABA_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  GP_slider_plot <-make_slider_score_plot(GABA_df$GABA_score_column/100, GABA_df$label)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_GABASynthesis.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  
  Genus_IPA_RAup <- c("Lactobacillus", 
                      "Clostridium", 
                      "Akkermansia", 
                      "Peptostreptococcus",
                      "Proteus")
  
  Genus_IPA_RAdown <- c("")
  
  
  total_genus_associated_with_IPA <- calculateProcessScoreGenus(inputFilePrefix, Genus_IPA_RAup, Genus_IPA_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  IPA_score_column <- total_genus_associated_with_IPA*100
  IPA_df <- data.frame(IPA_score_column) %>% 
    mutate(group=ifelse(IPA_score_column <33, "red", ifelse(IPA_score_column>=33 & IPA_score_column<65, "orange","green")),
           label=paste0(round(IPA_score_column,0), "%"), 
           IPA_col=paste0(round(IPA_score_column,2)))
  
  #print(IPA_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_IPAProduction.txt", sep = "")
  write.table(IPA_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  GP_slider_plot <-make_slider_score_plot(IPA_df$IPA_score_column/100, IPA_df$label)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_IPASynthesis.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  
  
  Anxiety_RAup <- c("Bacteroides",
                    "Eggerthella", "Ruminococcus") #removed Bifidobacterium --> moved to NEGATIVE association with Anxiety 
  
  Anxiety_RAdown <- c("Bifidobacterium", #NOTE the change 
                      "Faecalibacterium",
                      "Lachnospira",
                      "Sutterella",
                      "Dialister",
                      "Faecalibacterium") 
  #"Ruminococcus", moved to pos.assoc. with anxiety
  
  #Butyricicoccus low frequency, removed
  #"Enterococcus",             #Serotonin 
  #"Corynebacterium",          #Serotonin
  #"Roseburia",                #GABA
  #"Lactococcus",              #GABA
  #"Pediococcus",              #GABA
  #"Parabacteroides"          #GABA
  
  Stress_RAup <- c("Fusobacterium",
                   "Prevotella") #"Streptococcus" was removed - no reference
  
  
  Stress_RAdown <- c("Lactobacillus", "Bifidobacterium", "Streptococcus")
  
  #"Lactococcus", "Enterococcus", "Pediococcus", "Parabacteroides" are GABA producers; redundant 
  
  PoorSleep_RAup <- c("Blautia", "Bacteroides",  "Coprococcus",  "Streptococcus", "Veillonella")
  
  #Escherichia, Enterobacter, Klebsiella not in AGP, removed 
  #"Subdoligranulum" reclassified 
  
  PoorSleep_RAdown <-c("Bifidobacterium", 
                       "Anaerostipes", 
                       "Faecalibacterium", 
                       "Parabacteroides", 
                       "Phascolarctobacterium", 
                       "Prevotella", 
                       "Roseburia",
                       "Agathobacter", 
                       "Lachnospira", 
                       "Agathobaculum", #added to DB April 2024, SCFA 
                       "Gemmiger", #new
                       "Anaerobutyricum", #new, SCFA 
                       "Corynebacterium")    #new
  
  
  # "Ruminococcus", "Clostridium",  "Eubacterium", "Dialister", and "Alistipes" removed, controversial, need more data 
  
  # "Butyricicoccus" has a confirmed reference but its frequency is too low to be used for calculations (average Butyricoccus rel.abun in ref.population is ~0)
  
  
  
  
  MDD_RAup <- c("Eggerthella", "Veillonella", "Streptococcus") ##"Flavonifractor" 
  
  MDD_RAdown <- c("Faecalibacterium", "Coprococcus", "Dialister", "Bifidobacterium", "Lactobacillus", "Haemophilus", "Ruminococcus")
  
  
  
  Genus_OverallGutBrainAxis_RAup <- c(Genus_Serotonin_RAup, Genus_GABA_RAup, Genus_IPA_RAup, Anxiety_RAdown, Stress_RAdown, PoorSleep_RAdown, MDD_RAdown) #Bacteria associated with low Anxiety/Stress, improve GutBrainAxisScore. And bacteria negatively associated with PoorSleep/Insomnia also improve GutBrainAxisScore
  Genus_OverallGutBrainAxis_RAdown <- c(Anxiety_RAup, Stress_RAup, PoorSleep_RAup, MDD_RAup) #Bacteria associated with high Anxiety, decrease GutBrainAxisScore. Similarly, bacteria positively associated with PoorSleep/Insomnia decrease GutBrain AxisScore 
  
  
  
  total_genus_associated_with_OverallGutBrainAxis <- calculateProcessScoreGenus(inputFilePrefix, Genus_OverallGutBrainAxis_RAup, Genus_OverallGutBrainAxis_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  OverallGutBrainAxis_score_column <- total_genus_associated_with_OverallGutBrainAxis*100
  OverallGutBrainAxis_df <- data.frame(OverallGutBrainAxis_score_column) %>% 
    mutate(group=ifelse(OverallGutBrainAxis_score_column <33, "red", ifelse(OverallGutBrainAxis_score_column>=33 & OverallGutBrainAxis_score_column<65, "orange","green")),
           label=paste0(round(OverallGutBrainAxis_score_column,0), "%"), 
           GutBrainAxis_col=paste0(round(OverallGutBrainAxis_score_column,2)))
  
  OverallGutBrainAxis_label <- OverallGutBrainAxis_df$label
  OverallGutBrainAxis_group <- OverallGutBrainAxis_df$group
  
  #print(OverallGutBrainAxis_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_GutBrainAxis.txt", sep = "")
  write.table(OverallGutBrainAxis_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  make_rainbow_gauge_plot(OverallGutBrainAxis_df, OverallGutBrainAxis_score_column/100, OverallGutBrainAxis_label, OverallGutBrainAxis_group)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_OverallGutBrainAxis.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  make_rainbow_gauge_plot_short(OverallGutBrainAxis_df, OverallGutBrainAxis_score_column/100, OverallGutBrainAxis_label,OverallGutBrainAxis_group, "")
  plotname_GutProcess_short <- paste(output_dir,"/",inputFilePrefix,"_OverallGutBrainAxis_shortReport.png", sep = "")
  
  ggsave(plotname_GutProcess_short, plot = last_plot())
  
  
  
  #__________________Physical_Health___________________________
  
  
  
  Genus_MuscleStrength_RAup <- c("Lactobacillus", "Bifidobacterium", "Blautia", "Roseburia", "Faecalibacterium", "Anaerostipes", "Clostridium", "Eubacterium", "Coprococcus",  "Ruminococcus",
                                 "Veillonella", 
                                 "Gemmiger",       #new
                                 "Anaerobutyricum", #new
                                 "Agathobacter")   #new
  
  #"Agathobaculum" needs more data
  
  
  Genus_MuscleStrength_RAdown <- c("Bacteroides")
  
  
  total_genus_associated_with_MuscleStrength <- calculateProcessScoreGenus(inputFilePrefix, Genus_MuscleStrength_RAup, Genus_MuscleStrength_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  MuscleStrength_score_column <- total_genus_associated_with_MuscleStrength*100
  MuscleStrength_df <- data.frame(MuscleStrength_score_column) %>% 
    mutate(group=ifelse(MuscleStrength_score_column <33, "red", ifelse(MuscleStrength_score_column>=33 & MuscleStrength_score_column<65, "orange","green")),
           label=paste0(round(MuscleStrength_score_column,0), "%"), 
           MuscleStrength_col=paste0(round(MuscleStrength_score_column,2)))
  
  #print(MuscleStrength_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_MuscleStrength.txt", sep = "")
  write.table(MuscleStrength_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  GP_slider_plot <-make_slider_score_plot(MuscleStrength_df$MuscleStrength_score_column/100, MuscleStrength_df$label)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_MuscleStrength.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  Genus_AerobicEndurance_RAup <- c("Bifidobacterium", 
                                   "Faecalibacterium",
                                   "Collinsella", 
                                   "Mitsuokella", 
                                   "Blautia", 
                                   "Ruminococcus", 
                                   "Roseburia", 
                                   "Agathobacter", 
                                   "Anaerostipes", 
                                   "Coprococcus", 
                                   "Akkermansia")
  
  #Pseudobutyrivibrio not present in AGP, also removed family level data 
  
  Genus_AerobicEndurance_RAdown <- c("Eubacterium", "Clostridium", "Bacteroides")
  
  
  total_genus_associated_with_AerobicEndurance <- calculateProcessScoreGenus(inputFilePrefix, Genus_AerobicEndurance_RAup, Genus_AerobicEndurance_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  AerobicEndurance_score_column <- total_genus_associated_with_AerobicEndurance*100
  AerobicEndurance_df <- data.frame(AerobicEndurance_score_column) %>% 
    mutate(group=ifelse(AerobicEndurance_score_column <33, "red", ifelse(AerobicEndurance_score_column>=33 & AerobicEndurance_score_column<65, "orange","green")),
           label=paste0(round(AerobicEndurance_score_column,0), "%"), 
           AerobicEndurance_col=paste0(round(AerobicEndurance_score_column,2)))
  
  #print(AerobicEndurance_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_AerobicEndurance.txt", sep = "")
  write.table(AerobicEndurance_df, filename_GutProcess, quote = FALSE, row.names = FALSE,sep = "\t")
  
  GP_slider_plot <-make_slider_score_plot(AerobicEndurance_df$AerobicEndurance_score_column/100, AerobicEndurance_df$label)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_AerobicEndurance.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  
  
  
  Genus_OverallPhysicalHealth_RAup <- c(Genus_MuscleStrength_RAup, Genus_AerobicEndurance_RAup)
  Genus_OverallPhysicalHealth_RAdown <- c(Genus_MuscleStrength_RAdown, Genus_AerobicEndurance_RAdown)
  
  
  
  total_genus_associated_with_OverallPhysicalHealth <- calculateProcessScoreGenus(inputFilePrefix, Genus_OverallPhysicalHealth_RAup, Genus_OverallPhysicalHealth_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  OverallPhysicalHealth_score_column <- total_genus_associated_with_OverallPhysicalHealth*100
  OverallPhysicalHealth_df <- data.frame(OverallPhysicalHealth_score_column) %>% 
    mutate(group=ifelse(OverallPhysicalHealth_score_column <33, "red", ifelse(OverallPhysicalHealth_score_column>=33 & OverallPhysicalHealth_score_column<65, "orange","green")),
           label=paste0(round(OverallPhysicalHealth_score_column,0), "%"), 
           OverallPhysicalHealth_col=paste0(round(OverallPhysicalHealth_score_column,2)))
  
  OverallPhysicalHealth_label <- OverallPhysicalHealth_df$label
  OverallPhysicalHealth_group <- OverallPhysicalHealth_df$group
  
  #print(OverallPhysicalHealth_df)
  
  filename_GutProcess <- paste(output_dir,"/GP_",inputFilePrefix,"_OverallPhysicalHealth.txt", sep = "")
  write.table(OverallPhysicalHealth_df, filename_GutProcess, quote = FALSE, row.names = FALSE, sep = "\t")
  
  make_rainbow_gauge_plot(OverallPhysicalHealth_df, OverallPhysicalHealth_score_column/100, OverallPhysicalHealth_label, OverallPhysicalHealth_group)
  
  plotname_GutProcess <- paste(output_dir,"/",inputFilePrefix,"_OverallPhysicalHealth.png", sep = "")
  
  ggsave(plotname_GutProcess, plot = last_plot(), bg = "white")
  
  #read in all .txt files 
  
  
  Data.in <- lapply(list.files(path = output_dir, pattern = "GP_*", full.names = TRUE),read_tsv)  
  
  #combines all of the tables by column into one 
  Data.in <- do.call(cbind,Data.in)  %>%  select(ends_with("_col"))
  
  
  fileName_allGP <- paste(output_dir,"/", inputFilePrefix, "_all_GutPotential_Segments.txt", sep = "")
  write.table(Data.in, fileName_allGP, quote=FALSE, row.names = FALSE, sep = "\t")
  
  ### list all TXT files 
  
  files_to_remove_list <- list.files(path = output_dir, pattern = "GP_*", full.names = TRUE)
  file.remove(files_to_remove_list)
  
  
  ### UNCOMMENT IF USING HTML REPORTS 
  
  #read in all .png files and crop
  
  
  png_files_ls = list.files(path=output_dir, pattern="*.png", recursive = TRUE, full.names = TRUE)
  
  for (i in png_files_ls) {
    
    imported_png <- image_read(i)  
    image_write(image_trim(imported_png), i)  
  }
  
  
  
  
}





#_______________________checkSymptoms_GG2 _______________________________________________


## note change from March 25th, 2024 - print replaced with paste - console output not needed, results go into txt file

checkSymptoms_v2_GG2 <- function(inputFilePrefix, InputTable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus) {
  
  output_dir <- paste("./",inputFilePrefix,"/Symptoms/",sep="")
  
  dir.create(output_dir)
  
  sampleID <- paste(inputFilePrefix)  
  
  all_genus <- createGenusTable(inputFilePrefix, biomtable)
  genus_ref <- aggregate(RA ~ Genus, all_genus, sum)
  genus_ag <- aggregate(RA ~ Genus, genus_ref, sum)
  
  
  genus_ag$Percentage <- paste0(genus_ag$Genus,' [',100*round(genus_ag$RA,10),'%',']')
  
  ###NOTE: "up" is what INCREASES the likelihood of symptom (Positive association)
  ###      "down" is what DECREASES the likelihood of symptom (Negative association)
  
  ### Bloating
  
  
  Genus_Bloating_RAup <- c("Bacteroides", "Dorea", "Ruminococcus", "Clostridium") 
  Genus_Bloating_RAdown <- c("Bifidobacterium")  
  
  #Calculate score at Genus level  
  
  
  total_genus_associated_with_symptom_Bloating <- calculateSymptomScoreGenus(inputFilePrefix, Genus_Bloating_RAup, Genus_Bloating_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  Bloating_score_column <- total_genus_associated_with_symptom_Bloating
  Bloating_df <- data.frame(Bloating_score_column)
  
  
  ## Abdominal discomfort 
  
  Genus_AbdominalDiscomfort_RAup <- c("Bacteroides") 
  Genus_AbdominalDiscomfort_RAdown <- c("Faecalibacterium")
  
  
  
  total_genus_associated_with_symptom_AbdominalDiscomfort <- calculateSymptomScoreGenus(inputFilePrefix, Genus_AbdominalDiscomfort_RAup, Genus_AbdominalDiscomfort_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  AbdominalDiscomfort_score_column <- total_genus_associated_with_symptom_AbdominalDiscomfort
  AbdominalDiscomfort_df <- data.frame(AbdominalDiscomfort_score_column)
  
  
  ## Constipation Increase 
  
  Genus_Constipation_RAup <- c("Methanobrevibacter", "Alistipes", "Bacteroides") 
  Genus_Constipation_RAdown <- c("Bifidobacterium", "Lactobacillus")  
  
  #Calculate score at Genus level  
  
  
  total_genus_associated_with_symptom_Constipation <- calculateSymptomScoreGenus(inputFilePrefix, Genus_Constipation_RAup, Genus_Constipation_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  Constipation_score_column <- total_genus_associated_with_symptom_Constipation
  Constipation_df <- data.frame(Constipation_score_column)
  
  
  
  ## PoorSleep Increase 
  
  Genus_PoorSleep_RAup <- c("Blautia", "Bacteroides", "Coprococcus",  "Streptococcus", "Veillonella") 
  Genus_PoorSleep_RAdown <- c("Bifidobacterium", "Anaerostipes", "Faecalibacterium", "Parabacteroides", "Phascolarctobacterium", "Prevotella", 
                              "Roseburia","Agathobacter", "Lachnospira", "Agathobaculum", #added to DB April 2024, SCFA 
                              "Gemmiger","Anaerobutyricum", "Corynebacterium")    
  
  #Calculate score at Genus level  
  
  
  total_genus_associated_with_symptom_PoorSleep <- calculateSymptomScoreGenus(inputFilePrefix, Genus_PoorSleep_RAup, Genus_PoorSleep_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  PoorSleep_score_column <- total_genus_associated_with_symptom_PoorSleep
  PoorSleep_df <- data.frame(PoorSleep_score_column)
  
  
  
  ## Fatigue Increase 
  
  Genus_Fatigue_RAup <- c("Streptococcus") 
  Genus_Fatigue_RAdown <- c("Anaerostipes", "Bacteroides")
  #Calculate score at Genus level  
  
  
  total_genus_associated_with_symptom_Fatigue <- calculateSymptomScoreGenus(inputFilePrefix, Genus_Fatigue_RAup, Genus_Fatigue_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  Fatigue_score_column <- total_genus_associated_with_symptom_Fatigue
  Fatigue_df <- data.frame(Fatigue_score_column)
  
  
  
  ## Stress Increase 
  
  Genus_Stress_RAup <- c("Fusobacterium","Prevotella") 
  Genus_Stress_RAdown <- c("Lactobacillus", "Lactococcus", "Enterococcus", "Pediococcus", "Bifidobacterium", "Parabacteroides", "Streptococcus", 
                           "Veillonella") #Veillonella - HPA axis link 
  #Calculate score at Genus level  
  
  
  total_genus_associated_with_symptom_Stress <- calculateSymptomScoreGenus(inputFilePrefix, Genus_Stress_RAup, Genus_Stress_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  Stress_score_column <- total_genus_associated_with_symptom_Stress
  Stress_df <- data.frame(Stress_score_column)
  
  
  
  ## Anxiety Increase 
  
  Genus_Anxiety_RAup <- c("Bacteroides", "Eggerthella", "Ruminococcus") #removed Bifidobacterium --> moved to NEGATIVE association with Anxiety ) 
  Genus_Anxiety_RAdown <- c("Bifidobacterium", "Faecalibacterium", "Lachnospira", "Sutterella", "Dialister", "Faecalibacterium") 
  
  #Calculate score at Genus level  
  
  
  total_genus_associated_with_symptom_Anxiety <- calculateSymptomScoreGenus(inputFilePrefix, Genus_Anxiety_RAup, Genus_Anxiety_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  Anxiety_score_column <- total_genus_associated_with_symptom_Anxiety
  Anxiety_df <- data.frame(Anxiety_score_column)
  
  
  
  #____________DrySkin__________________
  
  ### Note some bacteria are new: from VIVO database 
  
  Genus_DrySkin_RAup <- c("Staphylococcus", "Clostridium", "Escherichia")
  
  
  Genus_DrySkin_RAdown <- c("Akkermansia", "Bifidobacterium", "Faecalibacterium", "Lactobacillus")
  
  
  
  
  total_genus_associated_with_symptom_DrySkin <- calculateSymptomScoreGenus(inputFilePrefix, Genus_DrySkin_RAup, Genus_DrySkin_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  DrySkin_score_column <- total_genus_associated_with_symptom_DrySkin
  DrySkin_df <- data.frame(DrySkin_score_column)
  
  
  
  #____PoorMucusProduction
  Genus_PoorMucusProduction_RAup <- c("")
  Genus_PoorMucusProduction_RAdown <- c("Lactobacillus", "Bifidobacterium", "Anaerostipes", "Akkermansia", "Alistipes", "Bacteroides", "Blautia", "Clostridium", "Coprococcus", "Dialister", "Eubacterium", "Faecalibacterium", "Parabacteroides", "Phascolarctobacterium", "Prevotella", "Roseburia","Ruminococcus", "Streptococcus", "Veillonella", "Agathobaculum", "Gemmiger", "Anaerobutyricum", "Agathobacter")   
  
  #"Butyricicoccus" removed due to low freq 
  
  total_genus_associated_with_symptom_PoorMucus <- calculateSymptomScoreGenus(inputFilePrefix, Genus_PoorMucusProduction_RAup, Genus_PoorMucusProduction_RAdown, InputTable = biomtable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus)
  
  PoorMucus_score_column <- total_genus_associated_with_symptom_PoorMucus
  PoorMucus_df <- data.frame(PoorMucus_score_column)
  
  
  ### Now combine scores 
  
  symptoms_df <- data.frame()
  symptoms_df = data.frame(matrix(nrow = 1, ncol = 0))
  
  symptoms_df$Name = inputFilePrefix
  symptoms_df$Bloating = Bloating_df$Bloating_score_column
  symptoms_df$Constipation = Constipation_df$Constipation_score_column
  symptoms_df$AbdominalDiscomfort <- AbdominalDiscomfort_df$AbdominalDiscomfort_score_column
  symptoms_df$Fatigue = Fatigue_df$Fatigue_score_column
  symptoms_df$PoorSleep = PoorSleep_df$PoorSleep_score_column
  symptoms_df$Anxiety = Anxiety_df$Anxiety_score_column
  symptoms_df$Stress = Stress_df$Stress_score_column
  symptoms_df$DrySkin = DrySkin_df$DrySkin_score_column
  symptoms_df$PoorMucus = PoorMucus_df$PoorMucus_score_column
  
  
  symptoms_df <- symptoms_df %>% 
    mutate(Bloating_symptom  = ifelse(Bloating <= 0.5, "Low", "High")) %>% ## bloating is output in combination with reduce food type list
    mutate(Digestion_symptom = ifelse(Constipation <= 0.5, "The bacteria within your gut are effective in helping you digest your meal with ease, linked to regular bowel movements.", "The bacteria within your gut are not effective in digesting your meal with ease, this can result in irregular bowel movements.")) %>% 
    mutate(Abdominal_pain_symptom  = ifelse(AbdominalDiscomfort <= 0.5, "Your current gut microbiome composition does not suggest an increased susceptibility to abdominal discomfort.", "Feelings of abdominal discomfort are linked with your current gut microbiome profile.")) %>% 
    mutate(Fatigue_symptom  = ifelse(Fatigue <= 0.5, "Based on the capabilities of your current gut microbiome you are less susceptible to feelings of tiredness and fatigue.", "Based on the capabilities of your current gut microbiome you are more susceptible to feelings of tiredness and fatigue.")) %>% 
    mutate(Sleep_symptom  = ifelse(PoorSleep <= 0.5, "Your current gut microbiome is not associated with sleeping impairment.", "Elements of your gut microbiome may be a contributing factor to problems connected to your sleep.")) %>% 
    mutate(Anxiety_symptom  = ifelse(Anxiety <= 0.5, "Your current gut microbiome composition is not associated with elevated feelings of anxiety.", "Frequent feelings of anxiety are linked with your current gut microbiome profile.")) %>% 
    mutate(Stress_symptom  = ifelse(Stress <= 0.5, "Your current gut microbiome is productive in helping you deal with stress.", "Your current gut microbiome is associated with increased levels of stress.")) %>% 
    mutate(Skin_symptom  = ifelse(DrySkin <= 0.5, "Elements of your gut microbiome may be protective against symptoms of mild skin conditions including itchy/dry skin.", "Elements of your gut microbiome may be contributing to certain skin conditions including itchy/dry skin.")) %>% 
    mutate(Gut_lining_cell_nutrients_symptom  = ifelse(PoorMucus <= 0.5, "Your microbiome’s ability to produce nutrients for intestinal cells is good.", "Your microbiome’s ability to produce nutrients for intestinal cells is low.")) %>% 
    mutate(Gut_lining_mucus_symptom  = ifelse(PoorMucus <= 0.5, "Your current gut microbiome is consistent with individuals who have good mucus production, which supports gut lining integrity.", "Your current gut microbiome is consistent with individuals who have poor mucus production, which might result in insufficient gut lining integrity.")) 
  
  
  #print(symptoms_df)
  
  
  symptomsDF_output_Filename <-  paste(output_dir, "/", inputFilePrefix, "_Symptoms.txt", sep = "")
  
  write.table(symptoms_df, symptomsDF_output_Filename, quote = FALSE, row.names = FALSE, sep = "\t")
  
  
}








#_________getBeneficialBacteria_GG2_________________________________________

getBeneficialBacteria_GG2_v3 <- function(inputFilePrefix, InputTable, language = "ENG", ReferenceQuartileTable = benefBac_reference) {
  
  
  output_dir <- paste("./",inputFilePrefix,"/BeneficialBacteria/", sep="")
  
  dir.create(output_dir)
  
  sampleID <- paste(inputFilePrefix)  
  
  all_genus <- createGenusTable(inputFilePrefix, biomtable) %>% 
    mutate(Genus = gsub(".*actobacillus", 'Lactobacillus', Genus)) %>%  #important: this is to show all Lactobacillus re-classified genera as Lactobacillus for the purposes of the plot; elsewhere keep GG2 re-classified genus names - this only for this plot
    mutate(Genus = gsub("Lacticaseibacillus", 'Lactobacillus', Genus)) #same as above
  
  genus_ref <- aggregate(RA ~ Genus, all_genus, sum)
  genus_ag <- aggregate(RA ~ Genus, genus_ref, sum)
  
  
  genus_ag$Percentage <- paste0(genus_ag$Genus,' [',100*round(genus_ag$RA,10),'%',']')
  
  
  
  beneficial_bacteria <- c("Akkermansia", 
                           "Alistipes",
                           "Bifidobacterium",
                           "Blautia",
                           "Eubacterium",
                           "Faecalibacterium",
                           "Lactobacillus",
                           "Roseburia",
                           "Ruminococcus")
  
  
  beneficial_bacteria_vec <- c() #Take the vector outside of the loop, so that whatever is output from a loop gets appended to it and the end result is one vector with multiple elements 
  
  q1_ben <- c()
  q2_ben <- c()
  q3_ben <- c()
  
  
  for(i in beneficial_bacteria) {
    
    
    sample_value_percent <-  (subset(genus_ag, Genus == i)$RA)*100
    sample_value_percent <- ifelse(length(sample_value_percent) == 0, 0, sample_value_percent)
    beneficial_bacteria_vec <- c(beneficial_bacteria_vec, sample_value_percent)
    
    q1_ben_val <- as.data.frame(ReferenceQuartileTable %>%  filter(Genus == i) %>% select(q25_percent))
    q1_ben <- c(q1_ben, q1_ben_val)
    
    q2_ben_val <- as.data.frame(ReferenceQuartileTable %>%  filter(Genus == i) %>% select(q50_percent))
    q2_ben <- c(q2_ben, q2_ben_val)
    
    q3_ben_val <- as.data.frame(ReferenceQuartileTable %>%  filter(Genus == i) %>% select(q75_percent))
    q3_ben <- c(q3_ben, q3_ben_val)
    
    
  }
  
  beneficial_bacteria_table <- data.frame(beneficial_bacteria, beneficial_bacteria_vec, unlist(q1_ben),unlist(q2_ben), unlist(q3_ben))  %>% 
    mutate(beneficial_bacteria_vec_percent =  round(beneficial_bacteria_vec,2)) %>% 
    mutate(q1_percent =  round(`unlist.q1_ben.`,2)) %>% 
    mutate(q2_percent =  round(`unlist.q2_ben.`,2)) %>%
    mutate(q3_percent =  round(`unlist.q3_ben.`,2)) %>% 
    plyr::rename(c("beneficial_bacteria_vec_percent" = "You"))  %>% 
    plyr::rename(c(q1_percent = "HealthyQ1")) %>% 
    plyr::rename(c(q2_percent = "Healthy")) %>%
    plyr::rename(c(q3_percent = "HealthyQ3")) %>% 
    plyr::rename(c(beneficial_bacteria = "BeneficialBacteria")) %>% 
    mutate(
      Status = case_when(
        You > 4 * Healthy ~ "Overgrowth",
        You < 0.5 * Healthy ~ "Low", TRUE ~ "Normal_range")) %>% #note an update! changed .default to TRUE #another update to match definition of "low" to that in Excel File 
    select(BeneficialBacteria, You, HealthyQ1, Healthy, HealthyQ3, Status)
  
  
  
  for(i in beneficial_bacteria) {
    beneficial_bacteria_table_gathered <- beneficial_bacteria_table %>% 
      select(BeneficialBacteria, You, Healthy) %>% 
      gather("Sample", "Percentage", -BeneficialBacteria) %>% 
      filter(BeneficialBacteria  == i)
    
    beneficial_bacteria_table_gathered <-  if (language == "ENG") { 
      beneficial_bacteria_table_gathered <- beneficial_bacteria_table_gathered
      
    } else {
      beneficial_bacteria_table_gathered <- beneficial_bacteria_table_gathered %>% 
        mutate(Sample = gsub("You", "Vous", Sample)) %>% 
        mutate(Sample = gsub("Healthy", "Sain", Sample))
    }  
    
    benefB_plot <- ggplot(beneficial_bacteria_table_gathered, aes(x=Sample, y=Percentage, color = Sample, fill = Sample)) +
      geom_bar(stat="identity") +
      scale_color_manual(values=c("#672C94", "#141EDC")) +
      scale_fill_manual(values=c("#672C94", "#141EDC")) +
      theme_minimal() +
      geom_text(
        aes(x = Sample, y = Percentage, label = paste0(Percentage,"%")),
        position = position_dodge(width = 1),
        vjust = -0.5, size = 8)  + 
      theme(legend.position="none") +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())  +
      theme(axis.text.x = element_text(color = "grey20", size = 22, vjust = 3.8),
            axis.title.x = element_text(color = "grey20", size = 22, face = "bold")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_y_continuous(limits = c(0, max(beneficial_bacteria_table_gathered$Percentage + 2))) +
      labs(x = i)
    
    
    plotname_BenefB <- paste(output_dir,"/",inputFilePrefix,i,".png", sep = "")
    
    ggsave(plotname_BenefB, plot = last_plot(), bg = "white")
    
    
    
  }
  
  
  tablename_BeneficialBacteria <- paste(output_dir,"/",inputFilePrefix,"_BeneficialBacteria.txt", sep = "")
  write.table(beneficial_bacteria_table, tablename_BeneficialBacteria, quote = FALSE, row.names = FALSE, sep = "\t")
  
  
}






getinflammatoryBacteria_GG2_v3 <- function(inputFilePrefix, InputTable, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus) {
  
  # output_dir <- paste("./",inputFilePrefix,"./inflammatoryBacteria/", sep="")
  output_dir <- paste("./",inputFilePrefix,"/inflammatoryBacteria/", sep="")
  
  dir.create(output_dir)
  
  sampleID <- paste(inputFilePrefix)  
  
  all_genus <- createGenusTable(inputFilePrefix, biomtable)
  
  genus_ref <- aggregate(RA ~ Genus, all_genus, sum)
  genus_ag <- aggregate(RA ~ Genus, genus_ref, sum)
  
  
  genus_ag$Percentage <- paste0(genus_ag$Genus,' [',100*round(genus_ag$RA,10),'%',']')
  
  
  inflammatory_bacteria <- c("Desulfovibrio", "Collinsella", "Fusobacterium")
  
  inflammatory_bacteria_vec <- c() #Take the vector outside of the loop, so that whatever is output from a loop gets appended to it and the end result is one vector with multiple elements 
  
  q1_ben <- c()
  q2_ben <- c()
  q3_ben <- c()
  for(i in inflammatory_bacteria) { 
    
    
    sample_value_percent <-  (subset(genus_ag, Genus == i)$RA)*100
    sample_value_percent <- ifelse(length(sample_value_percent) == 0, 0, sample_value_percent)
    inflammatory_bacteria_vec <- c(inflammatory_bacteria_vec, sample_value_percent)
    
    q1_ben_val <- as.data.frame(ReferenceQuartileTable %>%  filter(Genus == i) %>% select(q25_percent))
    q1_ben <- c(q1_ben, q1_ben_val)
    
    q2_ben_val <- as.data.frame(ReferenceQuartileTable %>%  filter(Genus == i) %>% select(q50_percent))
    q2_ben <- c(q2_ben, q2_ben_val)
    
    q3_ben_val <- as.data.frame(ReferenceQuartileTable %>%  filter(Genus == i) %>% select(q75_percent))
    q3_ben <- c(q3_ben, q3_ben_val)
    
    
  }
  inflammatory_bacteria_table <- data.frame(inflammatory_bacteria, inflammatory_bacteria_vec, unlist(q1_ben),unlist(q2_ben), unlist(q3_ben)) %>%
    mutate(inflammatory_bacteria_vec_percent =  round(inflammatory_bacteria_vec,2)) %>% 
    mutate(q1_percent =  round(`unlist.q1_ben.`,2)) %>% 
    mutate(q2_percent =  round(`unlist.q2_ben.`,2)) %>%
    mutate(q3_percent =  round(`unlist.q3_ben.`,2)) %>% 
    plyr::rename(c("inflammatory_bacteria_vec_percent" = "You"))  %>% 
    plyr::rename(c(q1_percent = "HealthyQ1")) %>% 
    plyr::rename(c(q2_percent = "HealthyMedian")) %>%
    plyr::rename(c(q3_percent = "HealthyQ3")) %>% 
    plyr::rename(c(inflammatory_bacteria = "inflammatoryBacteria")) %>% 
    mutate(
      Status = case_when(
        You > HealthyQ3 ~ "Above_normal",
        You == 0 ~ "Not_detected", TRUE ~ "Normal_range")) %>% #note change here as well: .default changed to TRUE
    select(inflammatoryBacteria, You, HealthyQ1, HealthyMedian, HealthyQ3, Status) 
  
  
  
  tablename_inflammatoryBacteria <- paste(output_dir,"/",inputFilePrefix,"_inflammatoryBacteria.txt", sep = "")
  write.table(inflammatory_bacteria_table, tablename_inflammatoryBacteria, quote = FALSE, row.names = FALSE, sep = "\t")
  
  
}







getProbioticBacteria_GG2_v3 <- function(inputFilePrefix, ReferenceQuartileTable = reference_valuesAGP_GG2_allsamples_genus) {
  
  
  output_dir <- paste("./",inputFilePrefix,"/ProbioticBacteria/", sep="")
  
  dir.create(output_dir)
  
  sampleID <- paste(inputFilePrefix)  
  all_genus <- createGenusTable(inputFilePrefix, biomtable)
  genus_ref <- aggregate(RA ~ Genus, all_genus, sum)
  genus_ag <- aggregate(RA ~ Genus, genus_ref, sum)
  genus_ag$Percentage <- paste0(genus_ag$Genus,' [',100*round(genus_ag$RA,10),'%',']')
  
  
  
  probiotic_bacteria <- c("Bacillus", 
                          "Bifidobacterium", 
                          "Limosilactobacillus", 
                          "Lactobacillus", 
                          "Lentilactobacillus", 
                          "Akkermansia")
  
  #Lactiplantibacillus and Alkalihalobacillus low frequency
  
  probiotic_bacteria_vec <- c() #Take the vector outside of the loop, so that whatever is output from a loop gets appended to it and the end result is one vector with multiple elements 
  
  q1_ben <- c()
  q2_ben <- c()
  q3_ben <- c()
  
  for(i in probiotic_bacteria) { 
    
    sample_value_percent <-  (subset(genus_ag, Genus == i)$RA)*100
    sample_value_percent <- ifelse(length(sample_value_percent) == 0, 0, sample_value_percent)
    probiotic_bacteria_vec <- c(probiotic_bacteria_vec, sample_value_percent)
    
    q1_ben_val <- as.data.frame(ReferenceQuartileTable %>%  filter(Genus == i) %>% select(q25_percent))
    q1_ben <- c(q1_ben, q1_ben_val)
    
    q2_ben_val <- as.data.frame(ReferenceQuartileTable %>%  filter(Genus == i) %>% select(q50_percent))
    q2_ben <- c(q2_ben, q2_ben_val)
    
    q3_ben_val <- as.data.frame(ReferenceQuartileTable %>%  filter(Genus == i) %>% select(q75_percent))
    q3_ben <- c(q3_ben, q3_ben_val)
    
    
  }
  
  
  
  probiotic_bacteria_table <- data.frame(probiotic_bacteria, probiotic_bacteria_vec, unlist(q1_ben),unlist(q2_ben), unlist(q3_ben)) %>%
    mutate(probiotic_bacteria_vec_percent =  round(probiotic_bacteria_vec,2)) %>% 
    mutate(q1_percent =  round(`unlist.q1_ben.`,2)) %>% 
    mutate(q2_percent =  round(`unlist.q2_ben.`,2)) %>%
    mutate(q3_percent =  round(`unlist.q3_ben.`,2)) %>% 
    plyr::rename(c("probiotic_bacteria_vec_percent" = "You"))  %>% 
    plyr::rename(c(q1_percent = "HealthyQ1")) %>% 
    plyr::rename(c(q2_percent = "Healthy")) %>%
    plyr::rename(c(q3_percent = "HealthyQ3")) %>% 
    plyr::rename(c(probiotic_bacteria = "ProbioticBacteria")) %>% 
    mutate(
      Recommendation = case_when(
        You == 0.00 ~ "Recommended",  #If customer has 0% of that genus 
        You < 0.5 * Healthy ~ "Recommended", #If customer % is less than a half of the median of the reference population
        TRUE ~ "Present at optimal level or higher")) %>% 
    select(ProbioticBacteria, You, HealthyQ1, Healthy, HealthyQ3, Recommendation) 
  
  
  probiotic_bacteria_table
  
  
  tablename_probioticBacteria <- paste(output_dir,"/",inputFilePrefix,"_ProbioticBacteria.txt", sep = "")
  write.table(probiotic_bacteria_table, tablename_probioticBacteria, quote = FALSE, row.names = FALSE, sep = "\t")
  
  
  
  
  
  
  
  
  
  
}








#____________________calculateGutScore_________________________________________


##this function remains the same as the function that was used for GG138 taxonomy
calculateGutScore <- function(inputFilePrefix, BatchCode) {
  output_dir <- paste("./",inputFilePrefix,sep="")
  
  dir.create(output_dir)
  
  
  
  sampleID <- paste(inputFilePrefix)
  GP_samplename <- paste(output_dir, "/GutPotential/",inputFilePrefix,"_all_GutPotential_Segments.txt", sep = "")
  
  
  All_GutPotential_Segments <- read_tsv(GP_samplename) %>% 
    select(GutBrainAxis_col, OverallMetabolism_col, OverallPhysicalHealth_col, OverallImmuneScore_col) %>% 
    mutate(MeanGutPotentialScore = mean(c_across(c(GutBrainAxis_col, OverallMetabolism_col, OverallPhysicalHealth_col, OverallImmuneScore_col)))/100) 
  
  # diversity <- read_tsv("./alpha-diversity.tsv", skip = 1, col_names = FALSE) %>%  # old: fails when alpha-diversity.tsv is a directory
  diversity <- read_tsv(get_alpha_diversity_path(), skip = 1, col_names = FALSE) %>% 
    plyr::rename(c("X1" = "SampleID")) %>% 
    plyr::rename(c("X2" = "shannon_entropy")) %>% 
    filter(SampleID == inputFilePrefix) %>% 
    mutate(shannon_entropy_percent = shannon_entropy/10) 
  
  benefBacteriaFileName <- paste(output_dir, "/BeneficialBacteria/", inputFilePrefix, "_BeneficialBacteria.txt", sep = "")
  beneficialBacteria <- read_tsv(benefBacteriaFileName) %>% 
    mutate(benefB_Score = ifelse(You >= HealthyQ1, 1, 0))
  
  
  combinedScores <- cbind(diversity, All_GutPotential_Segments)
  combinedScores$BeneficialBacteriaScore <- mean(beneficialBacteria$benefB_Score)
  
  
  diseasepredictionFileName <- paste("./Comparison_all_diseases_", BatchCode,".xlsx", sep = "")
  
  disease_prediction <- read_xlsx(diseasepredictionFileName) %>% 
    select(`OTU ID`, `Obesity Final Status`, `IBS Final Status`, `IBD Final Status`, `Diabetes Final Status`, `CVD Final Status`, `Arthritis Final Status`, `MASLD Final Status`, `Clostridium difficile infection Final Status`) %>%
    rename_with(~ str_remove(., " Final Status"), everything()) %>% 
    plyr::rename(c("OTU ID" = "SampleID")) %>% 
    mutate( across(
      .cols = everything(),
      ~str_replace( ., "healthy", "1" ))) %>% 
    mutate( across(
      .cols = everything(),
      ~str_replace( ., "low risk", "0.75" ))) %>% 
    mutate( across(
      .cols = everything(),
      ~str_replace( ., "moderate risk", "0.5" ))) %>% 
    mutate( across(
      .cols = everything(),
      ~str_replace( ., "high risk", "0" ))) %>% 
    filter(SampleID == inputFilePrefix) %>% 
    mutate_at(c(2:9), as.numeric)
  
  combinedScores$disease_prediction_mean_score <- rowMeans(disease_prediction[2:9], na.rm = TRUE)
  
  outputFileName <-  paste(output_dir,"/",inputFilePrefix,"_GutScore.txt", sep = "")
  
  combinedScores %>% 
    rowwise() %>% mutate(GutCheckScore = mean(c(shannon_entropy_percent, MeanGutPotentialScore, BeneficialBacteriaScore, disease_prediction_mean_score), na.rm = T)) %>% 
    write.table(outputFileName, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  
  
}











#______________________printInputTextboxesForReport_GG2____________________________________


printInputTextboxesForReport_GG2 <- function(inputFilePrefix, BatchCode, language = "ENG",
  # FoodMappingFile = "C:/Users/farid/OneDrive/Desktop/NucliqBiologics/Food/Canada/CANADA_food_categories_mapping_updated30June2025.csv") { # old hardcoded Windows path
  FoodMappingFile = "./CANADA_food_categories_mapping_updated.csv") {
  
  output_dir <- paste("./",inputFilePrefix,sep="")
  dir.create(output_dir)
  
  
  
  # diversity <- read_tsv("./alpha-diversity.tsv", skip = 1, col_names = FALSE) %>%  # old: fails when alpha-diversity.tsv is a directory
  diversity <- read_tsv(get_alpha_diversity_path(), skip = 1, col_names = FALSE) %>% 
    plyr::rename(c("X1" = "SampleID")) %>% 
    plyr::rename(c("X2" = "shannon_entropy")) %>% 
    filter(SampleID == inputFilePrefix) %>% 
    mutate(shannon_entropy_percent = shannon_entropy/10) %>% 
    mutate(shannon_entropy_round = round(shannon_entropy, 2)) %>% #NEW ADDITION 
    mutate(shannon_category = ifelse(shannon_entropy_round <4, "Low",
                                     ifelse(shannon_entropy_round >=4 & shannon_entropy_round<6, "Normal","High"))) 
  
  
  #Figure added Feb 22, 2024 and updated March 1st, 2024
  diversity_figure <- ggplot() + ylim(2.45, 2.75) + xlim(-1, 11) + 
    geom_rect(aes(xmin=0, 
                  xmax=10, ymin=2.598, ymax=2.603),  
              colour = "black", fill = "black") + 
    geom_rect(mapping=aes(xmin=4, 
                          xmax=6, ymin=2.598, ymax=2.603),  
              colour = "#5070ff", fill = "#5070ff") + 
    geom_point(aes(x=4.05, y = 2.6025), size = 7, colour = "#5070ff") + 
    geom_point(aes(x=5.95, y = 2.6025), size = 7, colour = "#5070ff") + 
    annotate('text', x = -0.5, y = 2.603, label = "0", color =  "black", size = 5) +
    annotate('text', x = 10.5, y = 2.603, label = "10", color =  "black", size = 5) +
    theme_void() +
    geom_point(aes(x=diversity$shannon_entropy_round, y = 2.625), size = 12, colour = "#00367f", fill = "#00367f", shape = 25) 
  
  
  diversity_figure <- if (language == "ENG") { 
    diversity_figure <- diversity_figure + 
      annotate('text', x = 0.5, y = 2.565, label = "Low Diversity", color =  "black", size = 3.5) +
      annotate('text', x = 9.5, y = 2.565, label = "High Diversity", color =  "black", size = 3.5) +
      annotate('text', x = 5, y = 2.565, label = "Normal range: 4 to 6", fontface = "bold", color =  "#5070ff", size = 3.5)  +
      annotate('text', x = diversity$shannon_entropy_round, y = 2.655, label = "You", color =  "#00367f", size = 5)
  } else {
    diversity_figure <- diversity_figure + 
      annotate('text', x = 0.5, y = 2.565, label = "Diversité faible", color =  "black", size = 3.5) +
      annotate('text', x = 9.5, y = 2.565, label = "Diversité élevée", color =  "black", size = 3.5) +
      annotate('text', x = 5, y = 2.565, label = "Plage de valeurs normales: 4 à 6 ", fontface = "bold", color =  "#5070ff", size = 3.5)  +
      annotate('text', x = diversity$shannon_entropy_round, y = 2.655, label = "Vous", color =  "#00367f", size = 5)
    
  }
  
  
  
  plotname_diversityfigure <- paste(output_dir,"/",inputFilePrefix,"_diversity_figure.png", sep = "")
  ggsave(plotname_diversityfigure, plot = diversity_figure, height = 1300, width  = 2100, units = "px")
  
  similarity_Filename <- paste("./Phylum_level_Similarity_GutCheck_", BatchCode, "Project.txt", sep="")
  similarity <- read_tsv(similarity_Filename) %>% 
    filter(SampleID == inputFilePrefix) %>% 
    mutate(Similarity_percent = Similarity_round*100) 
  
  GutScore_Filename <- paste(output_dir, "/", inputFilePrefix, "_GutScore.txt", sep = "")
  GutCheckScoreDF <- read_tsv(GutScore_Filename) %>%
    mutate(GutCheckScore_round = round(GutCheckScore, 2)) %>% 
    mutate(GutCheckScore_category = ifelse(GutCheckScore_round <= 0.40, "Poor",  ifelse(GutCheckScore_round>0.40 & GutCheckScore_round<=0.69, "Normal","Great"))) %>% 
    mutate(group=ifelse(GutCheckScore_round <=0.40, "red",  ifelse(GutCheckScore_round>0.40 & GutCheckScore_round<=0.69, "orange","green")),
           label=paste0("Overall\n Gutcheck\U2122 Score:\n", round(GutCheckScore,2)*100, "%"),  GutScore_Percentage = GutCheckScore_round*100)
  
  GutCheckScoreDF <- if (language == "ENG") { 
    GutCheckScoreDF <- GutCheckScoreDF
    
  } else {
    GutCheckScoreDF <- GutCheckScoreDF %>% 
      mutate(label = gsub("Overall\n Gutcheck\U2122 Score:\n", "Votre score\n de Gutcheck\U2122:\n", label)) 
    
  }
  
  ## Figure added Feb 22, 2024 and updated March 1st, 2024 
  
  GutCheckScoreFigure <- ggplot(GutCheckScoreDF, aes(fill = GutCheckScore_category, ymax = GutCheckScoreDF$GutCheckScore*2, ymin = 0, xmax = 1.3, xmin = 1)) +
    geom_rect(aes(ymax=2, ymin=0, xmax=1.3, xmin=1), fill = ifelse(GutCheckScoreDF$GutCheckScore_round <=0.40, "#FF7F7F",
                                                                   ifelse(GutCheckScoreDF$GutCheckScore_round>0.40 & GutCheckScoreDF$GutCheckScore_round<=0.69, "#fad8ac","#c0edbb")))  +  #updated Nov2024!
    geom_rect() + 
    coord_polar(theta = "y",  start=0) + xlim(c(0, 2.5)) + ylim(c(0,2)) +
    geom_text(aes(x = 0, y = 0, label = label, colour="black"), colour = "black", size=5.8)  + 
    #geom_text(aes(x=2.5, y=1, label= "GutCheckScore"), size=5.5) +
    #facet_wrap(~title, ncol = 5) +
    theme_void() +
    scale_fill_manual(name = NULL, values = c("Poor"="#ff1616", "Normal"="#e7c536", "Great"="#7dd956")) +
    scale_colour_manual(values = c("red"="#ff1616", "orange"="#e7c536", "green"="#7dd956")) +
    theme(strip.background = element_blank(), 
          strip.text.x = element_blank()) +
    guides(fill=FALSE) 
  
  plotname_GutCheckScoreFigure <- paste(output_dir,"/",inputFilePrefix,"_GutcheckScore_figure.png", sep = "")
  ggsave(plotname_GutCheckScoreFigure, plot = GutCheckScoreFigure, bg = "white")
  
  
  
  GutCheckScoreFigure_short <- ggplot(GutCheckScoreDF, aes(fill = GutCheckScore_category, ymax = GutCheckScoreDF$GutCheckScore*2, ymin = 0, xmax = 1.3, xmin = 1)) +
    geom_rect(aes(ymax=2, ymin=0, xmax=1.3, xmin=1), fill = ifelse(GutCheckScoreDF$GutCheckScore_round <0.40, "white",
                                                                   ifelse(GutCheckScoreDF$GutCheckScore_round>=0.40 & GutCheckScoreDF$GutCheckScore_round<0.69, "white","white")))  +
    geom_rect() + 
    coord_polar(theta = "y",  start=0) + xlim(c(0, 2.5)) + ylim(c(0,2)) +
    geom_text(aes(x = 0, y = 0, label = paste(GutScore_Percentage, "%", sep = "")), colour = "white", size=17, fontface = "bold.italic")  + 
    #geom_text(aes(x=2.5, y=1, label= "GutCheckScore"), size=5.5) +
    #facet_wrap(~title, ncol = 5) +
    theme_void() +
    scale_fill_manual(name = NULL, values = c("Poor"="#ff1616", "Normal"="#e7c536", "Great"="#7dd956")) +
    scale_colour_manual(values = c("red"="#ff1616", "orange"="#e7c536", "green"="#7dd956")) +
    theme(strip.background = element_blank(), 
          strip.text.x = element_blank()) +
    guides(fill=FALSE)
  #theme(plot.background = element_rect(fill = "#100871"))
  
  plotname_GutCheckScoreFigure_short <- paste(output_dir,"/",inputFilePrefix,"_GutcheckScore_figure_shortReport.png", sep = "")
  
  ggsave(plotname_GutCheckScoreFigure_short, plot = GutCheckScoreFigure_short)
  
  
  
  
  diseasepredictionFileName <- paste("./Comparison_all_diseases_", BatchCode,".xlsx", sep = "")
  
  diseasepredictionDF <- read_xlsx(diseasepredictionFileName) %>% 
    select(`OTU ID`, `Obesity Final Status`, `IBS Final Status`, `IBD Final Status`, `Diabetes Final Status`, `CVD Final Status`, `Arthritis Final Status`, `MASLD Final Status`, `Clostridium difficile infection Final Status`) %>%
    rename_with(~ str_remove(., " Final Status"), everything()) %>% 
    plyr::rename(c("OTU ID" = "SampleID")) %>% 
    filter(SampleID == inputFilePrefix)
  
  
  disease_high_risk <- diseasepredictionDF %>% 
    dplyr::rename(diabetes = Diabetes) %>%  #Change to lower case,these names will go in the middle of a sentence 
    dplyr::rename(arthritis = Arthritis) %>% 
    dplyr::rename(obesity = Obesity) %>% 
    select(where(~ any(. == "high risk"))) %>% 
    names()
  
  
  
  disease_moderate_risk <- diseasepredictionDF %>%  
    dplyr::rename(diabetes = Diabetes) %>%  #Change to lower case,these names will go in the middle of a sentence 
    dplyr::rename(arthritis = Arthritis) %>% 
    dplyr::rename(obesity = Obesity) %>%
    select(where(~ any(. == "moderate risk"))) %>% 
    names()
  
  disease_low_risk <- diseasepredictionDF %>% 
    dplyr::rename(diabetes = Diabetes) %>%  #Change to lower case,these names will go in the middle of a sentence 
    dplyr::rename(arthritis = Arthritis) %>% 
    dplyr::rename(obesity = Obesity) %>%
    select(where(~ any(. == "low risk"))) %>% 
    names()
  
  
  #changed print to paste to remove console output 
  disease_prediction_output <- if(length(disease_high_risk) > 0){  
    high_risk_text <- paste("Your gut bacteria correlate with a high risk of", paste(noquote(stri_replace_last(paste(disease_high_risk, collapse = ", "), fixed = ',', ' and')), sep = " ") )  #replace last comma with and
    paste(high_risk_text)                           
  } else if(length(disease_moderate_risk) > 0) {                                             
    moderate_risk_text <- paste("Your gut bacteria correlate with a moderate risk of", paste(noquote(stri_replace_last(paste(disease_moderate_risk, collapse = ", "), fixed = ',', ' and')), sep = " ") ) #replace last comma with and
    paste(moderate_risk_text)
  } else if (length(disease_low_risk) > 0) {
    low_risk_text <- paste("Your gut bacteria correlate with a low risk of", paste(noquote(stri_replace_last(paste(disease_low_risk, collapse = ", "), fixed = ',', ' and')), sep = " ") ) #replace last comma with and
    paste(low_risk_text)
  } else {
    paste("Your gut microbiome does not show similarities to any of the disease cohorts")  
  }
  
  
  
  
  
  benefBacteriaFileName <- paste(output_dir, "/BeneficialBacteria/", inputFilePrefix, "_BeneficialBacteria.txt", sep = "")
  
  number_of_low_abundance_benficial_bacteria <- read_tsv(benefBacteriaFileName) %>% 
    filter(You < 0.5*(Healthy)) %>% 
    nrow()
  
  
  number_of_overgrowing_beneficial_bacteria <- read_tsv(benefBacteriaFileName)  %>% 
    filter(You > 4*(Healthy)) %>% 
    nrow()
  
  
  
  GP_samplename <- paste(output_dir,"/GutPotential/",inputFilePrefix,"_all_GutPotential_Segments.txt", sep = "")
  All_GutPotential_SegmentsDF <- read_tsv(GP_samplename)  
  
  
  #changed print to paste to remove console output
  Protein_output <- if(All_GutPotential_SegmentsDF$Protein_col <= 25){  
    paste("poor")                           
  } else {                                             
    paste("good")
  } 
  
  Gluten_output <- if(All_GutPotential_SegmentsDF$Gluten_col <= 25){  
    paste("poor")                           
  } else {                                             
    paste("good")
  } 
  
  
  
  
  
  
  Symptoms_Filename <- paste(output_dir, "/Symptoms/", inputFilePrefix, "_Symptoms.txt", sep = "")
  Symptoms_DF <- read_tsv(Symptoms_Filename)
  
  
  #weight maintenance symptom logic is below
  
  obesity_df <- diseasepredictionDF %>% 
    select(SampleID, Obesity) 
  
  #changed print to paste to remove console output
  
  obesity_output <- if (obesity_df$Obesity == "high risk" | obesity_df$Obesity == "moderate risk") {
    paste("The bacteria within your gut are not supportive in helping you to maintain a healthy body weight.") 
  } else if(obesity_df$Obesity == "low risk") {
    paste("The bacteria within your gut are moderately supportive in helping you to maintain a healthy body weight.") 
  } else {
    paste("The bacteria within your gut are supportive in helping you to maintain a healthy body weight.")
  }
  
  
  
  
  foodTypeDF <- read_tsv(GP_samplename) #Load DF with metabolism elements
  names(foodTypeDF) <- str_replace(names(foodTypeDF), "_col", "")
  foodTypeDF <- foodTypeDF %>%  
    select(Protein, Lactose, Gluten, Fat, Carbohydrate) %>% 
    setNames(tolower(names(.)))
  
  below25foodtypes <- foodTypeDF %>% 
    mutate(across(.fns = ~ifelse(. <=25, "poor", .))) %>% #If any foodtypes are below or equal to 25%
    select(where(~ any(. == "poor"))) %>%                #convert the number value to 'poor'
    names()                                              #extract names of columns that are poor
  
  #changed print to paste to remove console output
  reduceFoodType <- if(length(below25foodtypes) > 0){  #if at least one is poor 
    below25foodtypes_list <- paste(below25foodtypes, collapse = ", ") #print its name
    below25foodtypes_list <- stri_replace_last(below25foodtypes_list, fixed = ',', ' and') #replace last comma with and
    below25food_output <- paste("Consuming high amounts of ",  below25foodtypes_list, "rich foods may induce abdominal discomfort", sep = " ")  
    paste(below25food_output)
  } else {                                             #if none are poor
    paste("Your microbiome contributes to sufficient metabolism of all macronutrients")
  }                                                     
  
  ##Previously would choose the macronutrient with the lowest value (even if above 25%)
  ## print(colnames(foodTypeDF)[which(foodTypeDF == min(foodTypeDF), arr.ind = TRUE)[ , 2]])
  #changed print to paste to remove console output
  bloating_output <- if(length(below25foodtypes) > 0){  #if at least one is poor 
    below25foodtypes_list <- paste(below25foodtypes, collapse = ", ") #print its name
    below25foodtypes_list <- stri_replace_last(below25foodtypes_list, fixed = ',', ' and') #replace last comma with and
    below25food_output <- paste("Symptoms of bloating/acid reflux after consuming high amounts of",  below25foodtypes_list, "rich foods are associated with your microbiome profile", sep = " ")  
    paste(below25food_output)
  } else {                                             #if none are poor
    paste("Symptoms of bloating/acid reflux are not related to your gut microbiome profile")
  }    
  
  #changed print to paste to remove console output
  Aerobic_output <- if(All_GutPotential_SegmentsDF$AerobicEndurance_col < 40){  #change as of Feb 9, 2024
    paste("not supportive") 
  } else if(All_GutPotential_SegmentsDF$AerobicEndurance_col >=  60) {
    paste("very supportive")
  } else {                                             
    paste("moderately supportive")
  }  
  
  taxonomy_Filename <- paste(output_dir, "/", inputFilePrefix, "_taxonomy.xlsx", sep = "")
  number_of_pathogens <- read_xlsx(taxonomy_Filename, sheet = "Species") %>% 
    filter(Species == "Helicobacter pylori" | Species == "Clostridium difficile" | Species == "Clostridium perfringens") %>% 
    tally()
  # B. fragilis removed 
  
  
  joint_pain <- diseasepredictionDF %>% select(SampleID, Arthritis) 
  
  #changed print to paste to remove console output
  joint_pain_output <- if (joint_pain$Arthritis == "high risk" | joint_pain$Arthritis == "moderate risk") {
    paste("Symptoms")  #changed March 25, 2024 to two categories 
  } else {
    paste("No symptoms")
  }
  
  
  probioticBacteriaFileName <- paste(output_dir, "/ProbioticBacteria/", inputFilePrefix, "_ProbioticBacteria.txt", sep = "")
  probiotic_recommendations <- read_tsv(probioticBacteriaFileName) %>% 
    filter(You < 0.5*Healthy) %>%  #If customer % is less than a half of the median of the reference population 
    select(ProbioticBacteria) 
  
  
  
  LowVitamins <- read_tsv(GP_samplename) 
  names(LowVitamins) <- str_replace(names(LowVitamins), "Vitamin", "Vitamin ")
  names(LowVitamins) <- str_replace(names(LowVitamins), "_col", "")
  
  LowVitamins <- LowVitamins %>% 
    mutate(across(.fns = ~ifelse(. <=25, "poor", "good"))) %>% 
    select(starts_with("Vitamin")) %>% 
    select(where(~ any(. == "poor"))) %>% 
    names()
  
  #changed print to paste to remove console output
  Vitamin_output_for_df <- if(length(LowVitamins) >0) {
    Vitamin_list <- paste(LowVitamins, collapse = ", ")
    Vitamin_list <- stri_replace_last(Vitamin_list, fixed = ',', ' and') # replace last comma with and
    Vitamin_output <- paste("Your gut microbiome profile indicates low synthesis of the following vitamins:", Vitamin_list, sep =  " ")
    paste(Vitamin_output) 
  } else {
    paste("Your microbiome profile shows strong levels of vitamin synthesis")
  }
  
  
  
  
  inflammatoryBacteriaFileName <- paste(output_dir, "/inflammatoryBacteria/", inputFilePrefix, "_inflammatoryBacteria.txt", sep = "")
  
  
  ###immune_score <- read_tsv(GP_samplename) %>% 
  ### select(OverallImmuneScore_col)  %>% 
  ### mutate(OverallImmuneScore_Dec = OverallImmuneScore_col/100) %>% 
  ###  select(-OverallImmuneScore_col)
  
  
  inflamm_diseases <- diseasepredictionDF %>% 
    select(SampleID, IBS, IBD, `Clostridium difficile infection`) %>%
    mutate( across(
      .cols = everything(),
      ~str_replace( ., "healthy", "1" ))) %>% 
    mutate( across(
      .cols = everything(),
      ~str_replace( ., "low risk", "0.75" ))) %>% 
    mutate( across(
      .cols = everything(),
      ~str_replace( ., "moderate risk", "0.5" ))) %>% 
    mutate( across(
      .cols = everything(),
      ~str_replace( ., "high risk", "0" ))) %>% 
    #filter(SampleID == inputFilePrefix) %>% already selected in diseasepredictionDF
    mutate_at(c(2:4), as.numeric)
  
  ## attention: column names have been shortened 
  inflammation_score_df <- inflamm_diseases 
  
  inflammatory_bacteria_DF <- read_tsv(inflammatoryBacteriaFileName) %>% 
    mutate(inflamBacteria_Score = ifelse(You >= HealthyQ3, 0, 1))
  
  inflammation_score_df$InflammatoryBacteriaScore <- mean(inflammatory_bacteria_DF$inflamBacteria_Score)
  
  
  inflammation_score_df$Inflammation_score <- rowMeans(inflammation_score_df[2:5], na.rm = TRUE)
  
  inflammation_output <- if(inflammation_score_df$Inflammation_score < 0.4){  #Make sure it is 0-1 scale - changed Dec 15, 2023 
    paste("poor") 
  } else if(inflammation_score_df$Inflammation_score >=  0.6) {
    paste("good")
  } else {                                             
    paste("moderately good")
  }  
  
  
  
  all_genera_filename <- paste(output_dir,"/",inputFilePrefix,"_all_genera.txt", sep = "")
  all_genera_df <- read_tsv(all_genera_filename)
  top_genus <- all_genera_df %>% 
    filter(Genus != "Other") %>% 
    arrange(desc(Percentage)) %>% 
    head(n=1)
  
  
  
  
  ##food 
  
  foodRec_Filename <- paste("./", "Food_recommendations_", BatchCode, 
                            ".xlsx", sep = "")
  
  #Make sure food list in foodRec and food_key match 
  food_key <- read_csv(FoodMappingFile)
  
  foodRec_DF <- read_excel(foodRec_Filename, sheet = inputFilePrefix) %>% 
    # plyr::rename(c("0" = "Food")) %>% already renamed in the file 
    plyr::rename(c("Amount" = "Recommendation")) %>% 
    plyr::rename(c("...1" = "N")) 
  
  
  foodRec_DF <- left_join(foodRec_DF, food_key, by = "Food") %>% 
    filter(Food != "Turkey tail mushroom" & 
             Food != "Lion's mane mushroom" & 
             Food !=  "White rice" & 
             Food !=  "White bread" & 
             Food !=  "Pork" & 
             Food !=  "Beef" & 
             Food !=  "Lamb" & 
             Food !=  "Egg yolk" & 
             Food !=  "Egg white" & 
             Food !=  "Chicken" & 
             Food !=  "Whole grain bread") %>% 
    filter(Category == "Fruit" | Category == "Vegetable" |  Category == "Legumes" |         Category == "Protein")
  #only fruit, veg, legumes and fish appear in roadmaps
  
  #filtering certain foods that should/should not be shown in roadmaps 
  #they will still appear in food tables 
  
  increase_foods <- foodRec_DF %>% filter(Recommendation == "Increased consumption")
  decrease_foods <- foodRec_DF %>% filter(Recommendation == "Decreased consumption")
  
  #loop to determine whether need to include moderate consumption foods to roadmaps
  food_list_for_roadmap <- if(nrow(increase_foods) >= 10){      
    foodRec_DF %>% 
      filter(Recommendation == "Increased consumption") %>%      
      sample_n(10) %>% 
      select(Food)                                              
  } else {                                                                           
    foodRec_DF %>% 
      filter(Recommendation == "Moderate consumption" | Recommendation == "Increased consumption") %>% 
      sample_n(10) %>% 
      select(Food) 
  }
  
  #if else logic added to check number of Decrease foods 
  
  two_fooods_summary_page <- if (nrow(decrease_foods) >= 2) {
    foodRec_DF %>%
      filter(Recommendation == "Decreased consumption") %>% 
      sample_n(2) %>% 
      select(Food) 
  } else {    
    foodRec_DF %>%
      filter(Recommendation == "Decreased consumption" | Recommendation == "Moderate consumption") %>%
      filter(!Food %in% food_list_for_roadmap$Food) %>% 
      sample_n(2) %>% 
      select(Food) 
  }
  
  
  
  
  combined_df <- data.frame()
  combined_df = data.frame(matrix(nrow = 1, ncol = 0)) 
  
  combined_df$Name = diversity$SampleID #Box1
  combined_df$GC_score_category = GutCheckScoreDF$GutCheckScore_category #Box2 and Box18
  combined_df$GC_score_value = round(GutCheckScoreDF$GutCheckScore, 2)
  combined_df$Diversity_category = diversity$shannon_category #Box3 and Box22 and Box19
  combined_df$Similarity_percent = similarity$Similarity_percent #Box4 and Box21
  combined_df$Disease_prediction <- disease_prediction_output #Box5 
  combined_df$Low_abundance_beneficial_bacteria <- paste(number_of_low_abundance_benficial_bacteria) #Box 6
  combined_df$Overgrowing_beneficial_bacteria <- paste(number_of_overgrowing_beneficial_bacteria)
  combined_df$Gluten_metabolism  <- paste(Gluten_output) #Box7 
  combined_df$VitaminSynthesis <- paste(Vitamin_output_for_df) #Box8
  combined_df$Weight_symptom <- paste(obesity_output) #Box9 and Box34
  combined_df$Joint_pain <- paste(joint_pain_output) #Box10
  combined_df$Anxiety_symptom <- Symptoms_DF$Anxiety_symptom
  combined_df$Aerobic_output <- paste(Aerobic_output) 
  combined_df$Inflammation_resistance <- paste(inflammation_output) 
  combined_df$Pathogens_number <- paste(number_of_pathogens) #Box14 
  combined_df$Food_to_moderate <- paste(noquote(paste(two_fooods_summary_page$Food, collapse = " and ")))
  combined_df$Food_for_roadmap <- paste(noquote(paste(food_list_for_roadmap$Food, collapse = ", "))) #update July 4th, 2024 
  combined_df$Reduce_food_type  <- paste(reduceFoodType, collapse=" and ") #Box16
  combined_df$Probiotic_recommendations <-  paste(noquote(paste(probiotic_recommendations$ProbioticBacteria, collapse = ", "))) #Box17
  combined_df$Diversity_score = round(diversity$shannon_entropy, 2) # will place this field next to diversity category 
  combined_df$Abdominal_pain_symptom = Symptoms_DF$Abdominal_pain_symptom
  combined_df$Digestion_symptom <- Symptoms_DF$Digestion_symptom 
  combined_df$Reduce_food_type  <- paste(reduceFoodType) 
  combined_df$Bloating_symptom_food_type <- paste(bloating_output)
  combined_df$Sleep_symptom <- Symptoms_DF$Sleep_symptom 
  combined_df$Stress_symptom <- Symptoms_DF$Stress_symptom
  combined_df$Fatigue_symptom <- Symptoms_DF$Fatigue_symptom
  combined_df$Gut_lining_cell_nutrients_symptom <- Symptoms_DF$Gut_lining_cell_nutrients_symptom
  combined_df$Gut_lining_mucus_symptom <- Symptoms_DF$Gut_lining_mucus_symptom
  combined_df$Joint_pain <- paste(joint_pain_output) 
  combined_df$Skin_symptom <- Symptoms_DF$Skin_symptom 
  # combined_df$Fat_metabolism_score <- All_GutPotential_SegmentsDF$Fat_col
  # combined_df$Carbohydrate_metabolism_score <- All_GutPotential_SegmentsDF$Carbohydrate_col
  # combined_df$Protein_metabolism_score <- All_GutPotential_SegmentsDF$Protein_col
  # combined_df$Gluten_metabolism_score <- All_GutPotential_SegmentsDF$Gluten_col
  # combined_df$Lactose_metabolism_score <- All_GutPotential_SegmentsDF$Lactose_col
  # combined_df$Vitamin_synthesis_score <- All_GutPotential_SegmentsDF$OverallVitaminSynthesis_col
  # combined_df$Serotonin_synthesis_score <- All_GutPotential_SegmentsDF$Serotonin_col
  # combined_df$GABA_synthesis_score <- All_GutPotential_SegmentsDF$GABA_col
  # combined_df$IPA_synthesis_score <- All_GutPotential_SegmentsDF$IPA_col
  # combined_df$Acetate_score <- All_GutPotential_SegmentsDF$Acetate_col
  # combined_df$Butyrate_score <- All_GutPotential_SegmentsDF$Butyrate_col
  # combined_df$Propionate_score <- All_GutPotential_SegmentsDF$Propionate_col
  combined_df$Top_genus <- paste(top_genus$Genus)
  combined_df$Top_genus_percentage <- paste(round(top_genus$Percentage, 0), "%", sep = "")
  
  
  All_GutPotential_SegmentsDF$Name <- inputFilePrefix #new March 18th, 2026
  combined_df <- combined_df %>% inner_join(All_GutPotential_SegmentsDF, by = "Name") #new March 18th, 2026
  
  print(combined_df)
  
  combinedDF_output_Filename <- paste(output_dir, "/", inputFilePrefix, "_GutCheck_Report_Input.tsv", sep = "")
  
  write.table(combined_df, combinedDF_output_Filename, quote = FALSE, row.names = FALSE, sep = "\t")
  
  short_report_DF_file_name <- paste(output_dir, "/", inputFilePrefix, "_GutCheck_reportInputs_Short.tsv", sep = "")
  
  read_tsv(combinedDF_output_Filename) %>% 
    select(Name, GC_score_category, Diversity_score, Similarity_percent, Top_genus, Top_genus_percentage, Acetate_col:VitaminK2_col, Probiotic_recommendations, Abdominal_pain_symptom:Skin_symptom) %>% 
    write.table(short_report_DF_file_name, quote = FALSE, row.names = FALSE, sep = "\t")
  
  
  ### This is to crop Diversity figure and Gutcheck score figure
  ### add this chunk of code to the print_input_text_for_reports function and uncomment if using HTML version of reports  
  
  function_crop <- function(inputFilePrefix) {
    
    output_dir <- paste("./",inputFilePrefix, sep="")
    
    dir.create(output_dir)  
    
    png_files_ls = list.files(path=output_dir, pattern=".*_figure.png|.*_shortReport.png", recursive = TRUE, full.names = TRUE)
    
    #print(png_files_ls)
    for (i in png_files_ls) {
      
      imported_png <- image_read(i)  
      image_write(image_trim(imported_png), i)  
      
    }
    
  }
  
  ## now run the function! 
  
  function_crop(inputFilePrefix)
  
  short_report_files_ls = list.files(path=output_dir, pattern=".*_shortReport.png", recursive = TRUE, full.names = TRUE)
  
  
  for (i in short_report_files_ls) {
    imported_png <- image_read(i) 
    
    final_img <- imported_png %>%
      image_background(color = "#130083") %>%
      image_flatten()
    
    image_write(final_img, i)  
    
    
  }
  
  
}





calculate_phylum_similarity <- function(Product = "Gutcheck", Batch) {
  
  ## Reference data  
  
  # ra_phyla_776 <- if (Product == "Gutcheck") {
  #   read_tsv("C:\\Users\\farid\\OneDrive/Desktop\\NucliqBiologics\\/Reference_datasets/AGP/rel-phyla-transposed-AGP-GG2.txt")
  # }  else if (Product == "GutcheckIndia") {
  #   read_tsv("C:\\Users\\farid\\OneDrive/Desktop\\NucliqBiologics\\/Reference_datasets/SouthAsia/SouthAsiaCohort/filteredRelAbunTables/rel-phyla-transposed-d2table-noblooms_SouthAsiasamples_healthy_min5s5f_noLowFreqSamples.txt")
  # }  else if (Product == "Vivo") {
  #   read_tsv("C:\\Users\\farid\\OneDrive/Desktop\\NucliqBiologics\\/Reference_datasets/Vivo/rel-phyla-transposed_AGP_VIVO.txt")
  # } else {
  #   "No product specified"
  # }

  ra_phyla_776 <- if (Product == "Gutcheck") {
  read_tsv(file.path(getwd(), "rel-phyla-transposed-AGP-GG2.txt"))
  } else if (Product == "GutcheckIndia") {
  read_tsv(file.path(getwd(), "rel-phyla-transposed-d2table-noblooms_SouthAsiasamples_healthy_min5s5f_noLowFreqSamples.txt"))
  } else if (Product == "Vivo") {
  read_tsv(file.path(getwd(), "rel-phyla-transposed_AGP_VIVO.txt"))
} else {
  "No product specified"
}
  
  
  mean_phyla_ra_776 <- as.data.frame(t(as.data.frame(colMeans(ra_phyla_776[2:length(ra_phyla_776)]))))
  
  # Batch data 
  
  
  ra_phyla_sample_data <- read_tsv("./rel-phyla.txt", skip = 1) %>% 
    pivot_longer(cols = -`#OTU ID`, names_to = "SampleID", values_to = "Value") %>%
    pivot_wider(names_from = `#OTU ID`, values_from = Value) 
  
  ra_phyla_sample_data <- if (Product == "Gutcheck") {
    
    # ra_phyla_sample_data %>% filter(grepl("NBIGC", SampleID, ignore.case = TRUE))
    ra_phyla_sample_data %>% filter(grepl("NBGIC|NBIGC|NBIVV", SampleID, ignore.case = TRUE))
  }   else if (Product == "Vivo") {
    
    ra_phyla_sample_data %>% filter(grepl("NBI[G,V]V", SampleID, ignore.case = TRUE))
  }   else if (Product == "GutcheckIndia") {
    
    ra_phyla_sample_data 
  } else {
    "No product specified"
  }
  
  ra_phyla_sample_data <- ra_phyla_sample_data %>% remove_rownames %>% column_to_rownames(var="SampleID")
  
  
  phyla_sample_vs_reference <- dplyr::bind_rows(mean_phyla_ra_776,ra_phyla_sample_data) %>% 
    replace(is.na(.), 0)
  
  
  
  distB_phyla <- vegdist(phyla_sample_vs_reference, method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
                         na.rm = TRUE) 
  
  # similarityFileName <- paste("Phylum_level_Similarity_", Product, "_", Batch, "Project.txt", sep = "")
  similarityFileName <- paste("Phylum_level_Similarity_GutCheck_", Batch, "Project.txt", sep = "")
  
  melt(as.matrix(distB_phyla), varnames = c("row", "col")) %>% 
    filter(row == "colMeans(ra_phyla_776[2:length(ra_phyla_776)])") %>% 
    mutate(Similarity = 1-value) %>% 
    mutate(Similarity_round = round(Similarity, 2)) %>% 
    filter(!str_detect(col, "colMeans")) %>% 
    select(col, Similarity_round) %>% 
    plyr::rename(c("col" = "SampleID")) %>% 
    write.table(similarityFileName, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  
}






translateInputTable <- function(inputFilePrefix) {
  
  output_dir <- paste("./",inputFilePrefix,sep="")
  dir.create(output_dir)
  
  combinedDF_output_Filename <- paste(output_dir, "/", inputFilePrefix, "_GutCheck_Report_Input.tsv", sep = "") 
  combinedDF_output_FilenameCopy <- paste(output_dir, "/", inputFilePrefix, "_GutCheck_Report_Input_originalEnglishversion.tsv", sep = "") 
  combined_df <- read_tsv(combinedDF_output_Filename)
  file.rename(combinedDF_output_Filename, combinedDF_output_FilenameCopy)
  
  ## Things to look out for in a dictionary: 
  ### For replacements, order matters - put longer inputs (e.g. sentences) first in the list, and words further down the list (these are replaced last)
  ### Escape the parentheses 
  ### Make sure punctuation marks and spaces are matching between input file and dictionary   
  
  ## For now use this this file but then move it to another directory
  #eng_fr_dict <- read_xlsx("C:/Users/farid/OneDrive/Desktop/NucliqBiologics/Dev/French_translations.xlsx", sheet = "Input_file_text_options") 
  eng_fr_dict <- read_xlsx(file.path(getwd(), "French_translations.xlsx"), sheet = "Input_file_text_options")
  eng_fr_dict 
  replacements <- c(eng_fr_dict$FR)
  names(replacements) <- c(eng_fr_dict$ENG)
  
  combined_df %>% 
    pivot_longer( everything(),  
                  names_to = "Headers", values_to = "Text", values_transform = list(Text = as.character)) %>% 
    mutate(FR_text = str_replace_all(Text, 
                                     pattern = replacements)) %>% 
    mutate(FR_text = str_remove_all(FR_text, 
                                    pattern = "spaceplaceholder")) %>% 
    select(-Text) %>% 
    pivot_wider(names_from = "Headers", values_from = "FR_text") %>% 
    write.table(combinedDF_output_Filename, sep = "\t", quote = FALSE, row.names = FALSE)
  
  
  
  
}








make_FR_subdir <- function(inputFilePrefix) {
  
  folder_path <- paste(getwd())
  
  FRpath <- paste(folder_path, "/FR/", sep="")
  dir.create(FRpath) 
  
  currentfolder <- paste(folder_path, "/", inputFilePrefix, sep="")    
  newlocation   <- paste(folder_path, "/FR/", sep="")  
  file.copy(from=currentfolder, to=newlocation, 
            overwrite = TRUE, recursive = TRUE, 
            copy.mode = TRUE)  
  
  
}




make_ENG_subdir <- function(inputFilePrefix) {
  
  folder_path <- paste(getwd())
  
  ENGpath <- paste(folder_path, "/ENG/", sep="")
  dir.create(ENGpath) 
  
  currentfolder <- paste(folder_path, "/", inputFilePrefix, sep="")    
  newlocation   <- paste(folder_path, "/ENG/", sep="")  
  file.copy(from=currentfolder, to=newlocation, recursive = TRUE,
            copy.mode = TRUE)  
  
  
  
}









arrange_folder_structure <- function(BatchCode, language = "ENG") {
  
  
  wd_path <- paste(getwd())
  
  
  folder_path <- if (language == "ENG") { 
    
    folder_path <- paste(wd_path, "/ENG", sep = "")
    
  } else {
    
    folder_path <- paste(wd_path, "/FR", sep = "")
    
  }
  
  
  # Create list of text files
  txt_files_ls = list.files(path=folder_path, pattern="*_GutCheck_Report_Input.tsv", recursive = TRUE) #recursive gets files from all subdirectories
  
  txt_files_ls
  
  txt_files_df <- lapply(txt_files_ls, function(x) {read_tsv(file = x)})
  # Combine them
  combined_df_with_txt_files <- do.call("rbind", lapply(txt_files_df, as.data.frame)) 
  combined_df_with_txt_files
  
  
  
  #Write output Excel file
  
  output_Excel_file_name <- paste(folder_path, "/Batch_", BatchCode, "_GutCheck_Report_Input.xlsx", sep = "")
  write_xlsx(combined_df_with_txt_files, output_Excel_file_name)
  
  
  ### same thing for Short Report
  
  # Create list of text files
  txt_files_ls = list.files(path=folder_path, pattern="*_GutCheck_reportInputs_Short.tsv", recursive = TRUE) #recursive gets files from all subdirectories
  
  txt_files_ls
  
  txt_files_df <- lapply(txt_files_ls, function(x) {read_tsv(file = x)})
  # Combine them
  combined_df_with_txt_files <- do.call("rbind", lapply(txt_files_df, as.data.frame)) 
  combined_df_with_txt_files
  
  
  
  #Write output Excel file
  
  output_Excel_file_name <- paste(folder_path, "/Batch_", BatchCode, "_GutCheck_Report_Input_SHORT.xlsx", sep = "")
  write_xlsx(combined_df_with_txt_files, output_Excel_file_name)
  
  
  
  #Next step: put all input files for pathogenic/inflammatory table into one folder 
  
  list_inflam_files <- list.files(path=folder_path, pattern="*inflammatoryBacteria", recursive = TRUE)
  list_taxonomy_files <- list.files(path=folder_path, pattern="*taxonomy.xlsx", recursive = TRUE)
  
  pathogen_dir_path <- paste(folder_path, "/pathogenic_bacteria_table/", sep = "") 
  
  dir.create(pathogen_dir_path) 
  for (f in list_inflam_files) {
    file.copy(f, pathogen_dir_path)
  }
  for (f in list_taxonomy_files) {
    file.copy(f, pathogen_dir_path)
  }
  
  
  beneficial_dir_path <- paste(folder_path, "/beneficial_bacteria_table/", sep = "")
  
  list_benef_files <- list.files(path=folder_path, pattern="*BeneficialBacteria.txt", recursive = TRUE)
  dir.create(beneficial_dir_path)
  for (f in list_benef_files) {
    file.copy(f, beneficial_dir_path)
  }
  
  
  probiotic_dir_path <- paste(folder_path, "/probiotic_bacteria_table/", sep = "")
  
  list_probioticbacteria_files <- list.files(path=folder_path, pattern="*ProbioticBacteria.txt", recursive = TRUE)
  dir.create(probiotic_dir_path)
  for (f in list_probioticbacteria_files) {
    file.copy(f, probiotic_dir_path)
  }
  
  
  ### Output disease predictions file with set column names 
  
  
  diseasepredictionFileName <- paste("./Comparison_all_diseases_", BatchCode,".xlsx", sep = "")
  diseasepredictionFileName2 <- paste(folder_path, "/Comparison_all_diseases_", BatchCode,".xlsx", sep = "")
  
  disease_prediction <- read_xlsx(diseasepredictionFileName) %>% 
    select(`OTU ID`, `Obesity Final Status`, `IBS Final Status`, `IBD Final Status`, `Diabetes Final Status`, `CVD Final Status`, `Arthritis Final Status`, `MASLD Final Status`, `Clostridium difficile infection Final Status`) %>%
    rename_with(~ str_remove(., " Final Status"), everything()) %>% 
    plyr::rename(c("OTU ID" = "SampleID")) 
  
  write_xlsx(disease_prediction, diseasepredictionFileName2)
  
  
  
  ### Clean up 
  
  wd_path <- paste(getwd())
  files_to_delete = list.files(path=wd_path, pattern="NBI", recursive = FALSE, full.names = FALSE, include.dirs = TRUE)
  
  folders_to_delete <- files_to_delete[ file.info(files_to_delete)$isdir ]
  
  for (directory_name in folders_to_delete) {
    
    check_pathENG <- paste("ENG/", directory_name, sep = "")
    check_pathFR <- paste("FR/", directory_name, sep = "")
    
    if (dir.exists(check_pathENG)) {
      print(paste(check_pathENG, "copy exists. Proceed to delete", directory_name, sep = " "))
      unlink(directory_name, recursive = TRUE) #WARNING this deletes entire dir
    } else if (dir.exists(check_pathFR)){
      print(paste(check_pathFR, "copy exists. Proceed to delete", directory_name, sep = " "))
      unlink(directory_name, recursive = TRUE) #WARNING this deletes entire dir
    } else {
      print("No copy found. Please double check")
    } 
  }
  
  
  
  
  
}






