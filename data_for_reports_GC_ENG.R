library(tidyverse)
library(pdftools)
library(readxl)
library(writexl)
library(magick)



args<-commandArgs(TRUE)

# Check the number of arguments and provide defaults if necessary
batch_id <- if (length(args) == 0) {
  stop("Provide batch ID.", call. = FALSE)
} else if (length(args) == 1) {
  args[1] 
}

biomtable <- read_tsv("./biom_table_with_taxonomy.txt")


ncolbiommunus1 <- ncol(biomtable)-1


biomtable[2:ncolbiommunus1] <- lapply(as.vector(biomtable[2:ncolbiommunus1]), function(x){x/sum(x, na.rm = TRUE)})


metadata_DF <- read_tsv("./metadata.txt")

####################################################
# ! ! ! This is a workaround --> 
#  later the preferred language should be taken directly from questionnaire 
#  this entire section will have to be removed 
####################################################

### Make all reportrs in ENG by default 
metadata_DF['Language'] <- 'ENG' ### Remove this line when language question has been added

df_GC_ENG <- metadata_DF %>%  filter(Language == "ENG") %>% filter(Sample_type == "Stool") %>% filter(Product == "Gutcheck")
sample_list_GC_ENG <- dput(df_GC_ENG$`sample-id`)


calculate_phylum_similarity("Gutcheck", batch_id)

for (i in sample_list_GC_ENG) {
  getSampleTaxonomicInfo_GG2_v3(i, language = "ENG")
  getGutPotentialScores_v3_GG2(i, Template = "Gutcheck")
  getinflammatoryBacteria_GG2_v3(i)
  getProbioticBacteria_GG2_v3(i)
  getBeneficialBacteria_GG2_v3(i)
  checkSymptoms_v2_GG2(i)
  calculateGutScore(i, batch_id)
  printInputTextboxesForReport_GG2(i, BatchCode = batch_id)
}



# run only after all samples finished running
for (i in sample_list_GC_ENG) {
  # Fix: Pass 'i' (the individual sample) instead of the whole list
  make_ENG_subdir(i)
}




# # run only after ENG directory was created and all files have been copied

arrange_folder_structure(batch_id, language = "ENG")

print("Analysis for batch completed. Check the files in the output folders.")
