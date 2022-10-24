

#LOAD UPON STARTUP
#Memory Limits
#Working Directory
#File Saving Path
#Varies Among Machines, Memory, and RAM


#Loading Packages & Increasing Memory size + Uploading to Directory
memory.limit(size=56000) #Setting Memory Limits
setwd("") #Setting the Working Directory
path <- "" #Setting File Pathways


#USED PROGRAMS AND PACKAGES

#-----------------------------------
BiocManager::install("GEOquery")
install.packages("BiocManager")
BiocManager::install()

install.packages("tidyverse")
install.packages("janitor")
install.packages("dplyr")
install.packages("gdata")
install.packages("lawstat")
install.packages("car")
BiocManager::install("GenomicRanges")
install.packages("utilities") 


#----------------------------------- Reload upon start up
library(tidyverse)
library(janitor)
library(readr)
library(BiocManager)
library(Biobase)
library(GEOquery)
library(gdata)
library(GenomicRanges)
library(data.table)
library(utils)
library(utilities)
library(stats)
library(lawstat)
library(car)
library(boot)
library(stats)
library(pairwiseComparisons)
#-----------------------------------


#DATA WRANGLING


#Attempting Variable Isolation from the NCBI Gene Expression Omnibus (Data repository) - raw/user data supplied

GDS596 <- getGEO('GDS596', destdir = ".") #Every other Normal Tissue

GDS3716 <- getGEO('GDS3716', destdir = ".") #Normal Breast Tissue, coming from a different data set, not sure that MAS5 protocols are exactly the same.

GDS1965 <- getGEO('GDS1965', destdir = ".") #Normal Melanocyte Tissue

GDS1096 <- getGEO('GDS1096', destdir = ".") #Normal Tissues (Lung) First Attempt
GDS3257 <- getGEO('GDS3257', destdir = ".") #Normal Tissues (Lung) Second Attempt

GDS505 <- getGEO('GDS505', destdir = ".") #First Normal Renal Tissue
GDS507 <- getGEO('GDS507', destdir = ".") #Second Normal Renal Tissue

GDS3057 <- getGEO('GDS3057', destdir = ".") #Normal Bone Marrow / Myeloid

GDS1375 <- getGEO('GDS1375', destdir = ".")#Normal Melanocyte Tissue 2



#Viewing the data and checking the first few rows/columns.

View(GDS596)
head(Meta(GDS596))

View(GDS3716)
head(Meta(GDS3716))

View(GDS1965)
head(Meta(GDS1965))

View(GDS1096)
head(Meta(GDS1096))

View(GDS3257)
head(Meta(GDS3257))

View(GDS505)
head(Meta(GDS505))

View(GDS507)
head(Meta(GDS507))

View(GDS3057)
head(Meta(GDS3057))

View(GDS1375)
head(Meta(GDS1375))


#Expression data from the data sets (LOG TRANSFORMATION OCCURS IN CODE FURTHER DOWN)

eset <- GDS2eSet(GDS596, do.log2 =FALSE) #Varying data for tissue types, Log2 transformation doesn't work.
View(eset)

eset2 <- GDS2eSet(GDS3716, do.log2 =FALSE) #Data for Breast tissue types
View(eset2)

eset3 <- GDS2eSet(GDS1965, do.log2 = FALSE) #Data for Melanocyte tissue types
View(eset3)

eset4 <- GDS2eSet(GDS1096, do.log2 = FALSE) # Data for Lung tissue types
View(eset4)

eset5 <- GDS2eSet(GDS3257, do.log2 = FALSE) # Data for Lung tissue types 2
View(eset5)

eset6 <- GDS2eSet(GDS505, do.log2 = FALSE) # Data for Renal tissue types 
View(eset6)

eset7 <- GDS2eSet(GDS507, do.log2 = FALSE) # Data for Renal tissue types 2
View(eset7)

eset8 <- GDS2eSet(GDS3057, do.log2 = FALSE) #Data for Bone Marrow / Myeloid
View(eset8)

eset9 <- GDS2eSet(GDS1375, do.log2 = FALSE)
View(eset9)


#Viewing the phenotype data tables from the GDS files

phenotype_data <- eset@phenoData@data #Normal Tissue samples
data_table_GDS596 <- GDS596@dataTable@table #Normal Tissue samples Table

breastphenotype <- eset2@phenoData@data #Normal Breast Tissue samples 
data_table_GDS3716 <- GDS3716@dataTable@table #Normal tissue, only isolate the reduction samples (average them as a reference for normal)

melanocytephenotype <- eset3@phenoData@data #Normal Melanocyte Tissue samples
data_table_GDS1965 <- GDS1965@dataTable@table #Normal Melanocyte Tissue Table

melanocytephenotype2 <- eset9@phenoData@data #Normal Melanocyte Tissue samples 2
data_table_GDS1375 <- GDS1375@dataTable@table #Normal Melanocyte Tissue Table 2

lungphenotype <- eset4@phenoData@data #Normal Lung Tissue samples
data_table_GDS1096 <- GDS1096@dataTable@table #Normal Lung Tissue Table

lungphenotype2 <- eset5@phenoData@data #Normal Lung Tissue samples
data_table_GDS3257 <- GDS3257@dataTable@table #Normal Lung Tissue Table 2

renalphenotype <- eset6@phenoData@data #Normal Renal Tissue samples
data_table_GDS505 <- GDS505@dataTable@table #Normal Renal Tissue Table 

renalphenotype2 <- eset7@phenoData@data #Normal Renal Tissue samples
data_table_GDS507 <- GDS507@dataTable@table #Normal Renal Tissue Table 2

myeloidphenotype <- eset8@phenoData@data #Normal bone marrow / Myeloid Tissue samples
data_table_GDS3057 <- GDS3057@dataTable@table 


#IF ESET DOESN'T WORK (TO LOG TRANSFORM) RUN THIS CODE!

#data_table_GDS596 <- data_table_GDS596 %>% mutate(across(GSM18927:GSM18964, log, 2))
#data_table_GDS3716 <- data_table_GDS3716 %>% mutate(across(GSM512539:GSM512565, log, 2)) #Log 2 transforming data without messing up the first two columns.
#data_table_GDS1965 <- data_table_GDS1965 %>% mutate(across(GSM102065:GSM102073, log, 2))
#data_table_GDS1096 <- data_table_GDS1096 %>% mutate(across(GSM44692:GSM44684, log, 2))
#data_table_GDS3257 <- data_table_GDS3257 %>% mutate(across(GSM254629:GSM254641, log, 2))
#data_table_GDS505 <- data_table_GDS505 %>% mutate(across(GSM11814:GSM12444, log, 2))
#data_table_GDS507 <- data_table_GDS507 %>% mutate(across(GSM11815:GSM12448, log, 2))
#data_table_GDS3057 <- data_table_GDS3057 %>% mutate(across(GSM239371:GSM239338, log, 2))
#data_table_GDS1375 <- data_table_GDS1375 %>% mutate(across(GSM71671:GSM71718, log, 2))


#Using Pivot longer to alter appearance of data_table to visualize what it would appear as in a database.
#data_table_GDS3716 %>% pivot_longer(cols = GSM512539:GSM512565, names_to= "sample", values_to="log2expression")

#Isolating the normal tissue samples from their respective phenotype tables
#Looking for specific tissue samples that match available tissue from NCI-60

sel_phenotype_data <- phenotype_data[str_detect(phenotype_data$tissue, "fetal-|lung|epithelial|kidney|breast|CD33+|bone marrow"),]

sel_phenotype_data1 <- sel_phenotype_data[3:dim(sel_phenotype_data)[1],]
dim(sel_phenotype_data)

sel_breast <- breastphenotype[str_detect(breastphenotype$specimen, "reduction mammoplasty"),]

sel_melanocyte <- melanocytephenotype[str_detect(melanocytephenotype$disease.state, "normal"),]
sel_melanocyte2 <- melanocytephenotype2[str_detect(melanocytephenotype2$disease.state, "normal"),]

sel_lung1 <- lungphenotype[str_detect(lungphenotype$tissue, "lung|epithelial|kidney|breast|CD33+|bone marrow"),] %>%
  filter(tissue == "bone marrow" | tissue == "breast" | tissue == "renal" | tissue == "lung" | tissue == "kidney")

sel_lung2 <- lungphenotype2 %>%
  filter(tissue == "normal" & individual == "never smoker" & disease.state == "stage I") %>% glimpse()

sel_renal1 <- renalphenotype[str_detect(renalphenotype$disease.state, "normal"),]
sel_renal2 <- renalphenotype2[str_detect(renalphenotype2$disease.state, "normal"),]

sel_myeloid <- myeloidphenotype %>%
  filter(disease.state == "normal" & cell.type == "bone marrow")

----------------------------------------------------------------------------------------------
  
  
#General Tissues
#Changing the current data table for each tissue to a tibble and isolating specific columns while also applying a log 2 transformation across the GDS accessions (action ,Re used to help in log 2 transformation and does not work properly without ,Re)


sphenotype <-  sel_phenotype_data1[,1]
f_GDS596 <-  as_tibble(data_table_GDS596) %>% select(c("ID_REF",sphenotype)) 
f_GDS596 <- f_GDS596 %>% mutate(across(starts_with("GSM"), log, 2)) 
f_GDS596 <- f_GDS596 %>% mutate(across(starts_with("GSM"), Re))


#Breast Tissue
sbreast <- sel_breast[,1]
f_GDS3716 <- as_tibble(data_table_GDS3716) %>% select(c("ID_REF",sbreast)) #There were warnings here
f_GDS3716 <- f_GDS3716 %>% mutate(across(starts_with("GSM"), log, 2))
f_GDS3716 <- f_GDS3716 %>% mutate(across(starts_with("GSM"), Re))


#Melanocyte Tissues
smelanocyte <- sel_melanocyte[,1]
smelanocyte2 <- sel_melanocyte2[,1]
f_GDS1965 <- as_tibble(data_table_GDS1965) %>% select(c("ID_REF",smelanocyte))
f_GDS1965 <- f_GDS1965 %>% mutate(across(starts_with("GSM"), log, 2))
f_GDS1965 <- f_GDS1965 %>% mutate(across(starts_with("GSM"), Re))


f_GDS1375 <- as_tibble(data_table_GDS1375) %>% select(c("ID_REF",smelanocyte2))
f_GDS1375 <- f_GDS1375 %>% mutate(across(starts_with("GSM"), log, 2))
f_GDS1375 <- f_GDS1375 %>% mutate(across(starts_with("GSM"), Re))


#General Tissues & Lung
slung1 <- sel_lung1[,1] #General Tissue
slung2 <- sel_lung2[,1]
f_GDS1096 <- as_tibble(data_table_GDS1096) %>% select(c("ID_REF",slung1))
f_GDS1096 <- f_GDS1096 %>% mutate(across(starts_with("GSM"), log, 2))
f_GDS1096 <- f_GDS1096 %>% mutate(across(starts_with("GSM"), Re))


f_GDS3257 <- as_tibble(data_table_GDS3257) %>% select(c("ID_REF",slung2))
f_GDS3257 <- f_GDS3257 %>% mutate(across(starts_with("GSM"), log, 2))
f_GDS3257 <- f_GDS3257 %>% mutate(across(starts_with("GSM"), Re))


#Renal Tissues
srenal <- sel_renal1[,1]
srenal2 <- sel_renal2[,1]
f_GDS505 <- as_tibble(data_table_GDS505) %>% select(c("ID_REF",srenal))
f_GDS505 <- f_GDS505 %>% mutate(across(starts_with("GSM"), log, 2))
f_GDS505 <- f_GDS505 %>% mutate(across(starts_with("GSM"), Re))


f_GDS507 <- as_tibble(data_table_GDS507) %>% select(c("ID_REF",srenal2))
f_GDS507 <- f_GDS507 %>% mutate(across(starts_with("GSM"), log, 2))
f_GDS507 <- f_GDS507 %>% mutate(across(starts_with("GSM"), Re))


#Rename columns for both 505 and 507 and do an rbind
rename_renal <- names(f_GDS505)
names(f_GDS507) <-rename_renal
all_normal_renal <- rbind(f_GDS505, f_GDS507)


#Myeloid Tissue
smyeloid <- sel_myeloid[,1]
f_GDS3057 <- as_tibble(data_table_GDS3057) %>% select(c("ID_REF",smyeloid))
f_GDS3057 <- f_GDS3057 %>% mutate(across(starts_with("GSM"), log, 2))
f_GDS3057 <- f_GDS3057 %>% mutate(across(starts_with("GSM"), Re))


#Changing the column names from ID_REF to probe_id
#Only need probe_id and NOT IDENTIFIER^
#Might have to use this version:  rename(c("ID_REF" = "probe_id")) instead of rename(c("probe_id" = "ID_REF"))
#Package updates change code arrangement procedure or it does not run correctly.

f_GDS596 <- f_GDS596 %>% rename(c("ID_REF" = "probe_id"))

f_GDS3716 <- f_GDS3716 %>%  rename(c("ID_REF" = "probe_id"))

f_GDS1965 <- f_GDS1965 %>%  rename(c("ID_REF" = "probe_id"))

f_GDS1375 <- f_GDS1375 %>%  rename(c("ID_REF" = "probe_id"))

f_GDS1096 <- f_GDS1096 %>%  rename(c("ID_REF" = "probe_id"))

f_GDS3257 <- f_GDS3257 %>%  rename(c("ID_REF" = "probe_id"))

f_GDS505 <- f_GDS505 %>%  rename(c("ID_REF" = "probe_id"))

f_GDS507 <- f_GDS507 %>%  rename(c("ID_REF" = "probe_id"))

all_normal_renal <- all_normal_renal %>%  rename(c("ID_REF" = "probe_id"))

f_GDS3057 <- f_GDS3057 %>%  rename(c("ID_REF" = "probe_id"))


#separate(col = #qname, into = c("chip","probe_id"), sep = ":") #used to separate punctuation in columns and rows
#Expression data as your first thing, want all data out of expression table. Do not want to do a "keep all."

#Reading the affy_U133AB into R and altering it accordingly (removing unnecessary punctuation).
#Affy_U133AB was downloaded from the UCSC Genome Table Browser as a table from the browser to obtain the entire genome track at once. 
#Necessary for the coordinates of the probes from the GRCh37/hg37 track. (Which probes go to which TADs).
#Essentially for mapping probe IDs to genomic coordinates

library(readr)
affy_U133AB_probe_coord <- read_delim("affy_U133AB_probe_coord", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE)
View(affy_U133AB_probe_coord)

affy_U133AB_probe_coord <- affy_U133AB_probe_coord %>% clean_names() 
affy_U133AB_probe_coord <- affy_U133AB_probe_coord %>% separate(col = number_q_name, into = c("chip","probe_id"), sep = ":") #separating singular columns into multiple and removing unnecessary punctuation
affy_U133AB_probe_coord <- affy_U133AB_probe_coord %>% mutate(probe_id = str_remove(probe_id,";")) %>% glimpse() 
View(affy_U133AB_probe_coord)

write.csv(affy_U133AB_probe_coord, file = "affy_U133AB_probe_coord_alt", row.names = FALSE)

#Processes of merging each altered GDS accession number (gene expressions from previous code) step by step instead of all at once to avoid numerous errors and issues (combining/merging on the probe_id column)

combined <- merge(affy_U133AB_probe_coord,f_GDS596, by.x="probe_id", by.y = "probe_id", all=T) 
combined_1 <- merge(combined,f_GDS3716, by = "probe_id")
combined_2 <- merge(combined_1,f_GDS1965, by = "probe_id")
combined_3 <- merge(combined_2, f_GDS1375, by = "probe_id")
combined_4 <- merge(combined_3, f_GDS1096, by = "probe_id")
combined_5 <- merge(combined_4, f_GDS3257, by = "probe_id")
combined_6 <- merge(combined_5, f_GDS3057, by = "probe_id")
combined_f <- merge(combined_6, all_normal_renal, by = "probe_id")

combined_f <- pivot_longer(as_tibble(combined_f), cols = starts_with("GSM"), names_to = "sample", values_to = "expression")

#Normal Expression Data
#compute Variance across Affy probes for each TAD, in NCI-60 the same thing is repeated and the variance is compared within a TAD between the two data sets.
#Across or within a tissue: here's the variability of the genes within a TAD, for each cancer cell line, here's the variability of the cancer genes within the same TAD, does the variability of the cancer cell line differ from the genes within the normal tissue? "If a structural rearrangement has altered the TAD within the cell line we will pick that up within the gene expression.
#Variability within the TAD."

#Information for sel_phenotype_data1
#adding the tissue column to the appropriate matching "cancerous version" into a uniformly named column for easier table merging

data_phenotype <- sel_phenotype_data1 %>%
  mutate(tissue = case_when(
    str_detect(tissue, "bronchial") ~ "skin",
    str_detect(tissue, "BM-CD") ~ "bone marrow",
    str_detect(tissue, "kidney") ~ "renal",
    str_detect(tissue, "bone") ~ "bone marrow",
    str_detect(tissue, "lung") ~ "lung",
    TRUE ~ NA_character_)) %>% select(c("sample","tissue"))

#Information for sel_breast
#May have to use rename("specimen" = "tissue") instead of rename("tissue" = "specimen")

data_breast <- sel_breast %>%
  rename("specimen" = "tissue") %>% mutate(tissue = dplyr::recode(tissue,"reduction mammoplasty" = "breast")) %>% select(c("sample", "tissue"))

#Information for sel_melanocyte & sel_melanocyte2
#May have to use rename("cell.type" = "tissue") instead of rename("tissue" = "cell.type")

data_melanocyte <- sel_melanocyte %>% 
  rename("cell.type" = "tissue") %>% mutate(tissue = dplyr::recode(tissue,
                                                                   "melanocyte" = "skin")) %>% select(c("sample", "tissue"))

data_melanocyte2 <- sel_melanocyte2 %>% mutate(sel_melanocyte2, tissue = "skin") %>% select(c("sample", "tissue"))


#Information for sel_lung1 & sel_lung2

data_lung1 <- sel_lung1 %>%
  mutate(tissue = dplyr::recode(tissue,
                                "bronchial epithelial cells" = "skin",
                                "BM-CD33+myeloid" = "bone marrow",
                                "kidney" = "renal")) %>% select(c("sample","tissue")) %>% filter(tissue == "bone marrow" | tissue == "breast" | tissue == "renal" | tissue == "lung")

data_lung2 <- sel_lung2 %>%
  mutate(tissue = dplyr::recode(tissue,
                                "normal" = "lung")) %>% select(c("sample","tissue"))

#Information for sel_renal1 & sel_renal2

data_renal1 <- sel_renal1 %>% mutate(sel_renal1, tissue = "renal") %>% select(c("sample", "tissue"))

data_renal2 <- sel_renal2 %>% mutate(sel_renal2, tissue = "renal") %>% select(c("sample", "tissue"))

#Information for sel_myeloid

#May have to use rename("cell.type" = "tissue") instead of rename("tissue" = "cell.type")

data_myeloid <- sel_myeloid %>%
  rename("cell.type" = "tissue") %>% select(c("sample", "tissue"))


#Combining the data_tissue tables into a single table after renaming columns and rows appropriately

Tissue_Sample_table <- rbind(data_phenotype, data_breast, data_melanocyte, data_melanocyte2, data_lung1,  data_lung2, data_renal1, data_renal2, data_myeloid)



#TADS should expand more of the genome and should be put first. 
#Build a table that has the following columns: TAD (coordinates) - do not need t_start or t_end or t_name 
#Output: TAD, Tissue, and Variance. 
#Compute variance for each normal tissue combined and for each cancer cell line individually
#Do normal(s) in one step and then do the NCI-60 in the next step.

#For each unique TAD ID
#Merge combined_f with Tissue_Sample_Table on "sample" - merge it so combined_f is first
#In for loop the next step is the code to group sub1 by tissue type and compute variance and assign that to a result and seeing the output.
#tidy the combined expression data so that columns are affy_probe_id, sample, expression value, do for both normal and NCI-60. 
#For normal, add tissue value. Based on Tissue_sample_table (merge on "sample" and make it a tibble again).



#The combined tables for Probe_id, Sample, Expression, and Tissue
combined_f <- combined_f %>%
  select("sample", "probe_id", "expression")

combined_tissue <- merge(combined_f, Tissue_Sample_table, by = "sample")
#combined_tissue <- combined_tissue %>%
#select("sample", "probe_id", "expression") #May have to use this "#" version based on package
combined_tissue <- as_tibble(combined_tissue)
combined_tissue <- unique(combined_tissue)

write.csv(combined_tissue, file = "combined_tissue", row.names = FALSE) #save the new combined table as an external file for modifications in another program or reloading later


#---------------------------------------------------------------------------------------------------------------------------------

#Loading Cancerous TAD Affy Expression data NCI60

library(readr) ###comes from CellMiner at discover.nci.nih.gov

RNA_Affy_HG_U133_A_B_MAS5 <- read_csv("nci60_RNA__Affy_HG_U133(A_B)_MAS5/output/RNA__Affy_HG_U133(A_B)_MAS5.csv")

RNA_Affy_HG_U133_A_B_MAS5 <- RNA_Affy_HG_U133_A_B_MAS5 %>% 
  rename("Identifier" = "probe_id") #rename("probe_id" = "Identifier")

RNA_Affy_HG_U133_A_B_MAS5 <- RNA_Affy_HG_U133_A_B_MAS5 %>% mutate_if(is.double,as.character)

RNA_Affy_HG_U133_A_B_MAS5 <- RNA_Affy_HG_U133_A_B_MAS5 %>% pivot_longer(cols = "BR_MCF7":"RE_UO-31", names_to= c("sample"), values_to = "expression") 

RNA_Affy_HG_U133_A_B_MAS5 <- RNA_Affy_HG_U133_A_B_MAS5 %>%
  filter(str_detect(sample, "BR_") | str_detect(sample, "LC_") | str_detect(sample, "RE_") | str_detect(sample, "ME_") | str_detect(sample, "LE_"))


RNA_Affy_HG_U133_A_B_MAS5 <-RNA_Affy_HG_U133_A_B_MAS5 %>% 
  mutate(tissue = case_when(
    str_detect(sample, "BR_") ~ "breast",
    str_detect(sample,"LC_") ~ "lung",
    str_detect(sample,"RE_") ~ "renal",
    str_detect(sample,"ME_") ~ "skin",
    str_detect(sample,"LE_") ~ "bone marrow", 
    TRUE ~ NA_character_))

RNA_Affy_HG_U133_A_B_MAS5 <- as_tibble(RNA_Affy_HG_U133_A_B_MAS5)

View(RNA_Affy_HG_U133_A_B_MAS5)

write.csv(RNA_Affy_HG_U133_A_B_MAS5, file = "RNA_Affy_HG_U133AB_Mas5_Alt", row.names = FALSE)


#-----------------------------------------------------------------------------------------------------------------------------

#Akdemir TADs Tables
#Loading in the Akdemir TAD data from Akdemir et al. (2020) after modifying their Supplementary Table 1 a small bit in Excel (renaming columns for easier import)

library(readr)
Akdemir_TADs_alt <- read_csv("Akdemir TADs_alt.csv")
View(Akdemir_TADs_alt)

Akdemir_TADs_alt <- Akdemir_TADs_alt %>% 
  rename("chromosome" = "t_name", "start" = "t_start", "end" = "t_end")
#rename("t_name" = "chromosome", "t_start" = "start", "t_end" = "end") #May have to use "#" out code dpending on loaded package, code is rearranged depending on package and is often fickle

Akdemir_TADs_alt <- Akdemir_TADs_alt %>%
  mutate(TAD = paste0(t_name,":", t_start, "-", t_end))
Akdemir_TADs_alt <- Akdemir_TADs_alt %>%
  select("t_name", "t_start", "t_end", "TAD")

Akdemir_TADs_alt <- as.data.table(Akdemir_TADs_alt)

setkey(affy_U133AB_probe_coord, t_name, t_start, t_end)
setkey(Akdemir_TADs_alt, t_name, t_start, t_end)

affy_U133AB_probe_coord <- as.data.table(affy_U133AB_probe_coord)
Akdemir_TADs_alt <- as.data.table(Akdemir_TADs_alt)

write.csv(Akdemir_TADs_alt, file = "Akdemir_TADS_Altered", row.names = FALSE)




#Creating an foverlaps table using the Akdemir TADs from Akdemir (2020) and the Affy_U133AB_Probe_Coord (UCSC Genome Table Browser) to create a "look-up" table that contains the TAD name and the start and stop genomic coordinates that is eventually combined with the cancerous and normal expression data below used to compute variance by TAD.

aux3  <- foverlaps( affy_U133AB_probe_coord, Akdemir_TADs_alt, 
                    by.x = c("t_name", "t_start",   "t_end"), 
                    by.y = c("t_name", "t_start",   "t_end"), 
                    type = "within", nomatch = NA)





#Computing Variance by TAD for Akdemir w/ normal expression data

variance_Akdemir <- as_tibble(data.table(tissue = NA_character_, variance = NA, tad= NA_character_, n_affy_probes= NA)) 
tad_names_Akdemir <- unique(Akdemir_TADs_alt$TAD)

for (m in 1:length(tad_names_Akdemir)){
  affy_probes5 <- aux3[TAD == tad_names_Akdemir[m], probe_id]
  sub5 <- combined_tissue %>% filter(probe_id %in% affy_probes5) #Replaced combined_f with combined_tissue
  result <- sub5 %>% group_by(tissue) %>% summarize(variance = var(expression, na.rm = TRUE))
  result$tad <- tad_names_Akdemir[m]
  result$n_affy_probes = length(affy_probes5)
  variance_Akdemir <- rbind(variance_Akdemir, result)
}


variance_Akdemir <- variance_Akdemir %>% filter(!is.na(tissue))

write.csv(variance_Akdemir, file = "Variance_Akdemir", row.names = FALSE)



#Computing Variance by TAD for Akdemir w/ NCI60 expression data

variance_Akdemir_affy <- as_tibble(data.table(sample = NA_character_, variance = NA, tad= NA_character_, n_affy_probes= NA))
tad_names_Akdemir_affy <- unique(Akdemir_TADs_alt$TAD)

for (n in 1:length(tad_names_Akdemir_affy)){
  affy_probes6 <- aux3[TAD == tad_names_Akdemir_affy[n], probe_id]
  sub6 <- RNA_Affy_HG_U133_A_B_MAS5 %>% filter(probe_id %in% affy_probes6) #grab expression values by affy_probe associated with TAD[i]
  result <- sub6 %>% group_by(sample) %>% summarize(variance = var(expression, na.rm = TRUE))
  result$tad <- tad_names_Akdemir_affy[n]
  result$n_affy_probes = length(affy_probes6)
  variance_Akdemir_affy <- rbind(variance_Akdemir_affy, result)
}

variance_Akdemir_affy <- variance_Akdemir_affy %>% filter(!is.na(sample))

variance_Akdemir_affy <- variance_Akdemir_affy %>%
  mutate(tissue = case_when(
    str_detect(sample, "BR_") ~ "breast",
    str_detect(sample,"LC_") ~ "lung",
    str_detect(sample,"RE_") ~ "renal",
    str_detect(sample,"ME_") ~ "skin",
    str_detect(sample,"LE_") ~ "bone marrow", 
    TRUE ~ NA_character_))

write.csv(variance_Akdemir_affy, file = "Variance_Akdemir_Affy", row.names = FALSE)





#------------------------------------------------------------------------------------------------------------

#Simple Visual Plots created to display the the cancerous and normal expression values separately

#ggplot original normal & NCI60 expression tables
#(2 tables)

Cancerous_expression <- RNA_Affy_HG_U133_A_B_MAS5
Normal_expression <- combined_tissue

Cancerous_expression <- Cancerous_expression %>% mutate(expression = as.numeric(expression))

Normal_expression <- Normal_expression %>% mutate(expression = as.numeric(expression()))

cancer_table <- ggplot(Cancerous_expression, aes(x=expression, color=tissue)) +
  geom_freqpoly()+
  scale_color_brewer()+
  theme_bw()
#facet_grid(~tissue, scales = 'free_x')

cancer_table + ggtitle("Plot of RNA_Affy Expression (Cancer)") +
  xlab("Expression") + ylab("Count")

ggplot(Cancerous_expression, aes(x=expression, color=tissue)) +
  geom_density()+
  scale_color_brewer()+
  theme_bw()

ggplot(Cancerous_expression, aes(x=expression, fill=tissue)) +
  geom_histogram(bins=100)+
  geom_freqpoly()+
  scale_fill_manual(values=c('gray', 'red', 'yellow', 'darkgreen', 'cyan'))+
  theme_bw()

#------------------------------------------------------------------------------------------------------------------------------------

normal_table <- ggplot(Normal_expression, aes(x=expression, color=tissue)) +
  geom_freqpoly()+
  scale_color_brewer()+
  theme_bw()
#facet_grid(~tissue, scales = 'free_x')

normal_table + ggtitle("Plot of combined_tissue Expression (Normal)") +
  xlab("Expression") + ylab("Count")


ggplot(Normal_expression, aes(x= expression, fill=tissue)) +
  geom_histogram(bins=100)+
  geom_freqpoly()+
  scale_fill_manual(values=c('gray', 'red', 'yellow', 'darkgreen', 'cyan'))+
  labs(x = "Log2 Expression") +
  theme_bw()


#------------------------------------------------------------------------------------------------------------------------------------

#Creating the table that computes the variance between the normal and cancerous expression for each TAD utilizing a Levene's equality of variances

#Testing the Levene's equality of variance by subsetting through both the cancerous and normal expression data with Akdemir TADs via loop

library(stats)
install.packages("pairwiseComparisons")
library(pairwiseComparisons)

#N_affy_probes is an indirect indicator of a gene. These are parts of a transcript, not whole genes.
Akdemir_table_L_all <- as_tibble(data.table(TAD = NA_character_, cell_line = NA_character_, tissue_names = NA_character_,
                                        n_affy_probes = NA, normal_sd = NA, cancer_sd = NA, f_value = NA, p_value = NA))

cell_line_list <- unique(RNA_Affy_HG_U133_A_B_MAS5$sample) #Obtaining the list of cell lines
unique_tad <- unique(variance_Akdemir_affy$tad) #Obtaining the unique TADS from the Affy table

for (n in 1:length(cell_line_list)){#length(cell_line_list)){ #Stepping through cell line by cell line
  print(n)
  subset <- RNA_Affy_HG_U133_A_B_MAS5 %>% filter(sample == cell_line_list[n]) #NCI-60 Affymetrix cancerous expression 
  sel_tissue <- unique(subset$tissue) #Tissues within the created subset
  for (m in 1:length(unique_tad)){ #Stepping into the TADs
    print(m) #What number am I on?
    affy_probes7 <- aux3[TAD == unique_tad[m], probe_id] #aux3 is the intersection between the Akdemir TADs and the Affymetrix Probe IDs
    #Picking up the TADs and the Probe_ID from the intersection
    if(length(affy_probes7) >=3) { #Checking to make sure there are enough Probe_IDs to compute variance by setting a value greater than                                      3
      sub_cell_line <- subset %>% filter(probe_id %in% affy_probes7) #Pulling cancerous expression data from the Probe_ID from the TAD
      sub_normal <- combined_tissue %>% filter(tissue == sel_tissue & probe_id %in% affy_probes7) #Pulling normal expression data from                                                                                                     the Probe_ID from the TAD
      cell_data_frame <- merge(sub_cell_line,sub_normal, by=c("probe_id"), all = TRUE) #Combining both Probe_ID info from cancerous and                                            normal expression into a single table
      cell_data_frame$expression.x <- as.numeric( cell_data_frame$expression.x) #Formatting the expression data as numeric (x=cancerous)
      cell_data_frame$expression.y <- as.numeric( cell_data_frame$expression.y) #Formatting the expression data as numeric (y=normal)
      normal_cell <- cell_data_frame %>% filter(!is.infinite(expression.y)) %>% select(expression.y) %>% filter(!is.na(expression.y))
      if (dim(normal_cell)[1] == 0) {next} #Filtering to make sure we do not possess infinite or NA's within normal & make sure if we do                                   not have something we pull back and exit the dim (if dim 1 = 0 we go to the next iteration of the loop)
      #Making sure we obtain results from expression data for the Probe_IDs
      cancer_cell <- cell_data_frame %>% filter(!is.infinite(expression.x)) %>% select(expression.x) %>% filter(!is.na(expression.x))
      if (dim(cancer_cell)[1] == 0) {next} #Filtering to make sure we do not possess infinite or NA's within cancerous & make sure if we                               do not have something we pull back and exit the dim (if dim 1 = 0 we go to the next iteration of the loop)
      #Making sure we obtain results from expression data for the Probe_IDs
      names(cancer_cell)[1] <-"expression.y" #Changed name on the cancer_cell so it is the same thing as the normal so we can paste into                                               the same array
      normal_cell$type <- "normal" #Added type to tell normal and cancerous apart (which values belong to which)
      cancer_cell$type <- "cancer"
      my.data <- as.data.frame(rbind(normal_cell, cancer_cell)) #Made a data frame that combines everything together
      my.data$type <- as.factor(my.data$type) #Changed the my.data to a factor for the Levene's analysis
      levene <- leveneTest(expression.y~type, my.data, center = mean) #Conducted Levene's test on normal and cancerous expression to                                                                            calculate different in variance (null: Variances are equal;                                                                               alternative: variances are not equal) greater robustness overall
      f_value <- levene$'F value'[1] #Grabbing the F-value from the Levene's result
      p_value <- levene$`Pr(>F)`[1] #Grabbing the p-value from the Levene's result
      normal_sd <- sd(normal_cell$expression.y)#Grabbing the standard deviation for the normal_cell variable
      cancer_sd <- sd(cancer_cell$expression.y)#Grabbing the standard deviation for the cancer_cell variable
      row_to_add <- data.table(TAD = unique_tad[m], cell_line = cell_line_list[n], tissue_names = sel_tissue, 
                               n_affy_probes = length(affy_probes7), normal_sd = normal_sd, cancer_sd = cancer_sd, f_value = f_value, p_value = p_value)
      Akdemir_table_L_all <- rbind(Akdemir_table_L_all,row_to_add)} else
      {next} #Combined everything into a single table called "Akdemir_table_L" or Akdemir Table Levene's
  }
}


#Creating various tables that isolate the Levene's Akdemir table "Akdemir_table_L" above depending on what we want to look for.


distinct(Akdemir_table_L_all) #Creating distinct (no duplicate) values

normal_sd_fivenumer_summary <- fivenum(Akdemir_table_L_all$normal_sd, na.rm = TRUE) #Outputting summary: min, quartiles, max for normal_sd (standard deviation) column
cancer_sd_fivenumer_summary <- fivenum(Akdemir_table_L_all$cancer_sd, na.rm = TRUE) #Outputting summary: min, quartiles, max for cancer_sd (standard deviation) column

#Creating the table used throughout the rest of the code as a "reference table" containing TADs, cell lines, tissues, number of affymetrix probes, and p-values.
Akdemir_table_L <- Akdemir_table_L_all %>% select("TAD", "cell_line", "tissue_names", "n_affy_probes", "f_value", "p_value")
Akdemir_table_L <- drop_na(Akdemir_table_L) #Removing NA variables
Akdemir_table_L <- Akdemir_table_L %>%
  mutate(adjusted_p = p.adjust(Akdemir_table_L$p_value, method = "BH")) #Creating an "adjusted_p" column from the p-value using the Benjamini-Hochenberg correction

summary(Akdemir_table_L) #Outputting information such as length, quartiles, and classes for each column
Akdemir_table_L_summary <- summary(Akdemir_table_L) #Taking the output and putting it into a variable
write.csv(Akdemir_table_L_summary, file = "Akdemir_table_L_summary", row.names = FALSE) #Exporting the output from the variable into a CSV


Akdemir_table_L_filtered <- Akdemir_table_L %>%
  filter(adjusted_p < 0.01) #Making a new table and filtering for p-values less than 0.01


#Creating a new table that properly separates columns and changes variables appropriately (characters turned to numeric, etc.)
Akdemir_table_L_t1  <- Akdemir_table_L  %>% separate(col = TAD, into = c("chr","range"), sep = ":", remove = FALSE) 
Akdemir_table_L_t1  <- Akdemir_table_L_t1  %>% separate(col = "range", into = c("t_start","t_stop"), sep = "-")
Akdemir_table_L_t1$t_start <- as.numeric(Akdemir_table_L_t1$t_start)
Akdemir_table_L_t1$t_stop <- as.numeric(Akdemir_table_L_t1$t_stop)

#Creating another table from the filtered version (p-value less than 0.01) and spreading a single column out.
Akdemir_table_L_filtered_t2 <- Akdemir_table_L_filtered  %>% separate(col = TAD, into = c("chr","t_start"), sep = ":", remove = FALSE)
Akdemir_table_L_filtered_t2  <- Akdemir_table_L_filtered_t2 %>% separate(col = t_start, into = c("t_start","t_stop"), sep = "-", remove = FALSE)

write.csv(Akdemir_table_L, file = "Akdemir_table_L", row.names = FALSE)


#----------------------------------------------------------------------------------------------------------------

#Filtering further from the table where adjusted p-values are less than 0.01 for the specific ENCODE cell line TADs

Akdemir_ENCODE_TAD <- Akdemir_table_L_filtered_t2 %>% filter(cell_line == "ME_SK-MEL-5" | cell_line == "LC_NCI-H460" | cell_line == "RE_ACHN" | cell_line == "LC_A549/ATCC" | cell_line == "BR_T-47D")


write.csv(Akdemir_ENCODE_TAD, file = "Akdemir_ENCODE_Filtered_TADs", row.names = FALSE)


#Tidyverse piping

Akdemir_Encode_TAD_count <- Akdemir_ENCODE_TAD %>% group_by(cell_line) %>% summarise(count = n())

write.csv(Akdemir_ENCODE_TAD_count, file = "Akdemir_ENCODE_Filtered_TADs_by_cell_line", row.names = FALSE)


#Akdemir_table_L_filtered_ENCODE <- Akdemir_table_L_filtered %>% filter(cell_line == "ME_SK-MEL-5" | cell_line == "LC_NCI-H460" | cell_line == "RE_ACHN" | cell_line == "LC_A549/ATCC" | cell_line == "BR_T-47D") %>% group_by(cell_line) %>% summarise(count=n())


#From the code towards the bottom of the project (see this later) filtering for chromosomes 10 & 19 from the Akdemir_ENCODE_TAD table

Akdemir_ENCODE_TAD_chr10_chr19 <- Akdemir_ENCODE_TAD %>% filter(chr == "chr10" | chr == "chr19")

write.csv(Akdemir_ENCODE_TAD_chr10_chr19, file = "Akdemir_ENCODE_Filtered_TADs_chr10_chr19", row.names = FALSE)


#---------------------------------------------------------------------------------------------------------------------------------------


#Removing excess memory consumers from the Environment and reducing visual anxiety


rm(affy_U133AB_probe_coord)
rm(all_normal_renal)
rm(aux)
rm(aux2)
rm(aux3)
rm(breastphenotype)
rm(cancer_cell)
rm(cell_data_frame)
rm(combined)
rm(combined_1)
rm(combined_2)
rm(combined_3)
rm(combined_4)
rm(combined_5)
rm(combined_6)
rm(combined_f)
rm(combined_tissue)
rm(data_breast)
rm(data_lung1)
rm(data_lung2)
rm(data_melanocyte)
rm(data_melanocyte2)
rm(data_myeloid)
rm(data_phenotype)
rm(data_renal1)
rm(data_renal2)
rm(data_table_GDS1096)
rm(data_table_GDS1375)
rm(data_table_GDS1965)
rm(data_table_GDS3057)
rm(data_table_GDS3257)
rm(data_table_GDS3716)
rm(data_table_GDS505)
rm(data_table_GDS507)
rm(data_table_GDS596)
rm(eset)
rm(eset2)
rm(eset3)
rm(eset4)
rm(eset5)
rm(eset6)
rm(eset7)
rm(eset8)
rm(eset9)
rm(f_GDS1096)
rm(f_GDS1375)
rm(f_GDS1965)
rm(f_GDS3057)
rm(f_GDS3257)
rm(f_GDS3716)
rm(f_GDS505)
rm(f_GDS507)
rm(f_GDS596)
rm(GDS1096)
rm(GDS1375)
rm(GDS1965)
rm(GDS3057)
rm(GDS3257)
rm(GDS3716)
rm(GDS505)
rm(GDS507)
rm(GDS596)

rm(lungphenotype)
rm(lungphenotype2)
rm(melanocytephenotype)
rm(melanocytephenotype2)
rm(myeloidphenotype)
rm(renalphenotype)
rm(renalphenotype2)
rm(phenotype_data)
rm(sel_breast)
rm(sel_lung1)
rm(sel_lung2)
rm(sel_melanocyte)
rm(sel_melanocyte2)
rm(sel_myeloid)
rm(sel_phenotype_data)
rm(sel_phenotype_data1)
rm(sel_renal1)
rm(sel_renal2)

rm(levene)
rm(my.data)
rm(normal_cell)
rm(result)
rm(RNA_Affy_HG_U133_A_B_MAS5)
rm(row_to_add)
rm(sub_cell_line)
rm(sub_cell_line_hESC)
rm(sub_cell_line_IMR)
rm(sub_normal)
rm(sub_normal_hESC)
rm(sub_normal_IMR)
rm(sub1)
rm(sub2)
rm(sub3)
rm(sub5)
rm(sub6)
rm(subset)
rm(Tissue_Sample_table)

rm(Akdemir_TADs_alt)
rm(IMR90_domains_hg19_bed)
rm(hESC_domains_hg19_bed)
rm(variance_Akdemir)
rm(variance_Akdemir_affy)
rm(variance_hESC)
rm(variance_hESCAffy)
rm(variance_IMR90)
rm(variance_IMRAffy)

rm(affy_probes)
rm(affy_probes1)
rm(affy_probes5)
rm(affy_probes6)
rm(affy_probes7)
rm(affy_probes8)
rm(affy_probes9)
rm(cell_line_list)
rm(cell_line_list_aCGH)
rm(cell_line_list_hESC)
rm(cell_line_list_IMR)
rm(f_value)
rm(i)
rm(j)
rm(k)
rm(l)
rm(m)
rm(n)
rm(p_value)
rm(rename_renal)
rm(sbreast)
rm(sel_tissue)
rm(slung1)
rm(slung2)
rm(smelanocyte)
rm(smelanocyte2)
rm(smyeloid)
rm(sphenotype)
rm(srenal)
rm(srenal2)
rm(tad_names_Akdemir)
rm(tad_names_Akdemir_affy)
rm(tad_names_hESC)
rm(tad_names_hESCAffy)
rm(tad_names_IMR90)
rm(tad_names_IMRAffy)
rm(unique_cell_line)
rm(unique_tad)
rm(unique_tad_Akdemir_L)
rm(unique_tad_hESC)
rm(unique_tad_IMR)


rm(Akdemir_table_L_filtered)
rm(Akdemir_table_L_filtered_t2)
rm(hESC_table_L_filtered)
rm(hESC_table_L_filtered_t2)
rm(IMR90_table_L_filtered)
rm(IMR90_table_L_filtered_t2)

#--------------------------------------------------------------------------------------------------------------------------



#Identify breakpoints in the genome from aCGH and map them to genomic coordinates and intersect - figure out which TADs have breakpoints and which do not. Count / classify TADs by "type of non-difference in variability vs those that have difference."
#Fisher's Test - rows 

#foverlaps with data table object
#Parse chromosome start / stop from unfiltered tables
#Use Combined file first due to resolution
#Breakpoints cell line by cell line basis (overlap based on cell line)


#Set all cell line columns as numeric

#Place this DNA_Combined_aCGH_log2.txt.txt file in your Working Directory.
#DNA_Combined_aCGH_log2.txt.txt contains DNA Copy number data from CellMiner (NCI-60) and is used to find correlation or causation between altered variance and breakpoints between TADs.

#Each column had to be manually entered and altered to be classified as "character" as Rstudio (even with readr) would not load the columns correctly. There may be human error to account for here.

library(readr)
DNA_Combined_aCGH_log2_txt <- read_delim("C:\\Users\\Vainc\\Desktop\\TAD Datasets\\nci60_DNA__Combined_aCGH_log2\\output\\DNA__Combined_aCGH_log2.txt.txt", 
                                         "\t", escape_double = FALSE, col_types = cols(`Chromosome` = col_character(), `BR_MCF7` = col_double(), 
                                                                                       `BR_MDA-MB-231` = col_double(), `BR_HS_578T` = col_double(),                                                                                    `BR_BT-549` = col_double(), 
                                                                                       `BR_T-47D` = col_double(), `CNS_SF-268` = col_double(), 
                                                                                       `CNS_SF-295` = col_double(), `CNS_SF-539` = col_double(), 
                                                                                       `CNS_SNB-19` = col_double(), `CNS_SNB-75` = col_double(),
                                                                                       `CNS_U251` = col_double(),
                                                                                       `CO_COLO 205` = col_double(), `CO_HCC-2998` = col_double(), 
                                                                                       `CO_HCT-116` = col_double(), `CO_HCT-15` = col_double(), 
                                                                                       `CO_HT29` = col_double(), `CO_KM12` = col_double(), `CO_SW-620` = col_double(),
                                                                                       `LE_CCRF-CEM` = col_double(), `LE_HL-60(TB)` = col_double(), 
                                                                                       `LE_K-562` = col_double(), `LE_MOLT-4` = col_double(), 
                                                                                       `LE_RPMI-8226` = col_double(), `LE_SR` = col_double(), `ME_LOX_IMVI` = col_double(), 
                                                                                       `ME_MALME-3M` = col_double(), `ME_M14` = col_double(), `ME_SK-MEL-2` = col_double(), `ME_SK-MEL-28` = col_double(),
                                                                                       `ME_SK-MEL-5` = col_double(), `ME_UACC-257` = col_double(), 
                                                                                       `ME_UACC-62` = col_double(), `ME_MDA-MB-435` = col_double(),                                                                                    `ME_MDA-N` = col_double(), `LC_A549/ATCC` = col_double(), 
                                                                                       `LC_EKVX` = col_double(), `LC_HOP-62` = col_double(), 
                                                                                       `LC_HOP-92` = col_double(), `LC_NCI-H226` = col_double(), 
                                                                                       `LC_NCI-H23` = col_double(), `LC_NCI-H322M` = col_double(), `LC_NCI-H460` = col_double(),
                                                                                       `LC_NCI-H522` = col_double(), `OV_IGROV1` = col_double(),
                                                                                       `OV_OVCAR-3` = col_double(), `OV_OVCAR-4` = col_double(),
                                                                                       `OV_OVCAR-5` = col_double(), `OV_OVCAR-8` = col_double(), 
                                                                                       `OV_SK-OV-3` = col_double(), `OV_NCI/ADR-RES` = col_double(),	`PR_PC-3` = col_double(),	`PR_DU-145` = col_double(), `RE_786-0` = col_double(),	`RE_A498` = col_double(),	`RE_ACHN` = col_double(), `RE_CAKI-1` = col_double(), `RE_RXF-393` = col_double(), `RE_SN12C` = col_double(), `RE_TK-10` = col_double(), `RE_UO-31` = col_double()), trim_ws = TRUE)

DNA_Combined_aCGH_log2_txt <- DNA_Combined_aCGH_log2_txt %>%
  select("Identifier", "Chromosome", "Start", "End", "BR_MCF7":"RE_UO-31") #Selecting specific columns from the uploaded DNA file

#DNA_Combined_aCGH_log2_txt <- drop_na(DNA_Combined_aCGH_log2_txt) #Depending on packages and their versions the "#" blocked out code may have to be re-enabled since different versions work with different arrangements of code.

#DNA_Combined_aCGH_log2_txt <- DNA_Combined_aCGH_log2_txt %>% rename("Chromosome" = "chr", "Start" = "t_start", "End" = "t_stop")
#DNA_Combined_aCGH_log2_txt <-DNA_Combined_aCGH_log2_txt %>% rename("chr" = "Chromosome", "t_start" = "Start", "t_stop" = "End")

DNA_Combined_aCGH_log2_txt <-DNA_Combined_aCGH_log2_txt %>% rename("Chromosome" = "chr", "Start" = "t_start", "End" = "t_stop") #Renaming columns for uniformity purposes
DNA_Combined_aCGH_log2_txt <- DNA_Combined_aCGH_log2_txt %>% mutate(chr = paste0("chr", chr)) #Placing "chr" in the chromosome column

#DNA_Combined_aCGH_log2_txt$`BR_BT-549` <- as.double(DNA_Combined_aCGH_log2_txt$`BR_BT-549`)

DNA_Combined_aCGH_log2_txt <- (as_tibble(DNA_Combined_aCGH_log2_txt) 
                               %>% pivot_longer(cols = "BR_MCF7":"RE_UO-31", names_to= c("sample"), values_to = "expression")) #Setting the table to a tibble and rearranging the table so the cancerous columns become the rows and are under a column called "sample" with their expression values being in the adjacent rows

DNA_Combined_aCGH_log2_txt <- DNA_Combined_aCGH_log2_txt %>% #Mutating the tissue names on a case_when basis to add a tissue names column with the extended name of the tissue for each cancer
  mutate(tissue_names = case_when(
    str_detect(sample, "BR_") ~ "breast",
    str_detect(sample,"LC_") ~ "lung",
    str_detect(sample,"RE_") ~ "renal",
    str_detect(sample,"ME_") ~ "skin",
    str_detect(sample,"LE_") ~ "bone marrow",
    str_detect(sample,"CNS_") ~ "central nervous system",
    str_detect(sample,"OV_") ~ "ovarian",
    str_detect(sample,"CO_") ~ "colon",
    str_detect(sample,"LC_") ~ "lung",
    str_detect(sample,"PR_") ~ "prostate",
    TRUE ~ NA_character_))



unique_cell_line <- unique(Akdemir_table_L$cell_line) #Finding unique cell line combinations (no duplicates)
DNA_Combined_aCGH_log2_txt <- DNA_Combined_aCGH_log2_txt %>% filter(sample %in% unique_cell_line) %>% glimpse()

#DNA_Combined_aCGH_log2_txt <- DNA_Combined_aCGH_log2_txt$chr %>% mutate(chr = paste0("chr", chr)) #Depending on package or version this code may have to be used instead of the currently enabled code for renaming or arrangement purposes
#DNA_Combined_aCGH_log2_txt$chr <- paste0("chr", DNA_Combined_aCGH_log2_txt$chr)

DNA_Combined_aCGH_log2_txt <- DNA_Combined_aCGH_log2_txt %>% unite("chr_position", t_start:t_stop, sep = "-", remove = FALSE) #Combining the chromosome position and genomic start and stop coordinates into a "TAD" column
DNA_Combined_aCGH_log2_txt <- DNA_Combined_aCGH_log2_txt %>% unite("chr_position", chr:chr_position, sep = ":", remove = FALSE) #Combining the chromosome position and genomic start and stop coordinates into a "TAD" column (second step invovles removing the colon)
DNA_Combined_aCGH_log2_txt %>% glimpse()

write.csv(DNA_Combined_aCGH_log2_txt, file = "DNA_Combined_aCGH_ALT", row.names = FALSE) #Exporting the altered DNA file to a cvs for any future alterations or imports


#Beginning steps for finding altered variance and breakpoints on the Akdemir table that includes the TADs and expression variance for the normal and cancerous expression data

Akdemir_table_L_t1 <- as.data.table(Akdemir_table_L_t1) #Setting as a data table
DNA_Combined_aCGH_log2_txt <- as.data.table(DNA_Combined_aCGH_log2_txt) #Setting as a data table

Akdemir_table_L_t1 <- drop_na(Akdemir_table_L_t1) #Removing NA's
DNA_Combined_aCGH_log2_txt <- drop_na(DNA_Combined_aCGH_log2_txt) #Removing NA's


#Adding a case_when column to the table that labels if the genomic coordinates follow a certain pattern

DNA_Combined_aCGH_log2_txt <- DNA_Combined_aCGH_log2_txt %>% mutate(test = case_when(t_start < t_stop ~ "less",
                                                                                     t_start > t_stop ~ "greater",
                                                                                     t_start == t_stop ~ "equal",
                                                                                     TRUE ~ NA_character_,))

#Filtering out the SNP from the DNA Table (those with same start and stop numbers) we do not want the SNP or start genomic coordinates that are greater than the stop to influence our data

DNA_Combined_aCGH_log2_txt <- DNA_Combined_aCGH_log2_txt %>% filter(test != "greater")
DNA_Combined_aCGH_log2_txt %>% glimpse()



#Creating the intersecting TADS w/ AffyProbes table

#SetKeys
setkey(Akdemir_table_L_t1, chr, t_start, t_stop) 	#set the data.table keys so that we can intersect genomic ranges
setkey(DNA_Combined_aCGH_log2_txt, chr, t_start, t_stop)

#Altered_variance cannot be summed unless converted from "character" class.
#Used 0.01 as it is more "conservative" than the 0.05.

#Establishing a column of "altered variance" by assigning binary outputs from adjusted p-values, p-values less than 0.01 have a breakpoint and are given "1", no altered variance is given a "0" when the adjusted p-value is greater than or equal to the 0.01 value. We gave a binary output for presence and absence of altered variance for summing purposes in future code.

Akdemir_table_L_t1 <- Akdemir_table_L_t1 %>% mutate(altered_variance = case_when(
  adjusted_p < 0.01 ~ "1",
  adjusted_p >= 0.01 ~ "0",
  TRUE ~ NA_character_)) %>% glimpse()
sum(Akdemir_table_L_t1$altered_variance)

write.csv(DNA_Combined_aCGH_log2_txt, file = "aCGH_DNA_alt", row.names = FALSE) #Exporting the altered DNA copy number file
write.csv(Akdemir_table_L_t1, file = "Akdemir_table_alt", row.names = FALSE) #Exporting the altered Akdemir_table_L_t1 file


#Create a table with only TADS and no extraneous data (Identifier 'general class of things' TAD table = TAD & aCGH = Identifier, Chromosome, t_start, t_stop) label of the region that we are overlapping.

#TAD_coords <- Akdemir_table_L_t1 %>% select(TAD:t_stop) %>% distinct()
#DNA_coords <- DNA_Combined_aCGH_log2_txt %>% select(Identifier, chr:t_stop) %>% distinct()
#setkey(TAD_coords, chr, t_start, t_stop)
#setkey(DNA_coords, chr, t_start, t_stop)

#lookup table Akdemir
#Akdemir_DNA_coords <- foverlaps(TAD_coords, DNA_coords, 
 #                               by.x = c("chr", "t_start", "t_stop"), 
  #                              by.y = c("chr", "t_start", "t_stop"), 
   #                             type = "any", nomatch = NA)


#---------------------------------------------------------------------------------------------------------------

#Foroverlaps


#This overlap worked
#Creating an overlap between the Akdemir TADs and the DNA Copy Number to find altered variance and breakpoint

aux_aCGH_Akdemir <- foverlaps(Akdemir_table_L_t1, DNA_Combined_aCGH_log2_txt, 
                              by.x = c("chr", "t_start", "t_stop"), 
                              by.y = c("chr", "t_start", "t_stop"), 
                              type = "any", nomatch = NA)

write.csv(aux_aCGH_Akdemir, file.path(path,"aux_aCGH_Akdemir_alt.csv"), row.names = FALSE) #Exporting the file for safety and future imports.



#Finding breakpoints to find correlation or causation of altered variance and breakpoints present within Akdemir TADs having both normal and cancerous expression values (variance)

cell_line_list_aCGH <- unique(DNA_Combined_aCGH_log2_txt$sample) #Creating a unique cell line list "sample" (cancerous cell lines)
unique_tad_Akdemir_L <- unique(Akdemir_table_L_t1$TAD) #Creating a unique Akdemir TAD list

for (n in 1:length(cell_line_list_aCGH)){ 
  print(paste0("cell line ", n))
  cell_line_acgh <- DNA_Combined_aCGH_log2_txt %>% filter(sample == cell_line_list_aCGH[n]) #take this subset and use the logic of the article pasted to create a variable indicating something about lagging. lag.ratio := c(NA, ...  ), use a while or an if statement to "if this is yes {next TAD} if it's no keep looking. Difference column is the if statement
  
  #Goal 1: What probes have breakpoint?
  #Goal 2: What is the difference between the probes?
  
  #Do we have a breakpoint within the TAD? 
  for (m in 1:length(unique_tad_Akdemir_L)){ #For each unique TAD
    print(m) #Printing "m" for positional reference in the loop
    tad_acgh <- Akdemir_DNA_coords %>% filter(TAD == unique_tad_Akdemir_L[m]) %>% select(Identifier) #Subset to pull aCGH probes contained by a TAD
    subset_acgh <- as.data.table(merge(tad_acgh,cell_line_acgh,by.x = "Identifier",by.y = "Identifier")) #Pulling aCGH data and combining tables by "Identifier"
    if(dim(subset_acgh)[1] <= 1){ #if no overlaps do <=1 #Checking to make sure there is more than one aCGH probe present within a TAD
      {next}} #Next iteration of the loop if there is not more than one aCGH probe per TAD
    setorder(subset_acgh, t_start) #Make sure the probes are in chromosomal order
    subset_acgh[, difference := expression - shift(expression, n=1L, fill = NA, type = "lag")] #"Windowing" through probe by probe and determining the difference between probe [i] and probe [i] +1 (shift)
    for (i in 2:dim(subset_acgh)[1]){ #Started at 2 since the first variable is "NA"
      if(subset_acgh[i,difference] >= abs(0.5)) { #If the difference in expression is greater than the absolute value of 0.5 then assign a breakpoint value of "1" or "present" if not head to the next iteration of the loop and assign no breakpoint or "0"
        Akdemir_table_L_t1[TAD==unique_tad_Akdemir_L[m] & cell_line ==cell_line_list_aCGH[n], breakpoint := 1]
        {break}
      } else
      {Akdemir_table_L_t1[TAD==unique_tad_Akdemir_L[m] & cell_line ==cell_line_list_aCGH[n], breakpoint := 0]}
    }
    
    #for (k in 1:dim(aCGH_probes)[1]) {
    #add reference code here, modify aCGH_probes
    # setorder(aCGH_probes, i.t_start)
    # aCGH_probes[, prev_expression := c(NA, expression[-.N] )]
    #aCGH_probes[, diff := expression - prev_expression]
    #for (i in 2:dim(aCGH_probes)[1]){
    #if(aCGH_probes[i, diff] >= abs(0.3)) { #abs(diff) < 0.3 = relatively unchanged
    # Akdemir_table_L_t1[TAD==unique_tad_Akdemir[m] & cell_line ==cell_line_list_aCGH[n], breakpoint := 1]
    #{break}} else
    # {Akdemir_table_L_t1[TAD==unique_tad_Akdemir[m] & cell_line ==cell_line_list_aCGH[n], breakpoint := 0]}
  } 
}



Akdemir_table_L_t1 <- Akdemir_table_L_t1 %>% mutate(altered_variance = as.numeric(altered_variance), breakpoint = as.numeric(breakpoint))
#write.csv(Akdemir_table_L_t1, file = "Akdemir_table_L_t1", row.names = FALSE)
write.csv(Akdemir_table_L_t1, file.path(path, "Akdemir_table_L_t1.csv"), row.names=FALSE)

#----------------------------------------------------------------------------------------------------------------------------------------------------

#FISHER'S MATRIX TEST (w/o UPSAMPLING)

#Testing a Fisher's Matrix with available data (not up sampled or altered)
Akdemir_test <- Akdemir_table_L_t1[!is.na(breakpoint),.N, by = .(altered_variance, breakpoint)]#filter before comma
head(Akdemir_test)
akdemir_fisher_table <-as.matrix(rbind(c(378,15388),c(0,75)))
akdemir_fisher_table_adj <- as.matrix(rbind(c(398, 15408), c(20, 95)))

#Variance - column, breakpoint - row
fisher.test(akdemir_fisher_table) #Trust this one, in line with the data from the ESC display. Odds ratio of altered variance associated with breakpoint is: (closer to 0 means the association is opposite or inverse 'protective', further from 0 means it is much more likely). Look at 95% confidence interval, interval does not span 1 then we can expect the p-value to be significant. This confidence interval doesn't contain no difference in the odds. 

fisher.test(akdemir_fisher_table_adj) #Do not trust as highly as the unaltered fisher.test above, not sure if correction is viable - viable correction needs a 0.5 Fisher should be *exact* (fudge factor of 20 to every single cell in the matrix so we had enough to run the Fisher)^

Akdemir_table_L_t1 <- as.data.table(Akdemir_table_L_t1) #setting as a data table
Akdemir_filter <- Akdemir_table_L_t1 %>% filter(altered_variance == 1 & breakpoint == 1) %>% glimpse() #Filtering by isolating the data that has both altered variance and breakpoint present

Akdemir_filter <- unique(Akdemir_filter$TAD) %>% glimpse() #Visualizing the data without eating too much memory
Akdemir_filter <- Akdemir_table_L_t1[,c("TAD", "chr", "t_start", "t_stop")] %>% glimpse() #Data table syntax to only pull out the columns mentioned: "TAD", "t_start", and "t_stop"

write.csv(Akdemir_test, file.path(path, "Akdemir_test_matrix.csv"), row.names=FALSE) #Exporting the Akdemir_test file for reference

-----------------------------------------------------------------------------------------------------------------------------------------------------
  
#Simple GG-Plots to visualize tables by specific variables
  
#This below (ggplot) matters *******
  
ggplot(Akdemir_table_L_t1, aes(x= tissue_names, fill=altered_variance)) +
  geom_col() +
  ggtitle("Tissue Types by Count of Altered Variance") +
  ylab("Count by n() of Altered Variance") +xlab("Tissue Names")

ggplot(Akdemir_table_L_t1, aes(x= tissue_names, y=altered_variance))+ 
  geom_col()+
  ggtitle("Tissue Types by Count of Altered Variance")+ 
  ylab("Count by n() of Altered Variance")+
  xlab("Tissue Names")

----------------------------------------------------------------------------------------------------------------------------------

  
Akdemir_summary <- Akdemir_table_L_t1 %>% group_by(tissue_names, altered_variance, breakpoint) %>% summarize(n()) #Creating an Akdemir table that groups by tissue, and altered variance and breakpoint
Akdemir_summary_alt <- drop_na(Akdemir_summary) #Removes NA's from the previous table
sum(Akdemir_summary_alt$`n()`) #Turns the table into count data
#lowercase s

#Uppercase "Summary" in the name
write.csv(Akdemir_summary_alt, file.path(path, "Akdemir_Summary_alt_NA_removed.csv"), row.names=FALSE) #Exporting the Akdeimr_summary_alt table for future uses  
  
  
#Fisher's Tests and Firth & Flic Analysis  
    
library(readr)
Akdemir_Summary_alt <- read_csv("C:/Users/Vainc/Desktop/TAD Datasets/Akdemir_Summary_alt_NA_removed_altered.csv") #Importing the updated/altered version of the Akdemir_Summary_alt_NA_removed file above after finding proportions and up scaling by 500,000
View(Akdemir_Summary_alt)

Akdemir_Matrix <- Akdemir_Summary_alt %>% group_by(altered_variance,breakpoint) %>% summarise(N=sum(upsampled_proportion)) #Grouping by altered variance and breakpoint as well as summarizing the count of the up sampled proportions
view(Akdemir_Matrix)

fisher.test(Akdemir_Matrix)

Akdemir_Matrix_adj <- as.matrix(rbind(c(11951, 485721), c(20, 2388))) #Manually creating the matrix from the summarized data in the Akdemir_Matrix (up scales each number by 20 so there is no value with 0)

#fisher.test(Akdemir_Matrix) #Doesn't work with categories containing 0 values
fisher.test(Akdemir_Matrix_adj) #Likelyhood of getting the data distributed in the proportions you have, how likely is it to come up with proportions "this" extreme or more extreme?  2x2 Odds Ratio - for the upper left hand corner a breakpoints is 2 1/2 times (2.6 times) more likely to have altered variance. What is the likelihood of this distribution of counts across these categories by chance? Does altered variance correspond to

write.csv(Akdemir_Matrix, file.path(path, "Akdemir Matrix_new.csv"), row.names = FALSE) #Exporting the Akdemir Matrix

#Small loop that expands the data into a firth_table that displays the tissue name and the presence or absence of altered variance by binary output
for (i in 1:dim(Akdemir_Summary_alt)[1]){
  replicates <- as.numeric(paste((Akdemir_Summary_alt[i,6])))
  xrow <- Akdemir_Summary_alt[i,1:3]
  x <- do.call("rbind", replicate(replicates,xrow, simplify = FALSE))
  if(i == 1) {
    firth_table <- x} else {firth_table <- rbind(firth_table,x)}
}

#Visualizing the table
firth_table %>% glimpse()

#Removing NA's from the table
firth_table <- na.omit(firth_table)

#---------------------------------------------------------------------------------------------------------------------------

#Logistic with up sampled data, Fisher test with 1 observation is alright (validified)
#LogistF package

library(logistf)
#Firth done to reconfirm Fisher's.

fit<-logistf(altered_variance ~ breakpoint+tissue_names, data=Akdemir_Summary_alt) #Assigning information into a variable called "fit" that associates altered variance to breakpoint & tissue names from the Akdemir_Summary_alt table created earlier
summary(fit)

firth_T_flic_F = logistf(altered_variance ~ breakpoint+tissue_names, data = firth_table, firth = TRUE, flic = FALSE) #Worked
firth_T_flic_F #maxit iterration errors (p-value might be off)

#firth_T_flic_T = logistf(altered_variance ~ breakpoint+tissue_names, data=firth_table, firth = TRUE, flic = TRUE)
#firth_T_flic_T

firth_F_flic_F = logistf(altered_variance ~ breakpoint+tissue_names, data=firth_table, firth = FALSE, flic = FALSE) #Worked
firth_F_flic_F #3 separate error messages, Maximum iterations, parameter plcontrol, need larger iteration # etc.

#firth_F_flic_T = logistf(fit,Akdemir_Summary_alt, firth = FALSE, flic = TRUE)
#firth_F_flic_T

summary(firth_T_flic_F) #Reconfirms the Fisher's Test the odds of obtaining the altered_variance + breakpoint (both being 1) is slim by chance.
summary(firth_F_flic_F)

write.csv(firth_table, file.path(path, "Firth_Table.csv"), row.names=FALSE)

#Y - altered_variance
#X - breakpoint

firth_T_flic_F_test = logistf(breakpoint ~ altered_variance+tissue_names, data = firth_table, firth = TRUE, flic = FALSE)
summary(firth_T_flic_F_test)

#Small ggplot to display the Akdemir_Summary_alt table that displays the up sampled proportions by tissue with the altered variance and breakpoints

ggplot(Akdemir_Summary_alt, aes(tissue_names, upsampled_proportion, fill = as.factor(altered_variance), group = as.factor(breakpoint))) + geom_col(position = "dodge") + scale_y_log10() + scale_fill_brewer(palette = "Paired") +labs(title="Tissue Type by Altered Variance with Upsampled Proportions") + xlab("Tissue Names") + ylab("Upsampled Proportion") + labs(fill = "Altered Variance")


#-------------------------------------------------------------------------------------------------------------------------------------

#IMPORTING ENCODE CELL LINES

#Loading the CANCER TAD Cell Lines

setwd("C:/Users/Vainc/Desktop/TAD Datasets")

#ENCODE and 4D Nucleus project TADs
#Had to rename columns to match with previous data tables

#A549
library(readr)
ENCODE3_A549_HindIII_hg19_genome_C_40000_iced_tads <- read_delim("ENCODE Cell Lines/ENCFF336WPU_A549.bed/ENCODE3-A549-HindIII__hg19__genome__C-40000-iced.tads.bed", 
                                                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                                                 trim_ws = TRUE)
View(ENCODE3_A549_HindIII_hg19_genome_C_40000_iced_tads) #Viewing the file for proper import
A549 <- ENCODE3_A549_HindIII_hg19_genome_C_40000_iced_tads %>% #Renaming the imported file columns
  rename("X1" = "chr", "X2" = "t_start", "X3" = "t_stop", "X4" = "TAD", "X5" = "?")

A549 <- as.data.table(A549) #Setting as a data table
View(A549)


#TD47
library(readr)
ENCODE3_T470_HindIII_hg19_genome_C_40000_iced_tads <- read_delim("ENCODE Cell Lines/ENCFF437EBV_TD47.bed/ENCODE3-T470-HindIII__hg19__genome__C-40000-iced.tads.bed", 
                                                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                                                 trim_ws = TRUE)
View(ENCODE3_T470_HindIII_hg19_genome_C_40000_iced_tads)
TD47 <- ENCODE3_T470_HindIII_hg19_genome_C_40000_iced_tads %>% 
  rename("X1" = "chr", "X2" = "t_start", "X3" = "t_stop", "X4" = "TAD", "X5" = "?")
TD47 <- as.data.table(TD47)
View(TD47)

#H460
library(readr)
ENCODE3_NCIH460_HindIII_hg19_genome_C_40000_iced_tads <- read_delim("ENCODE Cell Lines/ENCFF451MCF_H460.bed/ENCODE3-NCIH460-HindIII__hg19__genome__C-40000-iced.tads.bed", 
                                                                    "\t", escape_double = FALSE, col_names = FALSE, 
                                                                    trim_ws = TRUE)
View(ENCODE3_NCIH460_HindIII_hg19_genome_C_40000_iced_tads)
H460 <- ENCODE3_NCIH460_HindIII_hg19_genome_C_40000_iced_tads %>% 
  rename("X1" = "chr", "X2" = "t_start", "X3" = "t_stop", "X4" = "TAD", "X5" = "?")

H460 <- as.data.table(TD47)
View(H460)


#ACHN
library(readr)
ENCFF477VVJ_ACHN <- read_delim("ENCODE Cell Lines/ENCFF477VVJ_ACHN.bed/ENCFF477VVJ_ACHN.bed", 
                               "\t", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)
View(ENCFF477VVJ_ACHN)
ACHN <- ENCFF477VVJ_ACHN %>% 
  rename("X1" = "chr", "X2" = "t_start", "X3" = "t_stop", "X4" = "TAD", "X5" = "?")

ACHN <- as.data.table(ACHN)
View(ACHN)

#SK-MEL5
library(readr)
ENCODE3_SKMEL5_HindIII_hg19_genome_C_40000_iced_tads <- read_delim("ENCODE Cell Lines/ENCFF784LMI_SK-MEL-5.bed/ENCODE3-SKMEL5-HindIII__hg19__genome__C-40000-iced.tads.bed", 
                                                                   "\t", escape_double = FALSE, col_names = FALSE, 
                                                                   trim_ws = TRUE)
View(ENCODE3_SKMEL5_HindIII_hg19_genome_C_40000_iced_tads)
SKMEL5 <- ENCODE3_SKMEL5_HindIII_hg19_genome_C_40000_iced_tads %>% 
  rename("X1" = "chr", "X2" = "t_start", "X3" = "t_stop", "X4" = "TAD", "X5" = "?")

SKMEL5 <- as.data.table(SKMEL5)
View(SKMEL5)



#---------------------------------------------------------------------------------------------------------------------------------------


#Loading the TAD Datasets

#Use DataTable for the below cell lines, Tidyverse is easier to learn due to syntax but DataTable is more memory efficient. These cell lines contain cancerous TADs (they're a direct measurement of how the 3D structure has been altered within the cell lines) Compare with the TADs in Akdemir for each cancerous cell lines
#Going to use foverlaps and find out which ones are equivalent. 


Akdemir_table_L_t1_1 <- Akdemir_table_L_t1 %>%
  rename("TAD" = "AkdemirTAD") #Renaming TADs to differentiate between the Akdemir & ENCODE TADs
Akdemir_table_L_t1_1_lung <- Akdemir_table_L_t1_1 %>% 
  filter(tissue_names == "lung") #Filtering out only the lung tissue since A549 is the cancerous version of lung tissue
Akdemir_table_L_t1_1_lung <- as.data.table(Akdemir_table_L_t1_1_lung) #Setting as data table (not tibble)
View(Akdemir_table_L_t1_1_lung)

###Specifically A549
Akdemir_table_L_t1_1_lung_alt <- Akdemir_table_L_t1_1 %>%  #Searching only for the tissues that contain the "A549" string from the cell_line column
  filter(stringr::str_detect(cell_line, 'A549'))
Akdemir_table_L_t1_1_lung_alt <- as.data.table(Akdemir_table_L_t1_1_lung_alt) #Setting as data table
setkey(Akdemir_table_L_t1_1_lung_alt, "chr", "t_start", "t_stop") #Setting keys for the three specified columns for Akdemir & A549
setkey(A549, "chr", "t_start", "t_stop")

Akdemir_A549_alt <- foverlaps(Akdemir_table_L_t1_1_lung_alt, A549, #######  T-test using this table for difference in sizes between Akdemir and Encode (for box plot)
                              by.x = c("chr", "t_start", "t_stop"), 
                              by.y = c("chr", "t_start", "t_stop"), 
                              type = "any", nomatch = NA) %>% glimpse()


Akdemir_A549_alt <- Akdemir_A549_alt %>% mutate(Encode_size = abs(t_stop - t_start), Akdemir_size = abs(i.t_stop - i.t_start)) %>% glimpse() #Making a new column called "Encode_size" that is the absolute value of the start coordinate subtracted from the stop coordinate (would be negative otherwise)

Akdemir_A549_alt <- Akdemir_A549_alt %>% mutate(TAD_diff_size = Encode_size - Akdemir_size) #Creating a "TAD_diff_size" column from subtracting the Akdemir_size column from the Encode_size column to find the differences in sizes between the two

#Finding altered variance (present) when filtering both altered variance and breakpoint while excluding NA variables
Akdemir_TADs_A549_alt_altvarp <- Akdemir_A549_alt %>% filter(breakpoint == 1, altered_variance ==1, !is.na(TAD), !is.na(AkdemirTAD)) %>% glimpse()

#Finding altered variance (not present) when filtering both altered variance and breakpoint while excluding NA variables
Akdemir_TADs_A549_alt_altvarn <- Akdemir_A549_alt %>% filter(breakpoint ==1, altered_variance ==0, !is.na(TAD), !is.na(AkdemirTAD)) %>% glimpse()

#Conducting the Welch's T-test between the groups that have altered variance vs those than do not (comparing their TAD_diff_size) on a two-sided test
Akdemir_A549_diff_test <- t.test(Akdemir_TADs_A549_alt_altvarp$TAD_diff_size, Akdemir_TADs_A549_alt_altvarn$TAD_diff_size, alternative = "two.sided")

#Altered variance present TADs do not disrupt as much as those that do not have altered variance present. Did t-test on the cell lines because A549 suggested that there might be a difference. Things with altered variance, which should predict breakpoint, they tended to result in a size of TAD that was different from the ones from ENCODE vs Akdemir (a change in size). But it turns out that it's only true in the A549 (oddball) the other four cell lines do not show evidence. The disruption of the TAD of the other four cell lines doesn't appear to correspond to variable expression.

#Results of the Welch's T-test to find the differences in size between the Akdemir TADs that have altered variance vs those than do not
Akdemir_A549_diff_test 


#----------------------------------------------------------------------------------------------------------------------------

#Creating an Akdemir_A549 box plot 

Akdemir_A549_boxplot <- Akdemir_A549_alt %>%filter(!is.na(TAD), !is.na(AkdemirTAD), breakpoint ==1) %>% glimpse() #Removing NA variables while filtering for a present breakpoint

#Creating a boxplot utilizing altered variance and the difference in TAD sizes (ENCODE - Akdemir)
ggplot(Akdemir_A549_boxplot, aes(as.factor(altered_variance),TAD_diff_size)) + geom_boxplot()
A549_plot_f <- ggplot(Akdemir_A549_boxplot, aes(as.factor(altered_variance),TAD_diff_size)) + geom_boxplot()+
  labs(x="Presence or Absence of Altered Variance", y="Difference in TAD Sizes")
A549_plot_f + ggtitle("Plot of A549 vs Akdemir TADs") + #Adding the title and the x and y labels to the plot
  xlab("Presence or Absence of Altered Variance") + ylab("Difference in TAD Size")


#Separate Akdemir_A549 table that filters out TAD_diff_size less than one million and groups by altered variance and summarizing the missing TADS for both output by count
#Displays how many TADs do not overlap
#Making sure the 0 result obtain was not due to TADs with altered variance not overlapping (missing category)
#Tables also double check that the Akdemir TAD is not missing in the overlap, first table (g1) checks without breakpoint second table 
#(g2) includes breakpoint while checking for the Akdemir & ENCODE TADs (hence the double-negative of missing TAD and "false" )
#A549 only has two overlapping TADs which is also displayed in the Forest Plot below

Akdemir_A549_alt_g1 <- Akdemir_A549_alt %>% filter(TAD_diff_size <= 1000000) %>% group_by(alt.variance = altered_variance, missing_A549_TAD = is.na(TAD), missing_Akdemir_TAD = is.na(AkdemirTAD)) %>% summarize(N = n())  #Filter


Akdemir_A549_alt_g2 <- Akdemir_A549_alt %>% filter(TAD_diff_size <= 1000000) %>% group_by(alt.variance = altered_variance, breakpoint = breakpoint, missing_A549_TAD = is.na(TAD), missing_Akdemir_TAD = is.na(AkdemirTAD)) %>% summarize(N = n()) 

A549_test_g <- Akdemir_A549_alt %>% group_by(alt.variance = altered_variance, missing_A549_TAD = is.na(TAD), missing_Akdemir_TAD = is.na(AkdemirTAD)) %>% summarize(N = n()) #Created to test if the A549 TADs work properly convered into a matrix, is still missing one category

Akdemir_A549_alt_f <- Akdemir_A549_alt %>% filter(altered_variance == 1) #Filtering for presence of altered variance "1"
Pair_number <- seq(from = 1, to = dim(Akdemir_A549_alt_f)[1], 1) #Setting a sequencing from 1 to the dimensions of the Akdemir_A549_alt_f through all the TADs present
Akdemir_A549_alt_f <- Akdemir_A549_alt_f %>% mutate(Pair_number = Pair_number) %>% glimpse() #Adding a column called Pair_number to sort through within the Forest Plots using a random generator below

Akdemir_filter <- Akdemir_table_L_t1[,c("TAD", "chr", "t_start", "t_stop")] %>% glimpse() #Making a separate filter table for specific columns: TAD, Chromosome, Start & Stop coordinates to be used below


#Forest Plot to show position
#Objects are data.tables, must use proper syntax.

#A549_alt_f_forest_TAD <- Akdemir_A549_alt_f %>% select(TAD, chr, t_start, t_stop, Pair_number) %>% glimpse() #May have to use this version of code depending on table format or package update/version

A549_alt_f_forest_TAD <- Akdemir_A549_alt_f[,c("TAD", "chr", "t_start", "t_stop", "Pair_number")] %>%glimpse() #Creating a filtered table for the forest plot for A549 (ENCODE) 
names(A549_alt_f_forest_TAD) <- c("TAD","chr","start","stop","pair") #Renaming the selected columns so tables merge correctly


#A549_alt_f_forest_Akdemir <- Akdemir_A549_alt_f %>% select(AkdemirTAD, chr, i.t_start, i.t_stop,Pair_number) %>% glimpse()
A549_alt_f_forest_Akdemir <- Akdemir_A549_alt_f[,c("AkdemirTAD", "chr", "i.t_start", "i.t_stop", "Pair_number")] %>% glimpse() #Creating a filtered table for the forest plot for Akdemir TAD
names(A549_alt_f_forest_Akdemir) <- c("TAD","chr","start","stop","pair") #Renaming the selected columns so tables merge correctly

A549_f_forest_mini <- rbind(A549_alt_f_forest_TAD, A549_alt_f_forest_Akdemir) %>% glimpse() #Binding the Akdemir and ENCODE f_forest tables created above together to find the TADs that overlap
A549_f_forest_mini <- A549_f_forest_mini %>% mutate(Midpoint = (start + abs((stop-start)/2))) %>% glimpse() #Establishing a Midpoint between Akdemir and the A549 ENCODE TADs (the start coordinates plus the absolute value of the start coordinates subtracted from the stop coordinates divided by two)


A549_f_forest_mini <- A549_f_forest_mini %>% #Creating a new column on a case_when basis to add the ENCODE or Akdemir category to the TADs with lowercase "tad" string being ENCODE and the TADs starting with "chr" being Akdemir (so the two are not confused)
  mutate(class = case_when(
    str_detect(TAD, "tad") ~ "ENCODE",
    str_detect(TAD, "chr") ~ "Akdemir",
    TRUE ~ NA_character_)) %>% glimpse()

write.csv(A549_f_forest_mini, file.path(path, "A549_forest_mini.csv"), row.names=FALSE) #Writing the forest plot table as a csv for import and reloading purposes (mostly incase something breaks)



#A549
setkey(Akdemir_table_L_t1_1_lung, "chr", "t_start", "t_stop") #Using setkeys to set keys between the Akdemir table that only contains lung tissue and the A549 table loaded previously to do a foverlaps between the two
setkey(A549, "chr", "t_start", "t_stop")

Akdemir_A549 <- foverlaps(Akdemir_table_L_t1_1_lung, A549, #foverlaps to combine the data of the two tables that contain only the chromosome, start coordinates and stop coordinates while removing NA's present
                          by.x = c("chr", "t_start", "t_stop"), 
                          by.y = c("chr", "t_start", "t_stop"), 
                          type = "any", nomatch = NA) %>% glimpse()

#THIS PROCESS IS REPEATED FOR THE OTHER ENCODE CANCEROUS CELL LINES AND AKDEMIR TADS BELOW (AFTER FOREST PLOTS)
---------------------------------------------------------------------------------------------------------------------------------------

#CREATING FOREST PLOTS FOR THE NORMAL AND CANCEROUS TADS  
    
install.packages("grid")
install.packages("ggExtra")
install.packages("gridBase")
install.packages("gridExtra")
library(grid)
library(ggExtra)
library(gridBase)
library(gridExtra) 

#USE THIS CODE FOR CREATING FOREST PLOTS 


#This loop is created for making Forest Plots for the Akdemir and ENCODE cell lines (separately by cancerous cell line)
#It is a "Generalized" loop for the ENCODE cell lines, but the A549 has a different setup for the grid.arrange code section since it only displayed two overlapping TADs

unique(A549_f_forest_mini$pair) #Code to test the number of pairs (appears negative do not be startled)
A549_f_forest_mini$pair #Code to test number of pairs (confirms above pair number but for each "class")

pairs <- unique(A549_f_forest_mini$pair) #Creating a unique list for the pairs column (no duplicates)
#A549_f_forest_mini$pair %>% unique() for data.table syntax
plots <- vector(mode="list", length = 3) #creating a list called "plots" so that only three pairs are pulled when displaying the Forest Plot image (3 TADs)
randgen <- as.integer(runif(3, min = 1, max = length(pairs))) #Randomly generates plot numbers based on pair (can sometimes duplicate, measures taken to not duplicate)
for (i in 1:3){
  print(i)
  aux <- A549_f_forest_mini %>% filter(pair == randgen[i]) #Filtering pairs and applying the filtered pairs to the randgen found only from within the A549 forest plot table created above
  plot1 <- ggplot(data=aux,aes(x = pair,y = Midpoint, ymin = start, ymax = stop, shape = class)) + #Creating the Forest Plot using ggplot
    geom_pointrange(aes(col=class)) +
    labs(title = "Overlap Between TADS: Akdemir and ENCODE A549", x=paste(aux$chr[1]), y="Genomic Position in bp") +
    guides(shape =guide_legend(title=NULL), col=guide_legend(title=NULL))+
    geom_errorbar(aes(ymin=start, ymax=stop,col=class), width=0.5,cex=1)+ 
    theme(plot.title=element_text(size=16,face="bold"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_text(face="bold"),
          axis.title=element_text(size=12,face="bold"),
          strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+ 
    coord_flip()
  plots[[i]] <- plot1
}


plots[i] #Checks if the plot did its job (for at least one TAD and displays a picture of a TAD comparison between Akdemir and ENCODE)
grid.arrange(plots[[1]],plots[[2]],plots[[3]], name = "arrange", as.table = TRUE, nrow = 3, ncol=1) #Creates an arrangement of three TAD plots that show the size difference between Akdmier and ENCODE, for this case it's the A549 compared with Akdemir
#Will display "Error in plots[[3]] : subscript out of bounds when run for A549 since A549 only contains two overlapping TADs, not 3
#General loop above used for the rest of the ENCODE & Akdemir TAD size comparisons, should work fine without any special treatments compared to A549


#Updated Forest Plot code specifically for the A459 Cancerous Line. Only had two overlapping TAD pairs between Akdemir and ENCODE so the loop had to be altered as such.  

unique(A549_f_forest_mini$pair) #Code to test the number of pairs (appears negative do not be startled)
A549_f_forest_mini$pair #Code to test number of pairs (confirms above pair number but for each "class")

pairs <- unique(A549_f_forest_mini$pair)
#A549_f_forest_mini$pair %>% unique() for data.table syntax
plots <- vector(mode="list", length = 2) 
#randgen <- as.integer(runif(2, min = 1, max = length(pairs))) #Randomly generates plot numbers based on pair (can sometimes duplicate)
for (i in 1:2){
  print(i)
  aux <- A549_f_forest_mini %>% filter(pair == pairs[i])
  plot1 <- ggplot(data=aux,aes(x = pair,y = Midpoint, ymin = start, ymax = stop, shape = class)) +
    geom_pointrange(aes(col=class)) +
    labs(title = "Overlap Between TADS: Akdemir and ENCODE A549", x=paste(aux$chr[1]), y="Genomic Position in bp") +
    guides(shape =guide_legend(title=NULL), col=guide_legend(title=NULL))+
    geom_errorbar(aes(ymin=start, ymax=stop,col=class), width=0.5,cex=1)+ 
    theme(plot.title=element_text(size=16,face="bold"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_text(face="bold"),
          axis.title=element_text(size=12,face="bold"),
          strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+ 
    coord_flip()
  plots[[i]] <- plot1
}


plots[i]
grid.arrange(plots[[1]],plots[[2]], name = "arrange", as.table = TRUE, nrow = 2, ncol=1) 



----------------------------------------------------------------------------------------------------------------------------------------
  
  
#TD47
Akdemir_table_L_t1_1_breast = Akdemir_table_L_t1_1[tissue_names == "breast",] #WORKS IN CONSOLE NOT IN CODE SECTION!!!
#Akdemir_table_L_t1_1[tissue_names == "breast",]
#Akdemir_table_L_t1_1 %>% filter(tissue_names == "breast") Akdemir_table_L_t1_1[tissue_names == "breast",] #assignment is being weird
Akdemir_table_L_t1_1_breast <- as.data.table(Akdemir_table_L_t1_1_breast)
View(Akdemir_table_L_t1_1_breast)

setkey(Akdemir_table_L_t1_1_breast, "chr", "t_start", "t_stop")
setkey(TD47, "chr", "t_start", "t_stop")

Akdemir_TD47 <- foverlaps(Akdemir_table_L_t1_1_breast, TD47, 
                          by.x = c("chr", "t_start", "t_stop"), 
                          by.y = c("chr", "t_start", "t_stop"), 
                          type = "any", nomatch = NA) %>% glimpse()
View(Akdemir_TD47)

#Specifically TD47

Akdemir_table_L_t1_1_breast_alt <- Akdemir_table_L_t1_1 %>% 
  filter(stringr::str_detect(cell_line, 'T-47D'))
Akdemir_table_L_t1_1_breast_alt <- as.data.table(Akdemir_table_L_t1_1_breast_alt)
setkey(Akdemir_table_L_t1_1_breast_alt, "chr", "t_start", "t_stop")
setkey(TD47, "chr", "t_start", "t_stop")
Akdemir_TD47_alt <- foverlaps(Akdemir_table_L_t1_1_breast_alt, TD47, 
                              by.x = c("chr", "t_start", "t_stop"), 
                              by.y = c("chr", "t_start", "t_stop"), 
                              type = "any", nomatch = NA) %>% glimpse()


Akdemir_TD47_alt <- Akdemir_TD47_alt %>% mutate(Encode_size = abs(t_stop - t_start), Akdemir_size = abs(i.t_stop - i.t_start)) %>% glimpse()
Akdemir_TD47_alt <- Akdemir_TD47_alt %>% mutate(TAD_diff_size = Encode_size - Akdemir_size)
Akdemir_TADs_TD47_alt_altvarp <- Akdemir_TD47_alt %>% filter(breakpoint == 1, altered_variance ==1, !is.na(TAD), !is.na(AkdemirTAD)) %>% glimpse()
Akdemir_TADs_TD47_alt_altvarn <- Akdemir_TD47_alt %>% filter(breakpoint ==1, altered_variance ==0, !is.na(TAD), !is.na(AkdemirTAD)) %>% glimpse()
Akdemir_TD47_diff_test <- t.test(Akdemir_TADs_TD47_alt_altvarp$TAD_diff_size, Akdemir_TADs_TD47_alt_altvarn$TAD_diff_size, alternative = "two.sided") #Altered variance present TADs do not disrupt as much as those that do not have altered variance present.
Akdemir_TD47_diff_test


Akdemir_TD47_boxplot <- Akdemir_TD47_alt %>%filter(!is.na(TAD), !is.na(AkdemirTAD), breakpoint ==1) %>% glimpse()
ggplot(Akdemir_TD47_boxplot, aes(as.factor(altered_variance),TAD_diff_size)) + geom_boxplot()

TD47_plot_f <- ggplot(Akdemir_TD47_boxplot, aes(as.factor(altered_variance),TAD_diff_size)) + geom_boxplot()

TD47_plot_f + ggtitle("Plot of TD47 vs Akdemir TADs") +
  xlab("Presence or Absence of Altered Variance") + ylab("Difference in TAD Size")


Akdemir_TD47_alt_g <- Akdemir_TD47_alt %>% filter(TAD_diff_size <=1000000) %>% group_by(alt.variance = altered_variance, breakpoint = breakpoint, missing_TD47_TAD = is.na(TAD), missing_Akdemir_TAD = is.na(AkdemirTAD)) %>% summarize(N = n())


TD47_g_test <- Akdemir_TD47_alt %>% filter(TAD_diff_size <= 1000000) %>% group_by(alt.variance = altered_variance, missing_TD47_TAD = is.na(TAD), missing_Akdemir_TAD = is.na(AkdemirTAD)) %>% summarize(N = n())


Akdemir_TD47_alt_f <- Akdemir_TD47_alt %>% filter(altered_variance == 1, breakpoint == 1, !is.na(TAD), !is.na(AkdemirTAD)) #Returning "FALSE"
Pair_number <- seq(from = 1, to = dim(Akdemir_TD47_alt_f)[1], 1)
Akdemir_TD47_alt_f <- Akdemir_TD47_alt_f %>% mutate(Pair_number = Pair_number) %>% glimpse()
#Forest Plot to show position


#TD47_alt_f_forest_TAD <- Akdemir_TD47_alt_f %>% select(TAD, chr, t_start, t_stop, Pair_number) %>% glimpse()
TD47_alt_f_forest_TAD <- Akdemir_TD47_alt_f[,c("TAD", "chr", "t_start", "t_stop", "Pair_number")] %>%glimpse()
names(TD47_alt_f_forest_TAD) <- c("TAD","chr","start","stop","pair")


#TD47_alt_f_forest_Akdemir <- Akdemir_TD47_alt_f %>% select(AkdemirTAD, chr, i.t_start, i.t_stop,Pair_number) %>% glimpse()
TD47_alt_f_forest_Akdemir <- Akdemir_TD47_alt_f[,c("AkdemirTAD", "chr", "i.t_start", "i.t_stop", "Pair_number")] %>% glimpse()
names(TD47_alt_f_forest_Akdemir) <- c("TAD","chr","start","stop","pair")

TD47_f_forest_mini <- rbind(TD47_alt_f_forest_TAD, TD47_alt_f_forest_Akdemir) %>% glimpse()
TD47_f_forest_mini <- TD47_f_forest_mini %>% mutate(Midpoint = (start + abs((stop-start)/2))) %>% glimpse()

TD47_f_forest_mini <- TD47_f_forest_mini %>%
  mutate(class = case_when(
    str_detect(TAD, "tad") ~ "ENCODE",
    str_detect(TAD, "chr") ~ "Akdemir",
    TRUE ~ NA_character_)) %>% glimpse()

#Plots
pairs <- TD47_f_forest_mini$pair %>% unique()
plots <- vector(mode="list", length = 3)
randgen <- as.integer(runif(3, min = 1, max = length(pairs)))
for (i in 1:3){
  print(i)
  aux <- TD47_f_forest_mini %>% filter(pair == randgen[i])
  plot1 <- ggplot(data=aux,aes(x = pair,y = Midpoint, ymin = start, ymax = stop, shape = class)) +
    geom_pointrange(aes(col=class)) +
    labs(title = "Overlap Between TADS: Akdemir and ENCODE TD47", x=paste(aux$chr[1]), y="Genomic Position in bp") +
    guides(shape =guide_legend(title=NULL), col=guide_legend(title=NULL))+
    geom_errorbar(aes(ymin=start, ymax=stop,col=class), width=0.5,cex=1)+ 
    #facet_wrap(~call,strip.position="left",nrow=9,scales = "free_y") +
    theme(plot.title=element_text(size=16,face="bold"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_text(face="bold"),
          axis.title=element_text(size=12,face="bold"),
          strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+ 
    coord_flip()
  plots[[i]] <- plot1
}


#plots[i]
grid.arrange(plots[[1]],plots[[2]],plots[[3]], name = "arrange", as.table = TRUE, nrow = 3, ncol=1)


#--------------------------------------------------------------------------------------------------------------------------------------


#H460
setkey(Akdemir_table_L_t1_1_lung, "chr", "t_start", "t_stop")
setkey(H460, "chr", "t_start", "t_stop")

Akdemir_H460 <- foverlaps(Akdemir_table_L_t1_1_lung, H460, 
                          by.x = c("chr", "t_start", "t_stop"), 
                          by.y = c("chr", "t_start", "t_stop"), 
                          type = "any", nomatch = NA) %>% glimpse()

#Specifically H460
Akdemir_table_L_t1_1_lung_alt <- Akdemir_table_L_t1_1 %>% 
  filter(stringr::str_detect(cell_line, 'H460'))
Akdemir_table_L_t1_1_lung_alt <- as.data.table(Akdemir_table_L_t1_1_lung_alt)
setkey(Akdemir_table_L_t1_1_lung_alt, "chr", "t_start", "t_stop")
setkey(H460, "chr", "t_start", "t_stop")
Akdemir_H460_alt <- foverlaps(Akdemir_table_L_t1_1_lung_alt, H460, 
                              by.x = c("chr", "t_start", "t_stop"), 
                              by.y = c("chr", "t_start", "t_stop"), 
                              type = "any", nomatch = NA) %>% glimpse()


Akdemir_H460_alt <- Akdemir_H460_alt %>% mutate(Encode_size = abs(t_stop - t_start), Akdemir_size = abs(i.t_stop - i.t_start)) %>% glimpse()
Akdemir_H460_alt <- Akdemir_H460_alt %>% mutate(TAD_diff_size = Encode_size - Akdemir_size)
Akdemir_TADs_H460_alt_altvarp <- Akdemir_H460_alt %>% filter(breakpoint == 1, altered_variance ==1, !is.na(TAD), !is.na(AkdemirTAD)) %>% glimpse()
Akdemir_TADs_H460_alt_altvarn <- Akdemir_H460_alt %>% filter(breakpoint ==1, altered_variance ==0, !is.na(TAD), !is.na(AkdemirTAD)) %>% glimpse()
Akdemir_H460_diff_test <- t.test(Akdemir_TADs_H460_alt_altvarp$TAD_diff_size, Akdemir_TADs_H460_alt_altvarn$TAD_diff_size, alternative = "two.sided") #Altered variance present TADs do not disrupt as much as those that do not have altered variance present.

Akdemir_H460_diff_test

Akdemir_H460_boxplot <- Akdemir_H460_alt %>%filter(!is.na(TAD), !is.na(AkdemirTAD), breakpoint ==1) %>% glimpse()
ggplot(Akdemir_H460_boxplot, aes(as.factor(altered_variance),TAD_diff_size)) + geom_boxplot()
H460_plot_f <- ggplot(Akdemir_H460_boxplot, aes(as.factor(altered_variance),TAD_diff_size)) + geom_boxplot()

H460_plot_f + ggtitle("Plot of H460 vs Akdemir TADs") +
  xlab("Presence or Absence of Altered Variance") + ylab("Difference in TAD Size")

Akdemir_H460_alt_g <- Akdemir_H460_alt %>% filter(TAD_diff_size <= 1000000) %>% group_by(alt.variance = altered_variance, missing_H460_TAD = is.na(TAD), missing_Akdemir_TAD = is.na(AkdemirTAD)) %>% summarize(N = n())  #Filter

Akdemir_H460_alt_f <- Akdemir_H460_alt %>% filter(altered_variance == 1)
Pair_number <- seq(from = 1, to = dim(Akdemir_H460_alt_f)[1], 1)
#Pair_number <- seq(from = 1, to = 21, 1)
Akdemir_H460_alt_f <- Akdemir_H460_alt_f %>% mutate(Pair_number = Pair_number) %>% glimpse()


#Forest Plot to show position

#H460_alt_f_forest_TAD <- Akdemir_H460_alt_f %>% select(TAD, chr, t_start, t_stop, Pair_number) %>% glimpse()
H460_alt_f_forest_TAD<- Akdemir_H460_alt_f[,c("TAD", "chr", "t_start", "t_stop", "Pair_number")] %>%glimpse()
names(H460_alt_f_forest_TAD) <- c("TAD","chr","start","stop","pair")


#H460_alt_f_forest_Akdemir <- Akdemir_H460_alt_f %>% select(AkdemirTAD, chr, i.t_start, i.t_stop,Pair_number) %>% glimpse()
H460_alt_f_forest_Akdemir <- Akdemir_H460_alt_f[,c("AkdemirTAD", "chr", "i.t_start", "i.t_stop", "Pair_number")] %>%glimpse()
names(H460_alt_f_forest_Akdemir) <- c("TAD","chr","start","stop","pair")

H460_f_forest_mini <- rbind(H460_alt_f_forest_TAD, H460_alt_f_forest_Akdemir) %>% glimpse()
H460_f_forest_mini <- H460_f_forest_mini %>% mutate(Midpoint = (start + abs((stop-start)/2))) %>% glimpse()


H460_f_forest_mini <- H460_f_forest_mini %>%
  mutate(class = case_when(
    str_detect(TAD, "tad") ~ "ENCODE",
    str_detect(TAD, "chr") ~ "Akdemir",
    TRUE ~ NA_character_)) %>% glimpse()

#Plots
pairs <- H460_f_forest_mini$pair %>% unique()
plots <- vector(mode="list", length = 3)
randgen <- as.integer(runif(10, min = 1, max = length(pairs)))
for (i in 1:3){
  print(i)
  aux <- H460_f_forest_mini %>% filter(pair == randgen[i])
  plot1 <- ggplot(data=aux,aes(x = pair,y = Midpoint, ymin = start, ymax = stop, shape = class)) +
    geom_pointrange(aes(col=class)) +
    labs(title = "Overlap Between TADS: Akdemir and ENCODE H460", x=paste(aux$chr[1]), y="Genomic Position in bp") +
    guides(shape =guide_legend(title=NULL), col=guide_legend(title=NULL))+
    geom_errorbar(aes(ymin=start, ymax=stop,col=class), width=0.5,cex=1)+ 
    #facet_wrap(~call,strip.position="left",nrow=9,scales = "free_y") +
    theme(plot.title=element_text(size=16,face="bold"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_text(face="bold"),
          axis.title=element_text(size=12,face="bold"),
          strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+ 
    coord_flip()
  plots[[i]] <- plot1
}


plots[i]
grid.arrange(plots[[1]],plots[[2]],plots[[3]], name = "arrange", as.table = TRUE, nrow = 3, ncol=1)


#---------------------------------------------------------------------------------------------------------------------------------------

#ACHN
Akdemir_table_L_t1_1_renal <- Akdemir_table_L_t1_1 %>% 
  filter(tissue_names == "renal")
Akdemir_table_L_t1_1_renal <- as.data.table(Akdemir_table_L_t1_1_renal)
View(Akdemir_table_L_t1_1_renal)

setkey(Akdemir_table_L_t1_1_renal, "chr", "t_start", "t_stop")
setkey(ACHN, "chr", "t_start", "t_stop")

Akdemir_ACHN <- foverlaps(Akdemir_table_L_t1_1_renal, ACHN, 
                          by.x = c("chr", "t_start", "t_stop"), 
                          by.y = c("chr", "t_start", "t_stop"), 
                          type = "any", nomatch = NA) %>% glimpse()
View(Akdemir_ACHN)

#Specifically ACHN

Akdemir_table_L_t1_1_renal_alt <- Akdemir_table_L_t1_1 %>% 
  filter(stringr::str_detect(cell_line, 'ACHN'))
Akdemir_table_L_t1_1_renal_alt <- as.data.table(Akdemir_table_L_t1_1_renal_alt)
setkey(Akdemir_table_L_t1_1_renal_alt, "chr", "t_start", "t_stop")
setkey(ACHN, "chr", "t_start", "t_stop")
Akdemir_ACHN_alt <- foverlaps(Akdemir_table_L_t1_1_renal_alt, ACHN, 
                              by.x = c("chr", "t_start", "t_stop"), 
                              by.y = c("chr", "t_start", "t_stop"), 
                              type = "any", nomatch = NA) %>% glimpse()


Akdemir_ACHN_alt <- Akdemir_ACHN_alt %>% mutate(Encode_size = abs(t_stop - t_start), Akdemir_size = abs(i.t_stop - i.t_start)) %>% glimpse()
Akdemir_ACHN_alt <- Akdemir_ACHN_alt %>% mutate(TAD_diff_size = Encode_size - Akdemir_size)
Akdemir_TADs_ACHN_alt_altvarp <- Akdemir_ACHN_alt %>% filter(breakpoint == 1, altered_variance ==1, !is.na(TAD), !is.na(AkdemirTAD)) %>% glimpse()
Akdemir_TADs_ACHN_alt_altvarn <- Akdemir_ACHN_alt %>% filter(breakpoint ==1, altered_variance ==0, !is.na(TAD), !is.na(AkdemirTAD)) %>% glimpse()
Akdemir_ACHN_diff_test <- t.test(Akdemir_TADs_ACHN_alt_altvarp$TAD_diff_size, Akdemir_TADs_ACHN_alt_altvarn$TAD_diff_size, alternative = "two.sided") #Altered variance present TADs do not disrupt as much as those that do not have altered variance present.

Akdemir_ACHN_diff_test


Akdemir_ACHN_boxplot <- Akdemir_ACHN_alt %>%filter(!is.na(TAD), !is.na(AkdemirTAD), breakpoint ==1) %>% glimpse()
ggplot(Akdemir_ACHN_boxplot, aes(as.factor(altered_variance),TAD_diff_size)) + geom_boxplot()
ACHN_plot_f <- ggplot(Akdemir_ACHN_boxplot, aes(as.factor(altered_variance),TAD_diff_size)) + geom_boxplot()

ACHN_plot_f + ggtitle("Plot of ACHN vs Akdemir TADs") +
  xlab("Presence or Absence of Altered Variance") + ylab("Difference in TAD Size")

Akdemir_ACHN_alt_g <- Akdemir_ACHN_alt %>% filter(TAD_diff_size <= 1000000) %>% group_by(alt.variance = altered_variance, missing_ACHN_TAD = is.na(TAD), missing_Akdemir_TAD = is.na(AkdemirTAD)) %>% summarize(N = n())  #Filter

Akdemir_ACHN_alt_f <- Akdemir_ACHN_alt %>% filter(altered_variance == 1)
Pair_number <- seq(from = 1, to = dim(Akdemir_ACHN_alt_f)[1], 1)
#Pair_number <- seq(from = 1, to = 15, 1)
Akdemir_ACHN_alt_f <- Akdemir_ACHN_alt_f %>% mutate(Pair_number = Pair_number) %>% glimpse()


#Forest Plot to show position

#ACHN_alt_f_forest_TAD <- Akdemir_ACHN_alt_f %>% select(TAD, chr, t_start, t_stop, Pair_number) %>% glimpse()
ACHN_alt_f_forest_TAD <- Akdemir_ACHN_alt_f[,c("TAD", "chr", "t_start", "t_stop", "Pair_number")] %>% glimpse()
names(ACHN_alt_f_forest_TAD) <- c("TAD","chr","start","stop","pair")


#ACHN_alt_f_forest_Akdemir <- Akdemir_ACHN_alt_f %>% select(AkdemirTAD, chr, i.t_start, i.t_stop,Pair_number) %>% glimpse()
ACHN_alt_f_forest_Akdemir <- Akdemir_ACHN_alt_f[,c("AkdemirTAD", "chr", "i.t_start", "i.t_stop", "Pair_number")] %>%glimpse()
names(ACHN_alt_f_forest_Akdemir) <- c("TAD","chr","start","stop","pair")

ACHN_f_forest_mini <- rbind(ACHN_alt_f_forest_TAD, ACHN_alt_f_forest_Akdemir) %>% glimpse()
ACHN_f_forest_mini <- ACHN_f_forest_mini %>% mutate(Midpoint = (start + abs((stop-start)/2))) %>% glimpse()


ACHN_f_forest_mini <- ACHN_f_forest_mini %>%
  mutate(class = case_when(
    str_detect(TAD, "tad") ~ "ENCODE",
    str_detect(TAD, "chr") ~ "Akdemir",
    TRUE ~ NA_character_)) %>% glimpse()

#Plots
pairs <- ACHN_f_forest_mini$pair %>% unique()
plots <- vector(mode="list", length = 3)
randgen <- as.integer(runif(10, min = 1, max = length(pairs)))
for (i in 1:3){
  print(i)
  aux <- ACHN_f_forest_mini %>% filter(pair == randgen[i])
  plot1 <- ggplot(data=aux,aes(x = pair,y = Midpoint, ymin = start, ymax = stop, shape = class)) +
    geom_pointrange(aes(col=class)) +
    labs(title = "Overlap Between TADS: Akdemir and ENCODE ACHN", x=paste(aux$chr[1]), y="Genomic Position in bp") +
    guides(shape =guide_legend(title=NULL), col=guide_legend(title=NULL))+
    geom_errorbar(aes(ymin=start, ymax=stop,col=class), width=0.5,cex=1)+ 
    #facet_wrap(~call,strip.position="left",nrow=9,scales = "free_y") +
    theme(plot.title=element_text(size=16,face="bold"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_text(face="bold"),
          axis.title=element_text(size=12,face="bold"),
          strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+ 
    coord_flip()
  plots[[i]] <- plot1
}


plots[i]
grid.arrange(plots[[1]],plots[[2]],plots[[3]], name = "arrange", as.table = TRUE, nrow = 3, ncol=1)



#How many of the invariant TADs are still present, what is the overlap between the Encode and the Akdemir TADs? For the TADs with altered variance, what's the deal? Are they the same in the actual measurement for that cell line? The TADs with altered variance with a breakpoint present are not the same in the measured cell line and are larger than what the common (Akedmir "invariant" TAD in the region is).


#SKL-MEL5
Akdemir_table_L_t1_1_skin <- Akdemir_table_L_t1_1 %>% 
  filter(tissue_names == "skin")
Akdemir_table_L_t1_1_skin <- as.data.table(Akdemir_table_L_t1_1_skin)
View(Akdemir_table_L_t1_1_skin)

setkey(Akdemir_table_L_t1_1_skin, "chr", "t_start", "t_stop")
setkey(SKMEL5, "chr", "t_start", "t_stop")

Akdemir_SKMEL5 <- foverlaps(Akdemir_table_L_t1_1_skin, SKMEL5, 
                            by.x = c("chr", "t_start", "t_stop"), 
                            by.y = c("chr", "t_start", "t_stop"), 
                            type = "any", nomatch = NA) %>% glimpse()

View(Akdemir_SKMEL5)

#Specficially SKL-MEL5  SK-MEL-5

Akdemir_table_L_t1_1_skin_alt <- Akdemir_table_L_t1_1 %>% 
  filter(stringr::str_detect(cell_line, 'SK-MEL-5'))
Akdemir_table_L_t1_1_skin_alt <- as.data.table(Akdemir_table_L_t1_1_skin_alt)
setkey(Akdemir_table_L_t1_1_skin_alt, "chr", "t_start", "t_stop")
setkey(SKMEL5, "chr", "t_start", "t_stop")
Akdemir_SKMEL5_alt <- foverlaps(Akdemir_table_L_t1_1_skin_alt, SKMEL5, 
                                by.x = c("chr", "t_start", "t_stop"), 
                                by.y = c("chr", "t_start", "t_stop"), 
                                type = "any", nomatch = NA) %>% glimpse()




Akdemir_SKMEL5_alt <- Akdemir_SKMEL5_alt %>% mutate(Encode_size = abs(t_stop - t_start), Akdemir_size = abs(i.t_stop - i.t_start)) %>% glimpse()
Akdemir_SKMEL5_alt <- Akdemir_SKMEL5_alt %>% mutate(TAD_diff_size = Encode_size - Akdemir_size)
Akdemir_TADs_SKMEL5_alt_altvarp <- Akdemir_SKMEL5_alt %>% filter(breakpoint == 1, altered_variance ==1, !is.na(TAD), !is.na(AkdemirTAD)) %>% glimpse()
Akdemir_TADs_SKMEL5_alt_altvarn <- Akdemir_SKMEL5_alt %>% filter(breakpoint ==1, altered_variance ==0, !is.na(TAD), !is.na(AkdemirTAD)) %>% glimpse()
Akdemir_SKMEL5_diff_test <- t.test(Akdemir_TADs_SKMEL5_alt_altvarp$TAD_diff_size, Akdemir_TADs_SKMEL5_alt_altvarn$TAD_diff_size, alternative = "two.sided") #Altered variance present TADs do not disrupt as much as those that do not have altered variance present.

Akdemir_SKMEL5_diff_test


Akdemir_SKMEL5_boxplot <- Akdemir_SKMEL5_alt %>%filter(!is.na(TAD), !is.na(AkdemirTAD), breakpoint ==1) %>% glimpse()
ggplot(Akdemir_SKMEL5_boxplot, aes(as.factor(altered_variance),TAD_diff_size)) + geom_boxplot()
SKMEL5_plot_f <- ggplot(Akdemir_SKMEL5_boxplot, aes(as.factor(altered_variance),TAD_diff_size)) + geom_boxplot()

SKMEL5_plot_f + ggtitle("Plot of SKMEL5 vs Akdemir TADs") +
  xlab("Presence or Absence of Altered Variance") + ylab("Difference in TAD Size")


Akdemir_SKMEL5_alt_g <- Akdemir_SKMEL5_alt %>% group_by(alt.variance = altered_variance, missing_SKMEL5_TAD = is.na(TAD), missing_Akdemir_TAD = is.na(AkdemirTAD)) %>% summarize(N = n()) 

#Akdemir_SKMEL5_alt_g <- Akdemir_SKMEL5_alt %>% filter(TAD_diff_size <= 1000000) %>% group_by(alt.variance = altered_variance, missing_SKMEL5_TAD = is.na(TAD), missing_Akdemir_TAD = is.na(AkdemirTAD)) %>% summarize(N = n())  #Filter

Akdemir_SKMEL5_alt_f <- Akdemir_SKMEL5_alt %>% filter(altered_variance == 1)
Pair_number <- seq(from = 1, to = dim(Akdemir_SKMEL5_alt_f)[1], 1)
#Pair_number <- seq(from = 1, to = 21, 1)
Akdemir_SKMEL5_alt_f <- Akdemir_SKMEL5_alt_f %>% mutate(Pair_number = Pair_number) %>% glimpse()



#Forest Plot to show position

#SKMEL5_alt_f_forest_TAD <- Akdemir_SKMEL5_alt_f %>% select(TAD, chr, t_start, t_stop, Pair_number) %>% glimpse()
SKMEL5_alt_f_forest_TAD <- Akdemir_SKMEL5_alt_f[,c("TAD", "chr", "t_start", "t_stop", "Pair_number")] %>% glimpse()
names(SKMEL5_alt_f_forest_TAD) <- c("TAD","chr","start","stop","pair")


#SKMEL5_alt_f_forest_Akdemir <- Akdemir_SKMEL5_alt_f %>% select(AkdemirTAD, chr, i.t_start, i.t_stop,Pair_number) %>% glimpse()
SKMEL5_alt_f_forest_Akdemir <- Akdemir_SKMEL5_alt_f[,c("AkdemirTAD", "chr", "i.t_start", "i.t_stop", "Pair_number")] %>%glimpse()
names(SKMEL5_alt_f_forest_Akdemir) <- c("TAD","chr","start","stop","pair")

SKMEL5_f_forest_mini <- rbind(SKMEL5_alt_f_forest_TAD, SKMEL5_alt_f_forest_Akdemir) %>% glimpse()
SKMEL5_f_forest_mini <- SKMEL5_f_forest_mini %>% mutate(Midpoint = (start + abs((stop-start)/2))) %>% glimpse()


SKMEL5_f_forest_mini <- SKMEL5_f_forest_mini %>%
  mutate(class = case_when(
    str_detect(TAD, "tad") ~ "ENCODE",
    str_detect(TAD, "chr") ~ "Akdemir",
    TRUE ~ NA_character_)) %>% glimpse()


unique(SKMEL5_f_forest_mini$pair) #Code to test the number of pairs (appears negative do not be startled)
SKMEL5_f_forest_mini$pair


#Plots
pairs <- SKMEL5_f_forest_mini$pair %>% unique()
plots <- vector(mode="list", length = 3)
randgen <- as.integer(runif(10, min = 1, max = length(pairs)))
for (i in 1:3){
  print(i)
  aux <- SKMEL5_f_forest_mini %>% filter(pair == randgen[i])
  plot1 <- ggplot(data=aux,aes(x = pair,y = Midpoint, ymin = start, ymax = stop, shape = class)) +
    geom_pointrange(aes(col=class)) +
    labs(title = "Overlap Between TADS: Akdemir and ENCODE SKMEL5", x=paste(aux$chr[1]), y="Genomic Position in bp") +
    guides(shape =guide_legend(title=NULL), col=guide_legend(title=NULL))+
    geom_errorbar(aes(ymin=start, ymax=stop,col=class), width=0.5,cex=1)+ 
    #facet_wrap(~call,strip.position="left",nrow=9,scales = "free_y") +
    theme(plot.title=element_text(size=16,face="bold"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_text(face="bold"),
          axis.title=element_text(size=12,face="bold"),
          strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+ 
    coord_flip()
  plots[[i]] <- plot1
}


plots[i]
grid.arrange(plots[[1]],plots[[2]],plots[[3]], name = "arrange", as.table = TRUE, nrow = 3, ncol=1)


#-------------------------------------------------------------------------------------------------------------------------------------

#Tables for each cell line - to see the matrix numeric outputs for each category (checking code output for errors) make sure the Forest Plot with N=2 (for A549), not N=40 but only 2 overlapped with Akdemir TADs.
#Essentially taking the cell_line column from each table and isolating the specific string that contains the ENCODE cell lines used in the study (A549, TD47, H460, ACHN, and SKMEL5) then grouping them by altered variance and breakpoint and applying counts for each combination of breakpoint and altered variance into a matrix

#A549
akdemir_select_A549_lc <- Akdemir_table_L_t1[cell_line=="LC_A549/ATCC",.N, by=c("altered_variance","breakpoint")]

#TD47
akdemir_select_TD47_br <- Akdemir_table_L_t1[cell_line=="BR_T-47D",.N, by=c("altered_variance","breakpoint")]

#H460
akdemir_select_H460_lc <- Akdemir_table_L_t1[cell_line=="LC_NCI-H460",.N, by=c("altered_variance","breakpoint")]
#ACHN
akdemir_select_ACHN_re <- Akdemir_table_L_t1[cell_line=="RE_ACHN",.N, by=c("altered_variance","breakpoint")]
#SKMEL5
akdemir_select_SKMEL5_sk <- Akdemir_table_L_t1[cell_line=="ME_SK-MEL-5",.N, by=c("altered_variance","breakpoint")]


#-------------------------------------------------------------------------------------------------------------------------------------

#Last-minute analysis to find out which chromosomes or TADs are causing the most issues, used a cuttoff of 15 cell lines (how many cell lines the TAD may contain, could go up to 37 (max) but we set the bare minimum of cell lines present to 15 (slightly less than half since a majority displayed 15 cell lines)) then grouped by TAD and summarized the output

Akdemir_assay <- Akdemir_table_L_t1 %>% filter(altered_variance == 1) %>% group_by(TAD) %>% summarise(n()) %>% rename(c("n()" = "N")) #Filtered only for altered variance (present = 1) and summarized the output into count data
Akdemir_assay_filter <- Akdemir_assay %>% filter(N >= 15) %>% glimpse() #Takes the previous table and filtered for the number of cell lines at or greater than 15

Akdemir_assay_aux3 <- merge(Akdemir_assay_filter, aux3, by.x = "TAD", by.y = "TAD") #Combining the new table with cell line count greater than 15 with the "aux3" created from the foverlaps (refresh below) of the Affymetrix U133AB probe coordinates and the Akdemir TADs (essentially combining TADs with genomic coordinates)
#This table contains two sets of genomic coordinates, one for the Akdemir TADs and the other for the ENCODE TADs, but also contains the TAD and number of cell lines present (N) within the TAD

#aux3  <- foverlaps( affy_U133AB_probe_coord, Akdemir_TADs_alt, 
                    #by.x = c("t_name", "t_start",   "t_end"), 
                    #by.y = c("t_name", "t_start",   "t_end"), 
                    #type = "within", nomatch = NA)

Akdemir_assay_aux3_filtered <- Akdemir_assay_aux3 %>% filter(t_name == "chr10" | t_name == "chr19") #Filtering for the chromosomes that contain all the cell lines (37), manually noticed that both Chromsome 10 & 19 contained all 37 cell lines or better spoken: Chromosomes 10 & 19 contain 37 of the 60 NCI-60 cell lines (though we only used so many due to data availability.)

write.csv(Akdemir_assay_aux3, file = "Akdemir_assay_aux3_probe_id", row.names = FALSE) #Exporting the unfiltered file for use in DAVID (https://david.ncifcrf.gov/) to find out which processes are involved by probe_id

write.csv(Akdemir_assay_aux3_filtered, file = "Akdemir_assay_aux3_probe_id_filtered", row.names = FALSE) #Exporting the filtered file

#----------------------------------------------------------------------------------------------------------------------------------------

#Filtering further from the table where adjusted p-values are less than 0.01 for the specific ENCODE cell line TADs

Akdemir_ENCODE_TAD <- Akdemir_table_L_filtered_t2 %>% filter(cell_line == "ME_SK-MEL-5" | cell_line == "LC_NCI-H460" | cell_line == "RE_ACHN" | cell_line == "LC_A549/ATCC" | cell_line == "BR_T-47D")


write.csv(Akdemir_ENCODE_TAD, file = "Akdemir_ENCODE_Filtered_TADs", row.names = FALSE)


#Tidyverse piping for how many TADs are there per cell line by count?

Akdemir_Encode_TAD_count <- Akdemir_ENCODE_TAD %>% group_by(cell_line) %>% summarise(count = n())

write.csv(Akdemir_ENCODE_TAD_count, file = "Akdemir_ENCODE_Filtered_TADs_by_cell_line", row.names = FALSE)
