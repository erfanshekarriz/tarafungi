library(phyloseq)
library(tidyverse)



# INPUT DATA ####
inputdf <- read.csv("./data/tabular/TARA_V9_database.csv")
metadfraw <- read.csv("./data/tabular/TARA_sample_enviro.csv")

# FILTER FOR FUNGI ONLY
fungidf <- inputdf %>% select(md5sum, lineage) %>%
  filter(grepl("Fungi", lineage, ignore.case = TRUE))



# TAXONOMY TABLE
tax <- fungidf %>% select(md5sum, lineage) %>%
  mutate(lineage = str_replace_all(lineage, c(`\\+`=" ", 
                                              `_`=" ", 
                                              X=""))) %>%
  separate(lineage, sep = "\\|", into = c("Kingdom", "Supergroup", 
                                        "Group", "Phylum", 
                                        "Sub-phylum", "Class", 
                                        "Genus", "Species")) %>%
  mutate(Species = str_to_title(Species), 
         Genus = str_to_title(Genus), 
         Phylum = na_if(Phylum, "Fungi "),
         `Sub-phylum` = na_if(`Sub-phylum`, "Fungi "), 
         Class = na_if(Class, "Fungi "), 
         Genus = na_if(Class, "Fungi "),
         Species = na_if(Species, "Fungi  Sp.")) %>% 
  column_to_rownames(var = "md5sum") %>%
  mutate(maxrank = if_else(!is.na(Species), "Species", "Kingdom"),
         maxrank = if_else(is.na(Species), "Genus", maxrank), 
         maxrank = if_else(is.na(Genus), "Class", maxrank), 
         maxrank = if_else(is.na(Class), "Subphylum", maxrank), 
         maxrank = if_else(is.na(`Sub-phylum`), "Phylum", maxrank), 
         maxrank = if_else(is.na(Phylum), "Kingdom", maxrank),
         maxrank = factor(maxrank, levels = c("Kingdom", "Phylum", "Subphylum", "Class", "Genus", "Species"))) %>%
  # apply(., 2, str_trim) %>% 
  as.matrix()


# OTU TABLE
otu <- inputdf  %>%
  filter(grepl("Fungi", lineage, ignore.case = TRUE)) %>%
  select(md5sum, starts_with("Tara")) %>%
  column_to_rownames(var = "md5sum") %>%
  t() %>%
  as.matrix()


# METADATA 
meta <- metadfraw %>%
  filter(Sample.ID..TARA_barcode.. %in% rownames(otu)) %>%
  column_to_rownames(var = "Sample.ID..TARA_barcode..")



# MAKE PHYLOSEQ ###
OTU <- otu_table(otu, taxa_are_rows = FALSE)
TAX <- tax_table(tax)
MET <- sample_data(meta)



physeq <- phyloseq(OTU, TAX, MET)
saveRDS(physeq, "./data/rds/physeqTARAFungi.rds")




# GROUP-LEVEL PHYLOSEQ ####
ALLdf <- inputdf %>% select(md5sum, lineage) %>%
  filter(!grepl("Bacteria", lineage, ignore.case = TRUE), 
         !grepl("archae", lineage, ignore.case = TRUE)) %>%
  mutate(Superclade = str_extract(lineage, "\\|(.+?)\\|", group = 1), 
         Superclade = if_else(is.na(Superclade), "Unkown Eukaryote", Superclade), 
         Superclade = if_else(Superclade %in% c("Eukaryota_X", "Eukaryota_Mikro"),
                              "Unkown Eukaryote", Superclade), 
         Phylum = str_extract(lineage, "^.+?\\|.+?\\|(.+?)\\|", group = 1), 
         Phylum = if_else(grepl("X", Phylum), paste("Unkown", Superclade), Phylum), 
         Phylum = if_else(is.na(Phylum), paste("Unkown", Superclade), Phylum), 
         Group = if_else(Superclade == "Opisthokonta" & Phylum == "Fungi", "Fungi", Superclade), 
         Group = if_else(Superclade == "Opisthokonta" & Phylum == "Choanoflagellida", "Choanoflagellida", Group), 
         Group = if_else(Superclade == "Opisthokonta" & Phylum == "Metazoa", "Metazoa", Group), 
         Group = if_else(Superclade == "Opisthokonta" & Phylum == "Mesomycetozoa", "Mesomycetozoa", Group), 
         Group = if_else(Superclade == "Opisthokonta" & Phylum == "Opisthokonta_X", "Unkown Opisthokonta", Group), 
         Group = if_else(Superclade == "Opisthokonta" & is.na(Phylum),"Unkown Opisthokonta", Group))

unique(ALLdf$Superclade)
unique(ALLdf$Phylum)
unique(ALLdf$Group)


# TAXONOMY TABLE
taxALL <- ALLdf %>% 
  column_to_rownames(var = "md5sum") %>%
  as.matrix()


# OTU TABLE
otuALL <- inputdf  %>%
  filter(md5sum %in% rownames(taxALL)) %>%
  select(md5sum, starts_with("Tara")) %>%
  column_to_rownames(var = "md5sum") %>%
  t() %>%
  as.matrix()


# METADATA 
metaALL <- metadfraw %>%
  filter(Sample.ID..TARA_barcode.. %in% rownames(otuALL)) %>%
  column_to_rownames(var = "Sample.ID..TARA_barcode..")



# MAKE PHYLOSEQ ###
OTUALL <- otu_table(otuALL, taxa_are_rows = FALSE)
TAXALL <- tax_table(taxALL)
METALL <- sample_data(metaALL)



physeq <- phyloseq(OTUALL, TAXALL, METALL)
saveRDS(physeq, "./data/rds/physeqTARAALL.rds")


# CAZymes ####
input <- read.csv("./data/tabular/cazymes_abundancetable.csv") %>% 
  filter(Type == "MetaT") 

plant <- "GH1 GH2 GH3 GH4 GH5 GH8 GH9 GH11 GH12 GH15 GH16 GH17 GH26 GH27 GH28 GH29 GH36 GH39 GH43 GH44 GH48 GH51 GH53 GH55 GH67 GH74 GH78 GH93 GH94 GH95 GH115 GH117 GH121 PL1 PL2 PL6 PL7 PL9 PL11 PL15 PL22" %>% 
  strsplit(., " +")  %>% unlist()
animal <- "GH1 GH2 GH3 GH4 GH18 GH19 GH20 GH29 GH33 GH38 GH58 GH79 GH84 GH85 GH88 GH89 GH92 GH95 GH98 GH99 GH101 GH105 GH109 GH110 GH113 PL6 PL8 PL12 PL13 PL21" %>% 
  strsplit(., " +") %>% unlist()
pept <- "GH23 GH24 GH25 GH73 GH102 GH103 GH104 GH108" %>% 
  strsplit(., " +") %>% unlist()
fung <- "GH5 GH8 GH16 GH18 GH19 GH20 GH55 GH64 GH71 GH81" %>% 
  strsplit(., " +") %>% unlist()

cols <- c("Plant Cell Wall Carbohydrates", "Animal Carbohydrates", 
          "Peptidoglycan", "Fungal Carbohydrates")
rows <- unique(c(animal, plant, pept, fung))
taxkey <- matrix(nrow = length(rows), 
              ncol = 4,
              dimnames = list(rows, cols)) %>% 
  as.data.frame()

tax <- input %>% select(SeqID, cazy) %>% distinct() %>% 
  mutate(CAZYgroup = str_extract(cazy, "^(.+?)[_+]", group = 1), 
         CAZYgroup = if_else(is.na(CAZYgroup), cazy, CAZYgroup), 
         CAZYgroup = str_trim(CAZYgroup)) %>% 
  rename(CAZY = cazy) %>% 
  mutate(animal_carbs = if_else(CAZYgroup %in% animal, 1, 0), 
         plant_cell_wall = if_else(CAZYgroup %in% plant, 1, 0), 
         peptidoglycan = if_else(CAZYgroup %in% pept, 1, 0), 
         fungal_carbs = if_else(CAZYgroup %in% fung, 1, 0)) %>% 
  mutate(Target = if_else(animal_carbs == 1, "Animal Carbohydrates", NA), 
         Target = if_else(plant_cell_wall == 1, "Plant Cell Wall", Target), 
         Target = if_else(peptidoglycan == 1, "Peptidoglycan", Target), 
         Target = if_else(fungal_carbs == 1, "Fungal Carbohydrates", Target),
         Target = if_else(rowSums(.[4:7])>1, "Multitargeted", Target), 
         Target = if_else(rowSums(.[4:7])==0, "Other", Target), 
         ) %>% 
  distinct(SeqID, .keep_all = TRUE) %>% 
  column_to_rownames(var = "SeqID") %>% 
  as.matrix()



otu <- input %>% 
  select(Sample_Code, SeqID, Occurrence) %>% 
  distinct() %>% 
  # mutate(Occurrence = as.numeric(Occurrence)) %>% 
  pivot_wider(names_from = SeqID, 
              values_from = Occurrence, 
              values_fill = 0) %>% 
  column_to_rownames(var = "Sample_Code") %>% 
  as.matrix()
meta <- input %>% select(-c(SeqID, Occurrence, cazy)) %>% 
  distinct() %>% 
  column_to_rownames(var = "Sample_Code")
  

OTU <- otu_table(otu, taxa_are_rows = FALSE)
TAX <- tax_table(tax)
MET <- sample_data(meta)


physeq <- phyloseq(OTU, TAX, MET)
saveRDS(physeq, "./data/rds/physeqTARACazymes.rds")


# Peptidase ####
input <- read.csv("./data/tabular/protease_abundancetable.csv") %>% 
  filter(type == "MetaT") 

otu <- input %>% 
  select(Sample_Code, SeqID, Occurence) %>% 
  distinct() %>% 
  # mutate(Occurrence = as.numeric(Occurrence)) %>% 
  pivot_wider(names_from = SeqID, 
              values_from = Occurence, 
              values_fill = 0) %>% 
  column_to_rownames(var = "Sample_Code") %>% 
  as.matrix()
meta <- input %>% select(Sample_Code, type, Size, Filter, 
                         Stations, Depth) %>% 
  distinct() %>% 
  column_to_rownames(var = "Sample_Code") 
tax <- input %>% select(-c(Sample_Code, type, Size, Filter, 
                           Stations, Depth, Occurence)) %>% 
  distinct() %>% 
  column_to_rownames(var = "SeqID") %>% 
  mutate(class = str_replace(class, "Uncl.+$", "Unknown")) %>%
  as.matrix()



OTU <- otu_table(otu, taxa_are_rows = FALSE)
TAX <- tax_table(tax)
MET <- sample_data(meta)



physeq <- phyloseq(OTU, TAX, MET)
saveRDS(physeq, "./data/rds/physeqTARAPeptidase.rds")

