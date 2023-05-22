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


