library(tidyverse)
library(phyloseq)


# LOAD ####
physeq <- readRDS("./data/rds/physeqTARAFungi.rds")
metadata <- data.frame(sample_data(physeq))
features <- data.frame(features = colnames(metadata))

# DEPTH DISTRIBUTION
hist(metadata$`PAR..mol.quanta.m..2.day...for.a.period.of.8.days.around....`)


unique(metadata$OS.region...abbreviation..full.name..MRG....)

# TEMP 
