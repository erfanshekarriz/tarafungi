library(tidyverse)
library(phyloseq)


# LOAD ####
physeq <- readRDS("./data/rds/physeqTARAFungi.rds")
taxdf <- data.frame(tax_table(physeq)) %>% 
  mutate(maxrank = factor(maxrank, levels = c("Kingdom", "Phylum", "Subphylum", "Class", "Genus", "Species")))



# ASSIGN MAXIMUM RANK ####
taxdf %>%
  filter(maxrank != "Kingdom", 
         !is.na(Phylum)) %>%
  ggplot(aes(x = maxrank)) + 
  geom_bar(aes(fill = Phylum), position = "stack", color = NA) + 
  theme_minimal() + 
  theme(legend.position = "none", 
        legend.text = element_text(size = 5), 
        # axis.title.x = element_text(size = 5), 
        axis.title = element_blank(), 
        axis.text = element_text(size = 6),
        axis.text.x = element_text(angle = 55, hjust=1)) + 
  xlab("Maximum Classification Resolution") + 
  ylab("Count") + 
  scale_fill_brewer(palette = "Set1")

