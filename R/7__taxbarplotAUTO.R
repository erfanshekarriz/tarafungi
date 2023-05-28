library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(gridExtra)
library(cowplot) 


# INPUT PARAMETERS
fileoptions <- list.files("./data/rds/")
fileoptions
input <- "physeqTARACazymes.rds"

physeq <- readRDS(paste0("./data/rds/", input))
meta <- data.frame(sample_data(physeq))
tax <- data.frame(tax_table(physeq))
features <- data.frame(features = colnames(meta))
ranks <- data.frame(ranks = rank_names(physeq))

# SETUP FUNCTION ####
gettaxbarplot <- function(physeq, taxlevel, group, otherthresh){
  
  physeqGROUP <- merge_samples(physeq, group)
  
  otu <- data.frame(otu_table(physeq)) 
  colnames(otu) <- gsub("^X", "", colnames(otu))
  
  tax <- data.frame(tax_table(physeq)) %>%
    mutate(taxlevel = .data[[taxlevel]], 
           taxlevel = if_else(is.na(taxlevel), "Unknown", taxlevel),
           taxlevel = str_replace_all(taxlevel, "unknown", "Unknown"))
  rownames(tax) <- str_replace_all(rownames(tax), "-", "\\.")
  
  metadata <- data.frame(sample_data(physeq)) %>% 
    mutate(Group = .data[[group]])
  sample_data(physeq) <- metadata
  
  unique(tax$taxlevel)
  
  otutaxglom <- otu %>% 
    apply(1, function(x) x/sum(x)) %>%
    merge(tax, by = "row.names", all.x = TRUE) %>%
    mutate_at(rownames(otu), ~replace(., is.nan(.), 0)) %>% 
    column_to_rownames(var = "Row.names") %>% 
    group_by(taxlevel) %>%
    summarise(across(where(is.numeric), sum)) %>% 
    column_to_rownames(var = "taxlevel") %>%
    replace(is.na(.), 0) %>% 
    mutate(avgAbund = rowMeans(.)) %>%
    rownames_to_column(var = "taxlevel") %>%
    mutate(taxlevel = if_else(avgAbund > otherthresh, taxlevel, "Other")) %>%
    group_by(taxlevel) %>%
    summarize(across(where(is.numeric), sum)) %>%
    column_to_rownames(var = "taxlevel")
  
  colSums(otutaxglom) # check if columns add to 1 
  
  
  plotdf <- otutaxglom %>%
    t() %>% as.data.frame() %>%
    merge(metadata, by = "row.names") %>%
    pivot_longer(rownames(otutaxglom), 
                 names_to = "taxlevel", 
                 values_to = "relAbundance") 
  return(plotdf)
}
# PLOT CLASS ####
plotdf <- gettaxbarplot(readRDS(paste0("./data/rds/", input)), 
                        "peptidase.family", 
                        "Depth", 
                        0.03)


colourCount <- length(unique(plotdf$taxlevel))
getPalette <-  colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(colourCount)
names(cols) <- unique(plotdf$taxlevel)

plotdf %>%
  filter(Group %in% c("DCM", "SRF")) %>% 
  ggplot(aes(x=as.factor(Row.names), y=as.numeric(relAbundance), 
             fill=taxlevel, 
             group=taxlevel)) + 
  geom_bar(position="stack", stat="identity", width = 1) + 
  facet_grid(.~Group, scales = "free_x", space = "free") +
  scale_y_continuous(expand = c(0, 0), position = "left") + 
  theme_bw() + 
  theme(strip.background=element_blank(), 
        strip.text.x=element_text(face = "bold.italic", size = 7, vjust =-2), 
        strip.text.y=element_text(face = "bold", size = 7), 
        axis.title.x=element_blank(), 
        axis.ticks = element_line(linewidth = 0.3), 
        axis.title.y=element_text(size = 8, face ="italic", vjust = 4), 
        axis.text.y=element_text(size = 6, face ="italic"), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        legend.text=element_text(size = 6), 
        legend.background=element_blank(),
        # legend.title=element_text(face="bold", size=6, angle=90, vjust = 6, 
        #                           hjust = 0.2),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.2, 'cm'),
        panel.border=element_rect(colour="grey", size = 0), 
        panel.background=element_blank(), 
        plot.background=element_blank(), 
        plot.margin=unit(rep(0.5, 4), "cm"),
        panel.grid=element_blank(), 
        panel.spacing.y = unit(0.20, "lines"), 
        panel.spacing.x = unit(0.05, "lines"))  +
  ylab("Relative Abundance") +
  scale_fill_manual(values = cols) +
  guides(fill=guide_legend(title="Eukaryotic\nRelative\nAbundance", 
                           nrow = 4)) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = scales::percent)

