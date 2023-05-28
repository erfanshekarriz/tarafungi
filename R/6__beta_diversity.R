library(tidyverse)
library(vegan)
library(microbiome)
library(RColorBrewer)
library(ggsci)



# INPUT PARAMETERS
fileoptions <- list.files("./data/rds/")
fileoptions

input <- "physeqTARAFungi.rds"
physeq <- readRDS(paste0("./data/rds/", input))
meta <- data.frame(sample_data(physeq))
features <- data.frame(features = colnames(meta))

group <- "OS.region...abbreviation..full.name..MRG...."
colorgroup <- "Env.feature..abbreviation."
iterations <- 10
anosimperm <- 99
seednumber <- 123
distance <- "bray"

# UPLOAD DATA ####
physeq <- readRDS(paste0("./data/rds/", input))
physeqrel <- transform_sample_counts(physeq, function(x) x / sum(x) )

otudf <- data.frame(otu_table(physeq))

mindec <- min(as.matrix(otudf)[as.matrix(otudf) > 0])
decimals <- as.character(mindec) %>% str_extract("\\-(.+$)", group=1) %>% as.numeric() + 2
if (is.na(decimals)){otuINTEGER <- otudf} else {otuINTEGER <- ((otudf * 10^decimals) %>% ceiling())}


otudfrel <- data.frame(otu_table(physeqrel))
metadf <- data.frame(sample_data(physeq)) %>% mutate(Group = .data[[group]], 
                                                     ColorGroup = .data[[colorgroup]]) 
physeqclr <- microbiome::transform(physeq, "clr")
mindepth <- otuINTEGER %>% rowSums() %>% min() 



# SETUP THEME ####
customtheme <-
  theme_bw() + 
  theme(axis.title.x=element_text(face="italic", size=9, vjust=-3), 
        axis.title.y=element_text(face="italic", size=9, vjust=4), 
        axis.text=element_text(size = 6, face ="bold"), 
        panel.border = element_rect(linewidth = 0),
        legend.text=element_text(size = 8, face="bold"), 
        legend.background=element_rect(fill= ggplot2::alpha('grey', 0.2), 
                                       color = "grey90", linetype = 2, linewidth = 0),
        legend.key = element_rect(fill = "transparent"), 
        legend.key.size = unit(0.2, "cm"),
        legend.title=element_blank(),
        legend.position= c(0.5, 0.2),
        plot.margin=unit(c(0.5,1,1,1), "cm"),
        plot.background = element_blank(),
        panel.grid.minor=element_line(color="grey92"), 
        panel.grid.major=element_line(color="grey92"))


# FILTER LOW DEPTH SAMPLES (OPTIONAL) ####


# PERMANOVA ANOSIM ####
ano <- anosim(otudfrel, metadf$Group, distance = distance, permutations = anosimperm)
summary(ano)

# RAREFY ####
set.seed(seednumber)
rardist <- avgdist(otuINTEGER, dmethod = distance, iterations = iterations, 
                   sample = mindepth)


# RAREIFIED MDS ORDINATION ####
rarMDS <- cmdscale(rardist, 
                   k = 100,
                   eig = TRUE, 
                   add = TRUE) 

eigMDS <- rarMDS$eig %>% 
  data.frame(PC = 1:length(rarMDS$eig),
             Variance = .*100/sum(rarMDS$eig))
eigMDS %>% 
  ggplot(aes(x = PC, y = Variance)) + 
  geom_point(shape = 22, color = "blue", size = 4) + 
  geom_line(linetype = "dashed") + 
  theme_minimal()

kcutoff <- 4

rarMDS <- cmdscale(rardist, 
                   k = kcutoff,
                   eig = TRUE) 


rarMDSplot <- rarMDS$points %>%
  as.data.frame() %>% 
  merge(metadf, by = 'row.names') %>%
  column_to_rownames(var = "Row.names") 
colnames(rarMDSplot) <- str_replace(colnames(rarMDSplot), "^V(\\d)\\d?$", "PCA\\1")
  

  
corrMDSrar <- rarMDSplot %>% 
  select(where(is.numeric)) %>% 
  cor(method = c("spearman"), use="complete.obs") %>% 
  as.data.frame() %>% 
  select(starts_with("PCA", ignore.case = FALSE)) %>% 
  as.matrix() 

heatmap(x = corrMDSrar, symm = FALSE, col = terrain.colors(256), 
        margins =c(5,20),
        cexRow=0.5,
        cexCol=0.5)

getPalette <-  colorRampPalette((pal_frontiers("default")(10)))
cols <- getPalette(length(unique(rarMDSplot$Group)))


rarMDSplot %>%
  ggplot(aes(x=PCA1, y=PCA2)) + 
  geom_jitter(aes(fill = Group), color = "black", size=2.8, shape = 21,
              alpha = 0.6, width = 0) + 
  theme_bw() + 
  xlab(paste0("PCA1 (", round(rarMDS$eig[1],2), "%)")) + 
  ylab(paste0("PCA2 (", round(rarMDS$eig[2],2), "%)")) + 
  theme(axis.title.x=element_text(face="italic", size=9, vjust=-3), 
        axis.title.y=element_text(face="italic", size=9, vjust=4), 
        axis.text=element_text(size = 6, face ="bold"), 
        panel.border = element_rect(linewidth = 0),
        legend.text=element_text(size = 8, face="bold"), 
        legend.background=element_rect(fill= ggplot2::alpha('grey', 0.2), 
                                       color = "grey90", linetype = 2, linewidth = 0),
        legend.key = element_rect(fill = "transparent"), 
        legend.key.size = unit(0.2, "cm"),
        legend.title=element_blank(),
        legend.position= "bottom",
        plot.margin=unit(c(0.5,1,1,1), "cm"),
        plot.background = element_blank(),
        panel.grid.minor=element_line(color="grey92"), 
        panel.grid.major=element_line(color="grey92"))  + 
  geom_text(label = paste0("ANOSIM\n", 
                           "R: ", round(ano$statistic, 2), 
                           "\np-val: ", ano$signif), 
            # hjust = -1,
            # vjust = 2,
            x = min(rarMDSplot$PCA1) + diff(range(rarMDSplot$PCA1))*0.07,  
            y = min(rarMDSplot$PCA2)+ diff(range(rarMDSplot$PCA2))*0.93,
            size = 3, check_overlap = T, 
            fontface = "bold.italic", color = "red2") + 
  scale_fill_manual(values = cols) + 
  scale_color_manual(values = cols) + 
  guides(fill = guide_legend(ncol = 2, 
                             title = "Sampling Group"),
         shape = "none")  

# RAREIFIED NMDS ORDINATION ####
set.seed(seednumber)
rarNMDS <- metaMDS(rardist , k = 2, 
                   distance = "euclidean", 
                   engine="isoMDS") 

rarNMDSplot <- rarNMDS$points %>%
  as.data.frame() %>%
  rename(NMDS1=V1, NMDS2=V2) %>%
  merge(metadf, by = 'row.names') %>%
  column_to_rownames(var = "Row.names")

getPalette <-  colorRampPalette((pal_frontiers("default")(10)))
cols <- getPalette(length(unique(rarNMDSplot$ColorGroup)))



rarNMDSplot %>%
  ggplot(aes(x=NMDS1, y=NMDS2)) + 
  geom_jitter(aes(fill = ColorGroup), color = "black", size=2.8, shape = 21,
              alpha = 0.6, width = 0) + 
  customtheme + 
  geom_text(label = paste0("ANOSIM\n", 
                           "R: ", round(ano$statistic, 2), 
                           "\np-val: ", ano$signif), 
            hjust = 0,
            vjust = 1,
            x = min(rarNMDSplot$NMDS1) + diff(range(rarNMDSplot$NMDS1))*0.02,  
            y = min(rarNMDSplot$NMDS2)+ diff(range(rarNMDSplot$NMDS2))*0.98,
            size = 3, check_overlap = T, 
            fontface = "bold.italic", color = "red2") + 
  scale_fill_manual(values = cols) + 
  scale_color_manual(values = cols) + 
  guides(fill = guide_legend(ncol = 2, 
                             title = "Sampling Group"),
         shape = "none")  


# # NON-RARIFIED AITCHISON'S DISTANCE ####
pca  <- ordinate(physeqclr, method = "PCoA", distance = "euclidean")

plot_ordination(physeqclr, pca, type="samples", color=group) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")



# NON-RARIFIED NMDS ####
set.seed(12)

NMDS  <- ordinate(physeqclr,
                  method = "NMDS",
                  distance = "euclidean",
                  engine="isoMDS")

plot_ordination(physeqclr, NMDS, type="samples", color="biome") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

plot_ordination(physeqclr, NMDS, type="samples", color="Experiment.Region",
                shape = "biome") +
  theme_grey()


# RAREIFIED MDS ORDINATION ####
rarordplot <- cmdscale(rardist , k = 2) %>%
  as.data.frame() %>%
  rename(PC1=V1, PC2=V2) %>%
  merge(metadf, by = 'row.names') %>%
  column_to_rownames(var = "Row.names")

rarordplot %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color = Group)) +
  theme_minimal() +
  scale_color_brewer(palette = "Set2")





