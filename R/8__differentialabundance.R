library(tidyverse)
library(phyloseq)
library("MicrobiomeStat")
library(ggrepel)
library(RColorBrewer)
library(ggsci)

# INPUT PARAMETERS ####
fileoptions <- list.files("./data/rds/")
fileoptions
input <- "physeqTARAFungi.rds"

# Check your variables
physeq <- readRDS(paste0("./data/rds/", input))
meta <- data.frame(sample_data(physeq))
tax <- data.frame(tax_table(physeq))
features <- data.frame(features = colnames(meta))
ranks <- data.frame(ranks = rank_names(physeq))

# SET PARAMETERS
group <- "Env.feature..abbreviation."
variable1 <- "SRF"
variable2 <- "DCM"
covariates <-  c("OS.region...abbreviation..full.name..MRG....")

prevelancefilt <- 0.1
meanabundfilt <- 0 
zerohandling <- 'pseudo-count'
cores <- 4
alphapval <- 0.005
labelRank <- "Genus"
colorRank <- "Phylum"
lcfthresh <- 1
padjustcut <- 0.05


# UPLOAD DATA ####
physeq <- readRDS(paste0("./data/rds/", input))
physeqrel <- transform_sample_counts(physeq, function(x) x / sum(x) )

otudfrel <- data.frame(otu_table(physeqrel) %>% t())
metadf <- data.frame(sample_data(physeq)) %>% mutate(Group = .data[[group]]) 



# PERFORM LINDA ####
formulastr <- paste(paste0("~", group), paste(covariates, collapse = " + "), sep = " + ")
linda <- linda(feature.dat = otudfrel,
               meta.dat = metadf, 
               formula = formulastr,
               alpha = alphapval,
               feature.dat.type = "proportion",
               prev.filter = prevelancefilt, 
               mean.abund.filter = 0, 
               max.abund.filter = 0,
               is.winsor = TRUE,
               outlier.pct = 0.03,
               adaptive = TRUE,
               corr.cut = 0.1,
               zero.handling = "zerohandling",
               pseudo.cnt = 0.5,
               p.adj.method = "BH", 
               n.cores = cores)



# EXTRACT SIGNIFICANT RESULTS ####
variable1 <- linda$variables[[1]]
variable2 <- linda$variables[[2]]
variable1Fact <- str_extract(variable1, ".+\\.(.+$)", group = 1)
variable2Fact <- str_extract(variable2, ".+\\.(.+$)", group = 1)

resLinda <- rownames_to_column(linda$output[[variable1]])  %>% 
  dplyr::left_join(., rownames_to_column(tax), by="rowname") %>%
  mutate(Label = .data[[labelRank]]) %>%
  mutate(Group = if_else(log2FoldChange < 0, variable2Fact, variable1Fact))  




# merge results 
sigdf  <- resLinda %>%
  filter(padj<padjustcut) %>%
  mutate(Label=.data[[labelRank]]) %>%
  filter(!(Label %in% c("unknown", "uncultured")))



# PLOT VOLCANO PLOT ####

# This indicates the variable the log fold is relative to 
signifdf  <- resLinda %>%
  filter(padj<padjustcut & (abs(log2FoldChange) > lcfthresh)) %>%
  mutate(Label = .data[[labelRank]], 
         Color = .data[[colorRank]], 
         Color = if_else(is.na(Color), "Unkown", Color)) 

legendord <- signifdf %>%
  group_by(Color) %>%
  tally() %>%
  arrange(-n) %>%
  pull(Color)

signifdf <- signifdf %>% mutate(Color = factor(Color, levels = legendord))

unkowndf  <- resLinda %>%
  filter(padj<padjustcut & (abs(log2FoldChange) > lcfthresh)) %>%
  mutate(Label = .data[[labelRank]], 
         Color = .data[[colorRank]]) %>%
  filter((Label %in% c("unknown", "uncultured")) | (Color %in% c("unknown", "uncultured")))

nonsignifdf  <- resLinda %>%
  filter(!(rowname %in% signifdf$rowname)) %>%
  filter(!(rowname %in% unkowndf$rowname)) %>%
  mutate(Label=.data[[labelRank]], 
         Color = .data[[colorRank]])

set.seed(7)

colourCount <- length(unique(signifdf$Label))
getPalette <-  colorRampPalette(pal_frontiers("default")(9))
colors <- getPalette(colourCount)

nonsignifdf %>%
  ggplot(aes(x=log2FoldChange, y=-log(padj))) +
  annotate("rect", xmin = lcfthresh, xmax = Inf, 
           ymin = -Inf, ymax = Inf, 
           fill = pal_jco("default")(2)[2], alpha = 0.4) + 
  annotate("rect", xmin = -lcfthresh, xmax = -Inf, 
           ymin = -Inf, ymax = Inf, 
           fill = pal_jco("default")(2)[1], alpha = 0.4)  + 
  geom_boxplot(aes(fill = Group), size = 0, alpha = 0) + 
  geom_point(color = "grey90", size = 0.1, shape = 4) + 
  geom_point(data = unkowndf, color = "grey90", size = 0.1) + 
  geom_point(data = signifdf, aes(color = Color), size = 2) + 
  geom_text_repel(data = signifdf, aes(label = Label, color = Color), 
                  fontface = "bold.italic", size = 2) + 
  theme_bw() + 
  xlab("Log-fold Change (Log2)") + 
  ylab("- Log (p-value adjusted)") + 
  theme(strip.background = element_rect(colour="black",
                                        fill = "grey97", 
                                        linewidth = 0),
        strip.text = element_text(face="bold", size=12, color = "grey20"),
        axis.title.x=element_text(face="italic", size=14, vjust=-3), 
        axis.title.y=element_text(face="italic", size=14, vjust=5),
        axis.text.y=element_text(size = 12, face ="bold.italic", color="black"), 
        axis.text.x=element_text(size = 12, face ="bold", color="black"), 
        legend.text=element_text(size = 12, face="italic"), 
        legend.background=element_blank(),
        # legend.title=element_text(vjust = 2, face = "italic", size=8),
        # legend.title=element_blank(),
        legend.key.size = unit(0.3, 'cm'),
        legend.position = "right",
        panel.border=element_rect(colour="black", linewidth = 0), 
        panel.background=element_rect(fill="grey98"), 
        plot.background=element_blank(), 
        panel.grid = element_blank(),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.spacing.y = unit(0.4, "lines"), 
        panel.spacing.x = unit(0.05, "lines")) + 
  geom_vline(xintercept = lcfthresh, linetype = "dashed", 
             linewidth = 0.8, color = "grey") + 
  geom_vline(xintercept = -lcfthresh, linetype = "dashed", 
             linewidth = 0.8, color = "grey") + 
  geom_hline(yintercept = -log(0.05), linetype = "dashed", 
             linewidth = 0.8, color = "grey") + 
  scale_color_manual(values = colors) + 
  scale_fill_jco() + 
  guides(color = guide_legend(override.aes = list(size=2),
                              title = colorRank, 
                              ncol = 1), 
         fill = guide_legend(override.aes = list(alpha=1), 
                             title = group))



ggsave("./data/figures/svg/volcano.svg", 
       height = 9.8, 
       width = 29,
       units = "cm", 
       dpi = 200)


# SPECIES SPECIFIC CONSESUS PLOT ####
nonsignifdf %>%
  ggplot(aes(x=log2FoldChange, y=reorder(Label, +log2FoldChange))) +
  # geom_errorbarh(aes(xmin=log2FoldChange-lfcSE, 
  #                    xmax=log2FoldChange+lfcSE, 
  #                    color = Color), 
  #                linewidth=0.4, height=0.2) +
  annotate("rect", xmin = lcfthresh, xmax = Inf, 
           ymin = -Inf, ymax = Inf, 
           fill = pal_jco("default")(2)[2], alpha = 0.3) + 
  annotate("rect", xmin = -lcfthresh, xmax = -Inf, 
           ymin = -Inf, ymax = Inf, 
           fill = pal_jco("default")(2)[1], alpha = 0.3)  + 
  geom_point(aes(color = Group, 
                 shape = Group), 
             size = 1) +
  theme_bw() + 
  xlab("Log-fold Change (Log2)") + 
  scale_color_manual(values = colors) +
  theme(axis.title.x=element_text(face="bold", size=12, vjust=-3), 
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.text.y=element_text(size = 7, face ="bold.italic", color="black"), 
        axis.text.x=element_text(size = 8, face ="bold", color="black"), 
        legend.text=element_text(size = 7, face="bold"), 
        legend.background=element_blank(),
        legend.title=element_blank(),
        legend.position = "top",
        legend.key.size = unit(0.2, units = "cm"),
        panel.border=element_rect(colour="grey", linewidth = 1), 
        panel.background=element_rect(fill="grey98"), 
        plot.background=element_blank(), 
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.spacing.y = unit(0.4, "lines"), 
        panel.spacing.x = unit(0.05, "lines")) + 
  geom_vline(xintercept = lcfthresh, linetype = "dashed", 
             linewidth = 0.8, color = "grey") + 
  geom_vline(xintercept = -lcfthresh, linetype = "dashed", 
             linewidth = 0.8, color = "grey")  + 
  guides(shape = "none", 
         color = guide_legend(nrow = 2))

ggsave("./data/figures/png/consensus.png",
       width = 11,
       height = 13,
       units = "cm",
       dpi = 300 )

ggsave("./data/figures/svg/consensus.png",
       width = 11,
       height = 13,
       units = "cm",
       dpi = 300 )




