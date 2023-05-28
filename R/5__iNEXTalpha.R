library(tidyverse)
library(phyloseq)
library(iNEXT)
library(ggpubr)
library(corrplot)
library(ggsci)


# INPUT PARAMETERS
fileoptions <- list.files("./data/rds/")
fileoptions
input <- "physeqTARAFungi.rds"

physeq <- readRDS(paste0("./data/rds/", input))
meta <- data.frame(sample_data(physeq))
tax <- data.frame(tax_table(physeq))
features <- data.frame(features = colnames(meta))
ranks <- data.frame(ranks = rank_names(physeq))


SE <- FALSE # Whether to calculate the standard error (significantly increases run time)
group <- "Env.feature..abbreviation."
iterations <- 20
seednumber <- 123
coveragethresh <- 0.90
confidence <- 0.95
knots <- 20 # number of intervals between min and maximum seq samples
endpoint <- NULL # if null it will automatically do the maximum 


# 1) LOAD DATA ####
physeqGROUP <- merge_samples(physeq, group)

otu <- data.frame(otu_table(physeq))
otu <- otu[rowSums(otu[])>0,] %>% 
  round() # remove samples with no reads & round of decimals 
minotu <- min(rowSums(otu))
otuGROUP <- data.frame(otu_table(physeqGROUP)) 
otuGROUP <- otuGROUP[rowSums(otuGROUP[])>0,] %>% 
  round()
minotuGROUP <- min(rowSums(otuGROUP))

tax <- data.frame(tax_table(physeq))

sampledata <- data.frame(sample_data(physeq)) %>%
  as.data.frame() %>%
  mutate(SampleID = rownames(.), 
         Group = .data[[group]])  %>% 
  mutate(Group = str_replace(Group, " Dataset", ""))
sample_data(physeq) <- sampledata



# 2) RUN COMPUTATIONALLY HEAVY ANALYSIS & SAVE ####
set.seed(seednumber)
inextSAMP <- iNEXT(t(otu),
                   q=c(0,1,2),# hill number diversity order
                   datatype="abundance",
                   nboot=iterations,
                   se = SE, # Setting to TRUE increases computation time significantly
                   conf = confidence)
saveRDS(inextSAMP, paste0("./data/rds/rarcurveiNEXT_sample_", input))


set.seed(seednumber)
alphadiv <- estimateD(t(otu), 
                      datatype="abundance",
                      base="size", level=NULL, 
                      nboot=iterations, 
                      conf = confidence)

saveRDS(alphadiv, paste0("./data/rds/rarcurveiNEXT_estimateD_", input))

set.seed(seednumber)
inextGROUP <- iNEXT(t(otuGROUP),
                    q=c(0,1,2),# hill number diversity order
                    datatype="abundance",
                    nboot=iterations,
                    endpoint = NULL, 
                    se = SE,
                    knots = knots, 
                    conf = confidence)

saveRDS(inextGROUP, paste0("./data/rds/rarcurveiNEXT_group", group, "_", input))



# 3) PLOT SAMPLE BASED RAREFACTION ####
inextSAMP <- readRDS(paste0("./data/rds/rarcurveiNEXT_sample_", input)) # LOAD if necessary

plotdfSAMP <- inextSAMP$iNextEst$size_based %>% mutate(m = m/10000, 
                                                       qD = qD)  %>%
  filter(Order.q == 0) %>% 
  mutate(Group = sampledata[Assemblage, group]) %>%
  pivot_longer(cols = c(qD, SC))

observedSAMP <- plotdfSAMP %>% filter(Method == "Observed")
extrapolateSAMP <- plotdfSAMP %>% filter(Method == "Extrapolation")


plotdfSAMP %>%
  filter(Method == "Rarefaction")  %>% 
  ggplot(aes(x = m, y = value, group = Assemblage)) + 
  geom_line(color = "grey", size = 0.10, linetype="solid") + 
  geom_line(data = extrapolateSAMP, aes(x = m, y=value, group=Assemblage, color=Group), 
            size = 0.50, linetype="solid") + 
  geom_point(data = observedSAMP, aes(x = m, y = value), size =1, 
             color = "blue") + 
  facet_wrap(~name, scales = "free_y", ncol = 1) + 
  theme_bw() + 
  theme(axis.text =  element_text(size = 10), 
        axis.title.x = element_text(size = 12, vjust = -2), 
        axis.title.y = element_text(size = 12, vjust= 4), 
        plot.title = element_text(size = 12, face = "italic"),
        legend.text =   element_text(size = 9), 
        legend.background = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.border =element_rect(colour="black", size = 1.3), 
        panel.grid = element_line(color = "grey95"),        
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"), 
        legend.position = "right",
        legend.title = element_blank()) +
  expand_limits(x=c(0,22)) + 
  scale_color_brewer(palette = "Set1") + 
  xlab("Sequence Sampled (x10,000)") + 
  ylab("Richness (x1000)") + 
  ggtitle("Prokaryotic 16S Diversity")




# 4) PLOT POINT ESTIMATE ####
alphadiv <- readRDS(paste0("./data/rds/rarcurveiNEXT_estimateD_", input))

alphadiv %>%
  filter(Order.q == 0 ) %>% 
  ggplot(aes(x = SC)) +
  geom_histogram(fill = "blue", 
                 binwidth = 0.01) + theme_cleveland()


alphadivPLOT <- alphadiv %>%
  filter(SC >= coveragethresh) %>%
  merge(sampledata, by.x = "Assemblage", by.y = "row.names") 


groups <- unique(sampledata$Group)
my_comparisons <- list()
for (i in 1:(length(groups)-1)){
  comp1 <- groups[i]
  for (j in (i+1):(length(groups))){
    comp2 <- groups[j]
    innerlist <- c(comp1, comp2)
    my_comparisons <- append(my_comparisons, list(innerlist))
  }
}

alphadivPLOT %>%
  mutate(Order.q = str_replace(Order.q, "0", "q = 0"),
         Order.q = str_replace(Order.q, "1", "q = 1"), 
         Order.q = str_replace(Order.q, "2", "q = 2")) %>% 
  ggplot(aes(x = Group, y = qD, fill = Group)) + 
  geom_boxplot(outlier.shape = NA, color = "black",
               size=0.7, width=0.4) + 
  geom_point(aes(fill = Group), size = 0, shape = 22) +
  # geom_jitter() + 
  stat_boxplot(geom ='errorbar', width=0.2) +
  stat_compare_means(aes(label = ..p.signif..),
                     comparisons = my_comparisons,
                     size = 5,
                     method = "wilcox.test",
                     color = "red3",
                     vjust = 2, bracket.size = 0.5) +
  facet_wrap(~Order.q, scales = "free_y") + 
  xlab("Environment") + 
  ylab("Alpha Diversity (X100)") + 
  theme(strip.background=element_blank(), 
        strip.text.x=element_text(face = "bold.italic", size = 8, vjust = -1), 
        panel.border = element_rect(color="black", fill=NA, linewidth=1),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title.x = element_blank(), 
        legend.position = "bottom", 
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.title = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.y = element_text(size = 9, face = "italic", vjust = 5), 
        axis.text.y = element_text(size = 9, face = "italic"), 
        plot.margin=unit(rep(0.5, 4), "cm")) + 
  guides(fill = guide_legend(override.aes = list(size = 4), 
                             fill=NA)) + 
  scale_fill_jco() + 
  scale_color_jco()


ggsave(paste0("./data/figures/png/alphadiv.png"),
       height = 6.8,
       width = 14,
       units = "cm",
       dpi = 200)


ggsave(paste0("./data/figures/svg/alphadiv.svg"),
       height = 6.8,
       width = 14,
       units = "cm",
       dpi = 1000)




# 5) PLOT GROUP BASED RAREFACTION ####
inextGROUP <- readRDS(paste0("./data/rds/rarcurveiNEXT_group", group, "_", input)) # LOAD if necessary
plotdfGROUPED <- inextGROUP$iNextEst$size_based %>% mutate(m = m/10000, 
                                                           qD = qD/100)  %>%
  filter(Order.q == 0)
observedGROUPED <- plotdfGROUPED %>% filter(Method == "Observed")
extrapolateGROUPED <- plotdfGROUPED %>% filter(Method == "Extrapolation")


plotdfGROUPED %>%
  filter(Method == "Rarefaction") %>%
  ggplot(aes(x = m, y = qD, group = Assemblage, color = Assemblage)) + 
  geom_line(size = 1.5, linetype="solid") + 
  geom_line(data = extrapolateGROUPED, aes(x = m, y=qD, group=Assemblage),
            size = 1.5, linetype="dashed") +
  geom_point(data = observedGROUPED, aes(x = m, y = qD), size =1, 
             color = "blue") + 
  geom_text(data = observedGROUPED, aes(x = m, y = qD, label = Assemblage, color = Assemblage), 
            size = 3, hjust=-0.1, vjust=-2.5, fontface="bold") + 
  theme(plot.background = element_blank(),
        strip.background=element_blank(), 
        panel.border = element_rect(color="black", fill=NA, linewidth=1.2), 
        axis.title.x = element_text(size = 8, face = "italic", vjust = -3), 
        panel.grid = element_line(color = "grey95"),
        legend.position = "none", 
        legend.key.size = unit(0.4, 'cm'),
        legend.title = element_blank(),
        plot.title = element_text(size = 9, face = "bold.italic", vjust = -1), 
        axis.text.x = element_text(size = 8, face = "italic"), 
        axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 8, face = "italic"), 
        plot.margin=unit(rep(0.5, 4), "cm")) + 
  scale_color_jco() + 
  ggtitle("Fungal Richness X100 (q=0)") + 
  geom_vline(xintercept = min(observedGROUPED$m), linetype = "dashed", 
             color = "red", linewidth = 0.7, alpha = 0.5) + 
  xlab("Sequences Sampled (x10,000)") + 
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15)), 
                     position="left") + 
  ylab("Richness (q=1)") 



ggsave(paste0("./data/figures/png/alphadivGROUP.png"),
       height = 5.5,
       width = 5,
       units = "cm",
       dpi = 200)

ggsave(paste0("./data/figures/svg/alphadivGROUP.svg"),
       height = 5.5,
       width = 5,
       units = "cm",
       dpi = 1000)



# 6) PLOT METADATA CORRELATION ####
corrdata <- alphadivPLOT %>% 
  filter(SC > coveragethresh) %>% 
  select(c(qD, colnames(sampledata))) %>% 
  select(where(is.numeric)) 

COR <- cor(corrdata, use="pairwise.complete.obs")

p <- (cor.mtest(corrdata))$p
pfilt <- p[COR>0.5]

corrplot(COR, 
         type = "upper", 
         method = 'circle', 
         p.mat = p, 
         tl.cex = 0.4,
         tl.col = "grey50", 
         number.cex = 0.3,
         insig = "blank") # colorful number
