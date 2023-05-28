library(tidyverse)
library(phyloseq)
library(MicrobiomeStat)
library(ggrepel)
library(RColorBrewer)
library(ggsci)
library(metagMisc)
library(ALDEx2)
library(ANCOMBC)
library(Maaslin2)


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
fixed_covariates <-  c("OS.region...abbreviation..full.name..MRG....")
random_covariates <-  c("OS.region...abbreviation..full.name..MRG....")

prevelancefilt <- 0.01
meanabundfilt <- 0 
zerohandling <- 'pseudo-count'
cores <- 4
alphapval <- 0.005
labelRank <- "Genus"
colorRank <- "Phylum"
lcfthresh <- 1
padjustcut <- 0.05
rankglom <- NULL

aldex2.mcsamples <- 1000  # 128 recommened for ttest, 1000 for rigorous effect size calculation


# UPLOAD DATA ####
physeq <- readRDS(paste0("./data/rds/", input))


# Filter phyloseq samples & agglomerate taxa 
metadf <- data.frame(sample_data(physeq)) %>% 
  mutate(Group = .data[[group]]) %>% 
  filter(Group %in% c(variable1, variable2))

physeqfilt <- prune_samples(intersect(sample_names(physeq), rownames(metadf)),
                            physeq)

if (!is.null(rankglom)){
physeqfilt <- tax_glom(physeqfilt, 
                       taxrank=rankglom,
                       NArm=TRUE, 
                       bad_empty=c(NA, "", " ", "\t"))}

physeqfilt <- phyloseq_filter_prevalence(physeqfilt, 
                                         prev.trh = prevelancefilt, 
                                         abund.trh = NULL,
                                         threshold_condition = "OR", 
                                         abund.type = "total")


otufilt <- data.frame(otu_table(physeqfilt) %>% t())

# Get relative abundance physeq
physeqrelfilt <- transform_sample_counts(physeqfilt, function(x) x / sum(x) )
oturelfilt <- data.frame(otu_table(physeqrelfilt) %>% t())


# Rarefy physeq
physeqfiltRar <- rarefy_even_depth(physeqfilt, sample.size = min(sample_sums(physeqfilt)),
                                   rngseed = 123, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
otufiltRar <- data.frame(otu_table(physeqfiltRar) %>% t())

# Get taxa table
tax <- data.frame(tax_table(physeqrelfilt))


# ALDEX2 ####
cat("Running the ALDEX2 Model... \n")
start_time = Sys.time()
aldexCLR <- aldex.clr(
  reads = otufilt,
  conds = metadf[,group], 
  mc.samples = aldex2.mcsamples, 
  denom = "all",
  verbose = FALSE)
aldexTT <- aldex.ttest(
  aldexCLR, 
  paired.test = FALSE, 
  verbose = FALSE)
aldexEffect <- aldex.effect(aldexCLR, CI = TRUE, verbose = FALSE)
aldexOUTPUT <- data.frame(aldexTT, aldexEffect)
end_time = Sys.time()
duration = end_time - start_time
cat("Done! \nALDEX2 runtime ~", round(duration, digits = 2), 
    units(duration), "\n")



# ANCOM-BC2 ####
cat("Running the ANCOMBC2 Model... \n")
start_time = Sys.time()
ancombc2OUTPUT <- ancombc2(
  data = physeqfilt,
  fix_formula =  paste(group, sep = " + ", collapse = " + "), # for fixed effects
  # fix_formula =  paste(group, fixed_covariates, sep = " + ", collapse = " + "), # for fixed effects
  # rand_formula = "(timepoint | subject)", # for random effects 
  p_adj_method = "fdr", 
  prv_cut = 0, # no prev filtering necessary anymore 
  lib_cut = 0, 
  group = group, 
  struc_zero = TRUE, 
  neg_lb = TRUE, 
  iter_control = list(tol = 1e-2, max_iter = 100, verbose = FALSE),
  alpha = 0.05, 
  global = TRUE
)
end_time = Sys.time()
duration = end_time - start_time
cat("Done! \nANCOMBC2 runtime ~", round(duration, digits = 2), 
    units(duration), "\n")
ancombc2RES <- ancombc2OUTPUT$res



# MaAsLin2 Rarefied #### 
cat("Running the MaASLin2 Rarefied Model... \n")
start_time = Sys.time()
maaslin2rarFITDATA <- Maaslin2(
  input_data = otufiltRar,
  input_metadata = metadf,
  fixed_effects = c(group),
  # fixed_effects = c(group, covariates),
  # random_effects = c(random_covariates),
  standardize = TRUE, # Apply z-score so continuous metadata are on the same scale.
  output = "DAA output",
  transform = "AST",
  reference = paste(group, variable1, sep = ","),  
  normalization = "TSS",
  cores = cores,
  min_prevalence = 0 # prev filtering already done
)
end_time = Sys.time()
duration = end_time - start_time
cat("Done! \nMaASLin2rarefied runtime ~", round(duration, digits = 2), 
    units(duration), "\n")




# LINDA ####
cat("Running the LinDA  Model...")
start_time = Sys.time()
lindaRES <- linda(
  otufilt, 
  metadf, 
  formula = paste0("~", group), 
  alpha = 0.05, 
  # lib.cut = 1000, 
  #winsor.quan = 0.97,
  prev.filter = 0) # we already filtered )
end_time = Sys.time()
duration = end_time - start_time
cat("Done! \nLinDA runtime ~", round(duration, digits = 2), 
    units(duration), "\n")




# SUMMARY ####
if (is.null(rankglom)) {rankglomNAME <- "ASV"} else {rankglomNAME <- rankglom}

aldex2CLEAN <- rownames_to_column(aldexOUTPUT, var = rankglomNAME) %>% 
  dplyr::select(rankglomNAME, aldex2 = wi.eBH)
ancombc2CLEAN <- ancombc2RES %>% 
  dplyr::select(!!rankglomNAME := taxon, ancombc2 = starts_with(paste0("diff_", group)))
maaslin2rarCLEAN <- dplyr::select(maaslin2rarFITDATA$results, !!rankglomNAME := feature, maaslin2 = qval) %>%
  dplyr::mutate(!!rankglomNAME := str_remove(.data[[rankglomNAME]] , "X"))
lindaCLEAN <-rownames_to_column(lindaRES$output[[1]], var = rankglomNAME) %>% dplyr::select(rankglomNAME, linda = reject)


summaryOUTPUT <- full_join(aldex2CLEAN, ancombc2CLEAN) %>% 
  full_join(maaslin2rarCLEAN) %>%
  full_join(lindaCLEAN) %>%
  mutate(across(c(aldex2, maaslin2), ~ .x <= 0.05), 
         score = rowSums(across(c(aldex2, ancombc2, maaslin2, linda)))
  ) %>% 
  merge(tax, by.x = rankglomNAME, by.y = "row.names")

# How many genera were identified by each
summarise(summaryOUTPUT, across(where(is.logical), sum))
# How many genera were identified by all four
filter(summaryOUTPUT, score == 4)



# FILTER THE PREVIOUS RESULTS WITH HIGHLY CONFIDENT OUTPUTS

# VISUALIZATION ####










# EXTRA ####

#PLOT ALDEX2
# par(mfrow = c(1, 2))
# aldex.plot(
#   aldexOUTPUT,
#   type = "MA",
#   test = "welch",
#   xlab = "Log-ratio abundance",
#   ylab = "Difference",
#   cutoff = 0.05
# )
# aldex.plot(
#   aldexOUTPUT,
#   type = "MW",
#   test = "welch",
#   xlab = "Dispersion",
#   ylab = "Difference",
#   cutoff = 0.05
# )

# RUN DEPRECATED ANCOMBC1
# ancombcOUTPUT = ancombc(
#   data = physeqfilt, 
#   formula = group, 
#   p_adj_method = "fdr", 
#   prv_cut = 0, # no prev filtering necessary anymore 
#   lib_cut = 0, 
#   group = group, 
#   struc_zero = TRUE, 
#   neg_lb = TRUE, 
#   tol = 1e-5, 
#   max_iter = 100, 
#   conserve = TRUE, 
#   alpha = 0.05, 
#   global = TRUE
# )
# 
# ancombcRES <- ancombcOUTPUT$res$diff_abn


