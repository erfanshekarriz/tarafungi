library(tidyverse)
library(phyloseq)
library(scatterpie)

# FUNGI ####
# LOAD
physeqFung <- readRDS("./data/rds/physeqTARAFungi.rds")
metadataFung <- data.frame(sample_data(physeqFung))
features <- data.frame(features=colnames(metadataFung))

physeqPhylum <- tax_glom(physeqFung, taxrank = rank_names(physeqFung)[4])
taxdfPhylum <- data.frame(tax_table(physeqPhylum))
otuPhylum <- data.frame(otu_table(physeqPhylum)) %>%
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "OTUid") %>%
  mutate(OTUid = str_replace(OTUid, "X", "")) %>%
  mutate(Phylum = taxdfPhylum[OTUid, "Phylum"]) %>%
  column_to_rownames(var = "Phylum") %>%
  select(-OTUid) %>%
  apply(MARGIN = 2, FUN = function(x) x/sum(x)) %>%
  t() %>% as.data.frame() %>% 
  merge(metadataFung, by = "row.names") %>%
  mutate(SampleID = as.factor(Row.names)) %>%
  column_to_rownames(var = "Row.names")

 
otuPhylumlong <- otuPhylum %>%
  pivot_longer(cols = unique(taxdfPhylum$Phylum), 
               names_to = "Phylum", 
               values_to = "RelAbund") %>%
  mutate(Longitude = as.numeric(Longitude), 
         Latitude = as.numeric(Latitude)) %>% 
  rename(Region = `OS.region...abbreviation..full.name..MRG....`) 

otuPhylumlong%>%
  group_by(Phylum) %>% mutate(RelAbundPhylAvg = mean(RelAbund, na.rm = TRUE)) %>%
  mutate(Phylum = if_else(RelAbundPhylAvg < 0.02, "Other", Phylum), 
         Region = str_extract(Region, "\\((.+)\\)", group = 1)) %>%
  ggplot(aes(x = SampleID, y = RelAbund, fill = Phylum, group = Region)) + 
  geom_bar(position="fill", stat="identity",  width = 1) + 
  facet_grid(~Region, scales = "free_x", switch="y", space = "free") +
  scale_y_continuous(expand = c(0, 0), position = "left") + 
  themetaxabarplot + 
  ylab("Relative Abundance")  + 
  guides(fill=guide_legend(nrow=2)) + 
  scale_fill_brewer(palette = "Set1")
  




# GLOBAL MAP
otuPhylumRegion <- otuPhylumlong %>% 
  group_by(Phylum) %>% 
  mutate(RelAbundPhylAvg = mean(RelAbund, na.rm = TRUE)) %>%
  mutate(Phylum = if_else(RelAbundPhylAvg < 0.02, "Other", Phylum)) %>%
  select(Longitude, Latitude, RelAbund, Phylum) %>%
  mutate(Longitude = plyr::round_any(Longitude, 10) , 
         Latitude = plyr::round_any(Latitude, 10) , 
         Group = paste(Latitude, Longitude)) %>%
  group_by(Group, Phylum) %>%
  dplyr::summarise(meanRelAbund = mean(RelAbund, na.rm = TRUE), 
            long = as.numeric(Latitude), 
            lat = as.numeric(Longitude), 
            Phylum = Phylum, 
            region = as.factor(Group)) %>%
  pivot_wider(names_from = Phylum, values_from = meanRelAbund, 
              values_fn = mean)  %>%
  as.data.frame()

world <- map_data("world")
p <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill=NA, color="grey20") +
  coord_quickmap()


p +
  geom_scatterpie(aes(x = lat, y = long, group = region), 
                  data = otuPhylumRegion, 
                  cols = colnames(otuPhylumRegion)[5:8], 
                  pie_scale = 1.2) + 
  themeworldmap + 
  scale_fill_brewer(palette = "Set1") + 
  guides(fill = guide_legend(title = "Phylum", 
                             ncol = 1, 
                             title.position = "top"), 
         color = "none")


metadata$OS.region...abbreviation..full.name..MRG.... %>%  unique()


# ALL GROUPS ####
physeq <- readRDS("./data/rds/physeqTARAALL.rds")
metadata <- data.frame(sample_data(physeq))
features <- data.frame(features=colnames(metadata))
taxdf <- data.frame(tax_table(physeq))

otuClade <- data.frame(otu_table(physeq)) %>% t() %>% 
  apply(MARGIN = 2, FUN = function(x) x/sum(x)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "OTUid") %>%
  mutate(OTUid = str_replace(OTUid, "X", "")) %>%
  mutate(Clade = taxdfClade[OTUid, "Group"]) %>%
  pivot_longer(cols = sample_names(physeq), values_to = "RelAbund", names_to = "SampleID") %>%
  mutate(Clade = if_else(is.na(Clade), "Unkown Eukaryote", Clade)) %>% 
  group_by(Clade, SampleID) %>% 
  summarise(relAbund = sum(RelAbund, na.rm = TRUE)) %>% 
  pivot_wider(names_from = "SampleID", values_from = "relAbund", values_fill = 0) %>% 
  column_to_rownames(var = "Clade") %>% 
  t() %>% as.data.frame() %>% 
  merge(metadata, by = "row.names") %>%
  mutate(SampleID = as.factor(Row.names)) %>%
  column_to_rownames(var = "Row.names")




cols <- data.frame(row.names = length(unique(taxdf$Group)), 
                   Colors = RColorBrewer::brewer.pal(length(unique(taxdf$Group)), name = "Set2"))

otuCladelong <- otuClade %>%
  pivot_longer(cols = unique(taxdf$Group), 
               names_to = "Clade", 
               values_to = "RelAbund") %>%
  mutate(Longitude = as.numeric(Longitude), 
         Latitude = as.numeric(Latitude)) %>% 
  rename(Region = `OS.region...abbreviation..full.name..MRG....`) 



otuCladeProt <- otuCladelong%>%
  filter(!Clade %in% c("Metazoa", "Unkown Eukaryote")) %>% 
  group_by(Clade) %>% mutate(RelAbundPhylAvg = mean(RelAbund, na.rm = TRUE)) %>%
  mutate(Clade = if_else(RelAbundPhylAvg < 0.005, "Other", Clade), 
         Region = str_extract(Region, "\\((.+)\\)", group = 1))
otuCladeAll <- otuCladelong%>%
  group_by(Clade) %>% mutate(RelAbundPhylAvg = mean(RelAbund, na.rm = TRUE)) %>%
  mutate(Clade = if_else(RelAbundPhylAvg < 0.005, "Other", Clade), 
         Region = str_extract(Region, "\\((.+)\\)", group = 1))

otuCladePlot <- rbind(cbind(otuCladeProt, Groupings="Microbial"),
                      cbind(otuCladeAll, Groupings="All"))

otuCladePlot %>% 
  ggplot(aes(x = SampleID, y = RelAbund, fill = Clade, group = Region)) + 
  geom_bar(position="fill", stat="identity",  width = 1) + 
  facet_grid(Groupings~Region, scales = "free_x", space = "free") +
  scale_y_continuous(expand = c(0, 0), position = "left") + 
  themetaxabarplot + 
  ylab("Relative Abundance")  + 
  guides(fill=guide_legend(nrow=2)) + 
  scale_fill_brewer(palette = "Set1")





# GLOBAL MAP
otuCladeRegion <- otuCladeProt %>% 
  group_by(Clade) %>% 
  mutate(RelAbundPhylAvg = mean(RelAbund, na.rm = TRUE)) %>%
  mutate(Clade = if_else(RelAbundPhylAvg < 0.013, "Other", Clade)) %>%
  select(Longitude, Latitude, RelAbund, Clade) %>%
  mutate(Longitude = plyr::round_any(Longitude, 10) , 
         Latitude = plyr::round_any(Latitude, 10) , 
         Group = paste(Latitude, Longitude)) %>%
  group_by(Group, Clade) %>%
  dplyr::summarise(meanRelAbund = mean(RelAbund, na.rm = TRUE), 
                   long = as.numeric(Latitude), 
                   lat = as.numeric(Longitude), 
                   Clade = Clade, 
                   region = as.factor(Group)) %>%
  pivot_wider(names_from = Clade, values_from = meanRelAbund, 
              values_fn = mean)  %>%
  as.data.frame()

ncol(otuCladeRegion)-5


world <- map_data("world")
p <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill=NA, color="grey20") +
  coord_quickmap()


p +
  geom_scatterpie(aes(x = lat, y = long, group = region), 
                  data = otuCladeRegion, 
                  cols = colnames(otuCladeRegion)[5:length(colnames(otuCladeRegion))], 
                  pie_scale = 1.2) + 
  themeworldmap + 
  scale_fill_brewer(palette = "Dark2") + 
  guides(fill = guide_legend(title = "Clade", 
                             ncol = 1, 
                             title.position = "top"), 
         color = "none")



otuCladeRegion <- otuCladeAll %>% 
  group_by(Clade) %>% 
  mutate(RelAbundPhylAvg = mean(RelAbund, na.rm = TRUE)) %>%
  mutate(Clade = if_else(RelAbundPhylAvg < 0.17, "Other", Clade)) %>%
  select(Longitude, Latitude, RelAbund, Clade) %>%
  mutate(Longitude = plyr::round_any(Longitude, 10) , 
         Latitude = plyr::round_any(Latitude, 10) , 
         Group = paste(Latitude, Longitude)) %>%
  group_by(Group, Clade) %>%
  dplyr::summarise(meanRelAbund = mean(RelAbund, na.rm = TRUE), 
                   long = as.numeric(Latitude), 
                   lat = as.numeric(Longitude), 
                   Clade = Clade, 
                   region = as.factor(Group)) %>%
  pivot_wider(names_from = Clade, values_from = meanRelAbund, 
              values_fn = mean)  %>%
  as.data.frame()

ncol(otuCladeRegion)-5


world <- map_data("world")
p <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill=NA, color="grey20") +
  coord_quickmap()


p +
  geom_scatterpie(aes(x = lat, y = long, group = region), 
                  data = otuCladeRegion, 
                  cols = colnames(otuCladeRegion)[5:length(colnames(otuCladeRegion))], 
                  pie_scale = 1.2) + 
  themeworldmap + 
  scale_fill_brewer(palette = "Set2") + 
  guides(fill = guide_legend(title = "Clade", 
                             ncol = 1, 
                             title.position = "top"), 
         color = "none")


  


