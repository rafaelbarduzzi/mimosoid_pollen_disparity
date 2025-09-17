#
# Disparity analyses
#
# Script to analyse the morphological disparity of Mimosa and Stryphnodendron
# clade pollen grains
#
# 01-01-2023
#
# Rafael F. Barduzzi
#

# Libraries ===================================================================

library(ape) # 5.7-1
library(arules) # 1.7-6
library(Claddis) # 0.6.6
library(tidyverse) # 2.0.0
library(phytools) # 1.9.16
library(patchwork) # 1.1.2.9
library(paletteer) # 1.6.0
library(viridis) #  0.6.5
library(Rphylopars) # 0.3.10
library(proxy) # 0.4-27
library(dispRity) # 1.7.0
library(openxlsx) # 4.2.5.2
library(scales) # 1.2.1

# Data ========================================================================

dat_pollen <- read.csv("data/pollen_data.csv", na.strings = c(""))

tree <- read.tree("data/mimo_stryph_phylo.tre")

dat_eco <- read.csv("data/eco_data.csv")

# Prepare data ================================================================

## filter outgroups from data
dat_pollen <- dat_pollen %>%
  filter(!genus %in% c("Anadenanthera", "Vachellia"))

## prune phylogeny
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, 
                                      dat_pollen$name_phylogeny))

## get traits dataframe
traits <- cbind("taxon" = dat_pollen$name_phylogeny, dat_pollen[, 9:27])

## remove min and max traits - we will only work with mean values
traits <- traits %>%
  select(-contains("_min"), -contains("_max"))

## set seed for replicability
set.seed(7)

## create output directory
dir.create("output", showWarnings = F)

## Evaluate data completeness -------------------------------------------------

### traits with less than 60% completness (will be out of the analysis)
names(traits)[colMeans(!is.na(traits)) < 0.6]

### select traits with data completeness higher or equal to 60%
traits <- traits[, colMeans(!is.na(traits)) >= 0.6]

## check missing data after missing data filter

### overall
sum(is.na(traits)) / prod(dim(traits)) * 100

### traits
colSums(is.na(traits)) / nrow(traits) * 100

## Check continuous traits distribution and normality -------------------------

## get continuous traits
traits_cont <- c("shorter_diameter_mean", 
                 "longer_diameter_mean", 
                 "exine_thickness_mean")

for (x in colnames(traits)) {
  if (is.numeric(traits[[x]])) {
    plot(density(na.omit(traits[[x]])), main = x, xlab = x)
  }
}

shapiro.test(na.omit(traits[, c("shorter_diameter_mean")]))
shapiro.test(na.omit(traits[, c("longer_diameter_mean")]))
shapiro.test(na.omit(traits[, c("exine_thickness_mean")]))

## continuous data are not normally distributed

## log scale continuous data (Kriebel et al. 2017)
traits[, traits_cont] <- lapply(traits[, traits_cont], function(x) log(x))

## Cont data discretization ===================================================

### define the best number of breaks (character states)

# "FD" - In statistics, the Freedman Diaconis rule can be used to select the 
# width of the bins to be used in a histogram. Freedman-Diaconis rule is 
# designed roughly to minimize the integral of the squared difference between 
# the histogram (i.e., relative frequency density) and the density of the 
# theoretical probability distribution.

traits_breaks <- list(
  ### shorter_diameter_mean
  sd = hist(traits$shorter_diameter_mean,
            breaks = "FD",
            main = "shorter_diameter_mean")$counts,
  ### longer_diameter_mean
  ld = hist(traits$longer_diameter_mean,
            breaks = "FD",
            main = "longer_diameter_mean")$counts,
  ### exine_thickness_mean
  et = hist(traits$exine_thickness_mean,
            breaks = "FD",
            main = "exine_thickness_mean")$counts
)

### remove bins with 0 count
traits_breaks <- lapply(traits_breaks, function(x) x[x != 0])

### run the discretization function

traits_discretized <- arules::discretizeDF(
  traits,
  methods = list
  (shorter_diameter_mean = list(
    method = "interval",
    breaks = length(traits_breaks$sd),
    labels = c(0:(length(traits_breaks$sd)-1))),
    longer_diameter_mean = list(
      method = "interval",
      breaks = length(traits_breaks$ld),
      labels = c(0:(length(traits_breaks$ld)-1))),
    exine_thickness_mean = list(
      method = "interval",
      breaks = length(traits_breaks$et),
      labels = c(0:(length(traits_breaks$et)-1)))))

traits_discretized <- column_to_rownames(traits_discretized, var = "taxon")

as.matrix(traits_discretized) -> traits_discretized

# Pollen morphospace ==========================================================

## Ordinate matrix (PCoA) -----------------------------------------------------

## defining characters ordination
ord <- rep("unordered", 7)

## create Claddis cladistic matrix
cladistic_matrix <- build_cladistic_matrix(traits_discretized, ordering = ord)

## calculating distance matrix
dist <- calculate_morphological_distances(cladistic_matrix)

## pcoa
pcoa <- ape::pcoa(dist$distance_matrix, correction = "cailliez")

## calculating variance explained explained by each principal component
eigen <- pcoa$values$Corr_eig

# PC1
round(eigen[1] / sum(eigen[eigen > 0]) * 100, 2)
# PC2
round(eigen[2] / sum(eigen[eigen > 0]) * 100, 2)
# PC3
round(eigen[3] / sum(eigen[eigen > 0]) * 100, 2)

## Plot phylo-morphospace (Fig. 2) --------------------------------------------

### extract two first pcoa axis and filtering to taxa in tree

space_axis <- pcoa$vectors.cor[, c(1,2)]

space_axis <- space_axis[rownames(space_axis) %in% tree$tip.label, ]

# creating groups (of genera)
groups <- split(dat_pollen$name_phylogeny, dat_pollen$genus)

### creat plot using phytools

phylomorphospace <-
  phylomorphospace(
    tree_pruned,
    space_axis,
    xlab = "PC1",
    ylab = "PC2",
    label = "off"
  )

## extract configs

phylo_data <- data.frame(
  xstart = phylomorphospace$xx[phylomorphospace$edge[, 1]],
  ystart = phylomorphospace$yy[phylomorphospace$edge[, 1]],
  xstop = phylomorphospace$xx[phylomorphospace$edge[, 2]],
  ystop = phylomorphospace$yy[phylomorphospace$edge[, 2]],
  nodestart = phylomorphospace$edge[, 1],
  nodestop = phylomorphospace$edge[, 2]
)

data_vectors_genus <- as.data.frame(pcoa$vectors.cor[, c (1, 2)])

data_vectors_genus %>%
  add_column(as.factor(dat_pollen$genus), .before = T) -> data_vectors_genus
colnames(data_vectors_genus) <- c("Genus", "PC1", "PC2")

data_vectors_genus <-
  data_vectors_genus[row.names(data_vectors_genus) %in% tree$tip.label, ]

# check labels

phylomorphospace <-
  phylomorphospace(
    tree_pruned,
    space_axis,
    xlab = "PC1",
    ylab = "PC2",
    label = "horizontal"
  )

# plot

phylomorphospace <-
  ggplot() +
  geom_segment(
    data = phylo_data,
    aes(
      x = xstart,
      y = ystart,
      xend = xstop,
      yend = ystop
    ),
    linewidth = 0.3,
    color = "#4F5157"
  ) +
  geom_point(data = data_vectors_genus,
             aes(
               x = PC1,
               y = PC2,
               color = Genus,
               shape = Genus
             ),
             size = 3) +
  scale_color_manual(
    values = 
      paletteer_d("ggthemes::calc")[c(5,2,1,1,5,12,8,7,8,9,10,1,6)]) +
  scale_shape_manual(
    values = c(
      "Adenopodia" = 17,
      "Gwilymia" = 0,
      "Inga" = 4,
      "Lachesiodendron" = 4,
      "Marlimorimia" = 1,
      "Microlobius" = 2,
      "Mimosa" = 20,
      "Naiadendron" = 5,
      "Parapiptadenia" = 6,
      "Piptadenia" = 15,
      "Pityrocarpa" = 13,
      "Stryphnodendron" = 9,
      "Senegalia" = 4
    )
  ) +
  coord_fixed(
    ratio = 1,
    xlim = c(-5, 3),
    ylim = c(-4, 3.5),
    expand = TRUE,
    clip = "on"
  ) +
  xlab(paste(
    "PCo",1," ","(",round(pcoa$values$Rel_corr_eig[1] * 100, 2),"%",")",
    sep = ""
  )) + ylab(paste(
    "PCo",2," ","(",round(pcoa$values$Rel_corr_eig[2] * 100, 2),"%",")",
    sep = ""
  )) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# check

phylomorphospace

## Plot eco-morphospace (Fig. 5) ----------------------------------------------

### Preparing data ###

dat_eco <- subset(dat_eco, taxon %in% dat_pollen$name_phylogeny)

### Analysis ###

# filtering outgrops

dat_eco <- dat_eco %>% filter(!(
  grepl('Anadenanthera', .$taxon) |
    grepl('Vachellia', .$taxon)))

data_vectors_eco <- pcoa$vectors.cor[, c (1, 2)]

data_vectors_eco <- as.data.frame(data_vectors_eco) %>%
  add_column(as.factor(dat_eco$biome), .before = T) %>%
  add_column(as.factor(dat_pollen$genus), .before = T)
colnames(data_vectors_eco) <- c("Genus", "biome", "PC1", "PC2")

# checking taxa of each vegetation type #

groups_veg_type <- split(dat_eco$taxon, dat_eco$biome)

# finding the convex hulls (vegetation type) of the points being plotted...

# tropical rainforest
data_vectors_eco %>%
  filter(grepl("5", biome)) %>%
  select(PC1, PC2) %>%
  slice(chull(PC1, PC2)) -> hull_trf

# savanna
data_vectors_eco %>%
  filter(grepl("2", biome)) %>%
  select(PC1, PC2) %>%
  slice(chull(PC1, PC2)) -> hull_sav

# sdtf
data_vectors_eco %>%
  filter(grepl("3", biome)) %>%
  select(PC1, PC2) %>%
  slice(chull(PC1, PC2)) -> hull_sdtf

# flooded grassland
data_vectors_eco %>%
  filter(grepl("1", biome)) %>%
  select(PC1, PC2) %>%
  slice(chull(PC1, PC2)) -> hull_fg

# temperate
data_vectors_eco %>%
  filter(grepl("4", biome)) %>%
  select(PC1, PC2) %>%
  slice(chull(PC1, PC2)) -> hull_temp

# desert
data_vectors_eco %>%
  filter(grepl("0", biome)) %>%
  select(PC1, PC2) %>%
  slice(chull(PC1, PC2)) -> hull_des

# defining the scatterplot

morphospace <- ggplot(data_vectors_eco, mapping = aes(x = PC1, y = PC2)) +
  geom_point() +
  coord_fixed(
    ratio = 1,
    xlim = c(-5, 3),
    ylim = c(-4, 3.5),
    expand = TRUE,
    clip = "on"
  ) +
  xlab(
    paste("PCo",1," ","(",round(pcoa$values$Rel_corr_eig[1] * 100, 2),"%",")",
          sep = ""
    )) + ylab(
      paste("PCo",2," ","(",round(pcoa$values$Rel_corr_eig[2] * 100, 2),"%",")",
            sep = ""
      )) +
  theme_light() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

morphospace_eco <- 
  morphospace +
  geom_polygon(
    data = hull_trf,
    alpha = 0.3,
    inherit.aes = FALSE,
    mapping = aes(x = PC1, y = PC2, fill = "Tropical rainforest")
  ) +
  geom_polygon(
    data = hull_sav,
    alpha = 0.5,
    inherit.aes = FALSE,
    mapping = aes(x = PC1, y = PC2, fill = "Savanna")
  ) +
  geom_polygon(
    data = hull_sdtf,
    alpha = 0.4,
    inherit.aes = FALSE,
    mapping = aes(x = PC1, y = PC2, fill = "SDTF")
  ) +
  geom_polygon(
    data = hull_fg,
    alpha = 0.5,
    inherit.aes = FALSE,
    mapping = aes(x = PC1, y = PC2, fill = "Flooded grassland")
  ) +
  geom_polygon(
    data = hull_temp,
    alpha = 0.6,
    inherit.aes = FALSE,
    mapping = aes(x = PC1, y = PC2, fill = "Temperate")
  ) +
  geom_polygon(
    data = hull_des,
    alpha = 0.6,
    inherit.aes = FALSE,
    mapping = aes(x = PC1, y = PC2, fill = "Desert")
  ) +
  scale_fill_manual(
    values = c(
      "Tropical rainforest" = "#B85A0DFF",
      "Savanna" = "#3CB7CCFF",
      "SDTF" = "#FFD94AFF",
      "Flooded grassland" = "#FF7F0FFF",
      "Temperate" = "#39737CFF",
      "Desert" = "#32A251FF"
    )
  ) +
  labs(fill = "Biomes")

# checking

morphospace_eco

# Morphospace Occupancy Metrics ===============================================

## Defining groups ------------------------------------------------------------

# groups of genus (no outgroups) #

groups_genera <- split(dat_pollen$name_phylogeny, dat_pollen$genus)
groups_genera$Inga <- NULL
groups_genera$Senegalia <- NULL

# groups of principal clades #

groups_clades <- list(
  "Mimosa_clade" = filter(dat_pollen, genus %in% 
                            c("Mimosa", "Adenopodia", "Piptadenia")) %>%
    pull(name_phylogeny),
  "Stryphnodendron_clade" = filter(
    dat_pollen,!genus %in% 
      c("Mimosa", "Adenopodia", "Piptadenia", "Inga", "Senegalia", "Lachesiodendron")
  ) %>%
    pull(name_phylogeny)
)

# ecological groups #

groups_ecological <- list(
  "Desert" =  filter(dat_eco, grepl("0", biome)) %>%
    pull(taxon) %>%
    as.vector(),
  "Flooded_grassland" =  filter(dat_eco, grepl("1", biome)) %>%
    pull(taxon) %>%
    as.vector(),
  "Savanna" =  filter(dat_eco, grepl("2", biome)) %>%
    pull(taxon) %>%
    as.vector(),
  "SDTF" =  filter(dat_eco, grepl("3", biome)) %>%
    pull(taxon) %>%
    as.vector(),
  "Temperate" =  filter(dat_eco, grepl("4", biome)) %>%
    pull(taxon) %>%
    as.vector(),
  "Tropical_rainforest" =  filter(dat_eco, grepl("5", biome)) %>%
    pull(taxon) %>%
    as.vector()
)

# ecological lineages (based on biomes ancestral reconstruction delimitation) #

# Desert lineages (state 0) - no lineage delimitation was possible

# Flooded grassland lineages (state 1)

groups_flooded <- list (
  "Mimosa_Flooded" = dat_eco %>%
    filter(
      taxon %in% intersect(tree$tip.label[getDescendants(tree, getMRCA(tree, c(
        "Mimosa_dormiens", "Mimosa_pigra_dehiscens"
      )))], dat_pollen$name_phylogeny)
    ) %>%
    filter(grepl("1", biome)) %>%
    pull(taxon) %>%
    as.vector()
)

# savanna lineages (state 2)

groups_savanna <- list (
  "Mimosa_Savanna_1" = dat_eco %>%
    filter(
      taxon %in% intersect(tree$tip.label[getDescendants(tree, getMRCA(tree, c("Mimosa_filipes", "Mimosa_gatesiae")))], dat_pollen$name_phylogeny)
    ) %>%
    filter(grepl("2", biome)) %>%
    pull(taxon) %>%
    as.vector(),
  "Mimosa_Savanna_2" = dat_eco %>%
    filter(
      taxon %in% intersect(tree$tip.label[getDescendants(tree, getMRCA(tree, c(
        "Mimosa_minarum", "Mimosa_longepedunculata"
      )))], dat_pollen$name_phylogeny)
    ) %>%
    filter(grepl("2", biome)) %>%
    pull(taxon) %>%
    as.vector(),
  "Stryphnod_Savanna" = dat_eco %>%
    filter(
      taxon %in% intersect(tree$tip.label[getDescendants(tree, getMRCA(
        tree,
        c("Stryphnodendron_velutinum", "Stryphnodendron_adstringens")
      ))], dat_pollen$name_phylogeny)
    ) %>%
    filter(grepl("2", biome)) %>%
    pull(taxon) %>%
    as.vector()
)

# sdtf lineages (state 3)

groups_sdtf <- list (
  "Mimosa_SDTF_1" = dat_eco %>%
    filter(
      taxon %in% intersect(tree$tip.label[getDescendants(tree, getMRCA(tree, c(
        "Piptadenia_floribunda", "Mimosa_longepedunculata"
      )))], dat_pollen$name_phylogeny)
    ) %>%
    filter(grepl("3", biome)) %>%
    pull(taxon) %>%
    as.vector(),
  "Mimosa_SDTF_2" = dat_eco %>%
    filter(
      taxon %in% intersect(tree$tip.label[getDescendants(tree, getMRCA(tree, c("Mimosa_minarum", "Mimosa_cordistipula")))], dat_pollen$name_phylogeny)
    ) %>%
    filter(grepl("3", biome)) %>%
    pull(taxon) %>%
    as.vector(),
  "Mimosa_SDTF_3" = dat_eco %>%
    filter(
      taxon %in% intersect(tree$tip.label[getDescendants(tree, getMRCA(tree, c("Mimosa_woodii", "Mimosa_deamii")))], dat_pollen$name_phylogeny)
    ) %>%
    filter(grepl("3", biome)) %>%
    pull(taxon) %>%
    as.vector(),
  "Stryphnod_SDTF" = dat_eco %>%
    filter(
      taxon %in% intersect(tree$tip.label[getDescendants(tree, getMRCA(tree, c(
        "Pityrocarpa_brenanii", "Pityrocarpa_moniliformis"
      )))], dat_pollen$name_phylogeny)
    ) %>%
    filter(grepl("3", biome)) %>%
    pull(taxon) %>%
    as.vector()
)

# tempera lineages (state 4) - no lineage delimitation was possible

# tropical rainforest lineages (state 5)

groups_tropical <- list (
  "Mimosa_Tropical_1" = dat_eco %>%
    filter(
      taxon %in% intersect(tree$tip.label[getDescendants(tree, getMRCA(tree, c(
        "Piptadenia_robusta", "Piptadenia_adiantoides"
      )))], dat_pollen$name_phylogeny)
    ) %>%
    filter(grepl("5", biome)) %>%
    pull(taxon) %>%
    as.vector(),
  "Mimosa_Tropical_2" = dat_eco %>%
    filter(
      taxon %in% intersect(tree$tip.label[getDescendants(tree, getMRCA(tree, c(
        "Mimosa_lepidophora", "Mimosa_guilandinae"
      )))], dat_pollen$name_phylogeny)
    ) %>%
    filter(grepl("5", biome)) %>%
    pull(taxon) %>%
    as.vector(),
  "Mimosa_Tropical_3" = dat_eco %>%
    filter(
      taxon %in% intersect(tree$tip.label[getDescendants(tree, getMRCA(tree, c("Mimosa_gymnas", "Mimosa_involucrata")))], dat_pollen$name_phylogeny)
    ) %>%
    filter(grepl("5", biome)) %>%
    pull(taxon) %>%
    as.vector(),
  "Mimosa_Tropical_4" = dat_eco %>%
    filter(
      taxon %in% intersect(tree$tip.label[getDescendants(tree, getMRCA(tree, c("Mimosa_taimbensis", "Mimosa_urticaria")))], dat_pollen$name_phylogeny)
    ) %>%
    filter(grepl("5", biome)) %>%
    pull(taxon) %>%
    as.vector(),
  "Stryphnod_Tropical" = dat_eco %>%
    filter(
      taxon %in% intersect(tree$tip.label[getDescendants(tree, getMRCA(
        tree,
        c("Lachesiodendron_viridiflorum", "Microlobius_foetidus")
      ))], dat_pollen$name_phylogeny)
    ) %>%
    filter(grepl("5", biome)) %>%
    pull(taxon) %>%
    as.vector()
)

# removing subsets that will affect bootstrap/rarefaction output in step 6
# i.e. groups with less than three taxa

groups_genera <- groups_genera[sapply(groups_genera, length) >= 3]
groups_savanna <- groups_savanna[sapply(groups_savanna, length) >= 3]
groups_sdtf <- groups_sdtf[sapply(groups_sdtf, length) >= 3]
groups_tropical <- groups_tropical[sapply(groups_tropical, length) >= 3]

# flooded grassland ecological lineage does not have enough taxa
sapply(groups_flooded, length) >= 3

# create ecological lineages list
ecological_lineages <- c(groups_savanna, groups_sdtf, groups_tropical)

## Density - MPD --------------------------------------------------------------
### Mean Pairwise Distance - Pre-ordination metric

### Since dispRity do not deal with polymorphic matrices, I had to use Claddis

### MPD for genera
class(groups_genera) <- "taxonGroups"
mpd_genera <- calculate_MPD(dist, groups_genera)
class(groups_genera) <- "list"

### MPD for clades
class(groups_clades) <- "taxonGroups"
mpd_clades <- calculate_MPD(dist, groups_clades)
class(groups_clades) <- "list"


### MPD for ecological groups
class(groups_ecological) <- "taxonGroups"
mpd_ecological <- calculate_MPD(dist, groups_ecological)
class(groups_ecological) <- "list"


### MPD for ecological lineages
class(ecological_lineages) <- "taxonGroups"
mpd_ecolineages <- calculate_MPD(dist, ecological_lineages)
class(ecological_lineages) <- "list"

## Size - SV ------------------------------------------------------------------
# sum of variances - post-ordination metric

### This part of the script was based on the code of Casali et al. (2023)
### available at https://doi.org/10.5281/zenodo.7240006

bootstraps <- 1000
rarefaction <- "min" # "min" or FALSE

### SV - Genera
sub_genera <- custom.subsets(pcoa$vectors.cor, groups_genera)
bootstrapped_data <- boot.matrix(sub_genera, 
                                 bootstraps = bootstraps, 
                                 rarefaction = rarefaction)
genera_sv <- dispRity(bootstrapped_data, metric = c(sum, variances))
genera_sv_smry <- summary(genera_sv)
genera_sv_wilcox <- test.dispRity(genera_sv, 
                                  wilcox.test,
                                  correction = "holm")
genera_sv_bhatt <- test.dispRity(genera_sv, bhatt.coeff)

### SV - Clades
sub_clades <- custom.subsets(pcoa$vectors.cor, groups_clades)
bootstrapped_data <- boot.matrix(sub_clades, 
                                 bootstraps = bootstraps, 
                                 rarefaction = rarefaction)
clades_sv <- dispRity(bootstrapped_data, metric = c(sum, variances))
clades_sv_smry <- summary(clades_sv)
clades_sv_wilcox <- test.dispRity(clades_sv, 
                                  wilcox.test, 
                                  correction = "holm")
clades_sv_bhatt <- test.dispRity(clades_sv, 
                                 bhatt.coeff)

### SV - Ecological groups
sub_ecological <- custom.subsets(pcoa$vectors.cor, groups_ecological)
bootstrapped_data <- boot.matrix(sub_ecological, 
                                 bootstraps = bootstraps,
                                 rarefaction = rarefaction)
ecological_sv <- dispRity(bootstrapped_data, metric = c(sum, variances))
ecological_sv_smry <- summary(ecological_sv)
ecological_sv_wilcox <- test.dispRity(ecological_sv, 
                                      wilcox.test, 
                                      correction = "holm")
ecological_sv_bhatt <- test.dispRity(ecological_sv, bhatt.coeff)

### SV - Ecological lineages
sub_ecolineages <- custom.subsets(pcoa$vectors.cor, ecological_lineages)
bootstrapped_data <- boot.matrix(sub_ecolineages, 
                                 bootstraps = bootstraps,
                                 rarefaction = rarefaction)
ecolineages_sv <- dispRity(bootstrapped_data, metric = c(sum, variances))
ecolineages_sv_smry <- summary(ecolineages_sv)
ecolineages_sv_wilcox <- test.dispRity(ecolineages_sv,
                                       wilcox.test, 
                                       correction = "holm")
ecolineages_sv_bhatt <- test.dispRity(ecolineages_sv, bhatt.coeff)

## Summary and export ---------------------------------------------------------

labels_mpd <- c(rep("Genera", 7), (rep("Clades", 2)), (rep("Biomes", 6)), 
                (rep("Ecological lineages", 5)))

### MPD

summary_mpd <- data.frame(
  Subdivision = labels_mpd,
  Group = 
    names(c(mpd_genera, mpd_clades, mpd_ecological, mpd_ecolineages)),
  MPD = 
    round(unname(c(mpd_genera, mpd_clades, mpd_ecological, mpd_ecolineages)), 
          2))

labels <- c(rep("Genera", 7), (rep("Clades", 2)), (rep("Biomes", 6)), 
            (rep("Ecological lineages", 5)))


if (rarefaction == "min") { # rewrite if using rarefaction
  
  labels <- c(rep("Genera", 13), (rep("Clades", 3)), (rep("Biomes", 11)), 
              (rep("Ecological lineages", 9)))
  
}

### SV

genera_disp <- do.call(rbind, genera_sv_smry)
clade_disp <- do.call(rbind, clades_sv_smry)
ecological_disp <- do.call(rbind, ecological_sv_smry)
ecolineages_disp <- do.call(rbind, ecolineages_sv_smry)
colnames(genera_sv_smry)[8] <- "97.5%"
colnames(clades_sv_smry)[8] <- "97.5%"
colnames(ecological_sv_smry)[8] <- "97.5%"
colnames(ecolineages_sv_smry)[8] <- "97.5%"
info_disp <- rbind(genera_sv_smry, clades_sv_smry, ecological_sv_smry, 
                   ecolineages_sv_smry)
options(digits = 2)
summary_sv <- data.frame(labels, info_disp)
summary_sv <- summary_sv[-c(7:8)]
colnames(summary_sv) <-
  c("Subdivision",
    "Group",
    "Sample_size",
    "SV_Obs",
    "Median_BS",
    "0.025_BS",
    "0.975_BS")
summary_sv$Group[summary_sv$Group == "Mimosa_clade"] <-
  "Mimosa clade"
summary_sv$Group[summary_sv$Group == "Stryphnodendron_clade"] <-
  "Stryphnodendron clade"
summary_sv$Group[summary_sv$Group == "Flooded_grassland"] <-
  "Flooded grassland"
summary_sv$Group[summary_sv$Group == "Tropical_rainforest"] <-
  "Tropical rainforest"
summary_sv$Group[summary_sv$Group == "mimosa_sdtf_2"] <-
  "Mimosa SDTF 2"
summary_sv$Group[summary_sv$Group == "mimosa_sdtf_3"] <-
  "Mimosa SDTF 3"
summary_sv$Group[summary_sv$Group == "mimosa_sdtf_5"] <-
  "Mimosa SDTF 5"
summary_sv$Group[summary_sv$Group == "Stryphnod_SDTF"] <-
  "Stryphnodendron SDTF"
summary_sv$Group[summary_sv$Group == "mimosa_tropical_1"] <-
  "Mimosa Tropical 1"
summary_sv$Group[summary_sv$Group == "Stryphnod_Tropical"] <-
  "Stryphnodendron Tropical"

### SV statistics

statistic_W_sv <-
  c(
    unlist(genera_sv_wilcox)[c(1:21)],
    unlist(clades_sv_wilcox)[c(1)],
    unlist(ecological_sv_wilcox)[c(1:15)],
    unlist(ecolineages_sv_wilcox)[c(1:10)]
  )

p_value_sv <-
  c(
    unlist(genera_sv_wilcox)[c(22:42)],
    unlist(clades_sv_wilcox)[c(2)],
    unlist(ecological_sv_wilcox)[c(16:30)],
    unlist(ecolineages_sv_wilcox)[c(11:20)]
  )

b_coeff_sv <-
  c(
    unlist(genera_sv_bhatt),
    unlist(clades_sv_bhatt),
    unlist(ecological_sv_bhatt),
    unlist(ecolineages_sv_bhatt)
  )
gr_lab_sv <-
  c(
    rownames(genera_sv_bhatt),
    rownames(clades_sv_bhatt),
    rownames(ecological_sv_bhatt),
    rownames(ecolineages_sv_bhatt)
  )
pt_lab_sv <-
  c(rep("Genera", 21),
    "Clades",
    rep("Biomes", 15),
    rep("Ecological lineages", 10))

summary_stats_sv <-
  data.frame(pt_lab_sv,
             gr_lab_sv,
             round(statistic_W_sv, 2),
             round(p_value_sv, 2),
             round(b_coeff_sv, 2))

colnames(summary_stats_sv) <-
  c("Subdivision",
    "Groups",
    "Statistic W",
    "P-value",
    "Bhattacharyya Coefficient")

### saving workbooks

wb <- createWorkbook()
addWorksheet(wb, "MPD")
writeData(wb, 1, summary_mpd)
saveWorkbook(wb, "output/mean_pairwise_distance.xlsx", overwrite = TRUE)

wb <- createWorkbook()
addWorksheet(wb, "SV")
addWorksheet(wb, "Stats - SV")
writeData(wb, 1, summary_sv)
writeData(wb, 2, summary_stats_sv)
saveWorkbook(wb, "output/sum_of_variances.xlsx", overwrite = TRUE)

# Correlation between metrics =================================================

mpd_sv_corr <- cor.test(summary_mpd$MPD, summary_sv$SV_Obs[!is.na(summary_sv$SV_Obs)], 
                        method = "spearman",
                        exact = F,
                        alternative = "two.sided")
mpd_sv_corr$p.value
### highly positive correlation statistically significant between obsv. values

# SV Boxplots (Supplementary Fig. S3) =========================================

### This part of the script was based on the code of Casali et al. (2023)
### available at https://doi.org/10.5281/zenodo.7240006

### required functions

ci <- function(x)
{
  exit <- quantile(x, probs = c(0.025, 0.975))
  names(exit) <- c("ymin", "ymax")
  return(exit)
}

median <- function(x)
{
  exit <- quantile(x, probs = c(0.5))
  names(exit) <- c("y")
  return(exit)
}

#### Genera
F1 <- c(rep("Gwilymia", bootstraps))
F2 <- c(rep("Marlimorimia", bootstraps))
F3 <- c(rep("Mimosa", bootstraps))
F4 <- c(rep("Parapiptadenia", bootstraps))
F5 <- c(rep("Piptadenia", bootstraps))
F6 <- c(rep("Pityrocarpa", bootstraps))
F7 <- c(rep("Stryphnodendron", bootstraps))
genera_plot_disp <- vector(mode = "list", length = length(genera_sv))

sum_Gwilymia <- as.vector(genera_sv$disparity$Gwilymia[[2]])
sum_Marlimorimia <-
  as.vector(genera_sv$disparity$Marlimorimia[[2]])
sum_Mimosa <- as.vector(genera_sv$disparity$Mimosa[[2]])
sum_Parapiptadenia <-
  as.vector(genera_sv$disparity$Parapiptadenia[[2]])
sum_Piptadenia <- as.vector(genera_sv$disparity$Piptadenia[[2]])
sum_Pityrocarpa <-
  as.vector(genera_sv$disparity$Pityrocarpa[[2]])
sum_Stryphnodendron <-
  as.vector(genera_sv$disparity$Stryphnodendron[[2]])

SV_Gwilymia <- data.frame(F1, sum_Gwilymia)
colnames(SV_Gwilymia) <- c("Group", "Sum_of_variances")
SV_Marlimorimia <- data.frame(F2, sum_Marlimorimia)
colnames(SV_Marlimorimia) <- c("Group", "Sum_of_variances")
SV_Mimosa <- data.frame(F3, sum_Mimosa)
colnames(SV_Mimosa) <- c("Group", "Sum_of_variances")
SV_Parapiptadenia <- data.frame(F4, sum_Parapiptadenia)
colnames(SV_Parapiptadenia) <- c("Group", "Sum_of_variances")
SV_Piptadenia <- data.frame(F5, sum_Piptadenia)
colnames(SV_Piptadenia) <- c("Group", "Sum_of_variances")
SV_Pityrocarpa <- data.frame(F6, sum_Pityrocarpa)
colnames(SV_Pityrocarpa) <- c("Group", "Sum_of_variances")
SV_Stryphnodendron <- data.frame(F7, sum_Stryphnodendron)
colnames(SV_Stryphnodendron) <- c("Group", "Sum_of_variances")
genera_to_plot <-
  rbind(
    SV_Gwilymia,
    SV_Marlimorimia,
    SV_Mimosa,
    SV_Parapiptadenia,
    SV_Piptadenia,
    SV_Pityrocarpa,
    SV_Stryphnodendron
  )

genera_to_plot$Group <-
  factor(
    genera_to_plot$Group,
    levels = c(
      "Gwilymia",
      "Marlimorimia",
      "Mimosa",
      "Parapiptadenia",
      "Piptadenia",
      "Pityrocarpa",
      "Stryphnodendron"
    )
  )

plot_genera_sv <- ggplot(genera_to_plot,
                         aes(
                           x = Group,
                           y = Sum_of_variances,
                           colour = Group,
                           fill = Group
                         ),
                         alpha = 0.9) +
  geom_boxplot(alpha = 0.5, show.legend = F) + labs(x = "", y = "") +
  scale_fill_manual(values = as.character(paletteer_d(
    "ggthemes::excel_Office_2007_2010"))[1:7]
  ) +
  scale_colour_manual(values = as.character(paletteer_d(
    "ggthemes::excel_Office_2007_2010"))[1:7]
  ) +
  stat_summary(
    fun = ci,
    geom = "line",
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 2.5,
    show.legend = FALSE
  ) +
  theme_bw() + labs("") +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 10, face = "bold", color = "black",
                               angle = 40, 
                               vjust = 0.6),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(labels = c(
    "Gwilymia (7)",
    "Marlimorimia (4)",
    "Mimosa (135)",
    "Parapiptadenia (5)",
    "Piptadenia (18)",
    "Pityrocarpa (7)",
    "Stryphnod. (6)"
  )) +
  scale_y_continuous(labels = label_number(accuracy = 0.1 , decimal.mark =
                                             "."))  +
  guides(fill = guide_legend(nrow = 2)) +
  ggtitle("B") +
  theme(plot.title = element_text(
    size = 14,
    hjust = 0
  ))

#### Clades
C1 <- c(rep("Mimosa clade", bootstraps))
C2 <- c(rep("Stryphnodendron clade", bootstraps))
clades_plot_disp <- vector(mode = "list", length = length(clades_sv))

sum_Mimosa_clade <-
  as.vector(clades_sv$disparity$Mimosa_clade[[2]])
sum_Stryphnodendron_clade <-
  as.vector(clades_sv$disparity$Stryphnodendron_clade[[2]])

SV_Mimosa_clade <- data.frame(C1, sum_Mimosa_clade)
colnames(SV_Mimosa_clade) <- c("Group", "Sum_of_variances")
SV_Stryphnodendron_clade <- data.frame(C2, sum_Stryphnodendron_clade)
colnames(SV_Stryphnodendron_clade) <- c("Group", "Sum_of_variances")
clades_to_plot <- rbind(SV_Mimosa_clade, SV_Stryphnodendron_clade)

clades_to_plot$Group <-
  factor(clades_to_plot$Group,
         levels = c("Mimosa clade", "Stryphnodendron clade"))

plot_clades_sv <-
  ggplot(clades_to_plot,
         aes(
           x = Group,
           y = Sum_of_variances,
           colour = Group,
           fill = Group
         ),
         alpha = 0.9) +
  geom_boxplot(alpha = 0.5, show.legend = F) + labs(x = "", y = "Sum of variances") +
  scale_fill_manual(values = as.character(paletteer_d(
    "ggthemes::excel_Office_2007_2010"))[1:2]
  ) +
  scale_colour_manual(values = as.character(paletteer_d(
    "ggthemes::excel_Office_2007_2010"))[1:2]
  ) +
  stat_summary(
    fun = ci,
    geom = "line",
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 2.5,
    show.legend = FALSE
  ) +
  theme_bw() + labs("") +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 10, face = "bold", color = "black",
                               angle = 40,
                               vjust = 0.6),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(
    labels = c(
      "Mimosa\nclade (155)",
      "Stryphnod.\nclade (31)"
    )
  ) +
  scale_y_continuous(labels = label_number(accuracy = 0.1 , decimal.mark =
                                             ".")) +
  guides(fill = guide_legend(nrow = 1)) +
  ggtitle("A") +
  theme(plot.title = element_text(
    size = 14,
    hjust = 0,
    vjust = 1
  ))

#### Ecological groups
E1 <- c(rep("Desert", bootstraps))
E2 <- c(rep("Flooded grassland", bootstraps))
E3 <- c(rep("Savanna", bootstraps))
E4 <- c(rep("SDTF", bootstraps))
E5 <- c(rep("Temperate", bootstraps))
E6 <- c(rep("Tropical rainforest", bootstraps))
ecological_plot_disp <-
  vector(mode = "list", length = length(ecological_sv))

sum_Desert <- as.vector(ecological_sv$disparity$Desert[[2]])
sum_Flooded_grassland <-
  as.vector(ecological_sv$disparity$Flooded_grassland[[2]])
sum_Savanna <- as.vector(ecological_sv$disparity$Savanna[[2]])
sum_SDTF <- as.vector(ecological_sv$disparity$SDTF[[2]])
sum_Temperate <-
  as.vector(ecological_sv$disparity$Temperate[[2]])
sum_Tropical_rainforest <-
  as.vector(ecological_sv$disparity$Tropical_rainforest[[2]])

SV_Desert <- data.frame(E1, sum_Desert)
colnames(SV_Desert) <- c("Group", "Sum_of_variances")
SV_Flooded_grassland <- data.frame(E2, sum_Flooded_grassland)
colnames(SV_Flooded_grassland) <- c("Group", "Sum_of_variances")
SV_Savanna <- data.frame(E3, sum_Savanna)
colnames(SV_Savanna) <- c("Group", "Sum_of_variances")
SV_SDTF <- data.frame(E4, sum_SDTF)
colnames(SV_SDTF) <- c("Group", "Sum_of_variances")
SV_Temperate <- data.frame(E5, sum_Temperate)
colnames(SV_Temperate) <- c("Group", "Sum_of_variances")
SV_Tropical_rainforest <- data.frame(E6, sum_Tropical_rainforest)
colnames(SV_Tropical_rainforest) <- c("Group", "Sum_of_variances")
ecological_to_plot <-
  rbind(
    SV_Desert,
    SV_Flooded_grassland,
    SV_Savanna,
    SV_SDTF,
    SV_Temperate,
    SV_Tropical_rainforest
  )

ecological_to_plot$Group <-
  factor(
    ecological_to_plot$Group,
    levels = c(
      "Desert",
      "Flooded grassland",
      "Savanna",
      "SDTF",
      "Temperate",
      "Tropical rainforest"
    )
  )

plot_ecogroups_sv <-
  ggplot(
    ecological_to_plot,
    aes(
      x = Group,
      y = Sum_of_variances,
      colour = Group,
      fill = Group
    ),
    alpha = 0.9
  ) +
  geom_boxplot(alpha = 0.5, show.legend = F) + labs(x = "", y = "Sum of variances") +
  scale_fill_manual(values = as.character(paletteer_d(
    "ggthemes::excel_Office_2007_2010"))[1:6]
  ) +
  scale_colour_manual(values = as.character(paletteer_d(
    "ggthemes::excel_Office_2007_2010"))[1:6]
  ) +
  stat_summary(
    fun = ci,
    geom = "line",
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 2.5,
    show.legend = FALSE
  ) +
  theme_bw() + labs("") +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(size = 10, face = "bold", color = "black",
                               angle = 40,
                               vjust = 0.6),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(
    labels = c(
      "Desert (3)",
      "Flooded\ngrassland (9)",
      "Savanna (83)",
      "SDTF (66)",
      "Temperate (7)",
      "Tropical\nrainforest (59)"
    )
  ) +
  scale_y_continuous(labels = label_number(accuracy = 0.1 , decimal.mark =
                                             ".")) +
  guides(fill = guide_legend(nrow = 2)) +
  ggtitle("C") +
  theme(plot.title = element_text(
    size = 14,
    hjust = 0,
    vjust = 1
  ))

#### Ecological lineages
EL1 <- c(rep("Mimosa Savanna", bootstraps))
EL2 <- c(rep("Mimosa SDTF 1", bootstraps))
EL3 <- c(rep("Mimosa SDTF 2", bootstraps))
EL4 <- c(rep("Stryphnodendron SDTF", bootstraps))
EL5 <- c(rep("Stryphnodendron Tropical", bootstraps))
ecolineages_plot_disp <-
  vector(mode = "list", length = length(ecolineages_sv))

sum_Mimosa_Savanna_2 <-
  as.vector(ecolineages_sv$disparity$Mimosa_Savanna_2[[2]])
sum_Mimosa_SDTF_1 <-
  as.vector(ecolineages_sv$disparity$Mimosa_SDTF_1[[2]])
sum_Mimosa_SDTF_2 <-
  as.vector(ecolineages_sv$disparity$Mimosa_SDTF_2[[2]])
sum_Stryphnodendron_sdtf <-
  as.vector(ecolineages_sv$disparity$Stryphnod_SDTF[[2]])
sum_Stryphnodendron_tropical <-
  as.vector(ecolineages_sv$disparity$Stryphnod_Tropical[[2]])

SV_Mimosa_Savanna_2 <- data.frame(EL1, sum_Mimosa_Savanna_2)
colnames(SV_Mimosa_Savanna_2) <- c("Group", "Sum_of_variances")
SV_Mimosa_SDTF_1 <- data.frame(EL2, sum_Mimosa_SDTF_1)
colnames(SV_Mimosa_SDTF_1) <- c("Group", "Sum_of_variances")
SV_Mimosa_SDTF_2 <- data.frame(EL3, sum_Mimosa_SDTF_2)
colnames(SV_Mimosa_SDTF_2) <- c("Group", "Sum_of_variances")
SV_Stryphnodendron_sdtf <- data.frame(EL4, sum_Stryphnodendron_sdtf)
colnames(SV_Stryphnodendron_sdtf) <- c("Group", "Sum_of_variances")
SV_Stryphnodendron_tropical <-
  data.frame(EL5, sum_Stryphnodendron_tropical)
colnames(SV_Stryphnodendron_tropical) <- c("Group", "Sum_of_variances")
ecolineages_to_plot <-
  rbind(
    SV_Mimosa_Savanna_2,
    SV_Mimosa_SDTF_1,
    SV_Mimosa_SDTF_2,
    SV_Stryphnodendron_sdtf,
    SV_Stryphnodendron_tropical
  )

ecolineages_to_plot$Group <-
  factor(
    ecolineages_to_plot$Group,
    levels = c(
      "Mimosa Savanna",
      "Mimosa SDTF 1",
      "Mimosa SDTF 2",
      "Stryphnodendron SDTF",
      "Stryphnodendron Tropical"
    )
  )

plot_ecolineages_sv <-
  ggplot(
    ecolineages_to_plot,
    aes(
      x = Group,
      y = Sum_of_variances,
      colour = Group,
      fill = Group
    ),
    alpha = 0.9
  ) +
  geom_boxplot(alpha = 0.5, show.legend = F) + labs(x = "", y = "") +
  scale_fill_manual(values = as.character(paletteer_d(
    "ggthemes::excel_Office_2007_2010"))[1:5]
  ) +
  scale_colour_manual(values = as.character(paletteer_d(
    "ggthemes::excel_Office_2007_2010"))[1:5]
  ) +
  stat_summary(
    fun = ci,
    geom = "line",
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 2.5,
    show.legend = FALSE
  ) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 10, face = "bold", color = "black",
                               angle = 40,
                               vjust = 0.6),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(
    labels = c(
      "Mimosa\nSavanna (46)",
      "Mimosa\nSDTF 1 (33)",
      "Mimosa\nSDTF 2 (6)",
      "Stryphnod.\nSDTF (4)",
      "Stryphnod.\nTropical (36)"
    )
  ) +
  scale_y_continuous(labels = label_number(accuracy = 0.1 , decimal.mark =
                                             ".")) +
  guides(fill = guide_legend(nrow = 2)) +
  ggtitle("D") +
  theme(plot.title = element_text(
    size = 14,
    hjust = 0,
    vjust = 1
  ))

# SV Boxplots - rarefied (Fig. 3) =============================================

### This part of the script was based on the code of Casali et al. (2023)
### available at https://doi.org/10.5281/zenodo.7240006

### required functions

if (rarefaction == "min") {
  
  ci <- function(x)
  {
    exit <- quantile(x, probs = c(0.025, 0.975))
    names(exit) <- c("ymin", "ymax")
    return(exit)
  }
  
  median <- function(x)
  {
    exit <- quantile(x, probs = c(0.5))
    names(exit) <- c("y")
    return(exit)
  }
  
  #### Genera
  F1 <- c(rep("Gwilymia", bootstraps))
  F2 <- c(rep("Marlimorimia", bootstraps))
  F3 <- c(rep("Mimosa", bootstraps))
  F4 <- c(rep("Parapiptadenia", bootstraps))
  F5 <- c(rep("Piptadenia", bootstraps))
  F6 <- c(rep("Pityrocarpa", bootstraps))
  F7 <- c(rep("Stryphnodendron", bootstraps))
  genera_plot_disp <- vector(mode = "list", length = length(genera_sv))
  
  sum_Gwilymia <- as.vector(genera_sv$disparity$Gwilymia[[3]])
  sum_Marlimorimia <-
    as.vector(genera_sv$disparity$Marlimorimia[[2]])
  sum_Mimosa <- as.vector(genera_sv$disparity$Mimosa[[3]])
  sum_Parapiptadenia <-
    as.vector(genera_sv$disparity$Parapiptadenia[[3]])
  sum_Piptadenia <- as.vector(genera_sv$disparity$Piptadenia[[3]])
  sum_Pityrocarpa <-
    as.vector(genera_sv$disparity$Pityrocarpa[[3]])
  sum_Stryphnodendron <-
    as.vector(genera_sv$disparity$Stryphnodendron[[3]])
  
  SV_Gwilymia <- data.frame(F1, sum_Gwilymia)
  colnames(SV_Gwilymia) <- c("Group", "Sum_of_variances")
  SV_Marlimorimia <- data.frame(F2, sum_Marlimorimia)
  colnames(SV_Marlimorimia) <- c("Group", "Sum_of_variances")
  SV_Mimosa <- data.frame(F3, sum_Mimosa)
  colnames(SV_Mimosa) <- c("Group", "Sum_of_variances")
  SV_Parapiptadenia <- data.frame(F4, sum_Parapiptadenia)
  colnames(SV_Parapiptadenia) <- c("Group", "Sum_of_variances")
  SV_Piptadenia <- data.frame(F5, sum_Piptadenia)
  colnames(SV_Piptadenia) <- c("Group", "Sum_of_variances")
  SV_Pityrocarpa <- data.frame(F6, sum_Pityrocarpa)
  colnames(SV_Pityrocarpa) <- c("Group", "Sum_of_variances")
  SV_Stryphnodendron <- data.frame(F7, sum_Stryphnodendron)
  colnames(SV_Stryphnodendron) <- c("Group", "Sum_of_variances")
  genera_to_plot <-
    rbind(
      SV_Gwilymia,
      SV_Marlimorimia,
      SV_Mimosa,
      SV_Parapiptadenia,
      SV_Piptadenia,
      SV_Pityrocarpa,
      SV_Stryphnodendron
    )
  
  genera_to_plot$Group <-
    factor(
      genera_to_plot$Group,
      levels = c(
        "Gwilymia",
        "Marlimorimia",
        "Mimosa",
        "Parapiptadenia",
        "Piptadenia",
        "Pityrocarpa",
        "Stryphnodendron"
      )
    )
  
  plot_genera_sv_rare <- ggplot(genera_to_plot,
                                aes(
                                  x = Group,
                                  y = Sum_of_variances,
                                  colour = Group,
                                  fill = Group
                                ),
                                alpha = 0.9) +
    geom_boxplot(alpha = 0.5, show.legend = F) + labs(x = "", y = "") +
    scale_fill_manual(values = as.character(paletteer_d(
      "ggthemes::excel_Office_2007_2010"))[1:7]
    ) +
    scale_colour_manual(values = as.character(paletteer_d(
      "ggthemes::excel_Office_2007_2010"))[1:7]
    ) +
    stat_summary(
      fun = ci,
      geom = "line",
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      size = 2.5,
      show.legend = FALSE
    ) +
    theme_bw() + labs("") +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 10, face = "bold", color = "black",
                                 angle = 40, 
                                 vjust = 0.6),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_x_discrete(labels = c(
      "Gwilymia (7)",
      "Marlimorimia (4)",
      "Mimosa (135)",
      "Parapiptadenia (5)",
      "Piptadenia (18)",
      "Pityrocarpa (7)",
      "Stryphnod. (6)"
    )) +
    scale_y_continuous(labels = label_number(accuracy = 0.1 , decimal.mark =
                                               "."))  +
    guides(fill = guide_legend(nrow = 2)) +
    ggtitle("B") +
    theme(plot.title = element_text(
      size = 14,
      hjust = 0
    ))
  
  #### Clades
  C1 <- c(rep("Mimosa clade", bootstraps))
  C2 <- c(rep("Stryphnodendron clade", bootstraps))
  clades_plot_disp <- vector(mode = "list", length = length(clades_sv))
  
  sum_Mimosa_clade <-
    as.vector(clades_sv$disparity$Mimosa_clade[[3]])
  sum_Stryphnodendron_clade <-
    as.vector(clades_sv$disparity$Stryphnodendron_clade[[2]])
  
  SV_Mimosa_clade <- data.frame(C1, sum_Mimosa_clade)
  colnames(SV_Mimosa_clade) <- c("Group", "Sum_of_variances")
  SV_Stryphnodendron_clade <- data.frame(C2, sum_Stryphnodendron_clade)
  colnames(SV_Stryphnodendron_clade) <- c("Group", "Sum_of_variances")
  clades_to_plot <- rbind(SV_Mimosa_clade, SV_Stryphnodendron_clade)
  
  clades_to_plot$Group <-
    factor(clades_to_plot$Group,
           levels = c("Mimosa clade", "Stryphnodendron clade"))
  
  plot_clades_sv_rare <-
    ggplot(clades_to_plot,
           aes(
             x = Group,
             y = Sum_of_variances,
             colour = Group,
             fill = Group
           ),
           alpha = 0.9) +
    geom_boxplot(alpha = 0.5, show.legend = F) + labs(x = "", y = "Sum of variances") +
    scale_fill_manual(values = as.character(paletteer_d(
      "ggthemes::excel_Office_2007_2010"))[1:2]
    ) +
    scale_colour_manual(values = as.character(paletteer_d(
      "ggthemes::excel_Office_2007_2010"))[1:2]
    ) +
    stat_summary(
      fun = ci,
      geom = "line",
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      size = 2.5,
      show.legend = FALSE
    ) +
    theme_bw() + labs("") +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(size = 10, face = "bold", color = "black",
                                 angle = 40,
                                 vjust = 0.6),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
    ) +
    scale_x_discrete(
      labels = c(
        "Mimosa\nclade (155)",
        "Stryphnod.\nclade (31)"
      )
    ) +
    scale_y_continuous(labels = label_number(accuracy = 0.1 , decimal.mark =
                                               ".")) +
    guides(fill = guide_legend(nrow = 1)) +
    ggtitle("A") +
    theme(plot.title = element_text(
      size = 14,
      hjust = 0,
      vjust = 1
    ))
  
  #### Ecological groups
  E1 <- c(rep("Desert", bootstraps))
  E2 <- c(rep("Flooded grassland", bootstraps))
  E3 <- c(rep("Savanna", bootstraps))
  E4 <- c(rep("SDTF", bootstraps))
  E5 <- c(rep("Temperate", bootstraps))
  E6 <- c(rep("Tropical rainforest", bootstraps))
  ecological_plot_disp <-
    vector(mode = "list", length = length(ecological_sv))
  
  sum_Desert <- as.vector(ecological_sv$disparity$Desert[[2]])
  sum_Flooded_grassland <-
    as.vector(ecological_sv$disparity$Flooded_grassland[[3]])
  sum_Savanna <- as.vector(ecological_sv$disparity$Savanna[[3]])
  sum_SDTF <- as.vector(ecological_sv$disparity$SDTF[[3]])
  sum_Temperate <-
    as.vector(ecological_sv$disparity$Temperate[[3]])
  sum_Tropical_rainforest <-
    as.vector(ecological_sv$disparity$Tropical_rainforest[[3]])
  
  SV_Desert <- data.frame(E1, sum_Desert)
  colnames(SV_Desert) <- c("Group", "Sum_of_variances")
  SV_Flooded_grassland <- data.frame(E2, sum_Flooded_grassland)
  colnames(SV_Flooded_grassland) <- c("Group", "Sum_of_variances")
  SV_Savanna <- data.frame(E3, sum_Savanna)
  colnames(SV_Savanna) <- c("Group", "Sum_of_variances")
  SV_SDTF <- data.frame(E4, sum_SDTF)
  colnames(SV_SDTF) <- c("Group", "Sum_of_variances")
  SV_Temperate <- data.frame(E5, sum_Temperate)
  colnames(SV_Temperate) <- c("Group", "Sum_of_variances")
  SV_Tropical_rainforest <- data.frame(E6, sum_Tropical_rainforest)
  colnames(SV_Tropical_rainforest) <- c("Group", "Sum_of_variances")
  ecological_to_plot <-
    rbind(
      SV_Desert,
      SV_Flooded_grassland,
      SV_Savanna,
      SV_SDTF,
      SV_Temperate,
      SV_Tropical_rainforest
    )
  
  ecological_to_plot$Group <-
    factor(
      ecological_to_plot$Group,
      levels = c(
        "Desert",
        "Flooded grassland",
        "Savanna",
        "SDTF",
        "Temperate",
        "Tropical rainforest"
      )
    )
  
  plot_ecogroups_sv_rare <-
    ggplot(
      ecological_to_plot,
      aes(
        x = Group,
        y = Sum_of_variances,
        colour = Group,
        fill = Group
      ),
      alpha = 0.9
    ) +
    geom_boxplot(alpha = 0.5, show.legend = F) + labs(x = "", y = "Sum of variances") +
    scale_fill_manual(values = as.character(paletteer_d(
      "ggthemes::excel_Office_2007_2010"))[1:6]
    ) +
    scale_colour_manual(values = as.character(paletteer_d(
      "ggthemes::excel_Office_2007_2010"))[1:6]
    ) +
    stat_summary(
      fun = ci,
      geom = "line",
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      size = 2.5,
      show.legend = FALSE
    ) +
    theme_bw() + labs("") +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(size = 10, face = "bold", color = "black",
                                 angle = 40,
                                 vjust = 0.6),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_x_discrete(
      labels = c(
        "Desert (3)",
        "Flooded\ngrassland (9)",
        "Savanna (83)",
        "SDTF (66)",
        "Temperate (7)",
        "Tropical\nrainforest (59)"
      )
    ) +
    scale_y_continuous(labels = label_number(accuracy = 0.1 , decimal.mark =
                                               ".")) +
    guides(fill = guide_legend(nrow = 2)) +
    ggtitle("C") +
    theme(plot.title = element_text(
      size = 14,
      hjust = 0,
      vjust = 1
    ))
  
  #### Ecological lineages
  EL1 <- c(rep("Mimosa Savanna", bootstraps))
  EL2 <- c(rep("Mimosa SDTF 1", bootstraps))
  EL3 <- c(rep("Mimosa SDTF 2", bootstraps))
  EL4 <- c(rep("Stryphnodendron SDTF", bootstraps))
  EL5 <- c(rep("Stryphnodendron Tropical", bootstraps))
  ecolineages_plot_disp <-
    vector(mode = "list", length = length(ecolineages_sv))
  
  sum_Mimosa_Savanna_2 <-
    as.vector(ecolineages_sv$disparity$Mimosa_Savanna_2[[3]])
  sum_Mimosa_SDTF_1 <-
    as.vector(ecolineages_sv$disparity$Mimosa_SDTF_1[[3]])
  sum_Mimosa_SDTF_2 <-
    as.vector(ecolineages_sv$disparity$Mimosa_SDTF_2[[3]])
  sum_Stryphnodendron_sdtf <-
    as.vector(ecolineages_sv$disparity$Stryphnod_SDTF[[2]])
  sum_Stryphnodendron_tropical <-
    as.vector(ecolineages_sv$disparity$Stryphnod_Tropical[[3]])
  
  SV_Mimosa_Savanna_2 <- data.frame(EL1, sum_Mimosa_Savanna_2)
  colnames(SV_Mimosa_Savanna_2) <- c("Group", "Sum_of_variances")
  SV_Mimosa_SDTF_1 <- data.frame(EL2, sum_Mimosa_SDTF_1)
  colnames(SV_Mimosa_SDTF_1) <- c("Group", "Sum_of_variances")
  SV_Mimosa_SDTF_2 <- data.frame(EL3, sum_Mimosa_SDTF_2)
  colnames(SV_Mimosa_SDTF_2) <- c("Group", "Sum_of_variances")
  SV_Stryphnodendron_sdtf <- data.frame(EL4, sum_Stryphnodendron_sdtf)
  colnames(SV_Stryphnodendron_sdtf) <- c("Group", "Sum_of_variances")
  SV_Stryphnodendron_tropical <-
    data.frame(EL5, sum_Stryphnodendron_tropical)
  colnames(SV_Stryphnodendron_tropical) <- c("Group", "Sum_of_variances")
  ecolineages_to_plot <-
    rbind(
      SV_Mimosa_Savanna_2,
      SV_Mimosa_SDTF_1,
      SV_Mimosa_SDTF_2,
      SV_Stryphnodendron_sdtf,
      SV_Stryphnodendron_tropical
    )
  
  ecolineages_to_plot$Group <-
    factor(
      ecolineages_to_plot$Group,
      levels = c(
        "Mimosa Savanna",
        "Mimosa SDTF 1",
        "Mimosa SDTF 2",
        "Stryphnodendron SDTF",
        "Stryphnodendron Tropical"
      )
    )
  
  plot_ecolineages_sv_rare <-
    ggplot(
      ecolineages_to_plot,
      aes(
        x = Group,
        y = Sum_of_variances,
        colour = Group,
        fill = Group
      ),
      alpha = 0.9
    ) +
    geom_boxplot(alpha = 0.5, show.legend = F) + labs(x = "", y = "") +
    scale_fill_manual(values = as.character(paletteer_d(
      "ggthemes::excel_Office_2007_2010"))[1:5]
    ) +
    scale_colour_manual(values = as.character(paletteer_d(
      "ggthemes::excel_Office_2007_2010"))[1:5]
    ) +
    stat_summary(
      fun = ci,
      geom = "line",
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      size = 2.5,
      show.legend = FALSE
    ) +
    theme_bw() +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 10, face = "bold", color = "black",
                                 angle = 40,
                                 vjust = 0.6),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_x_discrete(
      labels = c(
        "Mimosa\nSavanna (46)",
        "Mimosa\nSDTF 1 (33)",
        "Mimosa\nSDTF 2 (6)",
        "Stryphnod.\nSDTF (4)",
        "Stryphnod.\nTropical (36)"
      )
    ) +
    scale_y_continuous(labels = label_number(accuracy = 0.1 , decimal.mark =
                                               ".")) +
    guides(fill = guide_legend(nrow = 2)) +
    ggtitle("D") +
    theme(plot.title = element_text(
      size = 14,
      hjust = 0,
      vjust = 1
    ))
  
}

# SV correlation matrices (Fig. 4) ============================================

library(reshape2)

summary(summary_stats_sv)

### Split the 'Groups' column into two new columns 'Taxa1' and 'Taxa2'

pairwise_dat <- transform(summary_stats_sv[, c(2,5)],
                          Taxon1 = sub(" : .*", "", Groups),
                          Taxon2 = sub(".* : ", "", Groups))
### Genus

create_bhatt_matrix <- function(pairwise_data) {
  
  # Step 1: Get unique taxa
  taxa <- unique(c(pairwise_data$Taxon1, pairwise_data$Taxon2))
  
  # Step 2: Replace underscores with spaces (no capitalization changes)
  taxa <- gsub("_", " ", taxa)
  
  # Step 3: Initialize an empty matrix
  bhatt_matrix <- matrix(NA, nrow = length(taxa), ncol = length(taxa),
                         dimnames = list(taxa, taxa))
  
  # Step 4: Fill the matrix with Bhattacharyya Coefficients
  for (i in 1:nrow(pairwise_data)) {
    taxon1 <- pairwise_data$Taxon1[i]
    taxon2 <- pairwise_data$Taxon2[i]
    coefficient <- pairwise_data$`Bhattacharyya Coefficient`[i]
    
    # Replace underscores with spaces (no capitalization changes)
    taxon1 <- gsub("_", " ", taxon1)
    taxon2 <- gsub("_", " ", taxon2)
    
    if (taxon1 %in% taxa && taxon2 %in% taxa) {
      bhatt_matrix[taxon1, taxon2] <- coefficient
      bhatt_matrix[taxon2, taxon1] <- coefficient
    }
  }
  
  # Step 5: Replace any remaining NA values with 0
  bhatt_matrix[is.na(bhatt_matrix)] <- 0
  
  # Step 6: Function to reorder the correlation matrix
  reorder_cormat <- function(cormat) {
    dd <- as.dist((1 - cormat) / 2)
    hc <- hclust(dd)
    cormat <- cormat[hc$order, hc$order]
    return(cormat)
  }
  
  # Step 7: Function to get the lower triangle of the matrix
  get_lower_tri <- function(cormat) {
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  
  # Step 8: Reorder the Bhattacharyya matrix
  bhatt_matrix <- reorder_cormat(bhatt_matrix)
  
  # Step 9: Get the lower triangular matrix
  lower_tri <- get_lower_tri(bhatt_matrix)
  
  # Step 10: Melt the matrix to a long format
  melted_matrix <- melt(lower_tri, na.rm = TRUE)
  
  # Step 11: Reverse the factor levels of Var1
  melted_matrix$Var1 <- factor(melted_matrix$Var1, 
                               levels = rev(levels(melted_matrix$Var1)))
  
  # Return the melted matrix
  return(melted_matrix)
}

bhatt_pairwise_genus <- summary_stats_sv[c(1:21), c(2,5)] %>%
  separate(Groups, into = c("Taxon1", "Taxon2"), sep = " : ")

melted_genus <- create_bhatt_matrix(bhatt_pairwise_genus)

plot_genus_corr <- ggplot(data = melted_genus, aes(x=Var1, y=Var2, 
                                                   fill=value)) + 
  geom_tile(color = "white")+
  scale_fill_viridis_c(option = "viridis", direction = 1, 
                       name="Bhattacaryya Coefficient", limits = c(0, 1)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, color = "black",
                                   size = 8, hjust = 1))+
  coord_fixed()

### Clades

bhatt_pairwise_clades <- summary_stats_sv[c(22), c(2,5)] %>%
  separate(Groups, into = c("Taxon1", "Taxon2"), sep = " : ")

melted_clades <- create_bhatt_matrix(bhatt_pairwise_clades)

plot_clades_corr <- ggplot(data = melted_clades, aes(x=Var1, y=Var2, 
                                                     fill=value)) + 
  geom_tile(color = "white")+
  scale_fill_viridis_c(option = "viridis", direction = 1, 
                       name="Bhattacaryya Coefficient", limits = c(0, 1)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, color = "black",
                                   size = 8, hjust = 1))+
  coord_fixed()

### Biomes

bhatt_pairwise_biomes <- summary_stats_sv[c(23:37), c(2,5)] %>%
  separate(Groups, into = c("Taxon1", "Taxon2"), sep = " : ")

melted_biomes <- create_bhatt_matrix(bhatt_pairwise_biomes)

plot_biomes_corr <- ggplot(data = melted_biomes, aes(x=Var1, y=Var2, 
                                                     fill=value)) + 
  geom_tile(color = "white")+
  scale_fill_viridis_c(option = "viridis", direction = 1, 
                       name="Bhattacaryya Coefficient", limits = c(0, 1)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, color = "black",
                                   size = 8, hjust = 1))+
  coord_fixed()

### Ecolineages

bhatt_pairwise_ecolineages <- summary_stats_sv[c(38:47), c(2,5)] %>%
  separate(Groups, into = c("Taxon1", "Taxon2"), sep = " : ")

melted_ecolineages <- create_bhatt_matrix(bhatt_pairwise_ecolineages)

plot_ecolineages_corr <- ggplot(data = melted_ecolineages, aes(x=Var1, y=Var2, 
                                                               fill=value)) + 
  geom_tile(color = "white")+
  scale_fill_viridis_c(option = "viridis", direction = 1, 
                       name="Bhattacaryya Coefficient", limits = c(0, 1)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, color = "black",
                                   size = 8, hjust = 1))+
  coord_fixed()

# Save plots ==================================================================

# Phylomorphospace

pdf("output/fig2_phylomorphospace.pdf")

phylomorphospace

dev.off()

# Morphospace + all vegetation types

pdf("output/fig5_morphospace-eco.pdf")

morphospace_eco

dev.off()

# SV

pdf("output/suppfig3_sum_of_variances.pdf", width = 10, height = 10)

(plot_clades_sv + plot_genera_sv) /
  (plot_ecogroups_sv + plot_ecolineages_sv)

dev.off()

# SV - rarefied

if (rarefaction == "min") {
  
  pdf("output/fig3_sum_of_variances_rarefied.pdf", width = 10, height = 10)
  
  sv_rarefied <- (plot_clades_sv_rare + plot_genera_sv_rare) / 
    (plot_ecogroups_sv_rare + plot_ecolineages_sv_rare)
  print(sv_rarefied)
  
  dev.off()
  
}

# Correlation matrices

pdf("output/fig4_corr_matrices.pdf", width = 9, height = 9)

plot_genus_corr + plot_biomes_corr + plot_ecolineages_corr + 
  plot_layout(nrow = 2, byrow = F) +
  plot_annotation(tag_levels = 'A')

dev.off()
