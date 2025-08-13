#
# Delimitation of different biome lineages
#
# Rafael F. Barduzzi
#
# 2024
#

# Libraries and data ==========================================================

require(phytools) #1.2-0
require(ape) #5.7-1
require(plyr)#1.8.8

# Data ========================================================================

eco_data <- read.csv("data/eco_data.csv",
                     na.strings = "?",
                     stringsAsFactors = F)

tree <- read.tree("data/mimo_stryph_phylo.tre")

pollen_matrix <- read.csv("data/pollen_data.csv")

tree <- force.ultrametric(tree, method = "extend")


# Prepare data ================================================================

# set seed for replication

set.seed(129)

## create output directory
dir.create("output", showWarnings = F)

# Pruning ecological data
traits <- eco_data[eco_data$taxon %in% tree$tip.label, ]
rownames(traits) <- traits$taxon

# Pruning tree
tree <- drop.tip(tree, setdiff(tree$tip.label, traits$taxon))

# Ancestral biomes estimate (Supplementary Fig. S4) ===========================

#The data.frame needs to be transformed into a named vector
#containing tip labels as names and character states as factors.
plyr::count(traits$biome)
cat.trait <- factor(
  traits$biome,
  levels = c(
    "0",
    "1",
    "2",
    "3",
    "4",
    "5",
    "1&2",
    "1&2&3",
    "1&2&3&5",
    "1&3",
    "1&5",
    "2&3",
    "2&3&5",
    "2&5",
    "3&4",
    "3&5",
    "missing"
  )
)

cat.trait[is.na(cat.trait)] <- "missing"
names(cat.trait) <- rownames(traits)

#The named vector, then, needs to be transformed into a binary matrix
#with tip labels as row names and character states as columns
bin.matrix <- to.matrix(cat.trait, levels(cat.trait))

#Dealing with polymorphisms
##When taxa with polymorphisms are present, we must assign a value
#equal to 1/(number of states) for each possible state and then remove
#the column that refers to the polymorphic state.

# Polymorphism 1&2
polymorph7 <- rownames(bin.matrix)[bin.matrix[, c(7)] == 1]
bin.matrix[polymorph7, c(2, 3)] <- 1 / 2

#Polymorphism 1&2&3
polymorph8 <- rownames(bin.matrix)[bin.matrix[, c(8)] == 1]
bin.matrix[polymorph8, c(2, 3, 4)] <- 1 / 3

#Polymorphism 1&2&3&5
polymorph9 <- rownames(bin.matrix)[bin.matrix[, c(9)] == 1]
bin.matrix[polymorph9, c(2, 3, 4, 6)] <- 1 / 4

# Polymorphism 1&3
polymorph10 <- rownames(bin.matrix)[bin.matrix[, c(10)] == 1]
bin.matrix[polymorph10, c(2, 4)] <- 1 / 2

#Polymorphism 1&5
polymorph11 <- rownames(bin.matrix)[bin.matrix[, c(11)] == 1]
bin.matrix[polymorph11, c(2, 6)] <- 1 / 2

#Polymorphism 2&3
polymorph12 <- rownames(bin.matrix)[bin.matrix[, c(12)] == 1]
bin.matrix[polymorph12, c(3, 4)] <- 1 / 2

#Polymorphism 2&3&5
polymorph13 <- rownames(bin.matrix)[bin.matrix[, c(13)] == 1]
bin.matrix[polymorph13, c(3, 4, 6)] <- 1 / 3

#Polymorphism 2&5
polymorph14 <- rownames(bin.matrix)[bin.matrix[, c(14)] == 1]
bin.matrix[polymorph14, c(3, 6)] <- 1 / 2

#Polymorphism 3&4
polymorph15 <- rownames(bin.matrix)[bin.matrix[, c(15)] == 1]
bin.matrix[polymorph15, c(4, 5)] <- 1 / 2

#Polymorphism 3&5
polymorph16 <- rownames(bin.matrix)[bin.matrix[, c(16)] == 1]
bin.matrix[polymorph16, c(4, 6)] <- 1 / 2

#Removing columns indicating polymorphism
bin.matrix <- bin.matrix[,-c(7, 8, 9, 10, 11, 12, 13, 14, 15, 16)]

#Which taxa show missing data?
missing.data <- names(which(bin.matrix[, "missing"] == 1))

#For missing data, we must assign a prior probability distribution on the tips
#that is flat across all possible states.
##Removing the column 'missing'
bin.matrix <- bin.matrix[,-c(ncol(bin.matrix))]
##Assigning 1/(number of states) for all possible states
bin.matrix[row.names(bin.matrix) %in% missing.data,] <-
  1 / ncol(bin.matrix)

#Changing state names to biomes
colnames(bin.matrix) <-
  c("Desert",
    "Flooded grassland",
    "Savanna",
    "SDTF",
    "Temperate",
    "Tropical rainforest")

#Running stochastic mapping (this can take a while)
trees <- make.simmap(tree, bin.matrix, model = "ER", nsim = 100)
obj <- summary(trees, plot = FALSE)

# putting the columns in the same order of the tip.labels, so we can plot polymophisms
bin.matrix <-
  bin.matrix[match(tree$tip.label, rownames(bin.matrix)), ]

#putting the columns of obj$ace in the same order of bin.matrix
obj$ace <- obj$ace[, colnames(bin.matrix)]

#plotting
cols <-
  setNames(palette()[1:length(colnames(bin.matrix))], colnames(bin.matrix))

par(lwd = 0.1)
plotTree(
  tree,
  type = "phylogram",
  fsize = 0.2,
  cex = 0.2,
  lwd = 1,
  ftype = "i"
)
nodelabels(pie = obj$ace,
           piecol = cols,
           cex = 0.2)
tiplabels(pie = bin.matrix,
          piecol = cols,
          cex = 0.2)
add.simmap.legend(
  colors = cols,
  x = -15,
  y = 10,
  prompt = FALSE,
  fsize = 0.5
)

#preparing for pollen data indication

pollen_presence <- data.frame(row.names = rownames(bin.matrix))
pollen_presence$species <- rownames(bin.matrix)
pollen_presence$pollen_data <-
  ifelse(pollen_presence$species %in% pollen_matrix$name_phylogeny,
         1,
         0)
pollen_presence <- pollen_presence[c(2)]

## Get clades and polen data indication

#delimitating clades based on previous plot

#mimosa savanna

nodes_sav <-
  c(getMRCA(tree, c("Mimosa_filipes", "Mimosa_gatesiae")),
    getMRCA(tree, c(
      "Mimosa_minarum", "Mimosa_longepedunculata"
    )))

labels_sav <- paste("Mimosa Savanna", 1:length(nodes_sav))

#mimosa sdtf

nodes_sdtf <-
  c(getMRCA(tree, c(
    "Piptadenia_floribunda", "Mimosa_longepedunculata"
  )),
  getMRCA(tree, c("Mimosa_minarum", "Mimosa_cordistipula")),
  getMRCA(tree, c("Mimosa_woodii", "Mimosa_deamii")))

labels_sdtf <- paste("Mimosa SDTF", 1:length(nodes_sdtf))

#mimosa tropical

nodes_tropical <-
  c(getMRCA(tree, c(
    "Piptadenia_robusta", "Piptadenia_adiantoides"
  )),
  getMRCA(tree, c(
    "Mimosa_lepidophora", "Mimosa_guilandinae"
  )),
  getMRCA(tree, c("Mimosa_gymnas", "Mimosa_involucrata")),
  getMRCA(tree, c("Mimosa_taimbensis", "Mimosa_urticaria")))

labels_tropical <- paste("Mimosa Tropical", 1:length(nodes_tropical))

#mimosa flooded

nodes_flooded <-
  c(getMRCA(tree, c(
    "Mimosa_dormiens", "Mimosa_pigra_dehiscens"
  )))

labels_flooded <- "Mimosa Flooded"

#stryphnod tropical

nodes_tropical2 <-
  c(getMRCA(
    tree,
    c("Marlimorimia_psilostachya", "Microlobius_foetidus")
  ))

labels_tropical2 <- "Stryphnodendron Tropical"

#stryphnod sdtf

nodes_sdtf2 <-
  c(getMRCA(tree, c(
    "Pityrocarpa_brenanii", "Pityrocarpa_moniliformis"
  )))

labels_sdtf2 <- "Stryphnodendron SDTF"

#stryphnod savanna

nodes_sav2 <-
  c(getMRCA(
    tree,
    c("Stryphnodendron_velutinum", "Stryphnodendron_adstringens")
  ))

labels_sav2 <- "Stryphnodendron Savanna"

# Plot / Save result ==========================================================

pdf("output/suppfig4_simmap_fan_delimitated.pdf")

par(lwd = 0.1)
plotTree(
  tree,
  type = "fan",
  fsize = 0.2,
  cex = 0.2,
  lwd = 1,
  ftype = "i",
  offset = 5
)
nodelabels(pie = obj$ace,
           piecol = cols,
           cex = 0.15)
tiplabels(pie = bin.matrix,
          piecol = cols,
          cex = 0.14)

#mimosa sdtf
for (i in 1:length(nodes_sdtf))
  arc.cladelabels(
    tree,
    labels_sdtf[i],
    nodes_sdtf[i],
    ln.offset = if (labels_sdtf[i] %in% c("Mimosa SDTF 1"))
      1.21
    else
      1.23,
    lab.offset = if (labels_sdtf[i] %in% c("Mimosa SDTF 1"))
      1.30
    else
      1.25,
    cex = if (labels_sdtf[i] %in% c("Mimosa SDTF 1"))
      0.3
    else
      0.25,
    lwd = 2,
    col = "#2297e6ff",
    orientation = "curved",
    mark.node = FALSE
  )
#mimosa savanna
for (i in 1:length(nodes_sav))
  arc.cladelabels(
    tree,
    labels_sav[i],
    nodes_sav[i],
    ln.offset = 1.22,
    lab.offset = 1.25,
    cex = 0.25,
    lwd = 2,
    col = "#61d04fff",
    orientation = "curved",
    mark.node = FALSE
  )
#mimosa tropical
for (i in 1:length(nodes_tropical))
  arc.cladelabels(
    tree,
    labels_tropical[i],
    nodes_tropical[i],
    ln.offset = 1.22,
    lab.offset = 1.25,
    cex = 0.25,
    lwd = 2,
    col = "#cd0bbcff",
    orientation = "curved",
    mark.node = FALSE
  )
#mimosa flooded
for (i in 1:length(nodes_flooded))
  arc.cladelabels(
    tree,
    labels_flooded[i],
    nodes_flooded[i],
    ln.offset = 1.22,
    lab.offset = 1.25,
    cex = 0.25,
    lwd = 2,
    col = "#df536bff",
    orientation = "curved",
    mark.node = FALSE
  )
#stryphnod tropical
for (i in 1:length(nodes_tropical2))
  arc.cladelabels(
    tree,
    labels_tropical2[i],
    nodes_tropical2[i],
    ln.offset = 1.21,
    lab.offset = 1.30,
    cex = 0.3,
    lwd = 2,
    col = "#cd0bbcff",
    orientation = "curved",
    mark.node = FALSE
  )
#stryphnod sdtf
for (i in 1:length(nodes_sdtf2))
  arc.cladelabels(
    tree,
    labels_sdtf2[i],
    nodes_sdtf2[i],
    ln.offset = 1.22,
    lab.offset = 1.25,
    cex = 0.25,
    lwd = 2,
    col = "#2297e6ff",
    orientation = "curved",
    mark.node = FALSE
  )
#stryphnod savanna
for (i in 1:length(nodes_sav2))
  arc.cladelabels(
    tree,
    labels_sav2[i],
    nodes_sav2[i],
    ln.offset = 1.22,
    lab.offset = 1.25,
    cex = 0.25,
    lwd = 2,
    col = "#61d04fff",
    orientation = "curved",
    mark.node = FALSE
  )

#pollen data
tiplabels(pie = pollen_presence,
          piecol = c("black", "white"),
          cex = 0.05)

add.simmap.legend(
  colors = cols,
  x = -35,
  y = 35,
  prompt = FALSE,
  fsize = 0.5
)

dev.off()

