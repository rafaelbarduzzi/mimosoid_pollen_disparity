#
# Supplementary analysis - Sum of Ranges metric
#
# Rafael F. Barduzzi
#
# 07-23-25
#

# Data ========================================================================

# This script requires the objects created by running "disparity_analysis.R"

# Morphospace Occupancy Metrics ===============================================

## Size - SR ------------------------------------------------------------------
## sum of ranges - post-ordination metric

### set seed for replicability
set.seed(7) 

### This part of the script was based on the code of Casali et al. (2023)
### available at https://doi.org/10.5281/zenodo.7240006

bootstraps <- 1000

### SR - Genera
sub_genera <- custom.subsets(pcoa$vectors.cor, groups_genera)
bootstrapped_data <- boot.matrix(sub_genera, 
                                 bootstraps = bootstraps, 
                                 rarefaction = "min")
genera_sr <- dispRity(bootstrapped_data, 
                      metric = c(sum, ranges))
genera_sr_smry <- summary(genera_sr)
genera_sr_wilcox <-test.dispRity(genera_sr, 
                                 wilcox.test, 
                                 correction = "holm")
genera_sr_bhatt <- test.dispRity(genera_sr, bhatt.coeff)

### SR - Clades
sub_clades <- custom.subsets(pcoa$vectors.cor, groups_clades)
bootstrapped_data <- boot.matrix(sub_clades, 
                                 bootstraps = bootstraps, 
                                 rarefaction = "min")
clades_sr <- dispRity(bootstrapped_data, 
                      metric = c(sum, ranges))
clades_sr_smry <- summary(clades_sr)
clades_sr_wilcox <- test.dispRity(clades_sr, 
                                  wilcox.test, 
                                  correction = "holm")
clades_sr_bhatt <- test.dispRity(clades_sr, bhatt.coeff)

### SR - Ecological groups
sub_ecological <- custom.subsets(pcoa$vectors.cor, groups_ecological)
bootstrapped_data <- boot.matrix(sub_ecological, 
                                 bootstraps = bootstraps,
                                 rarefaction = "min")
ecological_sr <- dispRity(bootstrapped_data, 
                          metric = c(sum, ranges))
ecological_sr_smry <- summary(ecological_sr)
ecological_sr_wilcox <- test.dispRity(ecological_sr, 
                                      wilcox.test, 
                                      correction = "holm")
ecological_sr_bhatt <- test.dispRity(ecological_sr, bhatt.coeff)

### SR - Ecological lineages
sub_ecolineages <- custom.subsets(pcoa$vectors.cor, ecological_lineages)
bootstrapped_data <- boot.matrix(sub_ecolineages, 
                                 bootstraps = bootstraps,
                                 rarefaction = "min")
ecolineages_sr <- dispRity(bootstrapped_data, metric = c(sum, ranges))
ecolineages_sr_smry <- summary(ecolineages_sr)
ecolineages_sr_wilcox <- test.dispRity(ecolineages_sr, 
                                       wilcox.test, 
                                       correction = "holm")
ecolineages_sr_bhatt <- test.dispRity(ecolineages_sr, bhatt.coeff)

## Summary and export ---------------------------------------------------------

labels <- c(rep("Genera", 13), (rep("Clades", 3)), (rep("Biomes", 11)), 
            (rep("Ecological lineages", 9)))

### SR

genera_disp <- do.call(rbind, genera_sr_smry)
clade_disp <- do.call(rbind, clades_sr_smry)
ecological_disp <- do.call(rbind, ecological_sr_smry)
ecolineages_disp <- do.call(rbind, ecolineages_sr_smry)
colnames(genera_sr_smry)[8] <- "97.5%"
colnames(clades_sr_smry)[8] <- "97.5%"
colnames(ecological_sr_smry)[8] <- "97.5%"
colnames(ecolineages_sr_smry)[8] <- "97.5%"
info_disp <- rbind(genera_sr_smry, clades_sr_smry, ecological_sr_smry, 
                   ecolineages_sr_smry)
options(digits = 2)
summary_sr <- data.frame(labels, info_disp)
summary_sr <- summary_sr[-c(7:8)]
colnames(summary_sr) <-
  c("Subdivision",
    "Group",
    "Sample_size",
    "SR_Obs",
    "Median_BS",
    "0.025_BS",
    "0.975_BS")
summary_sr$Group[summary_sr$Group == "Mimosa_clade"] <-
  "Mimosa clade"
summary_sr$Group[summary_sr$Group == "Stryphnodendron_clade"] <-
  "Stryphnodendron clade"
summary_sr$Group[summary_sr$Group == "Flooded_grassland"] <-
  "Flooded grassland"
summary_sr$Group[summary_sr$Group == "Tropical_rainforest"] <-
  "Tropical rainforest"
summary_sr$Group[summary_sr$Group == "mimosa_sdtf_2"] <-
  "Mimosa SDTF 2"
summary_sr$Group[summary_sr$Group == "mimosa_sdtf_3"] <-
  "Mimosa SDTF 3"
summary_sr$Group[summary_sr$Group == "mimosa_sdtf_5"] <-
  "Mimosa SDTF 5"
summary_sr$Group[summary_sr$Group == "Stryphnod_SDTF"] <-
  "Stryphnodendron SDTF"
summary_sr$Group[summary_sr$Group == "mimosa_tropical_1"] <-
  "Mimosa Tropical 1"
summary_sr$Group[summary_sr$Group == "Stryphnod_Tropical"] <-
  "Stryphnodendron Tropical"

### SR statistics

statistic_W_sr <-
  c(
    unlist(genera_sr_wilcox)[c(1:21)],
    unlist(clades_sr_wilcox)[c(1)],
    unlist(ecological_sr_wilcox)[c(1:15)],
    unlist(ecolineages_sr_wilcox)[c(1:10)]
  )

p_value_sr <-
  c(
    unlist(genera_sr_wilcox)[c(22:42)],
    unlist(clades_sr_wilcox)[c(2)],
    unlist(ecological_sr_wilcox)[c(16:30)],
    unlist(ecolineages_sr_wilcox)[c(11:20)]
  )

b_coeff_sr <-
  c(
    unlist(genera_sr_bhatt),
    unlist(clades_sr_bhatt),
    unlist(ecological_sr_bhatt),
    unlist(ecolineages_sr_bhatt)
  )
gr_lab_sr <-
  c(
    rownames(genera_sr_bhatt),
    rownames(clades_sr_bhatt),
    rownames(ecological_sr_bhatt),
    rownames(ecolineages_sr_bhatt)
  )
pt_lab_sr <-
  c(rep("Genera", 21),
    "Clades",
    rep("Biomes", 15),
    rep("Ecological lineages", 10))

summary_stats_sr <-
  data.frame(pt_lab_sr,
             gr_lab_sr,
             round(statistic_W_sr, 2),
             round(p_value_sr, 2),
             round(b_coeff_sr, 2))

colnames(summary_stats_sr) <-
  c("Subdivision",
    "Groups",
    "Statistic W",
    "P-value",
    "Bhattacharyya Coefficient")

### saving workbooks

wb <- createWorkbook()
addWorksheet(wb, "SR")
addWorksheet(wb, "Stats - SR")
writeData(wb, 1, summary_sr)
writeData(wb, 2, summary_stats_sr)
saveWorkbook(wb, "output/sum_of_ranges.xlsx", overwrite = TRUE)

# SR Boxplots (Supplementary Fig. S1) =========================================

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
genera_plot_disp <- vector(mode = "list", length = length(genera_sr))

sum_Gwilymia <- as.vector(genera_sr$disparity$Gwilymia[[2]])
sum_Marlimorimia <-
  as.vector(genera_sr$disparity$Marlimorimia[[2]])
sum_Mimosa <- as.vector(genera_sr$disparity$Mimosa[[2]])
sum_Parapiptadenia <-
  as.vector(genera_sr$disparity$Parapiptadenia[[2]])
sum_Piptadenia <- as.vector(genera_sr$disparity$Piptadenia[[2]])
sum_Pityrocarpa <-
  as.vector(genera_sr$disparity$Pityrocarpa[[2]])
sum_Stryphnodendron <-
  as.vector(genera_sr$disparity$Stryphnodendron[[2]])

SR_Gwilymia <- data.frame(F1, sum_Gwilymia)
colnames(SR_Gwilymia) <- c("Group", "Sum_of_ranges")
SR_Marlimorimia <- data.frame(F2, sum_Marlimorimia)
colnames(SR_Marlimorimia) <- c("Group", "Sum_of_ranges")
SR_Mimosa <- data.frame(F3, sum_Mimosa)
colnames(SR_Mimosa) <- c("Group", "Sum_of_ranges")
SR_Parapiptadenia <- data.frame(F4, sum_Parapiptadenia)
colnames(SR_Parapiptadenia) <- c("Group", "Sum_of_ranges")
SR_Piptadenia <- data.frame(F5, sum_Piptadenia)
colnames(SR_Piptadenia) <- c("Group", "Sum_of_ranges")
SR_Pityrocarpa <- data.frame(F6, sum_Pityrocarpa)
colnames(SR_Pityrocarpa) <- c("Group", "Sum_of_ranges")
SR_Stryphnodendron <- data.frame(F7, sum_Stryphnodendron)
colnames(SR_Stryphnodendron) <- c("Group", "Sum_of_ranges")
genera_to_plot <-
  rbind(
    SR_Gwilymia,
    SR_Marlimorimia,
    SR_Mimosa,
    SR_Parapiptadenia,
    SR_Piptadenia,
    SR_Pityrocarpa,
    SR_Stryphnodendron
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

plot_genera_sr <- ggplot(genera_to_plot,
                         aes(
                           x = Group,
                           y = Sum_of_ranges,
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
clades_plot_disp <- vector(mode = "list", length = length(clades_sr))

sum_Mimosa_clade <-
  as.vector(clades_sr$disparity$Mimosa_clade[[2]])
sum_Stryphnodendron_clade <-
  as.vector(clades_sr$disparity$Stryphnodendron_clade[[2]])

SR_Mimosa_clade <- data.frame(C1, sum_Mimosa_clade)
colnames(SR_Mimosa_clade) <- c("Group", "Sum_of_ranges")
SR_Stryphnodendron_clade <- data.frame(C2, sum_Stryphnodendron_clade)
colnames(SR_Stryphnodendron_clade) <- c("Group", "Sum_of_ranges")
clades_to_plot <- rbind(SR_Mimosa_clade, SR_Stryphnodendron_clade)

clades_to_plot$Group <-
  factor(clades_to_plot$Group,
         levels = c("Mimosa clade", "Stryphnodendron clade"))

plot_clades_sr <-
  ggplot(clades_to_plot,
         aes(
           x = Group,
           y = Sum_of_ranges,
           colour = Group,
           fill = Group
         ),
         alpha = 0.9) +
  geom_boxplot(alpha = 0.5, show.legend = F) + labs(x = "", y = "Sum of ranges") +
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
  vector(mode = "list", length = length(ecological_sr))

sum_Desert <- as.vector(ecological_sr$disparity$Desert[[2]])
sum_Flooded_grassland <-
  as.vector(ecological_sr$disparity$Flooded_grassland[[2]])
sum_Savanna <- as.vector(ecological_sr$disparity$Savanna[[2]])
sum_SDTF <- as.vector(ecological_sr$disparity$SDTF[[2]])
sum_Temperate <-
  as.vector(ecological_sr$disparity$Temperate[[2]])
sum_Tropical_rainforest <-
  as.vector(ecological_sr$disparity$Tropical_rainforest[[2]])

SR_Desert <- data.frame(E1, sum_Desert)
colnames(SR_Desert) <- c("Group", "Sum_of_ranges")
SR_Flooded_grassland <- data.frame(E2, sum_Flooded_grassland)
colnames(SR_Flooded_grassland) <- c("Group", "Sum_of_ranges")
SR_Savanna <- data.frame(E3, sum_Savanna)
colnames(SR_Savanna) <- c("Group", "Sum_of_ranges")
SR_SDTF <- data.frame(E4, sum_SDTF)
colnames(SR_SDTF) <- c("Group", "Sum_of_ranges")
SR_Temperate <- data.frame(E5, sum_Temperate)
colnames(SR_Temperate) <- c("Group", "Sum_of_ranges")
SR_Tropical_rainforest <- data.frame(E6, sum_Tropical_rainforest)
colnames(SR_Tropical_rainforest) <- c("Group", "Sum_of_ranges")
ecological_to_plot <-
  rbind(
    SR_Desert,
    SR_Flooded_grassland,
    SR_Savanna,
    SR_SDTF,
    SR_Temperate,
    SR_Tropical_rainforest
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

plot_ecogroups_sr <-
  ggplot(
    ecological_to_plot,
    aes(
      x = Group,
      y = Sum_of_ranges,
      colour = Group,
      fill = Group
    ),
    alpha = 0.9
  ) +
  geom_boxplot(alpha = 0.5, show.legend = F) + labs(x = "", y = "Sum of ranges") +
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
  vector(mode = "list", length = length(ecolineages_sr))

sum_Mimosa_Savanna_2 <-
  as.vector(ecolineages_sr$disparity$Mimosa_Savanna_2[[2]])
sum_Mimosa_SDTF_1 <-
  as.vector(ecolineages_sr$disparity$Mimosa_SDTF_1[[2]])
sum_Mimosa_SDTF_2 <-
  as.vector(ecolineages_sr$disparity$Mimosa_SDTF_2[[2]])
sum_Stryphnodendron_sdtf <-
  as.vector(ecolineages_sr$disparity$Stryphnod_SDTF[[2]])
sum_Stryphnodendron_tropical <-
  as.vector(ecolineages_sr$disparity$Stryphnod_Tropical[[2]])

SR_Mimosa_Savanna_2 <- data.frame(EL1, sum_Mimosa_Savanna_2)
colnames(SR_Mimosa_Savanna_2) <- c("Group", "Sum_of_ranges")
SR_Mimosa_SDTF_1 <- data.frame(EL2, sum_Mimosa_SDTF_1)
colnames(SR_Mimosa_SDTF_1) <- c("Group", "Sum_of_ranges")
SR_Mimosa_SDTF_2 <- data.frame(EL3, sum_Mimosa_SDTF_2)
colnames(SR_Mimosa_SDTF_2) <- c("Group", "Sum_of_ranges")
SR_Stryphnodendron_sdtf <- data.frame(EL4, sum_Stryphnodendron_sdtf)
colnames(SR_Stryphnodendron_sdtf) <- c("Group", "Sum_of_ranges")
SR_Stryphnodendron_tropical <-
  data.frame(EL5, sum_Stryphnodendron_tropical)
colnames(SR_Stryphnodendron_tropical) <- c("Group", "Sum_of_ranges")
ecolineages_to_plot <-
  rbind(
    SR_Mimosa_Savanna_2,
    SR_Mimosa_SDTF_1,
    SR_Mimosa_SDTF_2,
    SR_Stryphnodendron_sdtf,
    SR_Stryphnodendron_tropical
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

plot_ecolineages_sr <-
  ggplot(
    ecolineages_to_plot,
    aes(
      x = Group,
      y = Sum_of_ranges,
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

# SR Boxplots - rarefied (Supplementary Fig. S2) ==============================

#### Genera
F1 <- c(rep("Gwilymia", bootstraps))
F2 <- c(rep("Marlimorimia", bootstraps))
F3 <- c(rep("Mimosa", bootstraps))
F4 <- c(rep("Parapiptadenia", bootstraps))
F5 <- c(rep("Piptadenia", bootstraps))
F6 <- c(rep("Pityrocarpa", bootstraps))
F7 <- c(rep("Stryphnodendron", bootstraps))
genera_plot_disp <- vector(mode = "list", length = length(genera_sr))

sum_Gwilymia <- as.vector(genera_sr$disparity$Gwilymia[[3]])
sum_Marlimorimia <-
  as.vector(genera_sr$disparity$Marlimorimia[[2]])
sum_Mimosa <- as.vector(genera_sr$disparity$Mimosa[[3]])
sum_Parapiptadenia <-
  as.vector(genera_sr$disparity$Parapiptadenia[[3]])
sum_Piptadenia <- as.vector(genera_sr$disparity$Piptadenia[[3]])
sum_Pityrocarpa <-
  as.vector(genera_sr$disparity$Pityrocarpa[[3]])
sum_Stryphnodendron <-
  as.vector(genera_sr$disparity$Stryphnodendron[[3]])

SR_Gwilymia <- data.frame(F1, sum_Gwilymia)
colnames(SR_Gwilymia) <- c("Group", "Sum_of_ranges")
SR_Marlimorimia <- data.frame(F2, sum_Marlimorimia)
colnames(SR_Marlimorimia) <- c("Group", "Sum_of_ranges")
SR_Mimosa <- data.frame(F3, sum_Mimosa)
colnames(SR_Mimosa) <- c("Group", "Sum_of_ranges")
SR_Parapiptadenia <- data.frame(F4, sum_Parapiptadenia)
colnames(SR_Parapiptadenia) <- c("Group", "Sum_of_ranges")
SR_Piptadenia <- data.frame(F5, sum_Piptadenia)
colnames(SR_Piptadenia) <- c("Group", "Sum_of_ranges")
SR_Pityrocarpa <- data.frame(F6, sum_Pityrocarpa)
colnames(SR_Pityrocarpa) <- c("Group", "Sum_of_ranges")
SR_Stryphnodendron <- data.frame(F7, sum_Stryphnodendron)
colnames(SR_Stryphnodendron) <- c("Group", "Sum_of_ranges")
genera_to_plot <-
  rbind(
    SR_Gwilymia,
    SR_Marlimorimia,
    SR_Mimosa,
    SR_Parapiptadenia,
    SR_Piptadenia,
    SR_Pityrocarpa,
    SR_Stryphnodendron
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

plot_genera_sr_rare <- ggplot(genera_to_plot,
                              aes(
                                x = Group,
                                y = Sum_of_ranges,
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
clades_plot_disp <- vector(mode = "list", length = length(clades_sr))

sum_Mimosa_clade <-
  as.vector(clades_sr$disparity$Mimosa_clade[[3]])
sum_Stryphnodendron_clade <-
  as.vector(clades_sr$disparity$Stryphnodendron_clade[[2]])

SR_Mimosa_clade <- data.frame(C1, sum_Mimosa_clade)
colnames(SR_Mimosa_clade) <- c("Group", "Sum_of_ranges")
SR_Stryphnodendron_clade <- data.frame(C2, sum_Stryphnodendron_clade)
colnames(SR_Stryphnodendron_clade) <- c("Group", "Sum_of_ranges")
clades_to_plot <- rbind(SR_Mimosa_clade, SR_Stryphnodendron_clade)

clades_to_plot$Group <-
  factor(clades_to_plot$Group,
         levels = c("Mimosa clade", "Stryphnodendron clade"))

plot_clades_sr_rare <-
  ggplot(clades_to_plot,
         aes(
           x = Group,
           y = Sum_of_ranges,
           colour = Group,
           fill = Group
         ),
         alpha = 0.9) +
  geom_boxplot(alpha = 0.5, show.legend = F) + labs(x = "", y = "Sum of ranges") +
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
  vector(mode = "list", length = length(ecological_sr))

sum_Desert <- as.vector(ecological_sr$disparity$Desert[[2]])
sum_Flooded_grassland <-
  as.vector(ecological_sr$disparity$Flooded_grassland[[3]])
sum_Savanna <- as.vector(ecological_sr$disparity$Savanna[[3]])
sum_SDTF <- as.vector(ecological_sr$disparity$SDTF[[3]])
sum_Temperate <-
  as.vector(ecological_sr$disparity$Temperate[[3]])
sum_Tropical_rainforest <-
  as.vector(ecological_sr$disparity$Tropical_rainforest[[3]])

SR_Desert <- data.frame(E1, sum_Desert)
colnames(SR_Desert) <- c("Group", "Sum_of_ranges")
SR_Flooded_grassland <- data.frame(E2, sum_Flooded_grassland)
colnames(SR_Flooded_grassland) <- c("Group", "Sum_of_ranges")
SR_Savanna <- data.frame(E3, sum_Savanna)
colnames(SR_Savanna) <- c("Group", "Sum_of_ranges")
SR_SDTF <- data.frame(E4, sum_SDTF)
colnames(SR_SDTF) <- c("Group", "Sum_of_ranges")
SR_Temperate <- data.frame(E5, sum_Temperate)
colnames(SR_Temperate) <- c("Group", "Sum_of_ranges")
SR_Tropical_rainforest <- data.frame(E6, sum_Tropical_rainforest)
colnames(SR_Tropical_rainforest) <- c("Group", "Sum_of_ranges")
ecological_to_plot <-
  rbind(
    SR_Desert,
    SR_Flooded_grassland,
    SR_Savanna,
    SR_SDTF,
    SR_Temperate,
    SR_Tropical_rainforest
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

plot_ecogroups_sr_rare <-
  ggplot(
    ecological_to_plot,
    aes(
      x = Group,
      y = Sum_of_ranges,
      colour = Group,
      fill = Group
    ),
    alpha = 0.9
  ) +
  geom_boxplot(alpha = 0.5, show.legend = F) + labs(x = "", y = "Sum of ranges") +
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
  vector(mode = "list", length = length(ecolineages_sr))

sum_Mimosa_Savanna_2 <-
  as.vector(ecolineages_sr$disparity$Mimosa_Savanna_2[[3]])
sum_Mimosa_SDTF_1 <-
  as.vector(ecolineages_sr$disparity$Mimosa_SDTF_1[[3]])
sum_Mimosa_SDTF_2 <-
  as.vector(ecolineages_sr$disparity$Mimosa_SDTF_2[[3]])
sum_Stryphnodendron_sdtf <-
  as.vector(ecolineages_sr$disparity$Stryphnod_SDTF[[2]])
sum_Stryphnodendron_tropical <-
  as.vector(ecolineages_sr$disparity$Stryphnod_Tropical[[3]])

SR_Mimosa_Savanna_2 <- data.frame(EL1, sum_Mimosa_Savanna_2)
colnames(SR_Mimosa_Savanna_2) <- c("Group", "Sum_of_ranges")
SR_Mimosa_SDTF_1 <- data.frame(EL2, sum_Mimosa_SDTF_1)
colnames(SR_Mimosa_SDTF_1) <- c("Group", "Sum_of_ranges")
SR_Mimosa_SDTF_2 <- data.frame(EL3, sum_Mimosa_SDTF_2)
colnames(SR_Mimosa_SDTF_2) <- c("Group", "Sum_of_ranges")
SR_Stryphnodendron_sdtf <- data.frame(EL4, sum_Stryphnodendron_sdtf)
colnames(SR_Stryphnodendron_sdtf) <- c("Group", "Sum_of_ranges")
SR_Stryphnodendron_tropical <-
  data.frame(EL5, sum_Stryphnodendron_tropical)
colnames(SR_Stryphnodendron_tropical) <- c("Group", "Sum_of_ranges")
ecolineages_to_plot <-
  rbind(
    SR_Mimosa_Savanna_2,
    SR_Mimosa_SDTF_1,
    SR_Mimosa_SDTF_2,
    SR_Stryphnodendron_sdtf,
    SR_Stryphnodendron_tropical
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

plot_ecolineages_sr_rare <-
  ggplot(
    ecolineages_to_plot,
    aes(
      x = Group,
      y = Sum_of_ranges,
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

## Save plots -----------------------------------------------------------------

# SR

pdf("output/supfig2_sum_of_ranges.pdf", width = 10, height = 10)

(plot_clades_sr + plot_genera_sr) /
  (plot_ecogroups_sr + plot_ecolineages_sr)

dev.off()

# SR - rarefied

pdf("output/supfig1_sum_of_ranges_rarefied.pdf", width = 10, height = 10)

(plot_clades_sr_rare + plot_genera_sr_rare) /
  (plot_ecogroups_sr_rare + plot_ecolineages_sr_rare)

dev.off()
