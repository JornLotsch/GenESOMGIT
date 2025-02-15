#################################### Paths ########################################################################

pfad_o <- "/home/joern/Aktuell/GenerativeESOM/"
pfad_u1 <- "09Originale/"
pfad_r <- "08AnalyseProgramme/R/"
pfad_py <- "08AnalyseProgramme/Python/"
pfad_umx <- "04Umatrix/"
pfad_o0 <- "/home/joern/Aktuell/ProjectionsBiomed/"

#################################### Libraries ########################################################################

library(parallel)
library(doParallel)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggforce)
library(ggmosaic)
library(readxl)
library(cowplot)
library(missForest)
library(dplyr)

nProc <- detectCores() - 1 # round((detectCores() - 1) / 4) # /4 RF and MultiOmics data

nIter <- 1000
seed <- 42
list.of.seeds <- 1:nIter + seed - 1

source(paste0(pfad_o, pfad_r, "GenerateDataGMX.R"))

#################################### Functions ########################################################################

runAOVoneway <- function(x, y) {
  aov.res <- aov(x ~ as.factor(y))
  s <- summary(aov.res)
  pval <- s[[1]][["Pr(>F)"]][1]
  return(pval)
}

# https://stackoverflow.com/questions/31968549/holm-sidak-adjustment-using-r
Sidak <- function(vecP)
                  #
                  # This function corrects a vector of probabilities for multiple testing
                  # using the Bonferroni (1935) and Sidak (1967) corrections.
                  #
                  # References: Bonferroni (1935), Sidak (1967), Wright (1992).
                  #
                  # Bonferroni, C. E. 1935. Il calcolo delle assicurazioni su gruppi di teste.
                  # Pp. 13-60 in: Studi in onore del Professore Salvatore Ortu Carboni. Roma.
                  #
                  # Sidak, Z. 1967. Rectangular confidence regions for the means of multivariate
                  # normal distributions. Journal of the American Statistical Association 62:626-633.
                  #
                  # Wright, S. P. 1992. Adjusted P-values for simultaneous inference.
                  # Biometrics 48: 1005-1013.
                  #
#                  Pierre Legendre, May 2007
{
  k <- length(vecP)

  vecPB <- 0
  vecPS <- 0

  for (i in 1:k) {
    bonf <- vecP[i] * k
    if (bonf > 1) bonf <- 1
    vecPB <- c(vecPB, bonf)
    vecPS <- c(vecPS, (1 - (1 - vecP[i])^k))
  }
  #
  return(list(OriginalP = vecP, BonfP = vecPB[-1], SidakP = vecPS[-1]))
}

#################################### Get data  ########################################################################
Data_Behaviorial_and_lipid_profile_of_the_EAE_model_in_SJL_mice_and_effects_of_FTY720 <-
  data.frame(read_excel(paste0(pfad_o, pfad_u1, "Data Behaviorial and lipid profile of the EAE model in SJL mice and effects of FTY720-.xlsx")))
names(Data_Behaviorial_and_lipid_profile_of_the_EAE_model_in_SJL_mice_and_effects_of_FTY720) <-
  make.names(names(Data_Behaviorial_and_lipid_profile_of_the_EAE_model_in_SJL_mice_and_effects_of_FTY720))

Groups <- Data_Behaviorial_and_lipid_profile_of_the_EAE_model_in_SJL_mice_and_effects_of_FTY720$GROUP
MiceGroups <- Data_Behaviorial_and_lipid_profile_of_the_EAE_model_in_SJL_mice_and_effects_of_FTY720$MOUSE
dfMiceGroups <- cbind.data.frame(Groups, MiceGroups)

Mouse_lipids_transformed <- read.csv(paste0(pfad_o, pfad_r, "Mouse_lipids_transformed.csv"), row.names = 1)


#################################### Prepare reduced data set 1st 6 mice per group statistics  ########################################################################
pBonLimit05 <- 0.05 / ncol(Mouse_lipids_transformed)
pBonLimit01 <- 0.01 / ncol(Mouse_lipids_transformed)
n <- 6

MiceToAnalyze_6perGroup <-
  dfMiceGroups %>%
  group_by(Groups) %>%
  slice(1:n)
Mouse_lipids_transformed_6perGroup <- Mouse_lipids_transformed[rownames(Mouse_lipids_transformed) %in% MiceToAnalyze_6perGroup$MiceGroups, ]
set.seed(seed)
Mouse_lipids_transformed_imputed_6perGroup <- missForest(Mouse_lipids_transformed_6perGroup)$ximp

write.csv(cbind.data.frame(Group = MiceToAnalyze_6perGroup$Groups, Mouse_lipids_transformed_imputed_6perGroup), file = paste0(pfad_o, pfad_r, "Mouse_lipids_transformed_imputed_6perGroup.csv"))

# To Alfred
# WriteLRN(FileName = "Mouse_lipids_transformed_6perGroup.lrn",
#          Data = as.matrix(Mouse_lipids_transformed_imputed_6perGroup), Key = rownames(Mouse_lipids_transformed_imputed_6perGroup),
#          OutDirectory = paste0(pfad_o, pfad_r))
# WriteCLS(FileName = "Mouse_lipids_transformed_6perGroup.cls", Cls = MiceToAnalyze_6perGroup$Groups, Key = rownames(Mouse_lipids_transformed_6perGroup),
#          OutDirectory = paste0(pfad_o, pfad_r))

#################################### Generate data  ########################################################################
DensityRadius <- 1.355932203# 1.15 
GenPerData <- 2

set.seed(42)
MouseGenerated <- GeneratedDataGMX(
  Data = Mouse_lipids_transformed_imputed_6perGroup, Cls = MiceToAnalyze_6perGroup$Groups,
  DensityRadius = DensityRadius, GenPerData = GenPerData, PlotIt = 2
)
MouseGenerated$pGeneratedData + scale_color_viridis_d()

#################################### Do statistics  ########################################################################

pGeneratedGroups_Mouse_lipids <- data.frame(p_generated = apply(MouseGenerated$GeneratedData, 2, function(x) runAOVoneway(x = x, y = MouseGenerated$GeneratedClasses)))
pGeneratedGroups_Mouse_lipids$Variable <- rownames(pGeneratedGroups_Mouse_lipids)

#################################### Plot data and results  ########################################################################

pGeneratedPvalues <-
  ggplot(data = pGeneratedGroups_Mouse_lipids, aes(x = reorder(Variable, p_generated), y = -log10(p_generated))) +
  geom_bar(stat = "identity", position = "dodge", alpha = .5) + #, fill = "lightyellow4") +
  geom_hline(yintercept = c(-log10(0.01), -log10(0.05), -log10(pBonLimit01), -log10(pBonLimit05)), color = c("dodgerblue", "blue", "red", "salmon"), linetype = "dashed") +
  annotate("text",
    x = c(50), y = c(-log10(0.01), -log10(0.05), -log10(pBonLimit01), -log10(pBonLimit05)),
    color = c("dodgerblue", "blue", "red", "salmon"),
    label = c("Raw p = 0.01", "Raw p = 0.05", "Bonferroni p = 0.01", "Bonferroni p = 0.05"), vjust = -1
  ) +
  theme_light() +
  theme(
    legend.position = c(.8, .8),
    strip.background = element_rect(fill = "cornsilk"), strip.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(title = "Group ANOVA p-Values, ESOM-based generated data", x = "Lipid marker", y = "-log10(p)", fill = "Transformation") +
  scale_fill_colorblind()


OrigData <- cbind.data.frame(Group = MouseGenerated$OriginalClasses, MouseGenerated$OriginalData)
GeneratedData <- cbind.data.frame(Group = MouseGenerated$GeneratedClasses, MouseGenerated$GeneratedData)
OrigData_long <- reshape2::melt(OrigData, id.vars = "Group")
GeneratedData_long <- reshape2::melt(GeneratedData, id.vars = "Group")
dim(GeneratedData)

pGeneratedData <-
  ggplot() +
  geom_point(
    data = GeneratedData_long, aes(x = variable, y = value, fill = factor(Group), group = factor(Group)),
    alpha = .3, size = 1, position = position_dodge(width = .9), shape = 17, fill = "white"
  ) +
  geom_point(
    data = OrigData_long, aes(x = variable, y = value, color = factor(Group), fill = factor(Group)),
    alpha = 1, size = 1, shape = 19, position = position_dodge(width = .9)
  ) +
  theme_light() +
  theme(legend.position = c(.8, .9), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.direction = "horizontal") +
  labs(title = "Original data from 6 mice and generated data", x = NULL, y = "log concentration [ng/ml]", color = "Group", fill = "Group") +
  scale_color_colorblind() +
  scale_fill_colorblind()

yRange <- ggplot_build(pGeneratedData)$layout$panel_params[[1]]$y.range

pGeneratedData_PL_LPA16_0 <-
  ggplot() +
  geom_point(
    data = GeneratedData_long[GeneratedData_long$variable == "PL_LPA16_0", ], aes(x = variable, y = value, fill = factor(Group), group = factor(Group)),
    alpha = .3, size = 2, position = position_dodge(width = .9), shape = 17, fill = "white"
  ) +
  geom_point(
    data = OrigData_long[OrigData_long$variable == "PL_LPA16_0", ], aes(x = variable, y = value, color = factor(Group), fill = factor(Group)),
    alpha = 1, size = 2, shape = 19, position = position_dodge(width = .9), show.legend = F
  ) +
  theme_light() +
  theme(legend.position = c(.8, .9), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.direction = "horizontal") +
  labs(title = "PL_LPA16_0", x = NULL, y = "log concentration [ng/ml]", color = "Group", fill = "Group") +
  scale_color_colorblind() +
  scale_fill_colorblind() #+
#  ylim(yRange)


plot_grid(
  plot_grid(
    pGeneratedData,
    pGeneratedData_PL_LPA16_0,
    labels = LETTERS[1:2],
    nrow = 1, rel_widths = c(6.8, 1),
    align = "h", axis = "tb"
  ),
  plot_grid(
    pGeneratedPvalues,
    labels = LETTERS[3]
  ),
  ncol = 1,
  align = "v", axis = "lr"
)


#################################### QQplots generated versus orginal ########################################################################
quantiles <- seq(0, 1, 0.01)

quantilesMouse_lipids_transformed_imputed <- data.frame(apply(Mouse_lipids_transformed_imputed, 2, function(x) quantile(x, quantiles, na.rm = T)))
quantilesMouse_lipids_transformed_imputed_long <- reshape2::melt(quantilesMouse_lipids_transformed_imputed)
quantilesMouse_lipids_generated <- data.frame(apply(MouseGenerated$GeneratedData, 2, function(x) quantile(x, quantiles, na.rm = T)))
quantilesMouse_lipids_generated_long <- reshape2::melt(quantilesMouse_lipids_generated)

fun <- function(x, y) {
  ks.test(x, y)$p.value
}
ksPOrigGenerated <- data.frame(p = mapply(FUN = fun, quantilesMouse_lipids_transformed_imputed, quantilesMouse_lipids_generated))
ksPOrigGenerated$variable <- rownames(ksPOrigGenerated)
psych::describe(ksPOrigGenerated$p)
min(ksPOrigGenerated$p)

pKSp_generated <-
  ggplot(data = ksPOrigGenerated, aes(y = reorder(variable, p), x = p)) +
  geom_bar(stat = "identity", alpha = .5) +
  geom_vline(xintercept = 0.05, color = "salmon", linetype = "dashed") +
  theme_light() +
  labs(title = "KS tests", x = "p-value", y = NULL)


# dfQQ <- rbind.data.frame(cbind.data.frame(Set = "original", Mouse_lipids_transformed_long),
#                          cbind.data.frame(Set = "imputed", Mouse_lipids_transformed_imputed_long)
# )
dfQQ_generated <- cbind.data.frame(quantilesMouse_lipids_transformed_long, "value_generated" = quantilesMouse_lipids_generated_long$value)



pQQgeneration <-
  ggplot(data = dfQQ_generated, aes(x = value, y = value_generated)) +
  geom_point(color = "dodgerblue", alpha = .6) +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2, color = "salmon") +
  facet_wrap(Var2 ~ ., scales = "free") +
  theme_light() +
  theme(legend.position = c(.2, .2), strip.background = element_rect(fill = "cornsilk"), strip.text = element_text(colour = "black")) +
  labs(title = "QQ plots of generated versus original variable distribution", x = "Original", y = "Imputed")


plot_grid(
  pKSp_generated ,
  pQQgeneration,
  nrow = 1, rel_widths = c(1, 5),
  align = "h",  axis = "tb",
  labels = LETTERS[1:8]
)

