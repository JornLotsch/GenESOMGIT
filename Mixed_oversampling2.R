# TesData oversampling experiment

#################################### Paths ########################################################################

pfad_o <- "/home/joern/Aktuell/GenerativeESOM/"
pfad_u1 <- "09Originale/"
pfad_r <- "08AnalyseProgramme/R/"
pfad_py <- "08AnalyseProgramme/Python/"
pfad_umx <- "04Umatrix/"

OutFilename <- "UmxTest_oversampling"
OutDirectory <- paste0(pfad_o, pfad_umx)

source("/home/joern/Aktuell/GenerativeESOM/08AnalyseProgramme/R/HyperparametertuningRF.R")
source("/home/joern/Aktuell/GenerativeESOM/08AnalyseProgramme/R/generate_data1.R")

#################################### Libraries ########################################################################

library(parallel)
library(caret)
library(randomForest)
library(caTools)
library(pROC)
library(stringr)
library(matrixStats)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggstatsplot)
library(cowplot)

nProc <- detectCores() - 1 # 20#

nIter <- 100
seed <- 42
list.of.seeds <- 1:nIter + seed - 1

#################################### Functions ########################################################################

quantiles_95 <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

quantiles_100 <- function(x) {
  r <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

is.integer0 <- function(x) {
  is.integer(x) && length(x) == 0L
}

#################################### Get Test data ########################################################################

# Normal effect
n = 10
ma = seq(from = 1, to = 50, length.out = 50)
mb = seq(from = 1, to = 40, length.out = 50)

dfTest <- lapply(1:50, function(i) {
  set.seed(42)
  a <- rnorm(n, mean = ma[i], sd = 1)
  b <- rnorm(n, mean = mb[i], sd = 1)

  return(c(a, b))
})

dfTestall <- cbind.data.frame(Target = rep(1:2, each = n), do.call(cbind, dfTest))

apply(within(dfTestall, rm(Target)), 2, function(x) t.test(x)$p.value)

Test_DataCls <- dfTestall

# No effect

nVars <- 150
nCases <- 20
seed <- 42
list.of.seeds.nVars <- 1:nVars + seed - 1
Classes_noDiff <- rep(1:2, each = nCases / 2)

data_noDiff_list <- lapply(list.of.seeds.nVars, function(i) {
  set.seed(i)
  mean_i = sample(10:30, 1)
  SD_i_1 <- sample(seq(1, 3, 0.01), size = 1)
  Data_i_1 <- rnorm(n = nCases / 2, mean = mean_i, sd = SD_i_1)
  set.seed(i + 1000)
  SD_i_2 <- sample(seq(1, 3, 0.01), size = 1)
  SD_i_2 <- SD_i_1
  Data_i_2 <- rnorm(n = nCases / 2, mean = mean_i, sd = SD_i_2)
  Data_i <- c(Data_i_1, Data_i_2)
  return(list(Data_i = Data_i, SD_i_1 = SD_i_1, SD_i_2 = SD_i_2))
})

data_noDiff <- do.call(cbind, lapply(data_noDiff_list, "[[", 1))
data_noDiff_SD_1 <- unlist(lapply(data_noDiff_list, "[[", 2))
data_noDiff_SD_2 <- unlist(lapply(data_noDiff_list, "[[", 3))


pvals_t <- apply(data_noDiff, 2, function(x) t.test(x ~ Classes_noDiff)$p.value)
sum(pvals_t < 0.05)
nonSig05 <- which(pvals_t > 0.6)
length(nonSig05)

set.seed(seed)
data_noDiff_50 <- data_noDiff[, sample(nonSig05, 50)]
dim(data_noDiff_50)

Test_DataCls <- cbind.data.frame(Target = rep(1:2, each = nCases / 2), data_noDiff_50)
t.test(Test_DataCls$"40" ~ Test_DataCls$Target)

# Mouse

Test_DataCls <- read.csv(paste0(pfad_o, pfad_r, "Mouse_lipids_transformed_imputed_6perGroup.csv"), row.names = 1)
names(Test_DataCls)[1] <- "Target"

#################################### Standardize test data ########################################################################

Test_DataCls[,2:ncol(Test_DataCls)] <- apply(Test_DataCls[,2:ncol(Test_DataCls)],2,scale)  

#################################### U-matrix ########################################################################
umx_dfTestall <- Umatrix::iEsomTrain(Data = as.matrix(within(Test_DataCls, rm(Target))), Cls = Test_DataCls$Target)


Radius_umx_dfTestall <- Umatrix::calculate_Delauny_radius(Data = as.matrix(within(Test_DataCls, rm(Target))),
                                  BestMatches = umx_dfTestall$BestMatches, Columns = 50, Lines = 80, Toroid = F)

RadiusTesData <- Radius_umx_dfTestall$RadiusByEM
RadiusTesData <- 3.5 # NormalEffect scaled
RadiusTesData <- 6.6 # NoEffect
RadiusTesData <- 1.355932203 # Mouse

#################################### Split data set into training/test / validation  ########################################################################

TesDataTrainingTestValidation <- opdisDownsampling::opdisDownsampling(Data = within(Test_DataCls, rm(Target)), Cls = Test_DataCls$Target,
                                                                 Size = 0.8 * nrow(Test_DataCls), Seed = seed, nTrials = 10000, MaxCores = nProc)

Test_DataCls_TrainingTest <- Test_DataCls[rownames(Test_DataCls) %in% TesDataTrainingTestValidation$ReducedInstances,]
table(Test_DataCls_TrainingTest$Target)

Test_DataCls_Validation <- Test_DataCls[!rownames(Test_DataCls) %in% TesDataTrainingTestValidation$ReducedInstances,]
table(Test_DataCls_Validation$Target)

# Test_DataCls[,2:5] <- apply(Test_DataCls[,2:5],2,sample)
# Test_DataCls_TrainingTest[,2:5] <- apply(Test_DataCls_TrainingTest[,2:5],2,sample)


#################################### Evaluate variable importance  ########################################################################

DataSetSizes <- c("original", "engineered_0", "augmented_1_engineered", "augmented_5_engineered") #, "augmented_50_engineered")
GenPerData <- 1

Test_VarImps <- lapply(DataSetSizes, function(DatasetNr) {

  Imps_repeated <- pbmcapply::pbmclapply(list.of.seeds, function(x) {

    switch(DatasetNr,
           "original" = Test_DataCls_actual <- Test_DataCls,
           "reduced" = Test_DataCls_actual <- Test_DataCls_TrainingTest,
           "augmented_1" = {
      set.seed(x)
      Test_DataCls_generated <- generate_data1(Data = within(Test_DataCls, rm(Target)),
                                                      density_radius = RadiusTesData, gen_per_data = GenPerData, Cls = Test_DataCls$Target)

      Test_DataCls_actual <- rbind.data.frame(
               cbind.data.frame(Target = Test_DataCls_generated$original_classes, Test_DataCls_generated$original_data),
               cbind.data.frame(Target = Test_DataCls_generated$generated_classes, Test_DataCls_generated$generated_data)
             )

    },
           "augmented_5" = {
      set.seed(42)
      Test_DataCls_generated <- generate_data1(Data = within(Test_DataCls, rm(Target)),
                                                      density_radius = RadiusTesData, gen_per_data = 5 * GenPerData, Cls = Test_DataCls$Target)

      Test_DataCls_actual <- rbind.data.frame(
               cbind.data.frame(Target = Test_DataCls_generated$original_classes, Test_DataCls_generated$original_data),
               cbind.data.frame(Target = Test_DataCls_generated$generated_classes, Test_DataCls_generated$generated_data)
             )

    },
           "augmented_1_engineered" = {
      set.seed(x)
      Test_DataCls_engeneered <-
               cbind.data.frame(
                 Test_DataCls, apply(within(Test_DataCls, rm(Target)), 2, sample)
               )
      names(Test_DataCls_engeneered) <- c("Target", names(within(Test_DataCls, rm(Target))), paste0(names(within(Test_DataCls, rm(Target))), "_permuted"))
      set.seed(x)
      Test_DataCls_generated <- generate_data1(Data = within(Test_DataCls_engeneered, rm(Target)),
                                                      density_radius = RadiusTesData, gen_per_data = GenPerData, Cls = Test_DataCls_engeneered$Target)


      Test_DataCls_actual <- rbind.data.frame(
               cbind.data.frame(Target = Test_DataCls_generated$original_classes, Test_DataCls_generated$original_data),
               cbind.data.frame(Target = Test_DataCls_generated$generated_classes, Test_DataCls_generated$generated_data)
             )
    },

           "augmented_5_engineered" = {
      set.seed(x)
      Test_DataCls_engeneered <-
               cbind.data.frame(
                 Test_DataCls, apply(within(Test_DataCls, rm(Target)), 2, sample)
               )
      names(Test_DataCls_engeneered) <- c("Target", names(within(Test_DataCls, rm(Target))), paste0(names(within(Test_DataCls, rm(Target))), "_permuted"))
      set.seed(x)
      Test_DataCls_generated <- generate_data1(Data = within(Test_DataCls_engeneered, rm(Target)),
                                                      density_radius = RadiusTesData, gen_per_data = 5 * GenPerData, Cls = Test_DataCls_engeneered$Target)


      Test_DataCls_actual <- rbind.data.frame(
               cbind.data.frame(Target = Test_DataCls_generated$original_classes, Test_DataCls_generated$original_data),
               cbind.data.frame(Target = Test_DataCls_generated$generated_classes, Test_DataCls_generated$generated_data)
             )
    },
           "augmented_50_engineered" = {
      set.seed(x)
      Test_DataCls_engeneered <-
               cbind.data.frame(
                 Test_DataCls, apply(within(Test_DataCls, rm(Target)), 2, sample)
               )
      names(Test_DataCls_engeneered) <- c("Target", names(within(Test_DataCls, rm(Target))), paste0(names(within(Test_DataCls, rm(Target))), "_permuted"))
      set.seed(x)
      Test_DataCls_generated <- generate_data1(Data = within(Test_DataCls_engeneered, rm(Target)),
                                                      density_radius = RadiusTesData, gen_per_data = 50 * GenPerData, Cls = Test_DataCls_engeneered$Target)


      Test_DataCls_actual <- rbind.data.frame(
               cbind.data.frame(Target = Test_DataCls_generated$original_classes, Test_DataCls_generated$original_data),
               cbind.data.frame(Target = Test_DataCls_generated$generated_classes, Test_DataCls_generated$generated_data)
             )
    },

           "engineered_0" = {
      set.seed(x)
      Test_DataCls_actual <-
               cbind.data.frame(
                 Test_DataCls, apply(within(Test_DataCls, rm(Target)), 2, sample)
               )
      names(Test_DataCls_actual) <- c("Target", names(within(Test_DataCls, rm(Target))), paste0(names(within(Test_DataCls, rm(Target))), "_permuted"))
    }

    )

    set.seed(x)
    inTraining <- createDataPartition(Test_DataCls_actual$Target, p = .67, list = FALSE)
    training <- Test_DataCls_actual[inTraining,]

    Test_Boruta_actual <- Boruta::Boruta(Target ~ ., training, pValue = 0.0001, maxRuns = 100)
    # plot(Test_Boruta_actual)

    BorutaStats <- Boruta::attStats(Test_Boruta_actual)
    BorutaStats$Var <- rownames(BorutaStats)

    dfVarimp_Test_Boruta_actual_long <- reshape2::melt(Test_Boruta_actual$ImpHistory)
    dfVarimp_Test_Boruta_actual_long$ColorVar <- ifelse(dfVarimp_Test_Boruta_actual_long$Var2 %in% names(Test_DataCls_actual), 2, 1)
    dfVarimp_Test_Boruta_actual_long$ColorVar[grep("permuted", dfVarimp_Test_Boruta_actual_long$Var2)] <- 3

    return(list(
  Test_Boruta_actual = Test_Boruta_actual,
  BorutaStats = BorutaStats,
  dfVarimp_Test_Boruta_actual_long = dfVarimp_Test_Boruta_actual_long
    ))
  }, mc.cores = .5 * nProc)

  BorutaStats_all <- do.call(rbind.data.frame, lapply(Imps_repeated, function(x) x$BorutaStats))
  # boxplot(BorutaStats_all$medianImp ~ BorutaStats_all$Var, las = 2)

  ABCfeatureImportances_List <- list(
    ABCfeatureImportances = NA,
    AABCfeatureImportances = NA,
    AAABCfeatureImportances = NA,
    names_reducedFeatureSet_A = NA,
    names_reducedFeatureSet_AA = NA,
    names_reducedFeatureSet_AAA = NA,
    names_reducedFeatureSet_C = NA,
    dfVars = NA,
    limtFreq = NA,
    limtFreqRel = NA
  )

  if (!is.integer0(grep("permuted", BorutaStats_all$Var))) {
    BorutaStats_all_permuted_notrejected <- BorutaStats_all[grep("permuted", BorutaStats_all$Var),]
    BorutaStats_all_permuted_notrejected <- BorutaStats_all_permuted_notrejected[BorutaStats_all_permuted_notrejected$decision == "Confirmed",]

    BorutaStats_all_notrejected <- BorutaStats_all[-grep("permuted", BorutaStats_all$Var),]
    BorutaStats_all_notrejected <- BorutaStats_all_notrejected[BorutaStats_all_notrejected$decision == "Confirmed",]
  } else {
    BorutaStats_all_notrejected <- BorutaStats_all
    BorutaStats_all_notrejected <- BorutaStats_all_notrejected[BorutaStats_all_notrejected$decision == "Confirmed",]
  }

  SelectedTrue <- data.frame(table(BorutaStats_all_notrejected$Var))

  if (!is.integer0(grep("permuted", BorutaStats_all$Var))) {
    SelectedPermuted <- data.frame(table(BorutaStats_all_permuted_notrejected$Var))
    SelectedPermuted$Var2 <- gsub("_permuted", "", SelectedPermuted$Var1)
    SelectedPermuted <- SelectedPermuted[SelectedPermuted$Var2 %in% SelectedTrue$Var1,]
  }

  if (!is.integer0(grep("permuted", BorutaStats_all$Var))) {
    dfVars <- data.frame(Var = unique(BorutaStats_all$Var[-grep("permuted", BorutaStats_all$Var)]))
  } else {
    dfVars <- data.frame(Var = unique(BorutaStats_all$Var))
  }
  rownames(dfVars) <- dfVars$Var
  dfVars$SelectedTrue <- SelectedTrue$Freq[match(dfVars$Var, SelectedTrue$Var1)]
  if (!is.integer0(grep("permuted", BorutaStats_all$Var))) {
    dfVars$SelectedPermuted <- SelectedPermuted$Freq[match(dfVars$Var, SelectedPermuted$Var2)]
  } else {
    dfVars$SelectedPermuted <- NA
  }

  dfVars[is.na(dfVars)] <- 0
  dfVars$SelectedTrueCorr <- dfVars$SelectedTrue - dfVars$SelectedPermuted
  dfVars$SelectedTrueCorrRel <- (dfVars$SelectedTrue - dfVars$SelectedPermuted) * (dfVars$SelectedTrue + dfVars$SelectedPermuted)

  if (!is.integer0(grep("permuted", BorutaStats_all$Var))) {
    nBootstrap = 100000
    a <- dfVars$SelectedPermuted
    b <- dfVars$SelectedTrue
    set.seed(seed)
    a_sample <- sample(a, nBootstrap, replace = TRUE)
    set.seed(seed)
    b_sample <- sample(b, nBootstrap, replace = TRUE)
    FCs_ba_bootstrap <- b_sample - a_sample

    limtFreq <- quantile(FCs_ba_bootstrap, probs = 0.95) #max(dfVars$SelectedPermuted) #
    dfVars$Var[which(dfVars$SelectedTrueCorr > limtFreq)]

    nBootstrapRel = 100000
    a <- dfVars$SelectedPermuted
    b <- dfVars$SelectedTrue
    set.seed(seed)
    a_sample <- sample(a, nBootstrap, replace = TRUE)
    set.seed(seed)
    b_sample <- sample(b, nBootstrap, replace = TRUE)
    df_ba <- cbind.data.frame(b_sample, a_sample)
    FCs_ba_bootstrapRel <- apply(df_ba, 1, function(x)(x["a_sample"] - x["b_sample"]) * (x["a_sample"] + x["b_sample"]))
    limtFreqRel <- quantile(FCs_ba_bootstrapRel, probs = 0.95) #max(dfVars$SelectedPermuted) #
    
  } else {
    limtFreq <- NA
    limtFreqRel = NA
  }

  SelectFreq <- as.vector(dfVars$SelectedTrueCorr)
  names(SelectFreq) <- dfVars$Var
  SelectFreq <- na.omit(SelectFreq)
  attributes(SelectFreq)$na.action <- NULL

  ABCfeatureImportances = NA
  AABCfeatureImportances = NA
  AAABCfeatureImportances = NA
  names_reducedFeatureSet_A = NA
  names_reducedFeatureSet_AA = NA
  names_reducedFeatureSet_AAA = NA
  names_reducedFeatureSet_C = NA

  # if (length(SelectFreq[SelectFreq > 0]) > 0) {
  #   ABCfeatureImportances <- ABCanalysis(SelectFreq, PlotIt = T)
  #   AABCfeatureImportances <- ABCanalysis(SelectFreq[names(ABCfeatureImportances$Aind)], PlotIt = T)
  #   AAABCfeatureImportances <- ABCanalysis(SelectFreq[names(AABCfeatureImportances$Aind)], PlotIt = T)
  # 
  #   names_reducedFeatureSet_A <- names(ABCfeatureImportances$Aind)
  #   names_reducedFeatureSet_AA <- names(AABCfeatureImportances$Aind)
  #   names_reducedFeatureSet_AAA <- names(AAABCfeatureImportances$Aind)
  #   names_reducedFeatureSet_C <- names(ABCfeatureImportances$Cind)
  # }

  ABCfeatureImportances_List <- list(
      ABCfeatureImportances = ABCfeatureImportances,
      AABCfeatureImportances = AABCfeatureImportances,
      AAABCfeatureImportances = AAABCfeatureImportances,
      names_reducedFeatureSet_A = names_reducedFeatureSet_A,
      names_reducedFeatureSet_AA = names_reducedFeatureSet_AA,
      names_reducedFeatureSet_AAA = names_reducedFeatureSet_AAA,
      names_reducedFeatureSet_C = names_reducedFeatureSet_C,
      dfVars = dfVars,
      limtFreq = limtFreq,
      limtFreqRel = limtFreqRel
    )

  dfVarimp_Test_Boruta_actual_long_all <- do.call(rbind.data.frame, lapply(Imps_repeated, function(x) x$dfVarimp_Test_Boruta_actual_long))
  # boxplot(dfVarimp_Test_Boruta_actual_long_all$value ~ dfVarimp_Test_Boruta_actual_long_all$Var2)

  upperBorderOfNonImportance <- NA
  if (!is.integer0(grep("permuted", dfVarimp_Test_Boruta_actual_long_all$Var2)))
    upperBorderOfNonImportance <- quantile(BorutaStats_all$maxImp[grep("permuted", BorutaStats_all$Var)], prob = 1)

  pVarimp_Test_actual <-
  ggplot(data = dfVarimp_Test_Boruta_actual_long_all, aes(x = reorder(Var2, value), y = value, fill = factor(ColorVar))) +
  # geom_boxplot()+
  stat_summary(fun.data = quantiles_100, geom = "boxplot", alpha = 0.2, width = 0.5, position = "dodge") +
  scale_fill_manual(values = c("dodgerblue4", "chartreuse2", "salmon"), labels = c("Dummy", "True features", "Permuted features")) +
  labs(title = "Variable importances", y = "Importance [% decrease in accuracy]", x = NULL, fill = "Feature class") +
  theme_light() +
  theme(
    legend.position = c(.2, .8), legend.direction = "vertical",
    legend.background = element_rect(colour = "transparent", fill = ggplot2::alpha("white", 0.2)),
    strip.background = element_rect(fill = "cornsilk"), strip.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

  if (!is.na(upperBorderOfNonImportance)) {
    pVarimp_Test_actual <-
    pVarimp_Test_actual +
    geom_hline(yintercept = upperBorderOfNonImportance, linetype = "dashed", color = "red") +
    annotate("text", x = .5, y = 1.05 * upperBorderOfNonImportance, label = "Limit of alpha error inflation", color = "red", hjust = -.5)
  }

  if (!is.null(dim(ABCfeatureImportances_List$dfVars))) {
    pdfVarFreq2 <-
      ggplot(data = dfVars) +
      geom_bar(aes(y = reorder(Var, SelectedTrueCorr), x = SelectedTrue), fill = "dodgerblue", stat = "identity") +
      geom_bar(aes(y = reorder(Var, SelectedTrueCorr), x = -SelectedPermuted), fill = "cornsilk4", stat = "identity") +
      labs(title = "Variable selection frequency", x = "Times selected", y = NULL, fill = "Feature class") +
      theme_light()

    pdfVarFreq <-
      ggplot(data = dfVars) +
      geom_bar(aes(y = reorder(Var, SelectedTrueCorr), x = SelectedTrueCorr), fill = "dodgerblue", stat = "identity") +
      labs(title = "Variable selection frequency", x = "Times selected more than permuted copy", y = NULL, fill = "Feature class") +
      theme_light()

    if (!is.na(ABCfeatureImportances_List$limtFreq)) {
      pdfVarFreq <- pdfVarFreq + geom_vline(xintercept = ABCfeatureImportances_List$limtFreq, linetype = "dashed", color = "salmon")
    }
  }

  ##make returns
  return(list(
    importance_data_long = dfVarimp_Test_Boruta_actual_long_all,
  pVarimp_Test_actual = pVarimp_Test_actual,
  pdfVarFreq2 = pdfVarFreq2,
  pdfVarFreq = pdfVarFreq,
  ABCfeatureImportances_List = ABCfeatureImportances_List
  ))
})#, mc.cores = length(DataSetSizes))

names(Test_VarImps) <- DataSetSizes

#################################### Assemble plots of variable importance  ########################################################################

pTest_Varimps_all <-
  plot_grid(
    Test_VarImps$original$pVarimp_Test_actual + labs(title = "Normal data: Original"),
    Test_VarImps$engineered_0$pVarimp_Test_actual + labs(title = "Normal data: Engineered 0"),
    Test_VarImps$augmented_1_engineered$pVarimp_Test_actual + labs(title = "Normal data: Augmented 1, engineered"),
    Test_VarImps$augmented_5_engineered$pVarimp_Test_actual + labs(title = "Normal data: Augmented 5, engineered"),
    # Test_VarImps$augmented_50_engineered$pVarimp_Test_actual + labs(title = "Normal data: Augmented 50, engineered"),
# Test_VarImps$augmented_5$pVarimp_Test_actual + labs(title = "Normal data: Augmented 5"),
# Test_VarImps$augmented_engineered_1$pVarimp_Test_actual + labs(title = "Normal data: Augmented 1, engineered"),
# Test_VarImps$augmented_engineered_5$pVarimp_Test_actual + labs(title = "Normal data: Augmented 5, engineered"),
    labels = LETTERS[1:14],
    nrow = 1,
    align = "h", axis = "tb"
  )

print(pTest_Varimps_all)

ggsave(
  filename = paste0(pfad_o, pfad_r, paste0("pNormal_Varimps_aug100", ".svg")),
  plot = pTest_Varimps_all, width = 40, height = 10, limitsize = FALSE
)


pTest_Varfreqs_all <-
  plot_grid(
    Test_VarImps$original$pdfVarFreq + labs(title = "Normal data: Original"),
    Test_VarImps$engineered_0$pdfVarFreq + labs(title = "Normal data: Engineered 0"),
    Test_VarImps$augmented_1_engineered$pdfVarFreq + labs(title = "Normal data: Augmented 1, engineered"),
    Test_VarImps$augmented_5_engineered$pdfVarFreq + labs(title = "Normal data: Augmented 5, engineered"),
# Test_VarImps$augmented_5$pVarimp_Test_actual + labs(title = "Normal data: Augmented 5"),
# Test_VarImps$augmented_engineered_1$pVarimp_Test_actual + labs(title = "Normal data: Augmented 1, engineered"),
# Test_VarImps$augmented_engineered_5$pVarimp_Test_actual + labs(title = "Normal data: Augmented 5, engineered"),
    labels = LETTERS[1:14],
    nrow = 1,
    align = "h", axis = "tb"
  )

print(pTest_Varfreqs_all)

ggsave(
  filename = paste0(pfad_o, pfad_r, paste0("pnoEffect_Varfreqs_aug100", ".svg")),
  plot = pTest_Varfreqs_all, width = 20, height = 10, limitsize = FALSE
)

# pTest_Varfreqs_all_RF <-
#   plot_grid(
#     pBA_data_all,
#     Test_VarImps$original$pdfVarFreq + labs(title = "Normal data: Original"),
#     Test_VarImps$engineered_0$pdfVarFreq + labs(title = "Normal data: Engineered 0"),
#     Test_VarImps$augmented_1_engineered$pdfVarFreq + labs(title = "Normal data: Augmented 1, engineered"),
#     Test_VarImps$augmented_5_engineered$pdfVarFreq + labs(title = "Normal data: Augmented 5, engineered"),
#     # Test_VarImps$augmented_5$pVarimp_Test_actual + labs(title = "Normal data: Augmented 5"),
#     # Test_VarImps$augmented_engineered_1$pVarimp_Test_actual + labs(title = "Normal data: Augmented 1, engineered"),
#     # Test_VarImps$augmented_engineered_5$pVarimp_Test_actual + labs(title = "Normal data: Augmented 5, engineered"),
#     labels = LETTERS[1:14],
#     nrow = 1, rel_widths = c(1,1,1,1,1),
#     align = "h", axis = "tb"
#   )
# 
# print(pTest_Varfreqs_all_RF)
# 
# ggsave(
#   filename = paste0(pfad_o, pfad_r, paste0("pMouse_RF_Varfreqs_aug1000", ".svg")),
#   plot = pTest_Varfreqs_all_RF, width = 25, height = 10, limitsize = FALSE
# )
# 



Test_VarImps$original$ABCfeatureImportances_List$names_reducedFeatureSet_A
Test_VarImps$engineered_0$ABCfeatureImportances_List$names_reducedFeatureSet_A
Test_VarImps$augmented_1_engineered$ABCfeatureImportances_List$names_reducedFeatureSet_AA
Test_VarImps$augmented_5_engineered$ABCfeatureImportances_List$names_reducedFeatureSet_A
Test_VarImps$augmented_50_engineered$ABCfeatureImportances_List$names_reducedFeatureSet_A


Test_VarImps$original$ABCfeatureImportances_List$dfVars
names(xx) <- Test_VarImps$engineered_0$ABCfeatureImportances_List$dfVars$Var
barplot(xx, las = 2)

which(xx > 75)

# save.image(file = paste0(pfad_o, pfad_r, "Mouse_generativeAI_14_Sep_24.img"))
# load(file=paste0(pfad_o, pfad_r, "Mouse_generativeAI_14_Sep_24"))

