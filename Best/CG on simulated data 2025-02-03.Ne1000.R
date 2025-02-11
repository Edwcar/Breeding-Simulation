library(ASRgenomics)
library(asreml)
library(caret)
library(doParallel)
library(dplyr)

setwd("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in")

# It1

Data <- read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It1/Sim_Data_2025-02-03.txt",header = T)
G <- as.matrix(read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It1/Sim_G_Matrix_2025-02-03.txt",header = T))
Today <- Sys.Date()

Data <- Data %>%
  arrange(by = Genotype_ID)

G <- G + 0.01*diag(nrow(G))

colnames(G) <- gsub("^X", "", colnames(G))
rownames(G) <- colnames(G)

G<-G.inverse(G, sparseform = T)
G<-G$Ginv.sparse

Data$Genotype_ID<-as.factor(Data$Genotype_ID)
F1<-Data[Data$Gen=="F1",]
G0<-Data[Data$Gen=="G0",]

Genotypes<-unique(as.character(Data$Genotype_ID))
Genotypes_F1 <- unique(as.character(F1$Genotype_ID))
Genotypes_G0 <- unique(as.character(G0$Genotype_ID))

all_families <- unique(F1$Family)

# Calculate GEBVs on the whole dataset as pseudo phenotype
model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = Data,
                workspace = "500mb")

# Extract Genomic estimated breeidng values
GEBV_All <- summary(model, coef=T)$coef.random

F1$Family <- as.factor(F1$Family)
F1 <- F1 %>%
  group_by(Family) %>%
  mutate(H_WithinFam = (Pheno - mean(Pheno)))

model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = F1)

h2_val <- vpredict(model, h2~V1/(V1+V2))$Estimate

vc_val <- summary(model)$varcomp
Va_val <- vc_val[1,1]
Ve_val <- vc_val[2,1]

h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_G0, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running CG fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- G0[G0$Genotype_ID %in% Genotypes_G0[train_indices], ]
    
    # Create validation set
    validation_data <- G0[G0$Genotype_ID %in% Genotypes_G0[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- F1$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(F1$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~F1$Genotype_ID-1)
    
    model <- asreml(fixed = Pheno ~ 1,
                    random = ~vm(Genotype_ID, G),
                    residual = ~id(units),
                    maxiter = 50,
                    data = train_data)
    
    Converge<-model$converge
    
    
    # Variance components
    vc <- summary(model)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model, h2~V1/(V1+V2))$Estimate
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], F1$Pheno)
    
    # Prediction accuracy
    #Acc <- Pa/sqrt(h2_val)
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], F1$H_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], F1$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(F1$Family)) {
          fam_data <- F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}



CG1 <- bind_rows(results, .id = "Fold")

CG1$Scenario <- "CG1"
CG1$Date <- Today
CG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
CG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
CG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
CG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_CG1_", Today, ".txt")
write.table(CG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")


all_families <- unique(F1$Family)

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_F1, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG1 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- F1[F1$Genotype_ID %in% Genotypes_F1[train_indices], ]
    
    # Create validation set
    validation_data <- F1[F1$Genotype_ID %in% Genotypes_F1[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    # h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
   # vc <- summary(model_validation)$varcomp
   #  Va_val <- vc[1,1]
   #  Ve_val <- vc[2,1] 
   #  h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model_train, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    validation_data$Family <- as.factor(validation_data$Family)
    validation_data <- validation_data %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], validation_data$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(validation_data$Family)) {
          fam_data <- validation_data %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        #h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 =  Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG1 <- bind_rows(results, .id = "Fold")

WG1$Scenario <- "WG1"
WG1$Date <- as.character(Sys.Date())
WG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_WG1_", Today, ".txt")
write.table(WG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG2 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[2])
    
    # Create training set
    train_data <- Data[Data$Genotype_ID %in% Genotypes[train_indices], ]
    
    # Create validation set
    validation_data <- Data[Data$Genotype_ID %in% Genotypes[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model_train, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    #Acc <- Pa/sqrt(h2_val)
    
    Val_F1<-validation_data[validation_data$Gen=="F1",]
    Val_F1$Family <- as.factor(as.character(Val_F1$Family))
    
    # Remove F1 if only P0 is used
    Z_F1   <- model.matrix(~Val_F1$Genotype_ID-1)
    
    #model_val <- asreml(fixed = Pheno ~ 1,
    #                    random = ~vm(Genotype_ID, G),
    #                    residual = ~id(units),
    #                    maxiter = 50,
    #                    data = Val_F1)
  #  
   # 0.7 <- vpredict(model_val, h2~V1/(V1+V2))$Estimate
    
    Pa_F1<-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno)  
    
    Acc_F1 <- Pa_F1/sqrt(0.7)
    Acc_True_F1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Tbv)
    
    Val_F1$Family <- as.factor(Val_F1$Family)
    Val_F1 <- Val_F1 %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(Val_F1$Family)) {
          fam_data <- Val_F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        #h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1],
        Acc_F1 = Acc_F1,
        Acc_True_F1 =  Acc_True_F1[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG2 <- bind_rows(results, .id = "Fold")

WG2$Scenario <- "WG2"
WG2$Date <- as.character(Sys.Date())
WG2$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG2$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG2$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG2$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_WG2_", Today, ".txt")
write.table(WG2, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

#### It2 #####

Data <- read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It2/Sim_Data_2025-02-03.txt",header = T)
G <- as.matrix(read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It2/Sim_G_Matrix_2025-02-03.txt",header = T))
Today <- Sys.Date()

Data <- Data %>%
  arrange(by = Genotype_ID)

G <- G + 0.01*diag(nrow(G))

colnames(G) <- gsub("^X", "", colnames(G))
rownames(G) <- colnames(G)

G<-G.inverse(G, sparseform = T)
G<-G$Ginv.sparse

Data$Genotype_ID<-as.factor(Data$Genotype_ID)
F1<-Data[Data$Gen=="F1",]
G0<-Data[Data$Gen=="G0",]

all_families <- unique(F1$Family)

Genotypes<-unique(as.character(Data$Genotype_ID))
Genotypes_F1 <- unique(as.character(F1$Genotype_ID))
Genotypes_G0 <- unique(as.character(G0$Genotype_ID))

# Calculate GEBVs on the whole dataset as pseudo phenotype
model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = Data,
                workspace = "500mb")

# Extract Genomic estimated breeidng values
GEBV_All <- summary(model, coef=T)$coef.random

F1$Family <- as.factor(F1$Family)
F1 <- F1 %>%
  group_by(Family) %>%
  mutate(H_WithinFam = (Pheno - mean(Pheno)))

model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = F1)

h2_val <- vpredict(model, h2~V1/(V1+V2))$Estimate

vc_val <- summary(model)$varcomp
Va_val <- vc_val[1,1]
Ve_val <- vc_val[2,1]

h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_G0, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running CG fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- G0[G0$Genotype_ID %in% Genotypes_G0[train_indices], ]
    
    # Create validation set
    validation_data <- G0[G0$Genotype_ID %in% Genotypes_G0[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- F1$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(F1$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~F1$Genotype_ID-1)
    
    model <- asreml(fixed = Pheno ~ 1,
                    random = ~vm(Genotype_ID, G),
                    residual = ~id(units),
                    maxiter = 50,
                    data = train_data)
    
    Converge<-model$converge
    
    
    # Variance components
    vc <- summary(model)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model, h2~V1/(V1+V2))$Estimate
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], F1$Pheno)
    
    # Prediction accuracy
    Acc <- Pa/sqrt(h2_val)
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], F1$H_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], F1$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(F1$Family)) {
          fam_data <- F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        h2_val = h2_val,
        Pa = Pa[1],
        Acc = Acc[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


CG1 <- bind_rows(results, .id = "Fold")

CG1$Scenario <- "CG1"
CG1$Date <- Today
CG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
CG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
CG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
CG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_CG1_", Today, ".txt")
write.table(CG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")




results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_F1, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG1 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- F1[F1$Genotype_ID %in% Genotypes_F1[train_indices], ]
    
    # Create validation set
    validation_data <- F1[F1$Genotype_ID %in% Genotypes_F1[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model_train, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    #Acc <- Pa/sqrt(h2_val)
    
    validation_data$Family <- as.factor(validation_data$Family)
    validation_data <- validation_data %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], validation_data$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(validation_data$Family)) {
          fam_data <- validation_data %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
       # h2_val = h2_val,
        Pa = Pa[1],
       Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG1 <- bind_rows(results, .id = "Fold")

WG1$Scenario <- "WG1"
WG1$Date <- as.character(Sys.Date())
WG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_WG1_", Today, ".txt")
write.table(WG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG2 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- Data[Data$Genotype_ID %in% Genotypes[train_indices], ]
    
    # Create validation set
    validation_data <- Data[Data$Genotype_ID %in% Genotypes[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model_train, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    #Acc <- Pa/sqrt(h2_val)
    
    Val_F1<-validation_data[validation_data$Gen=="F1",]
    Val_F1$Family <- as.factor(as.character(Val_F1$Family))
    
    # Remove F1 if only P0 is used
    Z_F1   <- model.matrix(~Val_F1$Genotype_ID-1)
    
   # model_val <- asreml(fixed = Pheno ~ 1,
  #                      random = ~vm(Genotype_ID, G),
  #                      residual = ~id(units),
  #                      maxiter = 50,
  #                      data = Val_F1)
    
  #  0.7 <- vpredict(model_val, h2~V1/(V1+V2))$Estimate
    
    Pa_F1<-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno)  
    
    Acc_F1 <- Pa_F1/sqrt(0.7)
    Acc_True_F1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Tbv)
    
    Val_F1$Family <- as.factor(Val_F1$Family)
    Val_F1 <- Val_F1 %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(Val_F1$Family)) {
          fam_data <- Val_F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        #h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1],
        Acc_F1 = Acc_F1,
        Acc_True_F1 =  Acc_True_F1[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG2 <- bind_rows(results, .id = "Fold")

WG2$Scenario <- "WG2"
WG2$Date <- as.character(Sys.Date())
WG2$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG2$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG2$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG2$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim2_1000Ne_WG2_", Today, ".txt")
write.table(WG2, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

##### It3 #####

Data <- read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It3/Sim_Data_2025-02-03.txt",header = T)
G <- as.matrix(read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It3/Sim_G_Matrix_2025-02-03.txt",header = T))
Today <- Sys.Date()

Data <- Data %>%
  arrange(by = Genotype_ID)

G <- G + 0.01*diag(nrow(G))

colnames(G) <- gsub("^X", "", colnames(G))
rownames(G) <- colnames(G)

G<-G.inverse(G, sparseform = T)
G<-G$Ginv.sparse

Data$Genotype_ID<-as.factor(Data$Genotype_ID)
F1<-Data[Data$Gen=="F1",]
G0<-Data[Data$Gen=="G0",]

all_families <- unique(F1$Family)

Genotypes<-unique(as.character(Data$Genotype_ID))
Genotypes_F1 <- unique(as.character(F1$Genotype_ID))
Genotypes_G0 <- unique(as.character(G0$Genotype_ID))

# Calculate GEBVs on the whole dataset as pseudo phenotype
model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = Data,
                workspace = "500mb")

# Extract Genomic estimated breeidng values
GEBV_All <- summary(model, coef=T)$coef.random

F1$Family <- as.factor(F1$Family)
F1 <- F1 %>%
  group_by(Family) %>%
  mutate(H_WithinFam = (Pheno - mean(Pheno)))

model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = F1)

h2_val <- vpredict(model, h2~V1/(V1+V2))$Estimate

vc_val <- summary(model)$varcomp
Va_val <- vc_val[1,1]
Ve_val <- vc_val[2,1]

h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_G0, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running CG fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- G0[G0$Genotype_ID %in% Genotypes_G0[train_indices], ]
    
    # Create validation set
    validation_data <- G0[G0$Genotype_ID %in% Genotypes_G0[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- F1$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(F1$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~F1$Genotype_ID-1)
    
    model <- asreml(fixed = Pheno ~ 1,
                    random = ~vm(Genotype_ID, G),
                    residual = ~id(units),
                    maxiter = 50,
                    data = train_data)
    
    Converge<-model$converge
    
    
    # Variance components
    vc <- summary(model)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model, h2~V1/(V1+V2))$Estimate
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], F1$Pheno)
    
    # Prediction accuracy
    Acc <- Pa/sqrt(h2_val)
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], F1$H_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], F1$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(F1$Family)) {
          fam_data <- F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        h2_val = h2_val,
        Pa = Pa[1],
        Acc = Acc[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


CG1 <- bind_rows(results, .id = "Fold")

CG1$Scenario <- "CG1"
CG1$Date <- Today
CG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
CG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
CG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
CG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_CG1_", Today, ".txt")
write.table(CG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")


all_families <- unique(F1$Family)

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_F1, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG1 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- F1[F1$Genotype_ID %in% Genotypes_F1[train_indices], ]
    
    # Create validation set
    validation_data <- F1[F1$Genotype_ID %in% Genotypes_F1[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model_train, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    #Acc <- Pa/sqrt(h2_val)
    
    validation_data$Family <- as.factor(validation_data$Family)
    validation_data <- validation_data %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], validation_data$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(validation_data$Family)) {
          fam_data <- validation_data %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        #h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG1 <- bind_rows(results, .id = "Fold")

WG1$Scenario <- "WG1"
WG1$Date <- as.character(Sys.Date())
WG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_WG1_", Today, ".txt")
write.table(WG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG2 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- Data[Data$Genotype_ID %in% Genotypes[train_indices], ]
    
    # Create validation set
    validation_data <- Data[Data$Genotype_ID %in% Genotypes[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model_train, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    #Acc <- Pa/sqrt(h2_val)
    
    Val_F1<-validation_data[validation_data$Gen=="F1",]
    Val_F1$Family <- as.factor(as.character(Val_F1$Family))
    
    # Remove F1 if only P0 is used
    Z_F1   <- model.matrix(~Val_F1$Genotype_ID-1)
    
    #model_val <- asreml(fixed = Pheno ~ 1,
    #                    random = ~vm(Genotype_ID, G),
    #                    residual = ~id(units),
    #                    maxiter = 50,
    #                    data = Val_F1)
    
    # 0.7 <- vpredict(model_val, h2~V1/(V1+V2))$Estimate
    
    Pa_F1<-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno)  
    
    Acc_F1 <- Pa_F1/sqrt(0.7)
    Acc_True_F1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Tbv)
    
    Val_F1$Family <- as.factor(Val_F1$Family)
    Val_F1 <- Val_F1 %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(Val_F1$Family)) {
          fam_data <- Val_F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        #h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1],
        Acc_F1 = Acc_F1,
        Acc_True_F1 =  Acc_True_F1[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG2 <- bind_rows(results, .id = "Fold")

WG2$Scenario <- "WG2"
WG2$Date <- as.character(Sys.Date())
WG2$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG2$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG2$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG2$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim3_1000Ne_WG2_", Today, ".txt")
write.table(WG2, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")


#### It4 #####

Data <- read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It4/Sim_Data_2025-02-03.txt",header = T)
G <- as.matrix(read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It4/Sim_G_Matrix_2025-02-03.txt",header = T))
Today <- Sys.Date()

Data <- Data %>%
  arrange(by = Genotype_ID)

G <- G + 0.01*diag(nrow(G))

colnames(G) <- gsub("^X", "", colnames(G))
rownames(G) <- colnames(G)

G<-G.inverse(G, sparseform = T)
G<-G$Ginv.sparse

Data$Genotype_ID<-as.factor(Data$Genotype_ID)
F1<-Data[Data$Gen=="F1",]
G0<-Data[Data$Gen=="G0",]

Genotypes<-unique(as.character(Data$Genotype_ID))
Genotypes_F1 <- unique(as.character(F1$Genotype_ID))
Genotypes_G0 <- unique(as.character(G0$Genotype_ID))

# Calculate GEBVs on the whole dataset as pseudo phenotype
model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = Data,
                workspace = "500mb")

# Extract Genomic estimated breeidng values
GEBV_All <- summary(model, coef=T)$coef.random

F1$Family <- as.factor(F1$Family)
F1 <- F1 %>%
  group_by(Family) %>%
  mutate(H_WithinFam = (Pheno - mean(Pheno)))

model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = F1)

h2_val <- vpredict(model, h2~V1/(V1+V2))$Estimate

vc_val <- summary(model)$varcomp
Va_val <- vc_val[1,1]
Ve_val <- vc_val[2,1]

h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_G0, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running CG fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- G0[G0$Genotype_ID %in% Genotypes_G0[train_indices], ]
    
    # Create validation set
    validation_data <- G0[G0$Genotype_ID %in% Genotypes_G0[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- F1$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(F1$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~F1$Genotype_ID-1)
    
    model <- asreml(fixed = Pheno ~ 1,
                    random = ~vm(Genotype_ID, G),
                    residual = ~id(units),
                    maxiter = 50,
                    data = train_data)
    
    Converge<-model$converge
    
    
    # Variance components
    vc <- summary(model)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model, h2~V1/(V1+V2))$Estimate
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], F1$Pheno)
    
    # Prediction accuracy
    Acc <- Pa/sqrt(h2_val)
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], F1$H_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], F1$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(F1$Family)) {
          fam_data <- F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        h2_val = h2_val,
        Pa = Pa[1],
        Acc = Acc[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


CG1 <- bind_rows(results, .id = "Fold")

CG1$Scenario <- "CG1"
CG1$Date <- Today
CG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
CG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
CG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
CG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_CG1_", Today, ".txt")
write.table(CG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

all_families <- unique(F1$Family)


results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_F1, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG1 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- F1[F1$Genotype_ID %in% Genotypes_F1[train_indices], ]
    
    # Create validation set
    validation_data <- F1[F1$Genotype_ID %in% Genotypes_F1[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model_train, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    #Acc <- Pa/sqrt(h2_val)
    
    validation_data$Family <- as.factor(validation_data$Family)
    validation_data <- validation_data %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], validation_data$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(validation_data$Family)) {
          fam_data <- validation_data %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG1 <- bind_rows(results, .id = "Fold")

WG1$Scenario <- "WG1"
WG1$Date <- as.character(Sys.Date())
WG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_WG1_", Today, ".txt")
write.table(WG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG2 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- Data[Data$Genotype_ID %in% Genotypes[train_indices], ]
    
    # Create validation set
    validation_data <- Data[Data$Genotype_ID %in% Genotypes[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model_train, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    #Acc <- Pa/sqrt(h2_val)
    
    Val_F1<-validation_data[validation_data$Gen=="F1",]
    Val_F1$Family <- as.factor(as.character(Val_F1$Family))
    
    # Remove F1 if only P0 is used
    Z_F1   <- model.matrix(~Val_F1$Genotype_ID-1)
    
    #model_val <- asreml(fixed = Pheno ~ 1,
    #                    random = ~vm(Genotype_ID, G),
    #                    residual = ~id(units),
    #                    maxiter = 50,
    #                    data = Val_F1)
    
    # 0.7 <- vpredict(model_val, h2~V1/(V1+V2))$Estimate
    
    Pa_F1<-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno)  
    
    Acc_F1 <- Pa_F1/sqrt(0.7)
    Acc_True_F1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Tbv)
    
    Val_F1$Family <- as.factor(Val_F1$Family)
    Val_F1 <- Val_F1 %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(Val_F1$Family)) {
          fam_data <- Val_F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        # h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1],
        Acc_F1 = Acc_F1,
        Acc_True_F1 =  Acc_True_F1[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG2 <- bind_rows(results, .id = "Fold")

WG2$Scenario <- "WG2"
WG2$Date <- as.character(Sys.Date())
WG2$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG2$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG2$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG2$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim4_1000Ne_WG2_", Today, ".txt")
write.table(WG2, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")


#### It5  #####

Data <- read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It5/Sim_Data_2025-02-03.txt",header = T)
G <- as.matrix(read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It5/Sim_G_Matrix_2025-02-03.txt",header = T))
Today <- Sys.Date()

Data <- Data %>%
  arrange(by = Genotype_ID)

G <- G + 0.01*diag(nrow(G))

colnames(G) <- gsub("^X", "", colnames(G))
rownames(G) <- colnames(G)

G<-G.inverse(G, sparseform = T)
G<-G$Ginv.sparse

Data$Genotype_ID<-as.factor(Data$Genotype_ID)
F1<-Data[Data$Gen=="F1",]
G0<-Data[Data$Gen=="G0",]

Genotypes<-unique(as.character(Data$Genotype_ID))
Genotypes_F1 <- unique(as.character(F1$Genotype_ID))
Genotypes_G0 <- unique(as.character(G0$Genotype_ID))

# Calculate GEBVs on the whole dataset as pseudo phenotype
model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = Data,
                workspace = "500mb")

# Extract Genomic estimated breeidng values
GEBV_All <- summary(model, coef=T)$coef.random

F1$Family <- as.factor(F1$Family)
F1 <- F1 %>%
  group_by(Family) %>%
  mutate(H_WithinFam = (Pheno - mean(Pheno)))

model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = F1)

h2_val <- vpredict(model, h2~V1/(V1+V2))$Estimate

vc_val <- summary(model)$varcomp
Va_val <- vc_val[1,1]
Ve_val <- vc_val[2,1]

h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_G0, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running CG fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- G0[G0$Genotype_ID %in% Genotypes_G0[train_indices], ]
    
    # Create validation set
    validation_data <- G0[G0$Genotype_ID %in% Genotypes_G0[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- F1$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(F1$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~F1$Genotype_ID-1)
    
    model <- asreml(fixed = Pheno ~ 1,
                    random = ~vm(Genotype_ID, G),
                    residual = ~id(units),
                    maxiter = 50,
                    data = train_data)
    
    Converge<-model$converge
    
    
    # Variance components
    vc <- summary(model)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model, h2~V1/(V1+V2))$Estimate
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], F1$Pheno)
    
    # Prediction accuracy
    Acc <- Pa/sqrt(h2_val)
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], F1$H_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], F1$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(F1$Family)) {
          fam_data <- F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        h2_val = h2_val,
        Pa = Pa[1],
        Acc = Acc[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


CG1 <- bind_rows(results, .id = "Fold")

CG1$Scenario <- "CG1"
CG1$Date <- Today
CG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
CG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
CG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
CG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_CG1_", Today, ".txt")
write.table(CG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

all_families <- unique(F1$Family)


results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_F1, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG1 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- F1[F1$Genotype_ID %in% Genotypes_F1[train_indices], ]
    
    # Create validation set
    validation_data <- F1[F1$Genotype_ID %in% Genotypes_F1[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model_train, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    #Acc <- Pa/sqrt(h2_val)
    
    validation_data$Family <- as.factor(validation_data$Family)
    validation_data <- validation_data %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], validation_data$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(validation_data$Family)) {
          fam_data <- validation_data %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        #h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG1 <- bind_rows(results, .id = "Fold")

WG1$Scenario <- "WG1"
WG1$Date <- as.character(Sys.Date())
WG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_WG1_", Today, ".txt")
write.table(WG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG2 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- Data[Data$Genotype_ID %in% Genotypes[train_indices], ]
    
    # Create validation set
    validation_data <- Data[Data$Genotype_ID %in% Genotypes[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model_train, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    #Acc <- Pa/sqrt(h2_val)
    
    Val_F1<-validation_data[validation_data$Gen=="F1",]
    Val_F1$Family <- as.factor(as.character(Val_F1$Family))
    
    # Remove F1 if only P0 is used
    Z_F1   <- model.matrix(~Val_F1$Genotype_ID-1)
    
    #model_val <- asreml(fixed = Pheno ~ 1,
    #                    random = ~vm(Genotype_ID, G),
    #                    residual = ~id(units),
    #                    maxiter = 50,
    #                    data = Val_F1)
    
    # 0.7 <- vpredict(model_val, h2~V1/(V1+V2))$Estimate
    
    Pa_F1<-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno)  
    
    Acc_F1 <- Pa_F1/sqrt(0.7)
    Acc_True_F1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Tbv)
    
    Val_F1$Family <- as.factor(Val_F1$Family)
    Val_F1 <- Val_F1 %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(Val_F1$Family)) {
          fam_data <- Val_F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        #h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1],
        Acc_F1 = Acc_F1,
        Acc_True_F1 =  Acc_True_F1[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG2 <- bind_rows(results, .id = "Fold")

WG2$Scenario <- "WG2"
WG2$Date <- as.character(Sys.Date())
WG2$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG2$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG2$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG2$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim5_1000Ne_WG2_", Today, ".txt")
write.table(WG2, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")


##### It6 #### 

Data <- read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It6/Sim_Data_2025-02-03.txt",header = T)
G <- as.matrix(read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It6/Sim_G_Matrix_2025-02-03.txt",header = T))
Today <- Sys.Date()

Data <- Data %>%
  arrange(by = Genotype_ID)

G <- G + 0.01*diag(nrow(G))

colnames(G) <- gsub("^X", "", colnames(G))
rownames(G) <- colnames(G)

G<-G.inverse(G, sparseform = T)
G<-G$Ginv.sparse

Data$Genotype_ID<-as.factor(Data$Genotype_ID)
F1<-Data[Data$Gen=="F1",]
G0<-Data[Data$Gen=="G0",]

Genotypes<-unique(as.character(Data$Genotype_ID))
Genotypes_F1 <- unique(as.character(F1$Genotype_ID))
Genotypes_G0 <- unique(as.character(G0$Genotype_ID))

# Calculate GEBVs on the whole dataset as pseudo phenotype
model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = Data,
                workspace = "500mb")

# Extract Genomic estimated breeidng values
GEBV_All <- summary(model, coef=T)$coef.random

F1$Family <- as.factor(F1$Family)
F1 <- F1 %>%
  group_by(Family) %>%
  mutate(H_WithinFam = (Pheno - mean(Pheno)))

model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = F1)

h2_val <- vpredict(model, h2~V1/(V1+V2))$Estimate

vc_val <- summary(model)$varcomp
Va_val <- vc_val[1,1]
Ve_val <- vc_val[2,1]

h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_G0, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running CG fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- G0[G0$Genotype_ID %in% Genotypes_G0[train_indices], ]
    
    # Create validation set
    validation_data <- G0[G0$Genotype_ID %in% Genotypes_G0[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- F1$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(F1$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~F1$Genotype_ID-1)
    
    model <- asreml(fixed = Pheno ~ 1,
                    random = ~vm(Genotype_ID, G),
                    residual = ~id(units),
                    maxiter = 50,
                    data = train_data)
    
    Converge<-model$converge
    
    
    # Variance components
    vc <- summary(model)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model, h2~V1/(V1+V2))$Estimate
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], F1$Pheno)
    
    # Prediction accuracy
    Acc <- Pa/sqrt(h2_val)
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], F1$H_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], F1$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(F1$Family)) {
          fam_data <- F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        h2_val = h2_val,
        Pa = Pa[1],
        Acc = Acc[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


CG1 <- bind_rows(results, .id = "Fold")

CG1$Scenario <- "CG1"
CG1$Date <- Today
CG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
CG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
CG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
CG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_CG1_", Today, ".txt")
write.table(CG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

all_families <- unique(F1$Family)


results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_F1, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG1 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- F1[F1$Genotype_ID %in% Genotypes_F1[train_indices], ]
    
    # Create validation set
    validation_data <- F1[F1$Genotype_ID %in% Genotypes_F1[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model_train, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    # Acc <- Pa/sqrt(h2_val)
    
    validation_data$Family <- as.factor(validation_data$Family)
    validation_data <- validation_data %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], validation_data$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(validation_data$Family)) {
          fam_data <- validation_data %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        #h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG1 <- bind_rows(results, .id = "Fold")

WG1$Scenario <- "WG1"
WG1$Date <- as.character(Sys.Date())
WG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_WG1_", Today, ".txt")
write.table(WG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG2 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- Data[Data$Genotype_ID %in% Genotypes[train_indices], ]
    
    # Create validation set
    validation_data <- Data[Data$Genotype_ID %in% Genotypes[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    # Acc <- Pa/sqrt(h2_val)
    
    Val_F1<-validation_data[validation_data$Gen=="F1",]
    Val_F1$Family <- as.factor(as.character(Val_F1$Family))
    
    # Remove F1 if only P0 is used
    Z_F1   <- model.matrix(~Val_F1$Genotype_ID-1)
    
   # model_val <- asreml(fixed = Pheno ~ 1,
  #                      random = ~vm(Genotype_ID, G),
  #                      residual = ~id(units),
  #                      maxiter = 50,
  #                      data = Val_F1)
    
   # 0.7 <- vpredict(model_val, h2~V1/(V1+V2))$Estimate
    
    Pa_F1<-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno)  
    
    Acc_F1 <- Pa_F1/sqrt(0.7)
    Acc_True_F1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Tbv)
    
    Val_F1$Family <- as.factor(Val_F1$Family)
    Val_F1 <- Val_F1 %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(Val_F1$Family)) {
          fam_data <- Val_F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1],
        Acc_F1 = Acc_F1,
        Acc_True_F1 =  Acc_True_F1[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG2 <- bind_rows(results, .id = "Fold")

WG2$Scenario <- "WG2"
WG2$Date <- as.character(Sys.Date())
WG2$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG2$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG2$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG2$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim6_1000Ne_WG2_", Today, ".txt")
write.table(WG2, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")


#### It7 ####

Data <- read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It7/Sim_Data_2025-02-03.txt",header = T)
G <- as.matrix(read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It7/Sim_G_Matrix_2025-02-03.txt",header = T))
Today <- Sys.Date()

Data <- Data %>%
  arrange(by = Genotype_ID)

G <- G + 0.01*diag(nrow(G))

colnames(G) <- gsub("^X", "", colnames(G))
rownames(G) <- colnames(G)

G<-G.inverse(G, sparseform = T)
G<-G$Ginv.sparse

Data$Genotype_ID<-as.factor(Data$Genotype_ID)
F1<-Data[Data$Gen=="F1",]
G0<-Data[Data$Gen=="G0",]

Genotypes<-unique(as.character(Data$Genotype_ID))
Genotypes_F1 <- unique(as.character(F1$Genotype_ID))
Genotypes_G0 <- unique(as.character(G0$Genotype_ID))

# Calculate GEBVs on the whole dataset as pseudo phenotype
model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = Data,
                workspace = "500mb")

# Extract Genomic estimated breeidng values
GEBV_All <- summary(model, coef=T)$coef.random

F1$Family <- as.factor(F1$Family)
F1 <- F1 %>%
  group_by(Family) %>%
  mutate(H_WithinFam = (Pheno - mean(Pheno)))

model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = F1)

h2_val <- vpredict(model, h2~V1/(V1+V2))$Estimate

vc_val <- summary(model)$varcomp
Va_val <- vc_val[1,1]
Ve_val <- vc_val[2,1]

h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_G0, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running CG fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- G0[G0$Genotype_ID %in% Genotypes_G0[train_indices], ]
    
    # Create validation set
    validation_data <- G0[G0$Genotype_ID %in% Genotypes_G0[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- F1$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(F1$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~F1$Genotype_ID-1)
    
    model <- asreml(fixed = Pheno ~ 1,
                    random = ~vm(Genotype_ID, G),
                    residual = ~id(units),
                    maxiter = 50,
                    data = train_data)
    
    Converge<-model$converge
    
    
    # Variance components
    vc <- summary(model)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model, h2~V1/(V1+V2))$Estimate
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], F1$Pheno)
    
    # Prediction accuracy
    Acc <- Pa/sqrt(h2_val)
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], F1$H_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], F1$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(F1$Family)) {
          fam_data <- F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        h2_val = h2_val,
        Pa = Pa[1],
        Acc = Acc[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


CG1 <- bind_rows(results, .id = "Fold")

CG1$Scenario <- "CG1"
CG1$Date <- Today
CG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
CG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
CG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
CG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_CG1_", Today, ".txt")
write.table(CG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

all_families <- unique(F1$Family)


results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_F1, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG1 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- F1[F1$Genotype_ID %in% Genotypes_F1[train_indices], ]
    
    # Create validation set
    validation_data <- F1[F1$Genotype_ID %in% Genotypes_F1[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    # Acc <- Pa/sqrt(h2_val)
    
    validation_data$Family <- as.factor(validation_data$Family)
    validation_data <- validation_data %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], validation_data$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(validation_data$Family)) {
          fam_data <- validation_data %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        # h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG1 <- bind_rows(results, .id = "Fold")

WG1$Scenario <- "WG1"
WG1$Date <- as.character(Sys.Date())
WG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_WG1_", Today, ".txt")
write.table(WG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG2 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- Data[Data$Genotype_ID %in% Genotypes[train_indices], ]
    
    # Create validation set
    validation_data <- Data[Data$Genotype_ID %in% Genotypes[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    # Acc <- Pa/sqrt(h2_val)
    
    Val_F1<-validation_data[validation_data$Gen=="F1",]
    Val_F1$Family <- as.factor(as.character(Val_F1$Family))
    
    # Remove F1 if only P0 is used
    Z_F1   <- model.matrix(~Val_F1$Genotype_ID-1)
    
    #model_val <- asreml(fixed = Pheno ~ 1,
    #                    random = ~vm(Genotype_ID, G),
    #                    residual = ~id(units),
    #                    maxiter = 50,
    #                    data = Val_F1)
    
    #0.7 <- vpredict(model_val, h2~V1/(V1+V2))$Estimate
    
    Pa_F1<-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno)  
    
    Acc_F1 <- Pa_F1/sqrt(0.7)
    Acc_True_F1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Tbv)
    
    Val_F1$Family <- as.factor(Val_F1$Family)
    Val_F1 <- Val_F1 %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(Val_F1$Family)) {
          fam_data <- Val_F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        #h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1],
        Acc_F1 = Acc_F1,
        Acc_True_F1 =  Acc_True_F1[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG2 <- bind_rows(results, .id = "Fold")

WG2$Scenario <- "WG2"
WG2$Date <- as.character(Sys.Date())
WG2$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG2$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG2$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG2$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim7_1000Ne_WG2_", Today, ".txt")
write.table(WG2, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")


##### It8 #####

Data <- read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It8/Sim_Data_2025-02-03.txt",header = T)
G <- as.matrix(read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It8/Sim_G_Matrix_2025-02-03.txt",header = T))
Today <- Sys.Date()

Data <- Data %>%
  arrange(by = Genotype_ID)

G <- G + 0.01*diag(nrow(G))

colnames(G) <- gsub("^X", "", colnames(G))
rownames(G) <- colnames(G)

G<-G.inverse(G, sparseform = T)
G<-G$Ginv.sparse

Data$Genotype_ID<-as.factor(Data$Genotype_ID)
F1<-Data[Data$Gen=="F1",]
G0<-Data[Data$Gen=="G0",]

Genotypes<-unique(as.character(Data$Genotype_ID))
Genotypes_F1 <- unique(as.character(F1$Genotype_ID))
Genotypes_G0 <- unique(as.character(G0$Genotype_ID))

# Calculate GEBVs on the whole dataset as pseudo phenotype
model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = Data,
                workspace = "500mb")

# Extract Genomic estimated breeidng values
GEBV_All <- summary(model, coef=T)$coef.random

F1$Family <- as.factor(F1$Family)
F1 <- F1 %>%
  group_by(Family) %>%
  mutate(H_WithinFam = (Pheno - mean(Pheno)))

model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = F1)

h2_val <- vpredict(model, h2~V1/(V1+V2))$Estimate

vc_val <- summary(model)$varcomp
Va_val <- vc_val[1,1]
Ve_val <- vc_val[2,1]

h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_G0, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running CG fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- G0[G0$Genotype_ID %in% Genotypes_G0[train_indices], ]
    
    # Create validation set
    validation_data <- G0[G0$Genotype_ID %in% Genotypes_G0[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- F1$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(F1$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~F1$Genotype_ID-1)
    
    model <- asreml(fixed = Pheno ~ 1,
                    random = ~vm(Genotype_ID, G),
                    residual = ~id(units),
                    maxiter = 50,
                    data = train_data)
    
    Converge<-model$converge
    
    
    # Variance components
    vc <- summary(model)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model, h2~V1/(V1+V2))$Estimate
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], F1$Pheno)
    
    # Prediction accuracy
    Acc <- Pa/sqrt(h2_val)
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], F1$H_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], F1$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(F1$Family)) {
          fam_data <- F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        h2_val = h2_val,
        Pa = Pa[1],
        Acc = Acc[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


CG1 <- bind_rows(results, .id = "Fold")

CG1$Scenario <- "CG1"
CG1$Date <- Today
CG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
CG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
CG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
CG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_CG1_", Today, ".txt")
write.table(CG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")


all_families <- unique(F1$Family)

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_F1, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG1 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- F1[F1$Genotype_ID %in% Genotypes_F1[train_indices], ]
    
    # Create validation set
    validation_data <- F1[F1$Genotype_ID %in% Genotypes_F1[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    # Acc <- Pa/sqrt(h2_val)
    
    validation_data$Family <- as.factor(validation_data$Family)
    validation_data <- validation_data %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], validation_data$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(validation_data$Family)) {
          fam_data <- validation_data %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        #h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG1 <- bind_rows(results, .id = "Fold")

WG1$Scenario <- "WG1"
WG1$Date <- as.character(Sys.Date())
WG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_WG1_", Today, ".txt")
write.table(WG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG2 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- Data[Data$Genotype_ID %in% Genotypes[train_indices], ]
    
    # Create validation set
    validation_data <- Data[Data$Genotype_ID %in% Genotypes[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    # Acc <- Pa/sqrt(h2_val)
    
    Val_F1<-validation_data[validation_data$Gen=="F1",]
    Val_F1$Family <- as.factor(as.character(Val_F1$Family))
    
    # Remove F1 if only P0 is used
    Z_F1   <- model.matrix(~Val_F1$Genotype_ID-1)
    
    #model_val <- asreml(fixed = Pheno ~ 1,
    #                    random = ~vm(Genotype_ID, G),
    #                    residual = ~id(units),
    #                    maxiter = 50,
    #                    data = Val_F1)
    
   # 0.7 <- vpredict(model_val, h2~V1/(V1+V2))$Estimate
    
    Pa_F1<-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno)  
    
    Acc_F1 <- Pa_F1/sqrt(0.7)
    Acc_True_F1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Tbv)
    
    Val_F1$Family <- as.factor(Val_F1$Family)
    Val_F1 <- Val_F1 %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(Val_F1$Family)) {
          fam_data <- Val_F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        #h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1],
        Acc_F1 = Acc_F1,
        Acc_True_F1 =  Acc_True_F1[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG2 <- bind_rows(results, .id = "Fold")

WG2$Scenario <- "WG2"
WG2$Date <- as.character(Sys.Date())
WG2$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG2$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG2$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG2$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim8_1000Ne_WG2_", Today, ".txt")
write.table(WG2, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")


##### It9 #####

Data <- read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It9/Sim_Data_2025-02-03.txt",header = T)
G <- as.matrix(read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It9/Sim_G_Matrix_2025-02-03.txt",header = T))
Today <- Sys.Date()

Data <- Data %>%
  arrange(by = Genotype_ID)

G <- G + 0.01*diag(nrow(G))

colnames(G) <- gsub("^X", "", colnames(G))
rownames(G) <- colnames(G)

G<-G.inverse(G, sparseform = T)
G<-G$Ginv.sparse

Data$Genotype_ID<-as.factor(Data$Genotype_ID)
F1<-Data[Data$Gen=="F1",]
G0<-Data[Data$Gen=="G0",]

Genotypes<-unique(as.character(Data$Genotype_ID))
Genotypes_F1 <- unique(as.character(F1$Genotype_ID))
Genotypes_G0 <- unique(as.character(G0$Genotype_ID))

# Calculate GEBVs on the whole dataset as pseudo phenotype
model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = Data,
                workspace = "500mb")

# Extract Genomic estimated breeidng values
GEBV_All <- summary(model, coef=T)$coef.random

F1$Family <- as.factor(F1$Family)
F1 <- F1 %>%
  group_by(Family) %>%
  mutate(H_WithinFam = (Pheno - mean(Pheno)))

model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = F1)

h2_val <- vpredict(model, h2~V1/(V1+V2))$Estimate

vc_val <- summary(model)$varcomp
Va_val <- vc_val[1,1]
Ve_val <- vc_val[2,1]

h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_G0, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running CG fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- G0[G0$Genotype_ID %in% Genotypes_G0[train_indices], ]
    
    # Create validation set
    validation_data <- G0[G0$Genotype_ID %in% Genotypes_G0[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- F1$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(F1$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~F1$Genotype_ID-1)
    
    model <- asreml(fixed = Pheno ~ 1,
                    random = ~vm(Genotype_ID, G),
                    residual = ~id(units),
                    maxiter = 50,
                    data = train_data)
    
    Converge<-model$converge
    
    
    # Variance components
    vc <- summary(model)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model, h2~V1/(V1+V2))$Estimate
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], F1$Pheno)
    
    # Prediction accuracy
    Acc <- Pa/sqrt(h2_val)
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], F1$H_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], F1$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(F1$Family)) {
          fam_data <- F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        h2_val = h2_val,
        Pa = Pa[1],
        Acc = Acc[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


CG1 <- bind_rows(results, .id = "Fold")

CG1$Scenario <- "CG1"
CG1$Date <- Today
CG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
CG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
CG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
CG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_CG1_", Today, ".txt")
write.table(CG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")


all_families <- unique(F1$Family)

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_F1, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG1 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- F1[F1$Genotype_ID %in% Genotypes_F1[train_indices], ]
    
    # Create validation set
    validation_data <- F1[F1$Genotype_ID %in% Genotypes_F1[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    # Acc <- Pa/sqrt(h2_val)
    
    validation_data$Family <- as.factor(validation_data$Family)
    validation_data <- validation_data %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], validation_data$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(validation_data$Family)) {
          fam_data <- validation_data %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        #h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG1 <- bind_rows(results, .id = "Fold")

WG1$Scenario <- "WG1"
WG1$Date <- as.character(Sys.Date())
WG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_WG1_", Today, ".txt")
write.table(WG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG2 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- Data[Data$Genotype_ID %in% Genotypes[train_indices], ]
    
    # Create validation set
    validation_data <- Data[Data$Genotype_ID %in% Genotypes[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                          maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    Acc <- Pa/sqrt(0.7)
    
    Val_F1<-validation_data[validation_data$Gen=="F1",]
    Val_F1$Family <- as.factor(as.character(Val_F1$Family))
    
    # Remove F1 if only P0 is used
    Z_F1   <- model.matrix(~Val_F1$Genotype_ID-1)
    
    #model_val <- asreml(fixed = Pheno ~ 1,
    #                    random = ~vm(Genotype_ID, G),
    #                    residual = ~id(units),
    #                    maxiter = 50,
    #                    data = Val_F1)
    
  #  0.7 <- vpredict(model_val, h2~V1/(V1+V2))$Estimate
    
    Pa_F1<-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno)  
    
    Acc_F1 <- Pa_F1/sqrt(0.7)
    Acc_True_F1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Tbv)
    
    Val_F1$Family <- as.factor(Val_F1$Family)
    Val_F1 <- Val_F1 %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(Val_F1$Family)) {
          fam_data <- Val_F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        h2_val = h2_val,
        Pa = Pa[1],
        Acc = Acc[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1],
        Acc_F1 = Acc_F1,
        Acc_True_F1 =  Acc_True_F1[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG2 <- bind_rows(results, .id = "Fold")

WG2$Scenario <- "WG2"
WG2$Date <- as.character(Sys.Date())
WG2$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG2$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG2$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG2$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim9_1000Ne_WG2_", Today, ".txt")
write.table(WG2, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")


#### It10 ####

Data <- read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It10/Sim_Data_2025-02-03.txt",header = T)
G <- as.matrix(read.table("C:/Users/EDCA/OneDrive - Skogforsk/Documents/PhD Data/Simulations/Ne1000 No Burn-in/It10/Sim_G_Matrix_2025-02-03.txt",header = T))
Today <- Sys.Date()

Data <- Data %>%
  arrange(by = Genotype_ID)

G <- G + 0.01*diag(nrow(G))

colnames(G) <- gsub("^X", "", colnames(G))
rownames(G) <- colnames(G)

G<-G.inverse(G, sparseform = T)
G<-G$Ginv.sparse

Data$Genotype_ID<-as.factor(Data$Genotype_ID)
F1<-Data[Data$Gen=="F1",]
G0<-Data[Data$Gen=="G0",]

Genotypes<-unique(as.character(Data$Genotype_ID))
Genotypes_F1 <- unique(as.character(F1$Genotype_ID))
Genotypes_G0 <- unique(as.character(G0$Genotype_ID))

# Calculate GEBVs on the whole dataset as pseudo phenotype
model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = Data,
                workspace = "500mb")

# Extract Genomic estimated breeidng values
GEBV_All <- summary(model, coef=T)$coef.random

F1$Family <- as.factor(F1$Family)
F1 <- F1 %>%
  group_by(Family) %>%
  mutate(H_WithinFam = (Pheno - mean(Pheno)))

model <- asreml(fixed = Pheno ~ 1,
                random = ~vm(Genotype_ID, G),
                residual = ~id(units),
                maxiter = 50,
                data = F1)

h2_val <- vpredict(model, h2~V1/(V1+V2))$Estimate

vc_val <- summary(model)$varcomp
Va_val <- vc_val[1,1]
Ve_val <- vc_val[2,1]

h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_G0, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running CG fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- G0[G0$Genotype_ID %in% Genotypes_G0[train_indices], ]
    
    # Create validation set
    validation_data <- G0[G0$Genotype_ID %in% Genotypes_G0[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- F1$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(F1$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~F1$Genotype_ID-1)
    
    model <- asreml(fixed = Pheno ~ 1,
                    random = ~vm(Genotype_ID, G),
                    residual = ~id(units),
                    maxiter = 50,
                    data = train_data)
    
    Converge<-model$converge
    
    
    # Variance components
    vc <- summary(model)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model, h2~V1/(V1+V2))$Estimate
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], F1$Pheno)
    
    # Prediction accuracy
    Acc <- Pa/sqrt(h2_val)
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], F1$H_WithinFam)
    AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], F1$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(F1$Family)) {
          fam_data <- F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        h2_val = h2_val,
        Pa = Pa[1],
        Acc = Acc[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


CG1 <- bind_rows(results, .id = "Fold")

CG1$Scenario <- "CG1"
CG1$Date <- Today
CG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
CG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
CG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
CG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_CG1_", Today, ".txt")
write.table(CG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")


all_families <- unique(F1$Family)

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes_F1, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG1 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- F1[F1$Genotype_ID %in% Genotypes_F1[train_indices], ]
    
    # Create validation set
    validation_data <- F1[F1$Genotype_ID %in% Genotypes_F1[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    #Acc <- Pa/sqrt(h2_val)
    
    validation_data$Family <- as.factor(validation_data$Family)
    validation_data <- validation_data %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Zval%*%GEBV[,1], validation_data$Pheno_WithinFam)
     AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(validation_data$Family)) {
          fam_data <- validation_data %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        #h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG1 <- bind_rows(results, .id = "Fold")

WG1$Scenario <- "WG1"
WG1$Date <- as.character(Sys.Date())
WG1$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG1$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG1$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG1$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim1_1000Ne_WG1_", Today, ".txt")
write.table(WG1, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")

results <- list()

foreach (j = (1:20)) %do% {
  cv_folds <- createFolds(Genotypes, k = 5, list = TRUE, returnTrain = TRUE)
  # Run through each fold
  results_list <- vector("list", length = 5)
  
  foreach(i = (1:5)) %do% {
    # Remember i
    cat("Running WG2 fold number:",j , i, "\n")
    
    train_indices <- unlist(cv_folds[i])
    
    # Create training set
    train_data <- Data[Data$Genotype_ID %in% Genotypes[train_indices], ]
    
    # Create validation set
    validation_data <- Data[Data$Genotype_ID %in% Genotypes[-train_indices], ]
    
    Train_geno<-train_data$Genotype_ID
    Validation_geno <- validation_data$Genotype_ID
    
    Gen_train<-length(unique(train_data$Genotype_ID))
    Gen_val<-length(unique(validation_data$Genotype_ID))
    
    Ztrain <- model.matrix(~train_data$Genotype_ID-1)
    Zval <- model.matrix(~validation_data$Genotype_ID-1)
    
    model_train <- asreml(fixed = Pheno ~ 1,
                          random = ~vm(Genotype_ID, G),
                          residual = ~id(units),
                          maxiter = 50,
                          data = train_data)
    
    #model_validation <- asreml(fixed = Pheno ~ 1,
    #                           random = ~vm(Genotype_ID, G),
    #                           residual = ~id(units),
    #                           maxiter = 50,
    #                           data = validation_data)
    
    Converge<-model_train$converge
    
    # Variance components
    vc <- summary(model_train)$varcomp
    Va <- vc[1,1]
    Ve <- vc[2,1] 
    VaPerCh = vc[1,5]
    VePerCh = vc[2,5]
    
    # Heritability
    h2_train <- vpredict(model_train, h2~V1/(V1+V2))$Estimate
    #h2_val <- vpredict(model_validation, h2~V1/(V1+V2))$Estimate
    #vc <- summary(model_validation)$varcomp
    #Va_val <- vc[1,1]
    #Ve_val <- vc[2,1] 
    #h2_wfam = 0.5*Va_val / (0.5*Va_val + Ve_val)
    
    # GEBVs  
    GEBV <- summary(model, coef=T)$coef.random
    rownames(GEBV) <- Genotypes 
    
    # Predictive ability
    Pa<-cor(Zval%*%GEBV[,1], validation_data$Pheno)
    
    # Prediction accuracy
    # Acc <- Pa/sqrt(h2_val)
    
    Val_F1<-validation_data[validation_data$Gen=="F1",]
    Val_F1$Family <- as.factor(as.character(Val_F1$Family))
    
    # Remove F1 if only P0 is used
    Z_F1   <- model.matrix(~Val_F1$Genotype_ID-1)
    
   # model_val <- asreml(fixed = Pheno ~ 1,
   #                      random = ~vm(Genotype_ID, G),
   #                     residual = ~id(units),
   #                    maxiter = 50,
   #                     data = Val_F1)
    
    #0.7 <- vpredict(model_val, h2~V1/(V1+V2))$Estimate
    
    Pa_F1<-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno)  
    
    Acc_F1 <- Pa_F1/sqrt(0.7)
    Acc_True_F1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Tbv)
    
    Val_F1$Family <- as.factor(Val_F1$Family)
    Val_F1 <- Val_F1 %>%
      group_by(Family) %>%
      mutate(Pheno_WithinFam = (Pheno - mean(Pheno)))
    
    # Family PA 
    PaFam1 <-cor(Z_F1%*%GEBV[,1], Val_F1$Pheno_WithinFam)
     AccWfam1 <- PaFam1 / sqrt(h2_wfam)
    
    Acc.h2.70 <- Pa/sqrt(0.70)
    Acc.PsPh<-cor(Zval%*%GEBV[,1], Zval%*%GEBV_All[,1])
    Acc_True<-cor(Zval%*%GEBV[,1], validation_data$Tbv)
    
    PaFam <- setNames(
      sapply(all_families, function(fam) {
        if (fam %in% as.character(Val_F1$Family)) {
          fam_data <- Val_F1 %>% filter(as.character(Family) == fam)  # Ensure matching type
          Zfam <- model.matrix(~ fam_data$Genotype_ID - 1)
          cor(Zfam %*% GEBV[, 1], fam_data$Pheno)
        } else {
          NA  # Assign NA if the family is missing in this iteration
        }
      }),
      all_families  # Set all family names as column names
    )
    
    PaFam2<-mean(PaFam, na.rm = T)
    
    # Store results, ensuring PaFam values are named columns
    results[[paste0("Fold_", j, "_", i)]] <- c(
      list(
        Converge = Converge,
        Gen_train = Gen_train,
        Gen_val = Gen_val,
        Va = Va, 
        Ve = Ve,
        VaPerCh = VaPerCh,
        VePerCh = VePerCh,
        h2_train = h2_train,
        #h2_val = h2_val,
        Pa = Pa[1],
        Acc.h2.70 = Acc.h2.70[1],
        PaFam1 = PaFam1[1],
        PaFam2 = PaFam2,
        AccWfam1 = AccWfam1[1],
        Acc.PsPh = Acc.PsPh[1],
        Acc_True = Acc_True[1],
        Acc_F1 = Acc_F1,
        Acc_True_F1 =  Acc_True_F1[1]
      ),
      as.list(PaFam) # Convert named vector to list, creating separate columns
    )
  }
}


WG2 <- bind_rows(results, .id = "Fold")

WG2$Scenario <- "WG2"
WG2$Date <- as.character(Sys.Date())
WG2$Fixed <- as.character(paste((model[["formulae"]][["fixed"]]),collapse = ""))
WG2$Random <- model[["G.param"]][["vm(Genotype_name, G)"]][["vm(Genotype_name, G)"]][["facnam"]]
WG2$Residual <- as.character(paste((model[["formulae"]][["residual"]][[2]][]), collapse = ""))
WG2$Genotypes <- length(unique(G0$Genotype_name))

filename <- paste0("Sim10_1000Ne_WG2_", Today, ".txt")
write.table(WG2, file = filename, quote = F, col.names = T, row.names = F, sep = "\t")
