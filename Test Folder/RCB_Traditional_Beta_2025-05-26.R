# AlphaSim Simulation 

# Date of simulation
today_date <- paste(format(Sys.Date(), "%Y-%m-%d"))

#### Load Packages ####
library(AlphaSimR)
library(pedigree)
library(rrBLUP)
library(ASRgenomics)
library(caret)
library(doParallel)
library(dplyr)
library(pedigree)
library(gridExtra)
library(purrr)
library(ggplot2)
library(tidyr)
library(kinship2)
library(AGHmatrix)
library(breedR)
library(ggridges)
library(bWGR)

############################### Set initial values #############################

# Genetic Parameters
# Nr of chromosomes 
PaChr = 12 # Source: Nystedt et al., 2013 https://doi.org/10.1038/nature12211

# Chromosome length in bp 
ChrLen = ((19.6*10^9)/PaChr) # Source: Source: Nystedt et al., 2013 https://doi.org/10.1038/nature12211

# Chromosome length in morgans 
PG_Chr_M <- ((3556/100)/PaChr) # Source: Bernhardsson et. al. 2019. https://doi.org/10.1534/g3.118.200840

# Average age age at reproduction during burn-in
MatAge = 45         

# Mutation rate between generations
ConMR = (2.2*10^-9) * MatAge # Source: Nystedt et al., 2013 https://doi.org/10.1038/nature12211

# Number of segregating sites to simulate
SegSite = 10000

# Trait Parameters
# Need two if using multivariate as a way to integrate two measurements 
# Trait
Trait = "A"         # Make the trait completely additive for simplicity

# Nr of QTLs per chromosome
nQtl =  round(500/PaChr)   #Source: Hall et al., (2016) https://doi.org/10.1007/s11295-016-1073-0
# Assume growth trait

#Distribution of QTL effects
GAMMA = T 
GShape = 0.5
# Source: Hall et al., 2016 https://doi.org/10.1007/s11295-016-1073-0
#         Li et al., 2018 https://doi.org/10.1371/journal.pone.0208232

# Set mean value and varaince to mimik real data
# Height at 7 and 15 years

initMeanG = c(200,800)  
CV = 0.3 # The standard deviation ~ 30 % of the mean 
initVarG = (CV * initMeanG)^2  # Initial trait genetic variance

# Dominance and Inbreeding Depression
ID = -0.35 # Source:
meanDD = ID/(2*nQtl*PaChr) # 2 time the number of QTLs
initMeanD = c(meanDD*initMeanG[1],meanDD*initMeanG[2])
initVarD = initVarG*0.2

# Error variance
initVarE = 4*(CV * initMeanG)^2

# Additve genetic correlation between traits
CorA <- 0.9
# Phenotypic correlation between traits
CorP <- 0.8

# Heritability
h2_G0 = 0.2
h2_G1 = 0.2

# Rate of inbreeding deprssion
InbDepr = -0.5 # Source:


# Founder Parameters    
# Founder Population
nFounders = 1000   # Number of founders in the population
neFounders = 10   # Effective population size of founders

nGenBurnIn <- 100  # Number of descrete burn-in generations
nCrossBurnIn <- 100 # Number of crosses within each burn-in generation
nProgBurnIn <- 10

# Breeding Parameters    
# Size of the G0 population
nG0 = 500

# Crosses
nParents = 10   # Number of parents to be selected
nCross = 20     # Number of families to create from Parents
nProg = 10      # Number of progeny per cross

# Mating schemes
MatePlan = "RandCross"

# Phenotyping efforts
sites = 3
ramets = 4 # Ramets per site

# Number of replicates
n_reps = sites * ramets 

# Number of generations to simulate
# Define number of generations

# Genotyping Parameters
# SNP chip for GS
nSNP = 100     # Nr of SNPs per chromosome
minMAF = 0.01  # Minimum MAF for inclusion

# "Imaginary" SNp chip for estimating true inbreeding coefficients
nSNP2 = 300

############################## Functions #######################################
# Calculating MAF
calculate_maf <- function(column) {
  
  total_alleles <- 2 * length(column)
  
  minor_allele_count <- sum(column, na.rm = T)
  
  minor_allele_frequency <- minor_allele_count / total_alleles
  
  maf <- min(minor_allele_frequency, 1 - minor_allele_frequency)
  
  return(maf)
}


################################################################################


###### Initiate Founder Population ######
founderPop = runMacs2(nInd = nFounders,         
                      nChr = PaChr,
                      segSites = SegSite,
                      inbred   = FALSE, 
                      Ne = neFounders,
                      genLen = PG_Chr_M,
                      mutRate = ConMR,
                      bp = ChrLen,
                      ploidy = 2L)

SP = SimParam$new(founderPop)

# Global settings
SP$restrSegSites(overlap = FALSE) # Ensure QTLs and SNPs do not overlap
SP$setTrackPed(TRUE)

# Add trait
SP$addTraitAD(nQtlPerChr   = nQtl,
             mean         = initMeanG,
             meanDD       = initMeanD, 
             var          = initVarG,
             varDD        = initVarD,
             gamma        = T,
             shape        = GShape,
             corA         = matrix(c(1, CorA,
                                     CorA, 1), nrow=2))
#             corDD        = matrix(c(1, CorA,
#                              CorA, 1), nrow=2)) # Add AxA correlation between two height measurements

# Define heritabilities for both traits (this is REQUIRED before using varE)
SP$setVarE(h2 = c(0.2, 0.2))

# Set residual covariance matrix used by the setPheno function
SP$setVarE(varE = matrix(c(1, 0.8,
                           0.8, 1), nrow=2))

# Create new pop object
G0_pop = newPop(founderPop)

Acc_df <- data.frame(Trait1 = 0,
                     Trait2 = 0)

################################## Burn-in #####################################
VarA_list1<-list()
VarA_G0 <- varA(G0_pop)
VarA_list1 <- c(VarA_list1, VarA_G0[1,1])

VarD_list1<-list()
VarD_G0 <- varD(G0_pop)
VarD_list1 <- c(VarD_list1, VarD_G0[1,1])

MeanA_list1<-list()
MeanA_G0 <- meanG(G0_pop)
MeanA_list1 <- c(MeanA_list1, MeanA_G0[1])

VarA_list2<-list()
VarA_G0 <- varA(G0_pop)
VarA_list2 <- c(VarA_list2, VarA_G0[2,2])

VarD_list2<-list()
VarD_G0 <- varD(G0_pop)
VarD_list2 <- c(VarD_list2, VarD_G0[2,2])

MeanA_list2<-list()
MeanA_G0 <- meanG(G0_pop)
MeanA_list2 <- c(MeanA_list2, MeanA_G0[2])

for (i in 1:nGenBurnIn) {
  # Perform random crossing
  G0_pop <- randCross(
    pop = G0_pop,            # The population to cross
    nCrosses = nCrossBurnIn, # Total number of crosses
    nProgeny = nProgBurnIn,  # Number of progeny per cross
    ignoreSexes = TRUE)      # Ignore sexes in crossing
  VarA_G0 <- varA(G0_pop)
  VarA_list1 <- c(VarA_list1, VarA_G0[1,1])
  VarD_G0 <- varD(G0_pop)
  VarD_list1 <- c(VarD_list1, VarD_G0[1,1])
  MeanA_G0 <- meanG(G0_pop)
  MeanA_list1 <- c(MeanA_list1, MeanA_G0[1])
  
  VarA_G0 <- varA(G0_pop)
  VarA_list2 <- c(VarA_list2, VarA_G0[2,2])
  VarD_G0 <- varD(G0_pop)
  VarD_list2 <- c(VarD_list2, VarD_G0[2,2])
  MeanA_G0 <- meanG(G0_pop)
  MeanA_list2 <- c(MeanA_list2, MeanA_G0[2])
}


time <- seq_along(VarA_list) 

df <- data.frame(
  Time = rep(time, 3),  # Repeat time for both parameters
  Value = c(unlist(VarA_list1), unlist(MeanA_list1), unlist(VarD_list1)),  # Combine the values of both parameters
  Parameter = rep(c("VarA", "MeanG", "VarD"), each = length(time))  # Label each parameter
)

ggplot(df, aes(x = Time, y = Value, color = Parameter)) +
  geom_line(linewidth = 1) +         # Lines for both parameters
  geom_point(size = 1) + 
  scale_color_manual( values = c("VarA" = "darkgreen", "VarD" = "orange", "MeanG" = "darkred")  # Custom colors
  ) +# Points for both parameters
  labs(title = "",
       x = "Generation",
       y = "Value",
       color = "Parameter"
  ) +
  theme_minimal()

df <- data.frame(
  Time = rep(time, 3),  # Repeat time for both parameters
  Value = c(unlist(VarA_list2), unlist(MeanA_list2), unlist(VarD_list2)),  # Combine the values of both parameters
  Parameter = rep(c("VarA", "MeanG", "VarD"), each = length(time))  # Label each parameter
)

ggplot(df, aes(x = Time, y = Value, color = Parameter)) +
  geom_line(linewidth = 1) +         # Lines for both parameters
  geom_point(size = 1) + 
  scale_color_manual( values = c("VarA" = "darkgreen", "VarD" = "orange", "MeanG" = "darkred")  # Custom colors
  ) +# Points for both parameters
  labs(title = "",
       x = "Generation",
       y = "Value",
       color = "Parameter"
  ) +
  theme_minimal()

DDD<-dd(G0_pop)
genicVarD(G0_pop)
varD(G0_pop)
varA(G0_pop)
genicVarA(G0_pop)

###################### G0 Founder generation ######################################


# Select clones for first real G0 generation
G0_pop <- selectInd(G0_pop, nG0, use = "rand")


# Save pedigree
Pedigree_All <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G0_pop@id,])
Pedigree_All$mother <- 0 # Reset pedigree for Plus trees
Pedigree_All$father <- 0 # Reset pedigree for Plus trees


Kinship_matrix<- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                           dadid = Pedigree_All$father,
                           momid = Pedigree_All$mother)


#### Add SNP chips ####

# Randomly assigns eligible SNPs to SNP chip
# Use the latest population as reference
SP$addSnpChip(nSNP, minSnpFreq = minMAF,name = "GS", refPop = G0_pop)
SP$addSnpChip(nSNP2, minSnpFreq = minMAF,name = "Inb", refPop = G0_pop)

# Save MAF
G_new<-pullSnpGeno(G0_pop, snpChip = "Inb")
maf_values_df <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                            Gen = "G0")

# Add first phenotypes
G0_Pheno <- data.frame()

for (i in 1:1) {
  pheno <- setPheno(G0_pop,
                    varE = initVarE,
                    corE = matrix(c(1, CorP,
                                    CorP, 1), nrow = 2),
                    rep = 1,
                    simParam = NULL,
                    onlyPheno = T)
  
  G0_Pheno <- rbind(
    G0_Pheno,
    data.frame(ID = G0_pop@id, Rep = i, Pheno = pheno)
  )
}

G0_Pheno$Rep<-NULL

G0_Pheno_mean<- G0_Pheno %>%
  group_by(ID) %>%
  summarise(
    Trait1 = mean(Pheno.Trait1, na.rm = TRUE),
    Trait2 = mean(Pheno.Trait2, na.rm = TRUE)
  )

G0_Pheno_mean$ID <- as.character(G0_Pheno_mean$ID)
GLOBAL_Phenotypes <- G0_Pheno

# Reorder df_means based on the rownames of K
G0_Pheno_mean <- G0_Pheno_mean[match(rownames(Kinship_matrix), G0_Pheno_mean$ID), ]

y<-as.matrix(G0_Pheno_mean[,2:3])

#### Model ####
new_model<-mkr(y,Kinship_matrix)

# Save EBVs
EBVs<-as.data.frame(new_model$hat)
EBVs$ID <- G0_Pheno_mean$ID

EBVs <- EBVs[EBVs$ID %in% G0_pop@id,]
new_Pheno <- G0_Pheno_mean[G0_Pheno_mean$ID %in% G0_pop@id,]

G0_pop@ebv <- as.matrix(EBVs[,1:2])
G0_pop@pheno <- as.matrix(new_Pheno[,2:3])

##### Select parents ######
Parents <- selectWithinFam(G0_pop,
                           nInd = 1,
                           use = "ebv",
                           selectTop = TRUE)

# Heritability
h2<-var(bv(G0_pop))/varP(G0_pop)
h2_hat<-new_model$h2

True_h2_1 = h2[1,1]
True_h2_2 = h2[2,2]
Est_h2_1 = h2_hat[1]
Est_h2_2 = h2_hat[2]

# Additive genetic value   
MeanG <- meanG(G0_pop)  
VarA <- varA(G0_pop)
Cor_A = cov2cor(VarA)[2, 1]

# Dominance genetic value
VarD <- varD(G0_pop)

# Phenotypes
Phen.Mean <- meanP(G0_pop)  
VarP <- varP(G0_pop)
Cor_P = cov2cor(VarP)[2, 1]

# Prediction accuracy 
Acc <- cor(bv(G0_pop),G0_pop@ebv) 

# Estimate Bias
Bv <- bv(G0_pop)
Bias_1 <- mean(G0_pop@ebv[,1] - Bv[1])
Bias_2 <- mean(G0_pop@ebv[,2] - Bv[2])

# Estimate disversion 
bv1 <- Bv[,1]
bv2 <- Bv[,2]
ebv1 <- G0_pop@ebv[,1]
ebv2 <- G0_pop@ebv[,2]

# Fit regression model
bias_model1 <- lm(bv1 ~ ebv1)
bias_model2 <- lm(bv2 ~ ebv2)

# Optional: view slope too
slope1 <- coef(bias_model1)[2]
slope2 <- coef(bias_model2)[2]

# Selection Intensity
#Estimated Si
mean_EBV_p_1 <- mean(Parents@ebv[,1])
mean_EBV_p_2 <- mean(Parents@ebv[,2])

mean_EBV_1 <- mean(G0_pop@ebv[,1])
mean_EBV_2 <- mean(G0_pop@ebv[,2])

sd_EBV_1 <- sd(G0_pop@ebv[,1])
sd_EBV_2 <- sd(G0_pop@ebv[,2])

EBV_Intensity_1 = ((mean_EBV_p_1 - mean_EBV_1) / sd_EBV_1)
EBV_Intensity_2 = ((mean_EBV_p_2 - mean_EBV_2) / sd_EBV_2)

#True SI
bv_P<-bv(Parents)

mean_TBV_p_1 <- mean(bv_P[,1])
mean_TBV_p_2 <- mean(bv_P[,2])

mean_TBV_1 <- mean(Bv[,1])
mean_TBV_2 <- mean(Bv[,2])

sd_TBV_1 <- sd(Bv[,1])
sd_TBV_2 <- sd(Bv[,2])

TBV_Intensity_1 = ((mean_TBV_p_1 - mean_TBV_1) / sd_TBV_1)
TBV_Intensity_2 = ((mean_TBV_p_2 - mean_TBV_2) / sd_TBV_2)

GLOBAL_trait <- data.frame(Trait1_h2_est = Est_h2_1,
                        Trait1_h2 = True_h2_1,
                        Trait1_MeanG = MeanG[1],
                        Trait1_VarA = VarA[1,1],
                        Trait1_VarD = VarD[1,1],
                        Trait1_MeanP = Phen.Mean[1],
                        Trait1_VarP = VarP[1,1],
                        Trait1_Acc = Acc[1,1],
                        Trait1_WF_Acc = 0,
                        Trait1_Bias = Bias_1,
                        Trait1_Dispersion = slope1,
                        Trait1_SI_est = EBV_Intensity_1,
                        Trait1_SI = TBV_Intensity_1,
                        
                        Trait2_h2_est = Est_h2_2,
                        Trait2_h2 = True_h2_2,
                        Trait2_MeanG = MeanG[2],
                        Trait2_VarA = VarA[2,2],
                        Trait2_VarD = VarD[2,2],
                        Trait2_MeanP = Phen.Mean[2],
                        Trait2_VarP = VarP[2,2],
                        Trait2_Acc = Acc[2,2],
                        Trait2_WF_Acc = 0,
                        Trait2_Bias = Bias_2,
                        Trait2_Dispersion = slope2,
                        Trait2_SI_est = EBV_Intensity_2,
                        Trait2_SI = TBV_Intensity_2,
                        
                        CovA = VarA[1,2],
                        CorA = Cor_A,
                        CovP = VarP[1,2],
                        CorP = Cor_P,
                        
                        Gen = "G0")

Inbreeding<-mean(diag(Kinship_matrix))-1

# Plot regresson of inbreeding on phenotype

Inb<-diag(Kinship_matrix)-1

p1<-G0_pop@pheno[,1]

ID_model <- lm(p1 ~ Inb)
slope <- coef(ID_model)[2]


# Optional: view slope too
slope1 <- coef(bias_model1)[2]


GLOBAL_diversity <- data.frame(F.inb = Inbreeding,
                               Ns = 0,
                               Co = 0,
                               ID_est1 = 0,
                               ID_est2 = 0,
                               Gen = "G0")

### G1 ####

  # Generate progeny
  # Each cross gets a random numner of progenies, with minimum 1
  G_pops <- vector("list", nCross)
  for (i in 1:nCross) {
    Prog_number <- max(round(rnorm(1, mean = nProg, sd = 10)), 1)
    G_pops[[i]] <- randCross(
      pop = Parents, 
      nCrosses = 1, 
      nProgeny = Prog_number, 
      ignoreSexes = TRUE
    )
  }
  
  # Combine all families into one pop-object
    G1_pop <- do.call(c, G_pops)
    
  # Save New breeding values
  Bv <- bv(G1_pop)

  # Update pedigree
  Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G1_pop@id,])
  Pedigree_All <- rbind(Pedigree_All, Pedigree_New)
  
  A_Mat <- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                       dadid = Pedigree_All$father,
                       momid = Pedigree_All$mother)

  
  ####################### Calculate MAF Distributoin ###########################
  G_new<-pullSnpGeno(G1_pop, snpChip = "Inb")
  maf_values_new <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                               Gen = "G1")
  
  maf_values_df <- rbind(maf_values_df,maf_values_new) 
  
  ################ Assign phenotypic values ####################################
    New_Pheno <- data.frame()
    
  # Should I use the same initial error variance here?
  
  
    for (i in 1:10) {
      pheno <- setPheno(G1_pop,
                        varE = initVarE,
                        corE = matrix(c(1, CorP,
                                        CorP, 1), nrow = 2),
                        rep = 1,
                        simParam = NULL,
                        onlyPheno = T)
      New_Pheno <- rbind(
        New_Pheno,
        data.frame(ID = G1_pop@id, Rep = i, Pheno = pheno)
      )
    }
    
    New_Pheno$Rep<-NULL
    
    GLOBAL_Phenotypes<-rbind(GLOBAL_Phenotypes,New_Pheno)
    
    New_Pheno_mean<- GLOBAL_Phenotypes %>%
      group_by(ID) %>%
      summarise(
        Trait1 = mean(Pheno.Trait1, na.rm = TRUE),
        Trait2 = mean(Pheno.Trait2, na.rm = TRUE)
      )
    
    New_Pheno_mean$ID <- as.character(New_Pheno_mean$ID)
    dim(A_Mat)
    
    # Reorder df_means based on the rownames of K
    New_Pheno_mean <- New_Pheno_mean[match(rownames(A_Mat), New_Pheno_mean$ID), ]
    
    y<-as.matrix(New_Pheno_mean[,2:3])
    
    #### Model ####
    new_model<-mkr(y,A_Mat)

    # Save EBVs
    EBVs<-as.data.frame(new_model$hat)
    EBVs$ID <- New_Pheno_mean$ID
    
    EBVs <- EBVs[EBVs$ID %in% G1_pop@id,]
    new_Pheno <- New_Pheno_mean[New_Pheno_mean$ID %in% G1_pop@id,]
    
    G1_pop@ebv <- as.matrix(EBVs[,1:2])
    G1_pop@pheno <- as.matrix(new_Pheno[,2:3])
    
    ##### Select parents ######
    Parents <- selectWithinFam(G1_pop,
                               nInd = 1,
                               use = "ebv",
                               selectTop = TRUE)
    
    # Heritability
    h2<-var(bv(G1_pop))/varP(G1_pop)
    h2_hat<-new_model$h2
    
    True_h2_1 = h2[1,1]
    True_h2_2 = h2[2,2]
    Est_h2_1 = h2_hat[1]
    Est_h2_2 = h2_hat[2]
    
    # Additive genetic value   
    MeanG <- meanG(G1_pop)  
    VarA <- varA(G1_pop)
    Cor_A = cov2cor(VarA)[2, 1]
    
    # Dominance genetic value
    VarD <- varD(G1_pop)
    
    # Phenotypes
    Phen.Mean <- meanP(G1_pop)  
    VarP <- varP(G1_pop)
    Cor_P = cov2cor(VarP)[2, 1]
    
    # Prediction accuracy 
    Acc <- cor(bv(G1_pop),G1_pop@ebv) 
    
    # Within-family accuracy
    WF_df <- data.frame(ID = G1_pop@id,
                        Mum = G1_pop@mother,
                        Dad = G1_pop@father,
                        TBV1 = Bv[,1],
                        TBV2 = Bv[,2],
                        EBV1 = G1_pop@ebv[,1],
                        EBV2 = G1_pop@ebv[,2])
    WF_df$Family <- as.integer(interaction(WF_df$Mum, WF_df$Dad, drop = TRUE))
    
    summary_df <- WF_df %>%
      group_by(Family) %>%
      summarize(
        famSize = n(),
        accuracy1 = if(n() > 1) cor(TBV1, EBV1) else NA_real_,
        accuracy2 = if(n() > 1) cor(TBV2, EBV2) else NA_real_,
      ) %>%
      ungroup()
    
    
    WF_Acc1<-mean(summary_df$accuracy1, na.rm = T)
    WF_Acc2<-mean(summary_df$accuracy2, na.rm = T)
    
    T1_plot<-ggplot(summary_df, aes(x = Family, y = accuracy1, size = famSize)) +
      geom_point(alpha = 0.4, color = "red") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "G1",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    T2_plot<-ggplot(summary_df, aes(x = Family, y = accuracy2, size = famSize)) +
      geom_point(alpha = 0.4, color = "darkred") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    grid.arrange(T1_plot,T2_plot)

    # Prediction Bias
    Bv <- bv(G1_pop)
    Bias_1 <- mean(G1_pop@ebv[,1] - Bv[1])
    Bias_2 <- mean(G1_pop@ebv[,2] - Bv[2])
    
    # Estimate disversion 
    bv1 <- Bv[,1]
    bv2 <- Bv[,2]
    ebv1 <- G1_pop@ebv[,1]
    ebv2 <- G1_pop@ebv[,2]
    
    # Fit regression model
    bias_model1 <- lm(bv1 ~ ebv1)
    bias_model2 <- lm(bv2 ~ ebv2)
    
    # Optional: view slope too
    slope1 <- coef(bias_model1)[2]
    slope2 <- coef(bias_model2)[2]
    
    # Selection Intensity
    #Estimated Si
    mean_EBV_p_1 <- mean(Parents@ebv[,1])
    mean_EBV_p_2 <- mean(Parents@ebv[,2])
    
    mean_EBV_1 <- mean(G1_pop@ebv[,1])
    mean_EBV_2 <- mean(G1_pop@ebv[,2])
    
    sd_EBV_1 <- sd(G1_pop@ebv[,1])
    sd_EBV_2 <- sd(G1_pop@ebv[,2])
    
    EBV_Intensity_1 = ((mean_EBV_p_1 - mean_EBV_1) / sd_EBV_1)
    EBV_Intensity_2 = ((mean_EBV_p_2 - mean_EBV_2) / sd_EBV_2)
    
    #True SI
    bv_P<-bv(Parents)
    
    mean_TBV_p_1 <- mean(bv_P[,1])
    mean_TBV_p_2 <- mean(bv_P[,2])
    
    mean_TBV_1 <- mean(Bv[,1])
    mean_TBV_2 <- mean(Bv[,2])
    
    sd_TBV_1 <- sd(Bv[,1])
    sd_TBV_2 <- sd(Bv[,2])
    
    TBV_Intensity_1 = ((mean_TBV_p_1 - mean_TBV_1) / sd_TBV_1)
    TBV_Intensity_2 = ((mean_TBV_p_2 - mean_TBV_2) / sd_TBV_2)
    
    new_trait <- data.frame(Trait1_h2_est = Est_h2_1,
                               Trait1_h2 = True_h2_1,
                               Trait1_MeanG = MeanG[1],
                               Trait1_VarA = VarA[1,1],
                               Trait1_VarD = VarD[1,1],
                               Trait1_MeanP = Phen.Mean[1],
                               Trait1_VarP = VarP[1,1],
                               Trait1_Acc = Acc[1,1],
                               Trait1_WF_Acc = WF_Acc1,
                               Trait1_Bias = Bias_1,
                               Trait1_Dispersion = slope1,
                               Trait1_SI_est = EBV_Intensity_1,
                               Trait1_SI = TBV_Intensity_1,
                               
                               Trait2_h2_est = Est_h2_2,
                               Trait2_h2 = True_h2_2,
                               Trait2_MeanG = MeanG[2],
                               Trait2_VarA = VarA[2,2],
                               Trait2_VarD = VarD[2,2],
                               Trait2_MeanP = Phen.Mean[2],
                               Trait2_VarP = VarP[2,2],
                               Trait2_Acc = Acc[2,2],
                               Trait2_WF_Acc = WF_Acc2,
                               Trait2_Bias = Bias_2,
                               Trait2_Dispersion = slope2,
                               Trait2_SI_est = EBV_Intensity_2,
                               Trait2_SI = TBV_Intensity_2,
                               
                               CovA = VarA[1,2],
                               CorA = Cor_A,
                               CovP = VarP[1,2],
                               CorP = Cor_P,
                               
                               Gen = "G1")
    
    
    
    GLOBAL_trait <- rbind(GLOBAL_trait, new_trait)
    
    # Diversity metrics
    indices<-which((rownames(A_Mat) %in% G1_pop@id))
    A_Mat <-A_Mat[indices,indices]
    n <- nrow(A_Mat)
    
    group_coancestry <- sum(A_Mat) / (n^2)
    
    status_number <- 1 / (2*group_coancestry)
    
    Ns_new <- data.frame(Status_number = status_number,
                         Coancestry = group_coancestry,
                         Gen = "G1")

    
    Inbreeding<-mean(diag(A_Mat))-1
    
    GLOBAL_ID <- data.frame(Inb = (diag(A_Mat)-1),
                        Trait1 = G1_pop@pheno[,1],
                        Trait2 = G1_pop@pheno[,2],
                        Gen = "G1")
    
    Inb<-(diag(A_Mat)-1)
    
    Inb_df <- data.frame(I = Inb,
                         T1 = G1_pop@pheno[,1],
                         T2 = G1_pop@pheno[,2])
    
    model <- lm(T1 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID1_est<-slope/mean(Inb_df$T1)
    
    model <- lm(T2 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID2_est<-slope/mean(Inb_df$T2)
    
    new_diversity <- data.frame(F.inb = Inbreeding,
                                   Ns = status_number,
                                   Co = group_coancestry,
                                ID_est1 = 0,
                                ID_est2 = 0,
                                   Gen = "G1")
    
    GLOBAL_diversity <- rbind(GLOBAL_diversity, new_diversity)
    
#### G2 #####
  
    # Generate progeny
    # Each cross gets a random numner of progenies, with minimum 1
    G_pops <- vector("list", nCross)
    for (i in 1:nCross) {
      Prog_number <- max(round(rnorm(1, mean = nProg, sd = 10)), 1)
      G_pops[[i]] <- randCross(
        pop = Parents, 
        nCrosses = 1, 
        nProgeny = Prog_number, 
        ignoreSexes = TRUE
      )
    }
    
    # Combine all families into one pop-object
    G2_pop <- do.call(c, G_pops)
    
    
    Bv <- bv(G2_pop)
    
    # Update pedigree
    Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G2_pop@id,])
    Pedigree_All <- rbind(Pedigree_All, Pedigree_New)
    
    A_Mat <- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                       dadid = Pedigree_All$father,
                       momid = Pedigree_All$mother)
    
    
    ####################### Calculate MAF Distributoin ###########################
    G_new<-pullSnpGeno(G2_pop, snpChip = "Inb")
    maf_values_new <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                                 Gen = "G2")
    
    maf_values_df <- rbind(maf_values_df,maf_values_new) 
    
    ################ Assign phenotypic values ####################################
    New_Pheno <- data.frame()
    
    # Should I use the same initial error variance here?
    
    
    for (i in 1:10) {
      pheno <- setPheno(G2_pop,
                        varE = initVarE,
                        corE = matrix(c(1, CorP,
                                        CorP, 1), nrow = 2),
                        rep = 1,
                        simParam = NULL,
                        onlyPheno = T)
      New_Pheno <- rbind(
        New_Pheno,
        data.frame(ID = G2_pop@id, Rep = i, Pheno = pheno)
      )
    }
    
    New_Pheno$Rep<-NULL
    
    GLOBAL_Phenotypes<-rbind(GLOBAL_Phenotypes,New_Pheno)
    
    New_Pheno_mean<- GLOBAL_Phenotypes %>%
      group_by(ID) %>%
      summarise(
        Trait1 = mean(Pheno.Trait1, na.rm = TRUE),
        Trait2 = mean(Pheno.Trait2, na.rm = TRUE)
      )
    
    New_Pheno_mean$ID <- as.character(New_Pheno_mean$ID)
    dim(A_Mat)
    
    # Reorder df_means based on the rownames of K
    New_Pheno_mean <- New_Pheno_mean[match(rownames(A_Mat), New_Pheno_mean$ID), ]
    
    y<-as.matrix(New_Pheno_mean[,2:3])
    
    #### Model ####
    new_model<-mkr(y,A_Mat)
    
    # Save EBVs
    EBVs<-as.data.frame(new_model$hat)
    EBVs$ID <- New_Pheno_mean$ID
    
    EBVs <- EBVs[EBVs$ID %in% G2_pop@id,]
    new_Pheno <- New_Pheno_mean[New_Pheno_mean$ID %in% G2_pop@id,]
    
    G2_pop@ebv <- as.matrix(EBVs[,1:2])
    G2_pop@pheno <- as.matrix(new_Pheno[,2:3])
    
    ##### Select parents ######
    Parents <- selectWithinFam(G2_pop,
                               nInd = 1,
                               use = "ebv",
                               selectTop = TRUE)
    
    # Heritability
    h2<-var(bv(G2_pop))/varP(G2_pop)
    h2_hat<-new_model$h2
    
    True_h2_1 = h2[1,1]
    True_h2_2 = h2[2,2]
    Est_h2_1 = h2_hat[1]
    Est_h2_2 = h2_hat[2]
    
    # Additive genetic value   
    MeanG <- meanG(G2_pop)  
    VarA <- varA(G2_pop)
    Cor_A = cov2cor(VarA)[2, 1]
    
    # Dominance genetic value
    VarD <- varD(G2_pop)
    
    # Phenotypes
    Phen.Mean <- meanP(G2_pop)  
    VarP <- varP(G2_pop)
    Cor_P = cov2cor(VarP)[2, 1]
    
    # Prediction accuracy 
    Acc <- cor(bv(G2_pop),G2_pop@ebv) 
    
    # Within-family accuracy
    WF_df <- data.frame(ID = G2_pop@id,
                        Mum = G2_pop@mother,
                        Dad = G2_pop@father,
                        TBV1 = Bv[,1],
                        TBV2 = Bv[,2],
                        EBV1 = G2_pop@ebv[,1],
                        EBV2 = G2_pop@ebv[,2])
    WF_df$Family <- as.integer(interaction(WF_df$Mum, WF_df$Dad, drop = TRUE))
    
    summary_df <- WF_df %>%
      group_by(Family) %>%
      summarize(
        famSize = n(),
        accuracy1 = if(n() > 1) cor(TBV1, EBV1) else NA_real_,
        accuracy2 = if(n() > 1) cor(TBV2, EBV2) else NA_real_,
      ) %>%
      ungroup()
    
    
    WF_Acc1<-mean(summary_df$accuracy1, na.rm = T)
    WF_Acc2<-mean(summary_df$accuracy2, na.rm = T)
    
    T1_plot<-ggplot(summary_df, aes(x = Family, y = accuracy1, size = famSize)) +
      geom_point(alpha = 0.4, color = "red") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "G2",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    T2_plot<-ggplot(summary_df, aes(x = Family, y = accuracy2, size = famSize)) +
      geom_point(alpha = 0.4, color = "darkred") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    grid.arrange(T1_plot,T2_plot)
    
    # Prediction Bias
    Bias_1 <- mean(G2_pop@ebv[,1] - Bv[1])
    Bias_2 <- mean(G2_pop@ebv[,2] - Bv[2])
    
    # Estimate disversion 
    bv1 <- Bv[,1]
    bv2 <- Bv[,2]
    ebv1 <- G2_pop@ebv[,1]
    ebv2 <- G2_pop@ebv[,2]
    
    # Fit regression model
    bias_model1 <- lm(bv1 ~ ebv1)
    bias_model2 <- lm(bv2 ~ ebv2)
    
    # Optional: view slope too
    slope1 <- coef(bias_model1)[2]
    slope2 <- coef(bias_model2)[2]
    
    # Selection Intensity
    #Estimated Si
    mean_EBV_p_1 <- mean(Parents@ebv[,1])
    mean_EBV_p_2 <- mean(Parents@ebv[,2])
    
    mean_EBV_1 <- mean(G2_pop@ebv[,1])
    mean_EBV_2 <- mean(G2_pop@ebv[,2])
    
    sd_EBV_1 <- sd(G2_pop@ebv[,1])
    sd_EBV_2 <- sd(G2_pop@ebv[,2])
    
    EBV_Intensity_1 = ((mean_EBV_p_1 - mean_EBV_1) / sd_EBV_1)
    EBV_Intensity_2 = ((mean_EBV_p_2 - mean_EBV_2) / sd_EBV_2)
    
    #True SI
    bv_P<-bv(Parents)
    
    mean_TBV_p_1 <- mean(bv_P[,1])
    mean_TBV_p_2 <- mean(bv_P[,2])
    
    mean_TBV_1 <- mean(Bv[,1])
    mean_TBV_2 <- mean(Bv[,2])
    
    sd_TBV_1 <- sd(Bv[,1])
    sd_TBV_2 <- sd(Bv[,2])
    
    TBV_Intensity_1 = ((mean_TBV_p_1 - mean_TBV_1) / sd_TBV_1)
    TBV_Intensity_2 = ((mean_TBV_p_2 - mean_TBV_2) / sd_TBV_2)
    
    new_trait <- data.frame(Trait1_h2_est = Est_h2_1,
                            Trait1_h2 = True_h2_1,
                            Trait1_MeanG = MeanG[1],
                            Trait1_VarA = VarA[1,1],
                            Trait1_VarD = VarD[1,1],
                            Trait1_MeanP = Phen.Mean[1],
                            Trait1_VarP = VarP[1,1],
                            Trait1_Acc = Acc[1,1],
                            Trait1_WF_Acc = WF_Acc1,
                            Trait1_Bias = Bias_1,
                            Trait1_Dispersion = slope1,
                            Trait1_SI_est = EBV_Intensity_1,
                            Trait1_SI = TBV_Intensity_1,
                            
                            Trait2_h2_est = Est_h2_2,
                            Trait2_h2 = True_h2_2,
                            Trait2_MeanG = MeanG[2],
                            Trait2_VarA = VarA[2,2],
                            Trait2_VarD = VarD[2,2],
                            Trait2_MeanP = Phen.Mean[2],
                            Trait2_VarP = VarP[2,2],
                            Trait2_Acc = Acc[2,2],
                            Trait2_WF_Acc = WF_Acc2,
                            Trait2_Bias = Bias_2,
                            Trait2_Dispersion = slope2,
                            Trait2_SI_est = EBV_Intensity_2,
                            Trait2_SI = TBV_Intensity_2,
                            
                            CovA = VarA[1,2],
                            CorA = Cor_A,
                            CovP = VarP[1,2],
                            CorP = Cor_P,
                            
                            Gen = "G2")
    
    
    
    GLOBAL_trait <- rbind(GLOBAL_trait, new_trait)
    
    # Diversity metrics
    indices<-which((rownames(A_Mat) %in% G2_pop@id))
    A_Mat <-A_Mat[indices,indices]
    n <- nrow(A_Mat)
    
    group_coancestry <- sum(A_Mat) / (n^2)
    
    status_number <- 1 / (2*group_coancestry)
    
    Ns_new <- data.frame(Status_number = status_number,
                         Coancestry = group_coancestry,
                         Gen = "G2")
    
    new_ID <- data.frame(Inb = (diag(A_Mat)-1),
                        Trait1 = G2_pop@pheno[,1],
                        Trait2 = G2_pop@pheno[,2],
                        Gen = "G2")
    
    GLOBAL_ID <- rbind(GLOBAL_ID, new_ID)
    
    Inbreeding<-mean(diag(A_Mat))-1
    
    Inb<-(diag(A_Mat)-1)
    
    Inb_df <- data.frame(I = Inb,
                         T1 = G2_pop@pheno[,1],
                         T2 = G2_pop@pheno[,2])
    
    model <- lm(T1 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID1_est<-slope/mean(Inb_df$T1)
    
    model <- lm(T2 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID2_est<-slope/mean(Inb_df$T2)
    
    new_diversity <- data.frame(F.inb = Inbreeding,
                                Ns = status_number,
                                Co = group_coancestry,
                                ID_est1 = ID1_est,
                                ID_est2 = ID2_est,
                                Gen = "G2")
    
    GLOBAL_diversity <- rbind(GLOBAL_diversity, new_diversity)
    
####                                G3                      #####
  
    # Generate progeny
    # Each cross gets a random numner of progenies, with minimum 1
    G_pops <- vector("list", nCross)
    for (i in 1:nCross) {
      Prog_number <- max(round(rnorm(1, mean = nProg, sd = 10)), 1)
      G_pops[[i]] <- randCross(
        pop = Parents, 
        nCrosses = 1, 
        nProgeny = Prog_number, 
        ignoreSexes = TRUE
      )
    }
    
    # Combine all families into one pop-object
    G3_pop <- do.call(c, G_pops)
    
    # Save New breeding values
    Bv <- bv(G3_pop)
    
    # Update pedigree
    Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G3_pop@id,])
    Pedigree_All <- rbind(Pedigree_All, Pedigree_New)
    
    A_Mat <- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                       dadid = Pedigree_All$father,
                       momid = Pedigree_All$mother)
    
    
    ####################### Calculate MAF Distributoin ###########################
    G_new<-pullSnpGeno(G3_pop, snpChip = "Inb")
    maf_values_new <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                                 Gen = "G3")
    
    maf_values_df <- rbind(maf_values_df,maf_values_new) 
    
    ################ Assign phenotypic values ####################################
    New_Pheno <- data.frame()
    
    # Should I use the same initial error variance here?
    
    
    for (i in 1:10) {
      pheno <- setPheno(G3_pop,
                        varE = initVarE,
                        corE = matrix(c(1, CorP,
                                        CorP, 1), nrow = 2),
                        rep = 1,
                        simParam = NULL,
                        onlyPheno = T)
      New_Pheno <- rbind(
        New_Pheno,
        data.frame(ID = G3_pop@id, Rep = i, Pheno = pheno)
      )
    }
    
    New_Pheno$Rep<-NULL
    
    GLOBAL_Phenotypes<-rbind(GLOBAL_Phenotypes,New_Pheno)
    
    New_Pheno_mean<- GLOBAL_Phenotypes %>%
      group_by(ID) %>%
      summarise(
        Trait1 = mean(Pheno.Trait1, na.rm = TRUE),
        Trait2 = mean(Pheno.Trait2, na.rm = TRUE)
      )
    
    New_Pheno_mean$ID <- as.character(New_Pheno_mean$ID)
    dim(A_Mat)
    
    # Reorder df_means based on the rownames of K
    New_Pheno_mean <- New_Pheno_mean[match(rownames(A_Mat), New_Pheno_mean$ID), ]
    
    y<-as.matrix(New_Pheno_mean[,2:3])
    
    #### Model ####
    new_model<-mkr(y,A_Mat)
    
    # Save EBVs
    EBVs<-as.data.frame(new_model$hat)
    EBVs$ID <- New_Pheno_mean$ID
    
    EBVs <- EBVs[EBVs$ID %in% G3_pop@id,]
    new_Pheno <- New_Pheno_mean[New_Pheno_mean$ID %in% G3_pop@id,]
    
    G3_pop@ebv <- as.matrix(EBVs[,1:2])
    G3_pop@pheno <- as.matrix(new_Pheno[,2:3])
    
    ##### Select parents ######
    Parents <- selectWithinFam(G3_pop,
                               nInd = 1,
                               use = "ebv",
                               selectTop = TRUE)
    
    # Heritability
    h2<-var(bv(G3_pop))/varP(G3_pop)
    h2_hat<-new_model$h2
    
    True_h2_1 = h2[1,1]
    True_h2_2 = h2[2,2]
    Est_h2_1 = h2_hat[1]
    Est_h2_2 = h2_hat[2]
    
    # Additive genetic value   
    MeanG <- meanG(G3_pop)  
    VarA <- varA(G3_pop)
    Cor_A = cov2cor(VarA)[2, 1]
    
    # Dominance genetic value
    VarD <- varD(G3_pop)
    
    # Phenotypes
    Phen.Mean <- meanP(G3_pop)  
    VarP <- varP(G3_pop)
    Cor_P = cov2cor(VarP)[2, 1]
    
    # Prediction accuracy 
    Acc <- cor(bv(G3_pop),G3_pop@ebv) 
    
    # Within-family accuracy
    WF_df <- data.frame(ID = G3_pop@id,
                        Mum = G3_pop@mother,
                        Dad = G3_pop@father,
                        TBV1 = Bv[,1],
                        TBV2 = Bv[,2],
                        EBV1 = G3_pop@ebv[,1],
                        EBV2 = G3_pop@ebv[,2])
    WF_df$Family <- as.integer(interaction(WF_df$Mum, WF_df$Dad, drop = TRUE))
    
    summary_df <- WF_df %>%
      group_by(Family) %>%
      summarize(
        famSize = n(),
        accuracy1 = if(n() > 1) cor(TBV1, EBV1) else NA_real_,
        accuracy2 = if(n() > 1) cor(TBV2, EBV2) else NA_real_,
      ) %>%
      ungroup()
    
    
    WF_Acc1<-mean(summary_df$accuracy1, na.rm = T)
    WF_Acc2<-mean(summary_df$accuracy2, na.rm = T)
    
    T1_plot<-ggplot(summary_df, aes(x = Family, y = accuracy1, size = famSize)) +
      geom_point(alpha = 0.4, color = "red") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "G3",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    T2_plot<-ggplot(summary_df, aes(x = Family, y = accuracy2, size = famSize)) +
      geom_point(alpha = 0.4, color = "darkred") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    grid.arrange(T1_plot,T2_plot)
    
    # Prediction Bias
    Bv <- bv(G3_pop)
    Bias_1 <- mean(G3_pop@ebv[,1] - Bv[1])
    Bias_2 <- mean(G3_pop@ebv[,2] - Bv[2])
    
    # Estimate disversion 
    bv1 <- Bv[,1]
    bv2 <- Bv[,2]
    ebv1 <- G3_pop@ebv[,1]
    ebv2 <- G3_pop@ebv[,2]
    
    # Fit regression model
    bias_model1 <- lm(bv1 ~ ebv1)
    bias_model2 <- lm(bv2 ~ ebv2)
    
    # Optional: view slope too
    slope1 <- coef(bias_model1)[2]
    slope2 <- coef(bias_model2)[2]
    
    # Selection Intensity
    #Estimated Si
    mean_EBV_p_1 <- mean(Parents@ebv[,1])
    mean_EBV_p_2 <- mean(Parents@ebv[,2])
    
    mean_EBV_1 <- mean(G3_pop@ebv[,1])
    mean_EBV_2 <- mean(G3_pop@ebv[,2])
    
    sd_EBV_1 <- sd(G3_pop@ebv[,1])
    sd_EBV_2 <- sd(G3_pop@ebv[,2])
    
    EBV_Intensity_1 = ((mean_EBV_p_1 - mean_EBV_1) / sd_EBV_1)
    EBV_Intensity_2 = ((mean_EBV_p_2 - mean_EBV_2) / sd_EBV_2)
    
    #True SI
    bv_P<-bv(Parents)
    
    mean_TBV_p_1 <- mean(bv_P[,1])
    mean_TBV_p_2 <- mean(bv_P[,2])
    
    mean_TBV_1 <- mean(Bv[,1])
    mean_TBV_2 <- mean(Bv[,2])
    
    sd_TBV_1 <- sd(Bv[,1])
    sd_TBV_2 <- sd(Bv[,2])
    
    TBV_Intensity_1 = ((mean_TBV_p_1 - mean_TBV_1) / sd_TBV_1)
    TBV_Intensity_2 = ((mean_TBV_p_2 - mean_TBV_2) / sd_TBV_2)
    
    new_trait <- data.frame(Trait1_h2_est = Est_h2_1,
                            Trait1_h2 = True_h2_1,
                            Trait1_MeanG = MeanG[1],
                            Trait1_VarA = VarA[1,1],
                            Trait1_VarD = VarD[1,1],
                            Trait1_MeanP = Phen.Mean[1],
                            Trait1_VarP = VarP[1,1],
                            Trait1_Acc = Acc[1,1],
                            Trait1_WF_Acc = WF_Acc1,
                            Trait1_Bias = Bias_1,
                            Trait1_Dispersion = slope1,
                            Trait1_SI_est = EBV_Intensity_1,
                            Trait1_SI = TBV_Intensity_1,
                            
                            Trait2_h2_est = Est_h2_2,
                            Trait2_h2 = True_h2_2,
                            Trait2_MeanG = MeanG[2],
                            Trait2_VarA = VarA[2,2],
                            Trait2_VarD = VarD[2,2],
                            Trait2_MeanP = Phen.Mean[2],
                            Trait2_VarP = VarP[2,2],
                            Trait2_Acc = Acc[2,2],
                            Trait2_WF_Acc = WF_Acc2,
                            Trait2_Bias = Bias_2,
                            Trait2_Dispersion = slope2,
                            Trait2_SI_est = EBV_Intensity_2,
                            Trait2_SI = TBV_Intensity_2,
                            
                            CovA = VarA[1,2],
                            CorA = Cor_A,
                            CovP = VarP[1,2],
                            CorP = Cor_P,
                            
                            Gen = "G3")
    
    
    
    GLOBAL_trait <- rbind(GLOBAL_trait, new_trait)
    
    # Diversity metrics
    indices<-which((rownames(A_Mat) %in% G3_pop@id))
    A_Mat <-A_Mat[indices,indices]
    n <- nrow(A_Mat)
    
    group_coancestry <- sum(A_Mat) / (n^2)
    
    status_number <- 1 / (2*group_coancestry)
    
    Ns_new <- data.frame(Status_number = status_number,
                         Coancestry = group_coancestry,
                         Gen = "G3")
    
    new_ID <- data.frame(Inb = (diag(A_Mat)-1),
                         Trait1 = G3_pop@pheno[,1],
                         Trait2 = G3_pop@pheno[,2],
                         Gen = "G3")
    
    GLOBAL_ID <- rbind(GLOBAL_ID, new_ID)
    
    
    Inbreeding<-mean(diag(A_Mat))-1
    
    Inb<-(diag(A_Mat)-1)
    
    Inb_df <- data.frame(I = Inb,
                         T1 = G3_pop@pheno[,1],
                         T2 = G3_pop@pheno[,2])
    
    model <- lm(T1 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID1_est<-slope/mean(Inb_df$T1)
    
    model <- lm(T2 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID2_est<-slope/mean(Inb_df$T2)
    
    new_diversity <- data.frame(F.inb = Inbreeding,
                                Ns = status_number,
                                Co = group_coancestry,
                                ID_est1 = ID1_est,
                                ID_est2 = ID2_est,
                                Gen = "G3")
    
    GLOBAL_diversity <- rbind(GLOBAL_diversity, new_diversity)
    
### G4 ####
  
    # Generate progeny
    # Each cross gets a random numner of progenies, with minimum 1
    G_pops <- vector("list", nCross)
    for (i in 1:nCross) {
      Prog_number <- max(round(rnorm(1, mean = nProg, sd = 10)), 1)
      G_pops[[i]] <- randCross(
        pop = Parents, 
        nCrosses = 1, 
        nProgeny = Prog_number, 
        ignoreSexes = TRUE
      )
    }
    
    # Combine all families into one pop-object
    G4_pop <- do.call(c, G_pops)
    
    # Save New breeding values
    Bv <- bv(G4_pop)
    
    # Update pedigree
    Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G4_pop@id,])
    Pedigree_All <- rbind(Pedigree_All, Pedigree_New)
    
    A_Mat <- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                       dadid = Pedigree_All$father,
                       momid = Pedigree_All$mother)
    
    
    ####################### Calculate MAF Distributoin ###########################
    G_new<-pullSnpGeno(G4_pop, snpChip = "Inb")
    maf_values_new <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                                 Gen = "G4")
    
    maf_values_df <- rbind(maf_values_df,maf_values_new) 
    
    ################ Assign phenotypic values ####################################
    New_Pheno <- data.frame()
    
    # Should I use the same initial error variance here?
    
    
    for (i in 1:10) {
      pheno <- setPheno(G4_pop,
                        varE = initVarE,
                        corE = matrix(c(1, CorP,
                                        CorP, 1), nrow = 2),
                        rep = 1,
                        simParam = NULL,
                        onlyPheno = T)
      New_Pheno <- rbind(
        New_Pheno,
        data.frame(ID = G4_pop@id, Rep = i, Pheno = pheno)
      )
    }
    
    New_Pheno$Rep<-NULL
    
    GLOBAL_Phenotypes<-rbind(GLOBAL_Phenotypes,New_Pheno)
    
    New_Pheno_mean<- GLOBAL_Phenotypes %>%
      group_by(ID) %>%
      summarise(
        Trait1 = mean(Pheno.Trait1, na.rm = TRUE),
        Trait2 = mean(Pheno.Trait2, na.rm = TRUE)
      )
    
    New_Pheno_mean$ID <- as.character(New_Pheno_mean$ID)
    dim(A_Mat)
    
    # Reorder df_means based on the rownames of K
    New_Pheno_mean <- New_Pheno_mean[match(rownames(A_Mat), New_Pheno_mean$ID), ]
    
    y<-as.matrix(New_Pheno_mean[,2:3])
    
    #### Model ####
    new_model<-mkr(y,A_Mat)
    
    # Save EBVs
    EBVs<-as.data.frame(new_model$hat)
    EBVs$ID <- New_Pheno_mean$ID
    
    EBVs <- EBVs[EBVs$ID %in% G4_pop@id,]
    new_Pheno <- New_Pheno_mean[New_Pheno_mean$ID %in% G4_pop@id,]
    
    G4_pop@ebv <- as.matrix(EBVs[,1:2])
    G4_pop@pheno <- as.matrix(new_Pheno[,2:3])
    
    ##### Select parents ######
    Parents <- selectWithinFam(G4_pop,
                               nInd = 1,
                               use = "ebv",
                               selectTop = TRUE)
    
    # Heritability
    h2<-var(bv(G4_pop))/varP(G4_pop)
    h2_hat<-new_model$h2
    
    True_h2_1 = h2[1,1]
    True_h2_2 = h2[2,2]
    Est_h2_1 = h2_hat[1]
    Est_h2_2 = h2_hat[2]
    
    # Additive genetic value   
    MeanG <- meanG(G4_pop)  
    VarA <- varA(G4_pop)
    Cor_A = cov2cor(VarA)[2, 1]
    
    # Dominance genetic value
    VarD <- varD(G4_pop)
    
    # Phenotypes
    Phen.Mean <- meanP(G4_pop)  
    VarP <- varP(G4_pop)
    Cor_P = cov2cor(VarP)[2, 1]
    
    # Prediction accuracy 
    Acc <- cor(bv(G4_pop),G4_pop@ebv) 
    
    # Within-family accuracy
    WF_df <- data.frame(ID = G4_pop@id,
                        Mum = G4_pop@mother,
                        Dad = G4_pop@father,
                        TBV1 = Bv[,1],
                        TBV2 = Bv[,2],
                        EBV1 = G4_pop@ebv[,1],
                        EBV2 = G4_pop@ebv[,2])
    WF_df$Family <- as.integer(interaction(WF_df$Mum, WF_df$Dad, drop = TRUE))
    
    summary_df <- WF_df %>%
      group_by(Family) %>%
      summarize(
        famSize = n(),
        accuracy1 = if(n() > 1) cor(TBV1, EBV1) else NA_real_,
        accuracy2 = if(n() > 1) cor(TBV2, EBV2) else NA_real_,
      ) %>%
      ungroup()
    
    
    WF_Acc1<-mean(summary_df$accuracy1, na.rm = T)
    WF_Acc2<-mean(summary_df$accuracy2, na.rm = T)
    
    T1_plot<-ggplot(summary_df, aes(x = Family, y = accuracy1, size = famSize)) +
      geom_point(alpha = 0.4, color = "red") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "G4",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    T2_plot<-ggplot(summary_df, aes(x = Family, y = accuracy2, size = famSize)) +
      geom_point(alpha = 0.4, color = "darkred") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    grid.arrange(T1_plot,T2_plot)
    
    # Prediction Bias
    Bv <- bv(G4_pop)
    Bias_1 <- mean(G4_pop@ebv[,1] - Bv[1])
    Bias_2 <- mean(G4_pop@ebv[,2] - Bv[2])
    
    # Estimate disversion 
    bv1 <- Bv[,1]
    bv2 <- Bv[,2]
    ebv1 <- G4_pop@ebv[,1]
    ebv2 <- G4_pop@ebv[,2]
    
    # Fit regression model
    bias_model1 <- lm(bv1 ~ ebv1)
    bias_model2 <- lm(bv2 ~ ebv2)
    
    # Optional: view slope too
    slope1 <- coef(bias_model1)[2]
    slope2 <- coef(bias_model2)[2]
    
    # Selection Intensity
    #Estimated Si
    mean_EBV_p_1 <- mean(Parents@ebv[,1])
    mean_EBV_p_2 <- mean(Parents@ebv[,2])
    
    mean_EBV_1 <- mean(G4_pop@ebv[,1])
    mean_EBV_2 <- mean(G4_pop@ebv[,2])
    
    sd_EBV_1 <- sd(G4_pop@ebv[,1])
    sd_EBV_2 <- sd(G4_pop@ebv[,2])
    
    EBV_Intensity_1 = ((mean_EBV_p_1 - mean_EBV_1) / sd_EBV_1)
    EBV_Intensity_2 = ((mean_EBV_p_2 - mean_EBV_2) / sd_EBV_2)
    
    #True SI
    bv_P<-bv(Parents)
    
    mean_TBV_p_1 <- mean(bv_P[,1])
    mean_TBV_p_2 <- mean(bv_P[,2])
    
    mean_TBV_1 <- mean(Bv[,1])
    mean_TBV_2 <- mean(Bv[,2])
    
    sd_TBV_1 <- sd(Bv[,1])
    sd_TBV_2 <- sd(Bv[,2])
    
    TBV_Intensity_1 = ((mean_TBV_p_1 - mean_TBV_1) / sd_TBV_1)
    TBV_Intensity_2 = ((mean_TBV_p_2 - mean_TBV_2) / sd_TBV_2)
    
    new_trait <- data.frame(Trait1_h2_est = Est_h2_1,
                            Trait1_h2 = True_h2_1,
                            Trait1_MeanG = MeanG[1],
                            Trait1_VarA = VarA[1,1],
                            Trait1_VarD = VarD[1,1],
                            Trait1_MeanP = Phen.Mean[1],
                            Trait1_VarP = VarP[1,1],
                            Trait1_Acc = Acc[1,1],
                            Trait1_WF_Acc = WF_Acc1,
                            Trait1_Bias = Bias_1,
                            Trait1_Dispersion = slope1,
                            Trait1_SI_est = EBV_Intensity_1,
                            Trait1_SI = TBV_Intensity_1,
                            
                            Trait2_h2_est = Est_h2_2,
                            Trait2_h2 = True_h2_2,
                            Trait2_MeanG = MeanG[2],
                            Trait2_VarA = VarA[2,2],
                            Trait2_VarD = VarD[2,2],
                            Trait2_MeanP = Phen.Mean[2],
                            Trait2_VarP = VarP[2,2],
                            Trait2_Acc = Acc[2,2],
                            Trait2_WF_Acc = WF_Acc2,
                            Trait2_Bias = Bias_2,
                            Trait2_Dispersion = slope2,
                            Trait2_SI_est = EBV_Intensity_2,
                            Trait2_SI = TBV_Intensity_2,
                            
                            CovA = VarA[1,2],
                            CorA = Cor_A,
                            CovP = VarP[1,2],
                            CorP = Cor_P,
                            
                            Gen = "G4")
    
    
    
    GLOBAL_trait <- rbind(GLOBAL_trait, new_trait)
    
    # Diversity metrics
    indices<-which((rownames(A_Mat) %in% G4_pop@id))
    A_Mat <-A_Mat[indices,indices]
    n <- nrow(A_Mat)
    
    group_coancestry <- sum(A_Mat) / (n^2)
    
    status_number <- 1 / (2*group_coancestry)
    
    Ns_new <- data.frame(Status_number = status_number,
                         Coancestry = group_coancestry,
                         Gen = "G4")
    
    new_ID <- data.frame(Inb = (diag(A_Mat)-1),
                         Trait1 = G4_pop@pheno[,1],
                         Trait2 = G4_pop@pheno[,2],
                         Gen = "G4")
    
    GLOBAL_ID <- rbind(GLOBAL_ID, new_ID)
    
    Inbreeding<-mean(diag(A_Mat))-1
    
    Inb<-(diag(A_Mat)-1)
    
    Inb_df <- data.frame(I = Inb,
                         T1 = G4_pop@pheno[,1],
                         T2 = G4_pop@pheno[,2])
    
    model <- lm(T1 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID1_est<-slope/mean(Inb_df$T1)
    
    model <- lm(T2 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID2_est<-slope/mean(Inb_df$T2)
    
    new_diversity <- data.frame(F.inb = Inbreeding,
                                Ns = status_number,
                                Co = group_coancestry,
                                ID_est1 = ID1_est,
                                ID_est2 = ID2_est,
                                Gen = "G4")
    
    GLOBAL_diversity <- rbind(GLOBAL_diversity, new_diversity)
    
### G5 ####
    
    # Generate progeny
    # Each cross gets a random numner of progenies, with minimum 1
    G_pops <- vector("list", nCross)
    for (i in 1:nCross) {
      Prog_number <- max(round(rnorm(1, mean = nProg, sd = 10)), 1)
      G_pops[[i]] <- randCross(
        pop = Parents, 
        nCrosses = 1, 
        nProgeny = Prog_number, 
        ignoreSexes = TRUE
      )
    }
    
    # Combine all families into one pop-object
    G5_pop <- do.call(c, G_pops)
    
    # Save New breeding values
    Bv <- bv(G5_pop)
    
    # Update pedigree
    Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G5_pop@id,])
    Pedigree_All <- rbind(Pedigree_All, Pedigree_New)
    
    A_Mat <- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                       dadid = Pedigree_All$father,
                       momid = Pedigree_All$mother)
    
    
    ####################### Calculate MAF Distributoin ###########################
    G_new<-pullSnpGeno(G5_pop, snpChip = "Inb")
    maf_values_new <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                                 Gen = "G5")
    
    maf_values_df <- rbind(maf_values_df,maf_values_new) 
    
    ################ Assign phenotypic values ####################################
    New_Pheno <- data.frame()
    
    # Should I use the same initial error variance here?
    
    
    for (i in 1:10) {
      pheno <- setPheno(G5_pop,
                        varE = initVarE,
                        corE = matrix(c(1, CorP,
                                        CorP, 1), nrow = 2),
                        rep = 1,
                        simParam = NULL,
                        onlyPheno = T)
      New_Pheno <- rbind(
        New_Pheno,
        data.frame(ID = G5_pop@id, Rep = i, Pheno = pheno)
      )
    }
    
    New_Pheno$Rep<-NULL
    
    GLOBAL_Phenotypes<-rbind(GLOBAL_Phenotypes,New_Pheno)
    
    New_Pheno_mean<- GLOBAL_Phenotypes %>%
      group_by(ID) %>%
      summarise(
        Trait1 = mean(Pheno.Trait1, na.rm = TRUE),
        Trait2 = mean(Pheno.Trait2, na.rm = TRUE)
      )
    
    New_Pheno_mean$ID <- as.character(New_Pheno_mean$ID)
    dim(A_Mat)
    
    # Reorder df_means based on the rownames of K
    New_Pheno_mean <- New_Pheno_mean[match(rownames(A_Mat), New_Pheno_mean$ID), ]
    
    y<-as.matrix(New_Pheno_mean[,2:3])
    
    #### Model ####
    new_model<-mkr(y,A_Mat)
    
    # Save EBVs
    EBVs<-as.data.frame(new_model$hat)
    EBVs$ID <- New_Pheno_mean$ID
    
    EBVs <- EBVs[EBVs$ID %in% G5_pop@id,]
    new_Pheno <- New_Pheno_mean[New_Pheno_mean$ID %in% G5_pop@id,]
    
    G5_pop@ebv <- as.matrix(EBVs[,1:2])
    G5_pop@pheno <- as.matrix(new_Pheno[,2:3])
    
    ##### Select parents ######
    Parents <- selectWithinFam(G5_pop,
                               nInd = 1,
                               use = "ebv",
                               selectTop = TRUE)
    
    # Heritability
    h2<-var(bv(G5_pop))/varP(G5_pop)
    h2_hat<-new_model$h2
    
    True_h2_1 = h2[1,1]
    True_h2_2 = h2[2,2]
    Est_h2_1 = h2_hat[1]
    Est_h2_2 = h2_hat[2]
    
    # Additive genetic value   
    MeanG <- meanG(G5_pop)  
    VarA <- varA(G5_pop)
    Cor_A = cov2cor(VarA)[2, 1]
    
    # Dominance genetic value
    VarD <- varD(G5_pop)
    
    # Phenotypes
    Phen.Mean <- meanP(G5_pop)  
    VarP <- varP(G5_pop)
    Cor_P = cov2cor(VarP)[2, 1]
    
    # Prediction accuracy 
    Acc <- cor(bv(G5_pop),G5_pop@ebv) 
    
    # Within-family accuracy
    WF_df <- data.frame(ID = G5_pop@id,
                        Mum = G5_pop@mother,
                        Dad = G5_pop@father,
                        TBV1 = Bv[,1],
                        TBV2 = Bv[,2],
                        EBV1 = G5_pop@ebv[,1],
                        EBV2 = G5_pop@ebv[,2])
    WF_df$Family <- as.integer(interaction(WF_df$Mum, WF_df$Dad, drop = TRUE))
    
    summary_df <- WF_df %>%
      group_by(Family) %>%
      summarize(
        famSize = n(),
        accuracy1 = if(n() > 1) cor(TBV1, EBV1) else NA_real_,
        accuracy2 = if(n() > 1) cor(TBV2, EBV2) else NA_real_,
      ) %>%
      ungroup()
    
    
    WF_Acc1<-mean(summary_df$accuracy1, na.rm = T)
    WF_Acc2<-mean(summary_df$accuracy2, na.rm = T)
    
    T1_plot<-ggplot(summary_df, aes(x = Family, y = accuracy1, size = famSize)) +
      geom_point(alpha = 0.4, color = "red") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "G5",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    T2_plot<-ggplot(summary_df, aes(x = Family, y = accuracy2, size = famSize)) +
      geom_point(alpha = 0.4, color = "darkred") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    grid.arrange(T1_plot,T2_plot)
    
    # Prediction Bias
    Bv <- bv(G5_pop)
    Bias_1 <- mean(G5_pop@ebv[,1] - Bv[1])
    Bias_2 <- mean(G5_pop@ebv[,2] - Bv[2])
    
    # Estimate disversion 
    bv1 <- Bv[,1]
    bv2 <- Bv[,2]
    ebv1 <- G5_pop@ebv[,1]
    ebv2 <- G5_pop@ebv[,2]
    
    # Fit regression model
    bias_model1 <- lm(bv1 ~ ebv1)
    bias_model2 <- lm(bv2 ~ ebv2)
    
    # Optional: view slope too
    slope1 <- coef(bias_model1)[2]
    slope2 <- coef(bias_model2)[2]
    
    # Selection Intensity
    #Estimated Si
    mean_EBV_p_1 <- mean(Parents@ebv[,1])
    mean_EBV_p_2 <- mean(Parents@ebv[,2])
    
    mean_EBV_1 <- mean(G5_pop@ebv[,1])
    mean_EBV_2 <- mean(G5_pop@ebv[,2])
    
    sd_EBV_1 <- sd(G5_pop@ebv[,1])
    sd_EBV_2 <- sd(G5_pop@ebv[,2])
    
    EBV_Intensity_1 = ((mean_EBV_p_1 - mean_EBV_1) / sd_EBV_1)
    EBV_Intensity_2 = ((mean_EBV_p_2 - mean_EBV_2) / sd_EBV_2)
    
    #True SI
    bv_P<-bv(Parents)
    
    mean_TBV_p_1 <- mean(bv_P[,1])
    mean_TBV_p_2 <- mean(bv_P[,2])
    
    mean_TBV_1 <- mean(Bv[,1])
    mean_TBV_2 <- mean(Bv[,2])
    
    sd_TBV_1 <- sd(Bv[,1])
    sd_TBV_2 <- sd(Bv[,2])
    
    TBV_Intensity_1 = ((mean_TBV_p_1 - mean_TBV_1) / sd_TBV_1)
    TBV_Intensity_2 = ((mean_TBV_p_2 - mean_TBV_2) / sd_TBV_2)
    
    new_trait <- data.frame(Trait1_h2_est = Est_h2_1,
                            Trait1_h2 = True_h2_1,
                            Trait1_MeanG = MeanG[1],
                            Trait1_VarA = VarA[1,1],
                            Trait1_VarD = VarD[1,1],
                            Trait1_MeanP = Phen.Mean[1],
                            Trait1_VarP = VarP[1,1],
                            Trait1_Acc = Acc[1,1],
                            Trait1_WF_Acc = WF_Acc1,
                            Trait1_Bias = Bias_1,
                            Trait1_Dispersion = slope1,
                            Trait1_SI_est = EBV_Intensity_1,
                            Trait1_SI = TBV_Intensity_1,
                            
                            Trait2_h2_est = Est_h2_2,
                            Trait2_h2 = True_h2_2,
                            Trait2_MeanG = MeanG[2],
                            Trait2_VarA = VarA[2,2],
                            Trait2_VarD = VarD[2,2],
                            Trait2_MeanP = Phen.Mean[2],
                            Trait2_VarP = VarP[2,2],
                            Trait2_Acc = Acc[2,2],
                            Trait2_WF_Acc = WF_Acc2,
                            Trait2_Bias = Bias_2,
                            Trait2_Dispersion = slope2,
                            Trait2_SI_est = EBV_Intensity_2,
                            Trait2_SI = TBV_Intensity_2,
                            
                            CovA = VarA[1,2],
                            CorA = Cor_A,
                            CovP = VarP[1,2],
                            CorP = Cor_P,
                            
                            Gen = "G5")
    
    
    
    GLOBAL_trait <- rbind(GLOBAL_trait, new_trait)
    
    # Diversity metrics
    indices<-which((rownames(A_Mat) %in% G5_pop@id))
    A_Mat <-A_Mat[indices,indices]
    n <- nrow(A_Mat)
    
    group_coancestry <- sum(A_Mat) / (n^2)
    
    status_number <- 1 / (2*group_coancestry)
    
    Ns_new <- data.frame(Status_number = status_number,
                         Coancestry = group_coancestry,
                         Gen = "G5")
    
    new_ID <- data.frame(Inb = (diag(A_Mat)-1),
                         Trait1 = G5_pop@pheno[,1],
                         Trait2 = G5_pop@pheno[,2],
                         Gen = "G5")
    
    GLOBAL_ID <- rbind(GLOBAL_ID, new_ID)
    
    
    Inbreeding<-mean(diag(A_Mat))-1
    
    Inb<-(diag(A_Mat)-1)
    
    Inb_df <- data.frame(I = Inb,
                         T1 = G5_pop@pheno[,1],
                         T2 = G5_pop@pheno[,2])
    
    model <- lm(T1 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID1_est<-slope/mean(Inb_df$T1)
    
    model <- lm(T2 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID2_est<-slope/mean(Inb_df$T2)
    
    new_diversity <- data.frame(F.inb = Inbreeding,
                                Ns = status_number,
                                Co = group_coancestry,
                                ID_est1 = ID1_est,
                                ID_est2 = ID2_est,
                                Gen = "G5")
    
    GLOBAL_diversity <- rbind(GLOBAL_diversity, new_diversity)
    
#### G6 #####
    
    # Generate progeny
    # Each cross gets a random numner of progenies, with minimum 1
    G_pops <- vector("list", nCross)
    for (i in 1:nCross) {
      Prog_number <- max(round(rnorm(1, mean = nProg, sd = 10)), 1)
      G_pops[[i]] <- randCross(
        pop = Parents, 
        nCrosses = 1, 
        nProgeny = Prog_number, 
        ignoreSexes = TRUE
      )
    }
    
    # Combine all families into one pop-object
    G6_pop <- do.call(c, G_pops)
    
    # Save New breeding values
    Bv <- bv(G6_pop)
    
    # Update pedigree
    Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G6_pop@id,])
    Pedigree_All <- rbind(Pedigree_All, Pedigree_New)
    
    Pedigree_5G <- c(G1_pop@id,G2_pop@id,G3_pop@id,G4_pop@id,G5_pop@id,
                     G6_pop@id)
    Pedigree_5G <- Pedigree_All[row.names(Pedigree_All) %in% Pedigree_5G,]
    
    A_Mat <- 2*kinship(id = as.numeric(rownames(Pedigree_5G)),
                       dadid = Pedigree_5G$father,
                       momid = Pedigree_5G$mother)
    
    # Separate A-matrix for calculating co-ancestry
    A_Mat2 <- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                        dadid = Pedigree_All$father,
                        momid = Pedigree_All$mother)
    
    ####################### Calculate MAF Distributoin ###########################
    G_new<-pullSnpGeno(G6_pop, snpChip = "Inb")
    maf_values_new <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                                 Gen = "G6")
    
    maf_values_df <- rbind(maf_values_df,maf_values_new) 
    
    ################ Assign phenotypic values ####################################
    New_Pheno <- data.frame()
    
    # Should I use the same initial error variance here?
    
    
    for (i in 1:10) {
      pheno <- setPheno(G6_pop,
                        varE = initVarE,
                        corE = matrix(c(1, CorP,
                                        CorP, 1), nrow = 2),
                        rep = 1,
                        simParam = NULL,
                        onlyPheno = T)
      New_Pheno <- rbind(
        New_Pheno,
        data.frame(ID = G6_pop@id, Rep = i, Pheno = pheno)
      )
    }
    
    New_Pheno$Rep<-NULL
    
    GLOBAL_Phenotypes<-rbind(GLOBAL_Phenotypes,New_Pheno)
    Phenotypes_5G <- GLOBAL_Phenotypes[GLOBAL_Phenotypes$ID %in% rownames(Pedigree_5G),]
    
    New_Pheno_mean<- Phenotypes_5G %>%
      group_by(ID) %>%
      summarise(
        Trait1 = mean(Pheno.Trait1, na.rm = TRUE),
        Trait2 = mean(Pheno.Trait2, na.rm = TRUE)
      )
    
    New_Pheno_mean$ID <- as.character(New_Pheno_mean$ID)
    
    # Reorder df_means based on the rownames of K
    New_Pheno_mean <- New_Pheno_mean[match(rownames(A_Mat), New_Pheno_mean$ID), ]
    
    y<-as.matrix(New_Pheno_mean[,2:3])
    
    #### Model ####
    new_model<-mkr(y,A_Mat)
    
    # Save EBVs
    EBVs<-as.data.frame(new_model$hat)
    EBVs$ID <- New_Pheno_mean$ID
    
    EBVs <- EBVs[EBVs$ID %in% G6_pop@id,]
    new_Pheno <- New_Pheno_mean[New_Pheno_mean$ID %in% G6_pop@id,]
    
    G6_pop@ebv <- as.matrix(EBVs[,1:2])
    G6_pop@pheno <- as.matrix(new_Pheno[,2:3])
    
    ##### Select parents ######
    Parents <- selectWithinFam(G6_pop,
                               nInd = 1,
                               use = "ebv",
                               selectTop = TRUE)
    
    # Heritability
    h2<-var(bv(G6_pop))/varP(G6_pop)
    h2_hat<-new_model$h2
    
    True_h2_1 = h2[1,1]
    True_h2_2 = h2[2,2]
    Est_h2_1 = h2_hat[1]
    Est_h2_2 = h2_hat[2]
    
    # Additive genetic value   
    MeanG <- meanG(G6_pop)  
    VarA <- varA(G6_pop)
    Cor_A = cov2cor(VarA)[2, 1]
    
    # Dominance genetic value
    VarD <- varD(G6_pop)
    
    # Phenotypes
    Phen.Mean <- meanP(G6_pop)  
    VarP <- varP(G6_pop)
    Cor_P = cov2cor(VarP)[2, 1]
    
    # Prediction accuracy 
    Acc <- cor(bv(G6_pop),G6_pop@ebv) 
    
    # Prediction Bias
    Bv <- bv(G6_pop)
    Bias_1 <- mean(G6_pop@ebv[,1] - Bv[1])
    Bias_2 <- mean(G6_pop@ebv[,2] - Bv[2])
    
    # Estimate disversion 
    bv1 <- Bv[,1]
    bv2 <- Bv[,2]
    ebv1 <- G6_pop@ebv[,1]
    ebv2 <- G6_pop@ebv[,2]
    
    # Fit regression model
    bias_model1 <- lm(bv1 ~ ebv1)
    bias_model2 <- lm(bv2 ~ ebv2)
    
    # Optional: view slope too
    slope1 <- coef(bias_model1)[2]
    slope2 <- coef(bias_model2)[2]
    
    # Selection Intensity
    #Estimated Si
    mean_EBV_p_1 <- mean(Parents@ebv[,1])
    mean_EBV_p_2 <- mean(Parents@ebv[,2])
    
    mean_EBV_1 <- mean(G6_pop@ebv[,1])
    mean_EBV_2 <- mean(G6_pop@ebv[,2])
    
    sd_EBV_1 <- sd(G6_pop@ebv[,1])
    sd_EBV_2 <- sd(G6_pop@ebv[,2])
    
    EBV_Intensity_1 = ((mean_EBV_p_1 - mean_EBV_1) / sd_EBV_1)
    EBV_Intensity_2 = ((mean_EBV_p_2 - mean_EBV_2) / sd_EBV_2)
    
    #True SI
    bv_P<-bv(Parents)
    
    mean_TBV_p_1 <- mean(bv_P[,1])
    mean_TBV_p_2 <- mean(bv_P[,2])
    
    mean_TBV_1 <- mean(Bv[,1])
    mean_TBV_2 <- mean(Bv[,2])
    
    sd_TBV_1 <- sd(Bv[,1])
    sd_TBV_2 <- sd(Bv[,2])
    
    TBV_Intensity_1 = ((mean_TBV_p_1 - mean_TBV_1) / sd_TBV_1)
    TBV_Intensity_2 = ((mean_TBV_p_2 - mean_TBV_2) / sd_TBV_2)
    
    new_trait <- data.frame(Trait1_h2_est = Est_h2_1,
                            Trait1_h2 = True_h2_1,
                            Trait1_MeanG = MeanG[1],
                            Trait1_VarA = VarA[1,1],
                            Trait1_VarD = VarD[1,1],
                            Trait1_MeanP = Phen.Mean[1],
                            Trait1_VarP = VarP[1,1],
                            Trait1_Acc = Acc[1,1],
                            Trait1_WF_Acc = WF_Acc1,
                            Trait1_Bias = Bias_1,
                            Trait1_Dispersion = slope1,
                            Trait1_SI_est = EBV_Intensity_1,
                            Trait1_SI = TBV_Intensity_1,
                            
                            Trait2_h2_est = Est_h2_2,
                            Trait2_h2 = True_h2_2,
                            Trait2_MeanG = MeanG[2],
                            Trait2_VarA = VarA[2,2],
                            Trait2_VarD = VarD[2,2],
                            Trait2_MeanP = Phen.Mean[2],
                            Trait2_VarP = VarP[2,2],
                            Trait2_Acc = Acc[2,2],
                            Trait2_WF_Acc = WF_Acc2,
                            Trait2_Bias = Bias_2,
                            Trait2_Dispersion = slope2,
                            Trait2_SI_est = EBV_Intensity_2,
                            Trait2_SI = TBV_Intensity_2,
                            
                            CovA = VarA[1,2],
                            CorA = Cor_A,
                            CovP = VarP[1,2],
                            CorP = Cor_P,
                            
                            Gen = "G6")
    
    
    
    GLOBAL_trait <- rbind(GLOBAL_trait, new_trait)
    
    # Diversity metrics
    indices<-which((rownames(A_Mat2) %in% G6_pop@id))
    A_Mat2 <-A_Mat2[indices,indices]
    n <- nrow(A_Mat2)
    
    group_coancestry <- sum(A_Mat2) / (n^2)
    
    status_number <- 1 / (2*group_coancestry)
    
    Ns_new <- data.frame(Status_number = status_number,
                         Coancestry = group_coancestry,
                         Gen = "G6")
    
    new_ID <- data.frame(Inb = (diag(A_Mat2)-1),
                         Trait1 = G6_pop@pheno[,1],
                         Trait2 = G6_pop@pheno[,2],
                         Gen = "G6")
    
    GLOBAL_ID <- rbind(GLOBAL_ID, new_ID)
    
    Inbreeding<-mean(diag(A_Mat))-1
    
    Inb<-(diag(A_Mat2)-1)
    
    Inb_df <- data.frame(I = Inb,
                         T1 = G6_pop@pheno[,1],
                         T2 = G6_pop@pheno[,2])
    
    model <- lm(T1 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID1_est<-slope/mean(Inb_df$T1)
    
    model <- lm(T2 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID2_est<-slope/mean(Inb_df$T2)
    
    new_diversity <- data.frame(F.inb = Inbreeding,
                                Ns = status_number,
                                Co = group_coancestry,
                                ID_est1 = ID1_est,
                                ID_est2 = ID2_est,
                                Gen = "G6")
    
    GLOBAL_diversity <- rbind(GLOBAL_diversity, new_diversity)
    
### G7 ####    
    
    # Generate progeny
    # Each cross gets a random numner of progenies, with minimum 1
    G_pops <- vector("list", nCross)
    for (i in 1:nCross) {
      Prog_number <- max(round(rnorm(1, mean = nProg, sd = 10)), 1)
      G_pops[[i]] <- randCross(
        pop = Parents, 
        nCrosses = 1, 
        nProgeny = Prog_number, 
        ignoreSexes = TRUE
      )
    }
    
    # Combine all families into one pop-object
    G7_pop <- do.call(c, G_pops)
    
    # Save New breeding values
    Bv <- bv(G7_pop)
    
    # Update pedigree
    Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G7_pop@id,])
    Pedigree_All <- rbind(Pedigree_All, Pedigree_New)
    
    Pedigree_5G <- c(G2_pop@id,G3_pop@id,G4_pop@id,G5_pop@id,G6_pop@id,
                     G7_pop@id)
    Pedigree_5G <- Pedigree_All[row.names(Pedigree_All) %in% Pedigree_5G,]
    
    A_Mat <- 2*kinship(id = as.numeric(rownames(Pedigree_5G)),
                       dadid = Pedigree_5G$father,
                       momid = Pedigree_5G$mother)
    
    # Separate A-matrix for calculating co-ancestry
    A_Mat2 <- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                        dadid = Pedigree_All$father,
                        momid = Pedigree_All$mother)
    
    ####################### Calculate MAF Distributoin ###########################
    G_new<-pullSnpGeno(G7_pop, snpChip = "Inb")
    maf_values_new <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                                 Gen = "G7")
    
    maf_values_df <- rbind(maf_values_df,maf_values_new) 
    
    ################ Assign phenotypic values ####################################
    New_Pheno <- data.frame()
    
    # Should I use the same initial error variance here?
    
    
    for (i in 1:10) {
      pheno <- setPheno(G7_pop,
                        varE = initVarE,
                        corE = matrix(c(1, CorP,
                                        CorP, 1), nrow = 2),
                        rep = 1,
                        simParam = NULL,
                        onlyPheno = T)
      New_Pheno <- rbind(
        New_Pheno,
        data.frame(ID = G7_pop@id, Rep = i, Pheno = pheno)
      )
    }
    
    New_Pheno$Rep<-NULL
    
    GLOBAL_Phenotypes<-rbind(GLOBAL_Phenotypes,New_Pheno)
    Phenotypes_5G <- GLOBAL_Phenotypes[GLOBAL_Phenotypes$ID %in% rownames(Pedigree_5G),]
    
    New_Pheno_mean<- Phenotypes_5G %>%
      group_by(ID) %>%
      summarise(
        Trait1 = mean(Pheno.Trait1, na.rm = TRUE),
        Trait2 = mean(Pheno.Trait2, na.rm = TRUE)
      )
    
    New_Pheno_mean$ID <- as.character(New_Pheno_mean$ID)
    
    # Reorder df_means based on the rownames of K
    New_Pheno_mean <- New_Pheno_mean[match(rownames(A_Mat), New_Pheno_mean$ID), ]
    
    y<-as.matrix(New_Pheno_mean[,2:3])
    
    #### Model ####
    new_model<-mkr(y,A_Mat)
    
    # Save EBVs
    EBVs<-as.data.frame(new_model$hat)
    EBVs$ID <- New_Pheno_mean$ID
    
    EBVs <- EBVs[EBVs$ID %in% G7_pop@id,]
    new_Pheno <- New_Pheno_mean[New_Pheno_mean$ID %in% G7_pop@id,]
    
    G7_pop@ebv <- as.matrix(EBVs[,1:2])
    G7_pop@pheno <- as.matrix(new_Pheno[,2:3])
    
    ##### Select parents ######
    Parents <- selectWithinFam(G7_pop,
                               nInd = 1,
                               use = "ebv",
                               selectTop = TRUE)
    
    # Heritability
    h2<-var(bv(G7_pop))/varP(G7_pop)
    h2_hat<-new_model$h2
    
    True_h2_1 = h2[1,1]
    True_h2_2 = h2[2,2]
    Est_h2_1 = h2_hat[1]
    Est_h2_2 = h2_hat[2]
    
    # Additive genetic value   
    MeanG <- meanG(G7_pop)  
    VarA <- varA(G7_pop)
    Cor_A = cov2cor(VarA)[2, 1]
    
    # Dominance genetic value
    VarD <- varD(G7_pop)
    
    # Phenotypes
    Phen.Mean <- meanP(G7_pop)  
    VarP <- varP(G7_pop)
    Cor_P = cov2cor(VarP)[2, 1]
    
    # Prediction accuracy 
    Acc <- cor(bv(G7_pop),G7_pop@ebv) 
    
    # Within-family accuracy
    WF_df <- data.frame(ID = G7_pop@id,
                        Mum = G7_pop@mother,
                        Dad = G7_pop@father,
                        TBV1 = Bv[,1],
                        TBV2 = Bv[,2],
                        EBV1 = G7_pop@ebv[,1],
                        EBV2 = G7_pop@ebv[,2])
    WF_df$Family <- as.integer(interaction(WF_df$Mum, WF_df$Dad, drop = TRUE))
    
    summary_df <- WF_df %>%
      group_by(Family) %>%
      summarize(
        famSize = n(),
        accuracy1 = if(n() > 1) cor(TBV1, EBV1) else NA_real_,
        accuracy2 = if(n() > 1) cor(TBV2, EBV2) else NA_real_,
      ) %>%
      ungroup()
    
    
    WF_Acc1<-mean(summary_df$accuracy1, na.rm = T)
    WF_Acc2<-mean(summary_df$accuracy2, na.rm = T)
    
    T1_plot<-ggplot(summary_df, aes(x = Family, y = accuracy1, size = famSize)) +
      geom_point(alpha = 0.4, color = "red") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "G7",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    T2_plot<-ggplot(summary_df, aes(x = Family, y = accuracy2, size = famSize)) +
      geom_point(alpha = 0.4, color = "darkred") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    grid.arrange(T1_plot,T2_plot)
    
    # Prediction Bias
    Bv <- bv(G7_pop)
    Bias_1 <- mean(G7_pop@ebv[,1] - Bv[1])
    Bias_2 <- mean(G7_pop@ebv[,2] - Bv[2])
    
    # Estimate disversion 
    bv1 <- Bv[,1]
    bv2 <- Bv[,2]
    ebv1 <- G7_pop@ebv[,1]
    ebv2 <- G7_pop@ebv[,2]
    
    # Fit regression model
    bias_model1 <- lm(bv1 ~ ebv1)
    bias_model2 <- lm(bv2 ~ ebv2)
    
    # Optional: view slope too
    slope1 <- coef(bias_model1)[2]
    slope2 <- coef(bias_model2)[2]
    
    # Selection Intensity
    #Estimated Si
    mean_EBV_p_1 <- mean(Parents@ebv[,1])
    mean_EBV_p_2 <- mean(Parents@ebv[,2])
    
    mean_EBV_1 <- mean(G7_pop@ebv[,1])
    mean_EBV_2 <- mean(G7_pop@ebv[,2])
    
    sd_EBV_1 <- sd(G7_pop@ebv[,1])
    sd_EBV_2 <- sd(G7_pop@ebv[,2])
    
    EBV_Intensity_1 = ((mean_EBV_p_1 - mean_EBV_1) / sd_EBV_1)
    EBV_Intensity_2 = ((mean_EBV_p_2 - mean_EBV_2) / sd_EBV_2)
    
    #True SI
    bv_P<-bv(Parents)
    
    mean_TBV_p_1 <- mean(bv_P[,1])
    mean_TBV_p_2 <- mean(bv_P[,2])
    
    mean_TBV_1 <- mean(Bv[,1])
    mean_TBV_2 <- mean(Bv[,2])
    
    sd_TBV_1 <- sd(Bv[,1])
    sd_TBV_2 <- sd(Bv[,2])
    
    TBV_Intensity_1 = ((mean_TBV_p_1 - mean_TBV_1) / sd_TBV_1)
    TBV_Intensity_2 = ((mean_TBV_p_2 - mean_TBV_2) / sd_TBV_2)
    
    new_trait <- data.frame(Trait1_h2_est = Est_h2_1,
                            Trait1_h2 = True_h2_1,
                            Trait1_MeanG = MeanG[1],
                            Trait1_VarA = VarA[1,1],
                            Trait1_VarD = VarD[1,1],
                            Trait1_MeanP = Phen.Mean[1],
                            Trait1_VarP = VarP[1,1],
                            Trait1_Acc = Acc[1,1],
                            Trait1_WF_Acc = WF_Acc1,
                            Trait1_Bias = Bias_1,
                            Trait1_Dispersion = slope1,
                            Trait1_SI_est = EBV_Intensity_1,
                            Trait1_SI = TBV_Intensity_1,
                            
                            Trait2_h2_est = Est_h2_2,
                            Trait2_h2 = True_h2_2,
                            Trait2_MeanG = MeanG[2],
                            Trait2_VarA = VarA[2,2],
                            Trait2_VarD = VarD[2,2],
                            Trait2_MeanP = Phen.Mean[2],
                            Trait2_VarP = VarP[2,2],
                            Trait2_Acc = Acc[2,2],
                            Trait2_WF_Acc = WF_Acc2,
                            Trait2_Bias = Bias_2,
                            Trait2_Dispersion = slope2,
                            Trait2_SI_est = EBV_Intensity_2,
                            Trait2_SI = TBV_Intensity_2,
                            
                            CovA = VarA[1,2],
                            CorA = Cor_A,
                            CovP = VarP[1,2],
                            CorP = Cor_P,
                            
                            Gen = "G7")
    
    
    
    GLOBAL_trait <- rbind(GLOBAL_trait, new_trait)
    
    # Diversity metrics
    indices<-which((rownames(A_Mat2) %in% G7_pop@id))
    A_Mat2 <-A_Mat2[indices,indices]
    n <- nrow(A_Mat2)
    
    group_coancestry <- sum(A_Mat2) / (n^2)
    
    status_number <- 1 / (2*group_coancestry)
    
    Ns_new <- data.frame(Status_number = status_number,
                         Coancestry = group_coancestry,
                         Gen = "G7")
    
    new_ID <- data.frame(Inb = (diag(A_Mat2)-1),
                         Trait1 = G7_pop@pheno[,1],
                         Trait2 = G7_pop@pheno[,2],
                         Gen = "G7")
    
    GLOBAL_ID <- rbind(GLOBAL_ID, new_ID)
    
    Inbreeding<-mean(diag(A_Mat))-1
    
    Inb<-(diag(A_Mat2)-1)
    
    Inb_df <- data.frame(I = Inb,
                         T1 = G7_pop@pheno[,1],
                         T2 = G7_pop@pheno[,2])
    
    model <- lm(T1 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID1_est<-slope/mean(Inb_df$T1)
    
    model <- lm(T2 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID2_est<-slope/mean(Inb_df$T2)
    
    new_diversity <- data.frame(F.inb = Inbreeding,
                                Ns = status_number,
                                Co = group_coancestry,
                                ID_est1 = ID1_est,
                                ID_est2 = ID2_est,
                                Gen = "G7")
    
    GLOBAL_diversity <- rbind(GLOBAL_diversity, new_diversity)
    
#### G8 ####
    
    # Generate progeny
    # Each cross gets a random numner of progenies, with minimum 1
    G_pops <- vector("list", nCross)
    for (i in 1:nCross) {
      Prog_number <- max(round(rnorm(1, mean = nProg, sd = 10)), 1)
      G_pops[[i]] <- randCross(
        pop = Parents, 
        nCrosses = 1, 
        nProgeny = Prog_number, 
        ignoreSexes = TRUE
      )
    }
    
    # Combine all families into one pop-object
    G8_pop <- do.call(c, G_pops)
    
    # Save New breeding values
    Bv <- bv(G8_pop)
    
    # Update pedigree
    Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G8_pop@id,])
    Pedigree_All <- rbind(Pedigree_All, Pedigree_New)
    
    Pedigree_5G <- c(G3_pop@id,G4_pop@id,G5_pop@id,G6_pop@id,G7_pop@id,
                     G8_pop@id)
    Pedigree_5G <- Pedigree_All[row.names(Pedigree_All) %in% Pedigree_5G,]
    
    A_Mat <- 2*kinship(id = as.numeric(rownames(Pedigree_5G)),
                       dadid = Pedigree_5G$father,
                       momid = Pedigree_5G$mother)
    
    # Separate A-matrix for calculating co-ancestry
    A_Mat2 <- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                        dadid = Pedigree_All$father,
                        momid = Pedigree_All$mother)
    
    ####################### Calculate MAF Distributoin ###########################
    G_new<-pullSnpGeno(G8_pop, snpChip = "Inb")
    maf_values_new <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                                 Gen = "G8")
    
    maf_values_df <- rbind(maf_values_df,maf_values_new) 
    
    ################ Assign phenotypic values ####################################
    New_Pheno <- data.frame()
    
    # Should I use the same initial error variance here?
    
    
    for (i in 1:10) {
      pheno <- setPheno(G8_pop,
                        varE = initVarE,
                        corE = matrix(c(1, CorP,
                                        CorP, 1), nrow = 2),
                        rep = 1,
                        simParam = NULL,
                        onlyPheno = T)
      New_Pheno <- rbind(
        New_Pheno,
        data.frame(ID = G8_pop@id, Rep = i, Pheno = pheno)
      )
    }
    
    New_Pheno$Rep<-NULL
    
    GLOBAL_Phenotypes<-rbind(GLOBAL_Phenotypes,New_Pheno)
    Phenotypes_5G <- GLOBAL_Phenotypes[GLOBAL_Phenotypes$ID %in% rownames(Pedigree_5G),]
    
    New_Pheno_mean<- Phenotypes_5G %>%
      group_by(ID) %>%
      summarise(
        Trait1 = mean(Pheno.Trait1, na.rm = TRUE),
        Trait2 = mean(Pheno.Trait2, na.rm = TRUE)
      )
    
    New_Pheno_mean$ID <- as.character(New_Pheno_mean$ID)
    
    # Reorder df_means based on the rownames of K
    New_Pheno_mean <- New_Pheno_mean[match(rownames(A_Mat), New_Pheno_mean$ID), ]
    
    y<-as.matrix(New_Pheno_mean[,2:3])
    
    #### Model ####
    new_model<-mkr(y,A_Mat)
    
    # Save EBVs
    EBVs<-as.data.frame(new_model$hat)
    EBVs$ID <- New_Pheno_mean$ID
    
    EBVs <- EBVs[EBVs$ID %in% G8_pop@id,]
    new_Pheno <- New_Pheno_mean[New_Pheno_mean$ID %in% G8_pop@id,]
    
    G8_pop@ebv <- as.matrix(EBVs[,1:2])
    G8_pop@pheno <- as.matrix(new_Pheno[,2:3])
    
    ##### Select parents ######
    Parents <- selectWithinFam(G8_pop,
                               nInd = 1,
                               use = "ebv",
                               selectTop = TRUE)
    
    # Heritability
    h2<-var(bv(G8_pop))/varP(G8_pop)
    h2_hat<-new_model$h2
    
    True_h2_1 = h2[1,1]
    True_h2_2 = h2[2,2]
    Est_h2_1 = h2_hat[1]
    Est_h2_2 = h2_hat[2]
    
    # Additive genetic value   
    MeanG <- meanG(G8_pop)  
    VarA <- varA(G8_pop)
    Cor_A = cov2cor(VarA)[2, 1]
    
    # Dominance genetic value
    VarD <- varD(G8_pop)
    
    # Phenotypes
    Phen.Mean <- meanP(G8_pop)  
    VarP <- varP(G8_pop)
    Cor_P = cov2cor(VarP)[2, 1]
    
    # Prediction accuracy 
    Acc <- cor(bv(G8_pop),G8_pop@ebv) 
    
    # Within-family accuracy
    WF_df <- data.frame(ID = G8_pop@id,
                        Mum = G8_pop@mother,
                        Dad = G8_pop@father,
                        TBV1 = Bv[,1],
                        TBV2 = Bv[,2],
                        EBV1 = G8_pop@ebv[,1],
                        EBV2 = G8_pop@ebv[,2])
    WF_df$Family <- as.integer(interaction(WF_df$Mum, WF_df$Dad, drop = TRUE))
    
    summary_df <- WF_df %>%
      group_by(Family) %>%
      summarize(
        famSize = n(),
        accuracy1 = if(n() > 1) cor(TBV1, EBV1) else NA_real_,
        accuracy2 = if(n() > 1) cor(TBV2, EBV2) else NA_real_,
      ) %>%
      ungroup()
    
    
    WF_Acc1<-mean(summary_df$accuracy1, na.rm = T)
    WF_Acc2<-mean(summary_df$accuracy2, na.rm = T)
    
    T1_plot<-ggplot(summary_df, aes(x = Family, y = accuracy1, size = famSize)) +
      geom_point(alpha = 0.4, color = "red") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "G8",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    T2_plot<-ggplot(summary_df, aes(x = Family, y = accuracy2, size = famSize)) +
      geom_point(alpha = 0.4, color = "darkred") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    grid.arrange(T1_plot,T2_plot)
    
    # Prediction Bias
    Bv <- bv(G8_pop)
    Bias_1 <- mean(G8_pop@ebv[,1] - Bv[1])
    Bias_2 <- mean(G8_pop@ebv[,2] - Bv[2])
    
    # Estimate disversion 
    bv1 <- Bv[,1]
    bv2 <- Bv[,2]
    ebv1 <- G8_pop@ebv[,1]
    ebv2 <- G8_pop@ebv[,2]
    
    # Fit regression model
    bias_model1 <- lm(bv1 ~ ebv1)
    bias_model2 <- lm(bv2 ~ ebv2)
    
    # Optional: view slope too
    slope1 <- coef(bias_model1)[2]
    slope2 <- coef(bias_model2)[2]
    
    # Selection Intensity
    #Estimated Si
    mean_EBV_p_1 <- mean(Parents@ebv[,1])
    mean_EBV_p_2 <- mean(Parents@ebv[,2])
    
    mean_EBV_1 <- mean(G8_pop@ebv[,1])
    mean_EBV_2 <- mean(G8_pop@ebv[,2])
    
    sd_EBV_1 <- sd(G8_pop@ebv[,1])
    sd_EBV_2 <- sd(G8_pop@ebv[,2])
    
    EBV_Intensity_1 = ((mean_EBV_p_1 - mean_EBV_1) / sd_EBV_1)
    EBV_Intensity_2 = ((mean_EBV_p_2 - mean_EBV_2) / sd_EBV_2)
    
    #True SI
    bv_P<-bv(Parents)
    
    mean_TBV_p_1 <- mean(bv_P[,1])
    mean_TBV_p_2 <- mean(bv_P[,2])
    
    mean_TBV_1 <- mean(Bv[,1])
    mean_TBV_2 <- mean(Bv[,2])
    
    sd_TBV_1 <- sd(Bv[,1])
    sd_TBV_2 <- sd(Bv[,2])
    
    TBV_Intensity_1 = ((mean_TBV_p_1 - mean_TBV_1) / sd_TBV_1)
    TBV_Intensity_2 = ((mean_TBV_p_2 - mean_TBV_2) / sd_TBV_2)
    
    new_trait <- data.frame(Trait1_h2_est = Est_h2_1,
                            Trait1_h2 = True_h2_1,
                            Trait1_MeanG = MeanG[1],
                            Trait1_VarA = VarA[1,1],
                            Trait1_VarD = VarD[1,1],
                            Trait1_MeanP = Phen.Mean[1],
                            Trait1_VarP = VarP[1,1],
                            Trait1_Acc = Acc[1,1],
                            Trait1_WF_Acc = WF_Acc1,
                            Trait1_Bias = Bias_1,
                            Trait1_Dispersion = slope1,
                            Trait1_SI_est = EBV_Intensity_1,
                            Trait1_SI = TBV_Intensity_1,
                            
                            Trait2_h2_est = Est_h2_2,
                            Trait2_h2 = True_h2_2,
                            Trait2_MeanG = MeanG[2],
                            Trait2_VarA = VarA[2,2],
                            Trait2_VarD = VarD[2,2],
                            Trait2_MeanP = Phen.Mean[2],
                            Trait2_VarP = VarP[2,2],
                            Trait2_Acc = Acc[2,2],
                            Trait2_WF_Acc = WF_Acc2,
                            Trait2_Bias = Bias_2,
                            Trait2_Dispersion = slope2,
                            Trait2_SI_est = EBV_Intensity_2,
                            Trait2_SI = TBV_Intensity_2,
                            
                            CovA = VarA[1,2],
                            CorA = Cor_A,
                            CovP = VarP[1,2],
                            CorP = Cor_P,
                            
                            Gen = "G8")
    
    
    
    GLOBAL_trait <- rbind(GLOBAL_trait, new_trait)
    
    # Diversity metrics
    indices<-which((rownames(A_Mat2) %in% G8_pop@id))
    A_Mat2 <-A_Mat2[indices,indices]
    n <- nrow(A_Mat2)
    
    group_coancestry <- sum(A_Mat2) / (n^2)
    
    status_number <- 1 / (2*group_coancestry)
    
    Ns_new <- data.frame(Status_number = status_number,
                         Coancestry = group_coancestry,
                         Gen = "G8")
    
    new_ID <- data.frame(Inb = (diag(A_Mat2)-1),
                         Trait1 = G8_pop@pheno[,1],
                         Trait2 = G8_pop@pheno[,2],
                         Gen = "G8")
    
    GLOBAL_ID <- rbind(GLOBAL_ID, new_ID)
    
    Inbreeding<-mean(diag(A_Mat))-1
    
    Inb<-(diag(A_Mat2)-1)
    
    Inb_df <- data.frame(I = Inb,
                         T1 = G8_pop@pheno[,1],
                         T2 = G8_pop@pheno[,2])
    
    model <- lm(T1 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID1_est<-slope/mean(Inb_df$T1)
    
    model <- lm(T2 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID2_est<-slope/mean(Inb_df$T2)
    
    new_diversity <- data.frame(F.inb = Inbreeding,
                                Ns = status_number,
                                Co = group_coancestry,
                                ID_est1 = ID1_est,
                                ID_est2 = ID2_est,
                                Gen = "G8")
    
    GLOBAL_diversity <- rbind(GLOBAL_diversity, new_diversity)
    
#### G9 ####
    
    # Generate progeny
    # Each cross gets a random numner of progenies, with minimum 1
    G_pops <- vector("list", nCross)
    for (i in 1:nCross) {
      Prog_number <- max(round(rnorm(1, mean = nProg, sd = 10)), 1)
      G_pops[[i]] <- randCross(
        pop = Parents, 
        nCrosses = 1, 
        nProgeny = Prog_number, 
        ignoreSexes = TRUE
      )
    }
    
    # Combine all families into one pop-object
    G9_pop <- do.call(c, G_pops)
    
    # Save New breeding values
    Bv <- bv(G9_pop)
    
    # Update pedigree
    Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G9_pop@id,])
    Pedigree_All <- rbind(Pedigree_All, Pedigree_New)
    
    Pedigree_5G <- c(G4_pop@id,G5_pop@id,G6_pop@id,G7_pop@id,G8_pop@id,
                     G9_pop@id)
    Pedigree_5G <- Pedigree_All[row.names(Pedigree_All) %in% Pedigree_5G,]
    
    A_Mat <- 2*kinship(id = as.numeric(rownames(Pedigree_5G)),
                       dadid = Pedigree_5G$father,
                       momid = Pedigree_5G$mother)
    
    # Separate A-matrix for calculating co-ancestry
    A_Mat2 <- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                       dadid = Pedigree_All$father,
                       momid = Pedigree_All$mother)
    
    
    ####################### Calculate MAF Distributoin ###########################
    G_new<-pullSnpGeno(G9_pop, snpChip = "Inb")
    maf_values_new <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                                 Gen = "G9")
    
    maf_values_df <- rbind(maf_values_df,maf_values_new) 
    
    ################ Assign phenotypic values ####################################
    New_Pheno <- data.frame()
    
    # Should I use the same initial error variance here?
    
    
    for (i in 1:10) {
      pheno <- setPheno(G9_pop,
                        varE = initVarE,
                        corE = matrix(c(1, CorP,
                                        CorP, 1), nrow = 2),
                        rep = 1,
                        simParam = NULL,
                        onlyPheno = T)
      New_Pheno <- rbind(
        New_Pheno,
        data.frame(ID = G9_pop@id, Rep = i, Pheno = pheno)
      )
    }
    
    New_Pheno$Rep<-NULL
    
    GLOBAL_Phenotypes<-rbind(GLOBAL_Phenotypes,New_Pheno)
    Phenotypes_5G <- GLOBAL_Phenotypes[GLOBAL_Phenotypes$ID %in% rownames(Pedigree_5G),]
    
    New_Pheno_mean<- Phenotypes_5G %>%
      group_by(ID) %>%
      summarise(
        Trait1 = mean(Pheno.Trait1, na.rm = TRUE),
        Trait2 = mean(Pheno.Trait2, na.rm = TRUE)
      )
    
    New_Pheno_mean$ID <- as.character(New_Pheno_mean$ID)
    
    # Reorder df_means based on the rownames of K
    New_Pheno_mean <- New_Pheno_mean[match(rownames(A_Mat), New_Pheno_mean$ID), ]
    
    y<-as.matrix(New_Pheno_mean[,2:3])
    
    #### Model ####
    new_model<-mkr(y,A_Mat)
    
    # Save EBVs
    EBVs<-as.data.frame(new_model$hat)
    EBVs$ID <- New_Pheno_mean$ID
    
    EBVs <- EBVs[EBVs$ID %in% G9_pop@id,]
    new_Pheno <- New_Pheno_mean[New_Pheno_mean$ID %in% G9_pop@id,]
    
    G9_pop@ebv <- as.matrix(EBVs[,1:2])
    G9_pop@pheno <- as.matrix(new_Pheno[,2:3])
    
    ##### Select parents ######
    Parents <- selectWithinFam(G9_pop,
                               nInd = 1,
                               use = "ebv",
                               selectTop = TRUE)
    
    # Heritability
    h2<-var(bv(G9_pop))/varP(G9_pop)
    h2_hat<-new_model$h2
    
    True_h2_1 = h2[1,1]
    True_h2_2 = h2[2,2]
    Est_h2_1 = h2_hat[1]
    Est_h2_2 = h2_hat[2]
    
    # Additive genetic value   
    MeanG <- meanG(G9_pop)  
    VarA <- varA(G9_pop)
    Cor_A = cov2cor(VarA)[2, 1]
    
    # Dominance genetic value
    VarD <- varD(G9_pop)
    
    # Phenotypes
    Phen.Mean <- meanP(G9_pop)  
    VarP <- varP(G9_pop)
    Cor_P = cov2cor(VarP)[2, 1]
    
    # Prediction accuracy 
    Acc <- cor(bv(G9_pop),G9_pop@ebv) 
    
    # Within-family accuracy
    WF_df <- data.frame(ID = G9_pop@id,
                        Mum = G9_pop@mother,
                        Dad = G9_pop@father,
                        TBV1 = Bv[,1],
                        TBV2 = Bv[,2],
                        EBV1 = G9_pop@ebv[,1],
                        EBV2 = G9_pop@ebv[,2])
    WF_df$Family <- as.integer(interaction(WF_df$Mum, WF_df$Dad, drop = TRUE))
    
    summary_df <- WF_df %>%
      group_by(Family) %>%
      summarize(
        famSize = n(),
        accuracy1 = if(n() > 1) cor(TBV1, EBV1) else NA_real_,
        accuracy2 = if(n() > 1) cor(TBV2, EBV2) else NA_real_,
      ) %>%
      ungroup()
    
    
    WF_Acc1<-mean(summary_df$accuracy1, na.rm = T)
    WF_Acc2<-mean(summary_df$accuracy2, na.rm = T)
    
    T1_plot<-ggplot(summary_df, aes(x = Family, y = accuracy1, size = famSize)) +
      geom_point(alpha = 0.4, color = "red") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "G9",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    T2_plot<-ggplot(summary_df, aes(x = Family, y = accuracy2, size = famSize)) +
      geom_point(alpha = 0.4, color = "darkred") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    grid.arrange(T1_plot,T2_plot)
    
    # Prediction Bias
    Bv <- bv(G9_pop)
    Bias_1 <- mean(G9_pop@ebv[,1] - Bv[1])
    Bias_2 <- mean(G9_pop@ebv[,2] - Bv[2])
    
    # Estimate disversion 
    bv1 <- Bv[,1]
    bv2 <- Bv[,2]
    ebv1 <- G9_pop@ebv[,1]
    ebv2 <- G9_pop@ebv[,2]
    
    # Fit regression model
    bias_model1 <- lm(bv1 ~ ebv1)
    bias_model2 <- lm(bv2 ~ ebv2)
    
    # Optional: view slope too
    slope1 <- coef(bias_model1)[2]
    slope2 <- coef(bias_model2)[2]
    
    # Selection Intensity
    #Estimated Si
    mean_EBV_p_1 <- mean(Parents@ebv[,1])
    mean_EBV_p_2 <- mean(Parents@ebv[,2])
    
    mean_EBV_1 <- mean(G9_pop@ebv[,1])
    mean_EBV_2 <- mean(G9_pop@ebv[,2])
    
    sd_EBV_1 <- sd(G9_pop@ebv[,1])
    sd_EBV_2 <- sd(G9_pop@ebv[,2])
    
    EBV_Intensity_1 = ((mean_EBV_p_1 - mean_EBV_1) / sd_EBV_1)
    EBV_Intensity_2 = ((mean_EBV_p_2 - mean_EBV_2) / sd_EBV_2)
    
    #True SI
    bv_P<-bv(Parents)
    
    mean_TBV_p_1 <- mean(bv_P[,1])
    mean_TBV_p_2 <- mean(bv_P[,2])
    
    mean_TBV_1 <- mean(Bv[,1])
    mean_TBV_2 <- mean(Bv[,2])
    
    sd_TBV_1 <- sd(Bv[,1])
    sd_TBV_2 <- sd(Bv[,2])
    
    TBV_Intensity_1 = ((mean_TBV_p_1 - mean_TBV_1) / sd_TBV_1)
    TBV_Intensity_2 = ((mean_TBV_p_2 - mean_TBV_2) / sd_TBV_2)
    
    new_trait <- data.frame(Trait1_h2_est = Est_h2_1,
                            Trait1_h2 = True_h2_1,
                            Trait1_MeanG = MeanG[1],
                            Trait1_VarA = VarA[1,1],
                            Trait1_VarD = VarD[1,1],
                            Trait1_MeanP = Phen.Mean[1],
                            Trait1_VarP = VarP[1,1],
                            Trait1_Acc = Acc[1,1],
                            Trait1_WF_Acc = WF_Acc1,
                            Trait1_Bias = Bias_1,
                            Trait1_Dispersion = slope1,
                            Trait1_SI_est = EBV_Intensity_1,
                            Trait1_SI = TBV_Intensity_1,
                            
                            Trait2_h2_est = Est_h2_2,
                            Trait2_h2 = True_h2_2,
                            Trait2_MeanG = MeanG[2],
                            Trait2_VarA = VarA[2,2],
                            Trait2_VarD = VarD[2,2],
                            Trait2_MeanP = Phen.Mean[2],
                            Trait2_VarP = VarP[2,2],
                            Trait2_Acc = Acc[2,2],
                            Trait2_WF_Acc = WF_Acc2,
                            Trait2_Bias = Bias_2,
                            Trait2_Dispersion = slope2,
                            Trait2_SI_est = EBV_Intensity_2,
                            Trait2_SI = TBV_Intensity_2,
                            
                            CovA = VarA[1,2],
                            CorA = Cor_A,
                            CovP = VarP[1,2],
                            CorP = Cor_P,
                            
                            Gen = "G9")
    
    
    
    GLOBAL_trait <- rbind(GLOBAL_trait, new_trait)
    
    # Diversity metrics
    indices<-which((rownames(A_Mat2) %in% G9_pop@id))
    A_Mat2 <-A_Mat2[indices,indices]
    n <- nrow(A_Mat2)
    
    group_coancestry <- sum(A_Mat2) / (n^2)
    
    status_number <- 1 / (2*group_coancestry)
    
    Ns_new <- data.frame(Status_number = status_number,
                         Coancestry = group_coancestry,
                         Gen = "G9")
    
    
    Inbreeding<-mean(diag(A_Mat2))-1
    
    new_ID <- data.frame(Inb = (diag(A_Mat2)-1),
                         Trait1 = G9_pop@pheno[,1],
                         Trait2 = G9_pop@pheno[,2],
                         Gen = "G9")
    
    GLOBAL_ID <- rbind(GLOBAL_ID, new_ID)
    
    # Calculate the ID
    Inb<-(diag(A_Mat2)-1)
    
    Inb_df <- data.frame(I = Inb,
                         T1 = G9_pop@pheno[,1],
                         T2 = G9_pop@pheno[,2])
    
    model <- lm(T1 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID1_est<-slope/mean(Inb_df$T1)
    
    model <- lm(T2 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID2_est<-slope/mean(Inb_df$T2)
    
    new_diversity <- data.frame(F.inb = Inbreeding,
                                Ns = status_number,
                                Co = group_coancestry,
                                ID_est1 = ID1_est,
                                ID_est2 = ID2_est,
                                Gen = "G9")
    
    GLOBAL_diversity <- rbind(GLOBAL_diversity, new_diversity)
    
#### G10 ####
    
    # Generate progeny
    # Each cross gets a random numner of progenies, with minimum 1
    G_pops <- vector("list", nCross)
    for (i in 1:nCross) {
      Prog_number <- max(round(rnorm(1, mean = nProg, sd = 10)), 1)
      G_pops[[i]] <- randCross(
        pop = Parents, 
        nCrosses = 1, 
        nProgeny = Prog_number, 
        ignoreSexes = TRUE
      )
    }
    
    # Combine all families into one pop-object
    G10_pop <- do.call(c, G_pops)
    
    # Save New breeding values
    Bv <- bv(G10_pop)
    
    # Update pedigree
    Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G10_pop@id,])
    Pedigree_All <- rbind(Pedigree_All, Pedigree_New)
    
    Pedigree_5G <- c(G5_pop@id,G6_pop@id,G7_pop@id,G8_pop@id,G9_pop@id,
                     G10_pop@id)
    Pedigree_5G <- Pedigree_All[row.names(Pedigree_All) %in% Pedigree_5G,]
    
    A_Mat <- 2*kinship(id = as.numeric(rownames(Pedigree_5G)),
                       dadid = Pedigree_5G$father,
                       momid = Pedigree_5G$mother)
    
    # Separate A-matrix for calculating co-ancestry
    A_Mat2 <- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                        dadid = Pedigree_All$father,
                        momid = Pedigree_All$mother)
    
    ####################### Calculate MAF Distributoin ###########################
    G_new<-pullSnpGeno(G10_pop, snpChip = "Inb")
    maf_values_new <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                                 Gen = "G10")
    
    maf_values_df <- rbind(maf_values_df,maf_values_new) 
    
    ################ Assign phenotypic values ####################################
    New_Pheno <- data.frame()
    
    # Should I use the same initial error variance here?
    
    
    for (i in 1:10) {
      pheno <- setPheno(G10_pop,
                        varE = initVarE,
                        corE = matrix(c(1, CorP,
                                        CorP, 1), nrow = 2),
                        rep = 1,
                        simParam = NULL,
                        onlyPheno = T)
      New_Pheno <- rbind(
        New_Pheno,
        data.frame(ID = G10_pop@id, Rep = i, Pheno = pheno)
      )
    }
    
    New_Pheno$Rep<-NULL
    
    GLOBAL_Phenotypes<-rbind(GLOBAL_Phenotypes,New_Pheno)
    Phenotypes_5G <- GLOBAL_Phenotypes[GLOBAL_Phenotypes$ID %in% rownames(Pedigree_5G),]
    
    New_Pheno_mean<- Phenotypes_5G %>%
      group_by(ID) %>%
      summarise(
        Trait1 = mean(Pheno.Trait1, na.rm = TRUE),
        Trait2 = mean(Pheno.Trait2, na.rm = TRUE)
      )
    
    New_Pheno_mean$ID <- as.character(New_Pheno_mean$ID)
    
    # Reorder df_means based on the rownames of K
    New_Pheno_mean <- New_Pheno_mean[match(rownames(A_Mat), New_Pheno_mean$ID), ]
    
    y<-as.matrix(New_Pheno_mean[,2:3])
    
    #### Model ####
    new_model<-mkr(y,A_Mat)
    
    # Save EBVs
    EBVs<-as.data.frame(new_model$hat)
    EBVs$ID <- New_Pheno_mean$ID
    
    EBVs <- EBVs[EBVs$ID %in% G10_pop@id,]
    new_Pheno <- New_Pheno_mean[New_Pheno_mean$ID %in% G10_pop@id,]
    
    G10_pop@ebv <- as.matrix(EBVs[,1:2])
    G10_pop@pheno <- as.matrix(new_Pheno[,2:3])
    
    ##### Select parents ######
    Parents <- selectWithinFam(G10_pop,
                               nInd = 1,
                               use = "ebv",
                               selectTop = TRUE)
    
    # Heritability
    h2<-var(bv(G10_pop))/varP(G10_pop)
    h2_hat<-new_model$h2
    
    True_h2_1 = h2[1,1]
    True_h2_2 = h2[2,2]
    Est_h2_1 = h2_hat[1]
    Est_h2_2 = h2_hat[2]
    
    # Additive genetic value   
    MeanG <- meanG(G10_pop)  
    VarA <- varA(G10_pop)
    Cor_A = cov2cor(VarA)[2, 1]
    
    # Dominance genetic value
    VarD <- varD(G10_pop)
    
    # Phenotypes
    Phen.Mean <- meanP(G10_pop)  
    VarP <- varP(G10_pop)
    Cor_P = cov2cor(VarP)[2, 1]
    
    # Prediction accuracy 
    Acc <- cor(bv(G10_pop),G10_pop@ebv) 
    
    # Within-family accuracy
    WF_df <- data.frame(ID = G10_pop@id,
                        Mum = G10_pop@mother,
                        Dad = G10_pop@father,
                        TBV1 = Bv[,1],
                        TBV2 = Bv[,2],
                        EBV1 = G10_pop@ebv[,1],
                        EBV2 = G10_pop@ebv[,2])
    WF_df$Family <- as.integer(interaction(WF_df$Mum, WF_df$Dad, drop = TRUE))
    
    summary_df <- WF_df %>%
      group_by(Family) %>%
      summarize(
        famSize = n(),
        accuracy1 = if(n() > 1) cor(TBV1, EBV1) else NA_real_,
        accuracy2 = if(n() > 1) cor(TBV2, EBV2) else NA_real_,
      ) %>%
      ungroup()
    
    
    WF_Acc1<-mean(summary_df$accuracy1, na.rm = T)
    WF_Acc2<-mean(summary_df$accuracy2, na.rm = T)
    
    T1_plot<-ggplot(summary_df, aes(x = Family, y = accuracy1, size = famSize)) +
      geom_point(alpha = 0.4, color = "red") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "G10",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    T2_plot<-ggplot(summary_df, aes(x = Family, y = accuracy2, size = famSize)) +
      geom_point(alpha = 0.4, color = "darkred") +
      scale_size_continuous(name = "Family size") +
      labs(
        title = "",
        x = "Family ID",
        y = "Accuracy"
      ) +
      theme_minimal()
    
    grid.arrange(T1_plot,T2_plot)
    
    # Prediction Bias
    Bv <- bv(G10_pop)
    Bias_1 <- mean(G10_pop@ebv[,1] - Bv[1])
    Bias_2 <- mean(G10_pop@ebv[,2] - Bv[2])
    
    # Estimate disversion 
    bv1 <- Bv[,1]
    bv2 <- Bv[,2]
    ebv1 <- G10_pop@ebv[,1]
    ebv2 <- G10_pop@ebv[,2]
    
    # Fit regression model
    bias_model1 <- lm(bv1 ~ ebv1)
    bias_model2 <- lm(bv2 ~ ebv2)
    
    # Optional: view slope too
    slope1 <- coef(bias_model1)[2]
    slope2 <- coef(bias_model2)[2]
    
    # Selection Intensity
    #Estimated Si
    mean_EBV_p_1 <- mean(Parents@ebv[,1])
    mean_EBV_p_2 <- mean(Parents@ebv[,2])
    
    mean_EBV_1 <- mean(G10_pop@ebv[,1])
    mean_EBV_2 <- mean(G10_pop@ebv[,2])
    
    sd_EBV_1 <- sd(G10_pop@ebv[,1])
    sd_EBV_2 <- sd(G10_pop@ebv[,2])
    
    EBV_Intensity_1 = ((mean_EBV_p_1 - mean_EBV_1) / sd_EBV_1)
    EBV_Intensity_2 = ((mean_EBV_p_2 - mean_EBV_2) / sd_EBV_2)
    
    #True SI
    bv_P<-bv(Parents)
    
    mean_TBV_p_1 <- mean(bv_P[,1])
    mean_TBV_p_2 <- mean(bv_P[,2])
    
    mean_TBV_1 <- mean(Bv[,1])
    mean_TBV_2 <- mean(Bv[,2])
    
    sd_TBV_1 <- sd(Bv[,1])
    sd_TBV_2 <- sd(Bv[,2])
    
    TBV_Intensity_1 = ((mean_TBV_p_1 - mean_TBV_1) / sd_TBV_1)
    TBV_Intensity_2 = ((mean_TBV_p_2 - mean_TBV_2) / sd_TBV_2)
    
    new_trait <- data.frame(Trait1_h2_est = Est_h2_1,
                            Trait1_h2 = True_h2_1,
                            Trait1_MeanG = MeanG[1],
                            Trait1_VarA = VarA[1,1],
                            Trait1_VarD = VarD[1,1],
                            Trait1_MeanP = Phen.Mean[1],
                            Trait1_VarP = VarP[1,1],
                            Trait1_Acc = Acc[1,1],
                            Trait1_WF_Acc = WF_Acc1,
                            
                            Trait1_Bias = Bias_1,
                            Trait1_Dispersion = slope1,
                            Trait1_SI_est = EBV_Intensity_1,
                            Trait1_SI = TBV_Intensity_1,
                            
                            Trait2_h2_est = Est_h2_2,
                            Trait2_h2 = True_h2_2,
                            Trait2_MeanG = MeanG[2],
                            Trait2_VarA = VarA[2,2],
                            Trait2_VarD = VarD[2,2],
                            Trait2_MeanP = Phen.Mean[2],
                            Trait2_VarP = VarP[2,2],
                            Trait2_Acc = Acc[2,2],
                            Trait2_WF_Acc = WF_Acc2,
                            
                            Trait2_Bias = Bias_2,
                            Trait2_Dispersion = slope2,
                            Trait2_SI_est = EBV_Intensity_2,
                            Trait2_SI = TBV_Intensity_2,
                            
                            CovA = VarA[1,2],
                            CorA = Cor_A,
                            CovP = VarP[1,2],
                            CorP = Cor_P,
                            
                            Gen = "G10")
    
    
    
    GLOBAL_trait <- rbind(GLOBAL_trait, new_trait)
    
    # Diversity metrics
    indices<-which((rownames(A_Mat2) %in% G10_pop@id))
    A_Mat2 <-A_Mat2[indices,indices]
    n <- nrow(A_Mat2)
    
    group_coancestry <- sum(A_Mat2) / (n^2)
    
    status_number <- 1 / (2*group_coancestry)
    
    Ns_new <- data.frame(Status_number = status_number,
                         Coancestry = group_coancestry,
                         Gen = "G10")
    
    
    Inbreeding<-mean(diag(A_Mat2))-1
    
    # Calculate ID
    new_ID <- data.frame(Inb = (diag(A_Mat2)-1),
                         Trait1 = G10_pop@pheno[,1],
                         Trait2 = G10_pop@pheno[,2],
                         Gen = "G10")
    
    GLOBAL_ID <- rbind(GLOBAL_ID, new_ID)
    
    Inb<-(diag(A_Mat2)-1)
    
    Inb_df <- data.frame(I = Inb,
                         T1 = G10_pop@pheno[,1],
                         T2 = G10_pop@pheno[,2])
    
    model <- lm(T1 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID1_est<-slope/mean(Inb_df$T1)
    
    model <- lm(T2 ~ I, data = Inb_df)
    slope <- round(coef(model)[2], 2)
    ID2_est<-slope/mean(Inb_df$T2)
    
    new_diversity <- data.frame(F.inb = Inbreeding,
                                Ns = status_number,
                                Co = group_coancestry,
                                ID_est1 = ID1_est,
                                ID_est2 = ID2_est,
                                Gen = "G10")
    
    GLOBAL_diversity <- rbind(GLOBAL_diversity, new_diversity)    
    
  
######################### Visualize results ####################################

# Additive trait correlation
plot(x = c(1:11), y = GLOBAL_trait$CorA, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "Additive correlation")

# Phenotypic trait correlation
plot(x = c(1:11), y = GLOBAL_trait$CorP, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "Phenotypic correlation")

# Heritability
varRanges = range(c(GLOBAL_trait$Trait1_h2_est, GLOBAL_trait$Trait2_h2_est))
plot(x = c(1:11), y = GLOBAL_trait$Trait1_h2_est, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "h2_est", ylim = varRanges) +
  lines(x = c(1:11), y = GLOBAL_trait$Trait2_h2_est, type = "l", col = "purple", lty = 2, lwd = 3) 
legend(x = "topleft", legend = c("1", "2"), title = "Trait",
       lwd = 3, lty = c(1, 2), col = c("black", "purple"))

# True Heritability
varRanges = range(c(GLOBAL_trait$Trait1_h2, GLOBAL_trait$Trait2_h2))
plot(x = c(1:11), y = GLOBAL_trait$Trait1_h2, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "h2", ylim = varRanges) +
  lines(x = c(1:11), y = GLOBAL_trait$Trait2_h2, type = "l", col = "purple", lty = 2, lwd = 3) 
legend(x = "topleft", legend = c("1", "2"), title = "Trait",
       lwd = 3, lty = c(1, 2), col = c("black", "purple"))

# Bias
varRanges = range(c(GLOBAL_trait$Trait1_Bias, GLOBAL_trait$Trait2_Bias))
plot(x = c(1:11), y = GLOBAL_trait$Trait1_Bias, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "Bias", ylim = varRanges) +
  lines(x = c(1:11), y = GLOBAL_trait$Trait2_Bias, type = "l", col = "purple", lty = 2, lwd = 3) 
legend(x = "topleft", legend = c("1", "2"), title = "Trait",
       lwd = 3, lty = c(1, 2), col = c("black", "purple"))

# Additive variance
varRanges = range(c(GLOBAL_trait$Trait1_VarA, GLOBAL_trait$Trait2_VarA))
plot(x = c(1:11), y = GLOBAL_trait$Trait1_VarA, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "VarA", ylim = varRanges) +
  lines(x = c(1:11), y = GLOBAL_trait$Trait2_VarA, type = "l", col = "purple", lty = 2, lwd = 3) 
legend(x = "topleft", legend = c("1", "2"), title = "Trait",
       lwd = 3, lty = c(1, 2), col = c("black", "purple"))

# Genetic Value
varRanges = range(c(GLOBAL_trait$Trait1_MeanG, GLOBAL_trait$Trait2_MeanG))
plot(x = c(1:11), y = GLOBAL_trait$Trait1_MeanG, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "MeanG", ylim = varRanges) +
  lines(x = c(1:11), y = GLOBAL_trait$Trait2_MeanG, type = "l", col = "purple", lty = 2, lwd = 3) 
legend(x = "topleft", legend = c("1", "2"), title = "Trait",
       lwd = 3, lty = c(1, 2), col = c("black", "purple"))

# Dominance variance
varRanges = range(c(GLOBAL_trait$Trait1_VarD, GLOBAL_trait$Trait2_VarD))
plot(x = c(1:11), y = GLOBAL_trait$Trait1_VarD, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "VarD", ylim = varRanges) +
  lines(x = c(1:11), y = GLOBAL_trait$Trait2_VarD, type = "l", col = "purple", lty = 2, lwd = 3) 
legend(x = "topleft", legend = c("1", "2"), title = "Trait",
       lwd = 3, lty = c(1, 2), col = c("black", "purple"))

# Phenotypic variance
varRanges = range(c(GLOBAL_trait$Trait1_VarP, GLOBAL_trait$Trait2_VarP))
plot(x = c(1:11), y = GLOBAL_trait$Trait1_VarP, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "VarP", ylim = varRanges) +
  lines(x = c(1:11), y = GLOBAL_trait$Trait2_VarP, type = "l", col = "purple", lty = 2, lwd = 3) 
legend(x = "topleft", legend = c("1", "2"), title = "Trait",
       lwd = 3, lty = c(1, 2), col = c("black", "purple"))

# Phenotypic mean
varRanges = range(c(GLOBAL_trait$Trait1_MeanP, GLOBAL_trait$Trait2_MeanP))
plot(x = c(1:11), y = GLOBAL_trait$Trait1_MeanP, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "MeanP", ylim = varRanges) +
  lines(x = c(1:11), y = GLOBAL_trait$Trait2_MeanP, type = "l", col = "purple", lty = 2, lwd = 3) 
legend(x = "topleft", legend = c("1", "2"), title = "Trait",
       lwd = 3, lty = c(1, 2), col = c("black", "purple"))

# Prediction Accuracy
varRanges = range(c(GLOBAL_trait$Trait1_Acc, GLOBAL_trait$Trait2_Acc))
plot(x = c(1:11), y = GLOBAL_trait$Trait1_Acc, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "Acc", ylim = varRanges) +
  lines(x = c(1:11), y = GLOBAL_trait$Trait2_Acc, type = "l", col = "purple", lty = 2, lwd = 3) 
legend(x = "topleft", legend = c("1", "2"), title = "Trait",
       lwd = 3, lty = c(1, 2), col = c("black", "purple"))

varRanges = range(c(GLOBAL_trait$Trait1_WF_Acc, GLOBAL_trait$Trait2_WF_Acc))
plot(x = c(1:11), y = GLOBAL_trait$Trait1_WF_Acc, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "Within-Family Acc", ylim = varRanges) +
  lines(x = c(1:11), y = GLOBAL_trait$Trait2_WF_Acc, type = "l", col = "purple", lty = 2, lwd = 3) 
legend(x = "topleft", legend = c("1", "2"), title = "Trait",
       lwd = 3, lty = c(1, 2), col = c("black", "purple"))

# Dispersion
varRanges = range(c(GLOBAL_trait$Trait1_Dispersion, GLOBAL_trait$Trait2_Dispersion))
plot(x = c(1:11), y = GLOBAL_trait$Trait1_Dispersion, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "Dispersion", ylim = varRanges) +
  lines(x = c(1:11), y = GLOBAL_trait$Trait2_Dispersion, type = "l", col = "purple", lty = 2, lwd = 3) 
legend(x = "topleft", legend = c("1", "2"), title = "Trait",
       lwd = 3, lty = c(1, 2), col = c("black", "purple"))

# Selection intensities
varRanges = range(c(GLOBAL_trait$Trait1_SI_est, GLOBAL_trait$Trait2_SI_est))
plot(x = c(1:11), y = GLOBAL_trait$Trait1_SI_est, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "SI_est", ylim = varRanges) +
  lines(x = c(1:11), y = GLOBAL_trait$Trait2_SI_est, type = "l", col = "purple", lty = 2, lwd = 3) 
legend(x = "topleft", legend = c("1", "2"), title = "Trait",
       lwd = 3, lty = c(1, 2), col = c("black", "purple"))

# Selection intensity true
varRanges = range(c(GLOBAL_trait$Trait1_SI, GLOBAL_trait$Trait2_SI))
plot(x = c(1:11), y = GLOBAL_trait$Trait1_SI, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "SI", ylim = varRanges) +
  lines(x = c(1:11), y = GLOBAL_trait$Trait2_SI, type = "l", col = "purple", lty = 2, lwd = 3) 
legend(x = "topleft", legend = c("1", "2"), title = "Trait",
       lwd = 3, lty = c(1, 2), col = c("black", "purple"))


# Inbreeding
plot(x = c(1:11), y = GLOBAL_diversity$F.inb, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "Mean inbreeding")

# Status number
plot(x = c(1:11), y = GLOBAL_diversity$Ns, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "Status Number")

# Coancestry
plot(x = c(1:11), y = GLOBAL_diversity$Co, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "Group Coancestry")

# Estimated Inbreeding Depression
varRanges = range(c(GLOBAL_diversity$ID_est1, GLOBAL_diversity$ID_est2))
plot(x = c(1:11), y = GLOBAL_diversity$ID_est1, type = "l", col = "black", lwd = 3,
     xlab = "Generation", ylab = "ID", ylim = varRanges) +
  lines(x = c(1:11), y = GLOBAL_diversity$ID_est2, type = "l", col = "purple", lty = 2, lwd = 3) 
legend(x = "topleft", legend = c("1", "2"), title = "Trait",
       lwd = 3, lty = c(1, 2), col = c("black", "purple"))



# ID across all generations
GLOBAL_ID$Gen <- as.factor(GLOBAL_ID$Gen)
model <- lm(Trait2 ~ Inb + Gen, data = GLOBAL_ID)
slope_global <- coef(model)[2]

plot(GLOBAL_ID$Inb, GLOBAL_ID$Trait2)

maf_values_df$Gen <- as.factor(maf_values_df$Gen)
library(ggridges)
ggplot(maf_values_df, aes(x = Frequency, y = Gen)) +
  geom_density_ridges(fill = "skyblue", alpha = 0.8, scale = 1.1, rel_min_height = 0) +
  xlim(-0.1, 0.6) +
  xlab("Minor Allele Frequency") +
  ylab("Generation") +
  ggtitle("All Frequencies") +
  theme_minimal() +
  theme(legend.position = "none")

maf_values_df_0.01 <- maf_values_df[maf_values_df$Frequency>0.01,]

ggplot(maf_values_df_0.01, aes(x = Frequency, y = Gen)) +
  geom_density_ridges(fill = "deepskyblue", alpha = 0.8, scale = 1.1, rel_min_height = 0) +
  xlim(-0.1, 0.6) +
  xlab("Minor Allele Frequency") +
  ylab("Generation") +
  ggtitle("Min MAF 0.01") +
  theme_minimal() +
  theme(legend.position = "none")

"-----------------------------------End----------------------------------------"




################################################################################
# Record Simulated Parameters
Sim <- data.frame(Chromosomes = PaChr,
                  QTLs = nQtl,
                  SNPs = nSNP,
                  nSNP2 = nSNP2,
                  Chr.Size.bp = ChrLen,
                  Chr.Size.mo = PG_Chr_M,
                  Mut.Rate = ConMR,
                  FoundersNr = nFounders,
                  Founders.ne = neFounders,
                  BurnInGen = nGenBurnIn,
                  BurnInCross = nCrossBurnIn,
                  BurnInProg = nProgBurnIn,
                  BreedingGens = num_generations,
                  ParentsNr = nParents,
                  FamiliesNr = nCross,
                  ProgenyNr = nProg,
                  Trait = Trait,
                  Trait.mean = initMeanG,
                  GenVariance = initVarG,
                  EnvVariance = initVarE,
                  GeneticCorrelation = CorA, 
                  PhenotypicCorrelation = CorP, 
                  h2_G0 = h2_G0,
                  h2_G1 = h2_G1,
                  InbreedingDepression = InbDepr,
                  MatePlan = MatePlan,
                  Gamma = GAMMA,
                  GammaShape = GShape,
                  SegSite = SegSite)

# Simulation Parameters
filename <- paste0("Sim_Parameters_Ne1000_MultiGen", today_date, ".txt")
write.table(Sim, filename, quote = F, col.names = T, row.names = F, sep = "\t")