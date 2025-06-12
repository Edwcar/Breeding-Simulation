# AlphaSim Simulation 

# Scenario B1
# Blind Trauncation Selection GS selection 

# Ideas
# Add genetic gain relative to the founder population in percentage,
# for each new generation

# Compare estimates with AsremlR

# Only one phenotype  from a particular generation is available 
# at selection within that Generation

# Rolling front
# How can we do within-family selection with rolling front?
# First pick the best 50 families across all generations, 
# and then select within them?

# Add mechanism for removing measurements that decrease the height

# Should cross-check with 

options(scipen = 999)  # Avoid scientific notation globally

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
library(ggridges)
library(bWGR)

library(pheatmap)

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
SegSite = 1000

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


# Error variance
initVarE = 4*(CV * initMeanG)^2 

# Additve genetic correlation between traits
CorA <- 0.9
# Phenotypic correlation between traits
CorP <- 0.8

# Heritability
h2 = 0.2

# Founder Parameters    
# Founder Population
nFounders = 100   # Number of founders in the population
neFounders = 20   # Effective population size of founders

nGenBurnIn <- 100  # Number of descrete burn-in generations
nCrossBurnIn <- 100 # Number of crosses within each burn-in generation
nProgBurnIn <- 2

# Breeding Parameters    
# Size of the G0 population
nG0 = 200

# Crosses
nParents = 10   # Number of parents to be selected
nCross = 10     # Number of families to create from Parents
nProg = 20      # Number of progeny per cross
sdFam = 5

# Mating schemes
MatePlan = "RandCross"

# Phenotyping efforts
Ramets <- 12

# Number of generations to simulate
# Define number of generations

# Genotyping Parameters
# SNP chip for GS
nSNP = 100    # Nr of SNPs per chromosome
minMAF = 0.01  # Minimum MAF for inclusion


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
SP$addTraitA(nQtlPerChr   = nQtl,
             mean         = initMeanG,
             var          = initVarG,
             gamma        = T,
             shape        = GShape,
             corA         = matrix(c(1, CorA,
                                     CorA, 1), nrow=2))

# Specify default values for setPheno
# Define heritabilities for both traits (this is REQUIRED before using varE)
SP$setVarE(h2 = c(h2, h2))

# Set residual covariance matrix used by the setPheno function
E <- matrix(c(1, CorP,
              CorP, 1), nrow=2)
SP$setCorE(E)

# Create new pop object
G0_pop = newPop(founderPop)

Acc_df <- data.frame(Trait1 = 0,
                     Trait2 = 0)

################################## Burn-in #####################################
VarA_list1<-list()
VarA_G0 <- varA(G0_pop)
VarA_list1 <- c(VarA_list1, VarA_G0[1,1])

MeanA_list1<-list()
MeanA_G0 <- meanG(G0_pop)
MeanA_list1 <- c(MeanA_list1, MeanA_G0[1])

VarA_list2<-list()
VarA_G0 <- varA(G0_pop)
VarA_list2 <- c(VarA_list2, VarA_G0[2,2])

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
  
  MeanA_G0 <- meanG(G0_pop)
  MeanA_list1 <- c(MeanA_list1, MeanA_G0[1])
  
  VarA_G0 <- varA(G0_pop)
  VarA_list2 <- c(VarA_list2, VarA_G0[2,2])
  
  MeanA_G0 <- meanG(G0_pop)
  MeanA_list2 <- c(MeanA_list2, MeanA_G0[2])
}


time <- seq_along(VarA_list1) 

df <- data.frame(
  Time = rep(time, 2),  # Repeat time for both parameters
  Value = c(unlist(VarA_list1), unlist(MeanA_list1)),  # Combine the values of both parameters
  Parameter = rep(c("VarA", "MeanG"), each = length(time))  # Label each parameter
)

ggplot(df, aes(x = Time, y = Value, color = Parameter)) +
  geom_line(linewidth = 1) +         # Lines for both parameters
  geom_point(size = 1) + 
  scale_color_manual( values = c("VarA" = "darkgreen",  "MeanG" = "darkred")  # Custom colors
  ) +# Points for both parameters
  labs(title = "",
       x = "Generation",
       y = "Value",
       color = "Parameter"
  ) +
  theme_minimal()

df <- data.frame(
  Time = rep(time, 2),  # Repeat time for both parameters
  Value = c(unlist(VarA_list2), unlist(MeanA_list2)),  # Combine the values of both parameters
  Parameter = rep(c("VarA", "MeanG"), each = length(time))  # Label each parameter
)

ggplot(df, aes(x = Time, y = Value, color = Parameter)) +
  geom_line(linewidth = 1) +         # Lines for both parameters
  geom_point(size = 1) + 
  scale_color_manual( values = c("VarA" = "darkgreen",  "MeanG" = "darkred")  # Custom colors
  ) +# Points for both parameters
  labs(title = "",
       x = "Generation",
       y = "Value", 
       color = "Parameter"
  ) +
  theme_minimal()


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

# Save MAF
G_new<-pullSnpGeno(G0_pop)
maf_values_df <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                            Gen = "G0")

M_G0<-pullSnpGeno(G0_pop, snpChip = 1)-1
G = A.mat(M_G0, impute.method = "mean", min.MAF = 0.01)


# Add first phenotypes
G0_Pheno <- data.frame()

for (i in 1:1) {
  pheno <- setPheno(G0_pop,
                    rep = 1,
                    corE = E,
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

# Estimate inbreeding from A-matrix
Inb <- as.data.frame(diag(Kinship_matrix)-1)

# Match the order of the pop object
Inb <- Inb[match(rownames(Inb), G0_pop@id), ]

# Reorder df_means based on the rownames of K
G0_Pheno_mean <- G0_Pheno_mean[match((G0_pop@id), G0_Pheno_mean$ID), ]

G0_Pheno_mean$Trait1_F <- G0_Pheno_mean$Trait1 + G0_Pheno_mean$Trait1 * (Inb * -0.35)
G0_Pheno_mean$Trait2_F <- G0_Pheno_mean$Trait2 + G0_Pheno_mean$Trait2 * (Inb * -0.35)

y<-as.matrix(G0_Pheno_mean[,4:5])

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
                           Trait1_MeanP = Phen.Mean[1],
                           Trait1_VarP = VarP[1,1],
                           Trait1_Acc_Blind = NA,
                           Trait1_Acc = Acc[1,1],
                           Trait1_WF_Acc = NA,
                           Trait1_Bias = Bias_1,
                           Trait1_Dispersion = slope1,
                           Trait1_SI_est = EBV_Intensity_1,
                           Trait1_SI = TBV_Intensity_1,
                           
                           Trait2_h2_est = Est_h2_2,
                           Trait2_h2 = True_h2_2,
                           Trait2_MeanG = MeanG[2],
                           Trait2_VarA = VarA[2,2],
                           Trait2_MeanP = Phen.Mean[2],
                           Trait2_VarP = VarP[2,2],
                           Trait2_Acc_Blind = NA,
                           Trait2_Acc = Acc[2,2],
                           Trait2_WF_Acc = NA,
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
  Prog_number <- max(round(rnorm(1, mean = nProg, sd = sdFam)), 1)
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
G_new<-pullSnpGeno(G1_pop)
maf_values_new <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                             Gen = "G1")

maf_values_df <- rbind(maf_values_df,maf_values_new) 

M_G<-pullSnpGeno(G1_pop, snpChip = 1)-1
M <- rbind(M_G0, M_G)
G = A.mat(M, impute.method = "mean", min.MAF = 0.01)


################ Assign phenotypic values ####################################
New_Pheno <- data.frame()

# Should I use the same initial error variance here?


for (i in 1:Ramets) {
  pheno <- setPheno(G1_pop,
                    rep = 1,
                    corE = E,
                    simParam = NULL,
                    onlyPheno = T)
  New_Pheno <- rbind(New_Pheno,
                     data.frame(ID = G1_pop@id, Rep = i, Pheno = pheno))
}

New_Pheno$Rep<-NULL

New_Pheno_mean<- New_Pheno %>%
  group_by(ID) %>%
  summarise(
    Trait1 = mean(Pheno.Trait1, na.rm = TRUE),
    Trait2 = mean(Pheno.Trait2, na.rm = TRUE)
  )

New_Pheno_mean$ID <- as.character(New_Pheno_mean$ID)

# Reorder based on the rownames of K
New_Pheno_mean <- New_Pheno_mean[match((G1_pop@id), New_Pheno_mean$ID), ]

A_Mat_s <- A_Mat[G1_pop@id, G1_pop@id]

# Estimate inbreeding from A-matrix
Inb <- as.data.frame(diag(A_Mat_s)-1)

# Match the order of the pop object
Inb <- Inb[match(rownames(Inb), G1_pop@id), ]

New_Pheno_mean$Trait1_F <- New_Pheno_mean$Trait1 + New_Pheno_mean$Trait1 * (Inb * -0.35)
New_Pheno_mean$Trait2_F <- New_Pheno_mean$Trait2 + New_Pheno_mean$Trait2 * (Inb * -0.35)

GLOBAL_Phenotypes <- rbind(G0_Pheno_mean, New_Pheno_mean)

GLOBAL_Phenotypes <- GLOBAL_Phenotypes[match(rownames(A_Mat), GLOBAL_Phenotypes$ID), ]


y<-as.matrix(GLOBAL_Phenotypes[,4:5])

#### Model ####
new_model<-mkr(y,A_Mat)

# Save EBVs
EBVs<-as.data.frame(new_model$hat)
EBVs$ID <- rownames(A_Mat)

EBVs <- EBVs[EBVs$ID %in% G1_pop@id,]
EBVs <- EBVs[match(EBVs$ID, G1_pop@id), ]
G1_pop@ebv <- as.matrix(EBVs[,1:2])

New_Pheno_mean <- New_Pheno_mean[match(New_Pheno_mean$ID, G1_pop@id), ]
G1_pop@pheno <- as.matrix(New_Pheno_mean[,4:5])

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

# Phenotypes
Phen.Mean <- meanP(G1_pop)  
VarP <- varP(G1_pop)
Cor_P = cov2cor(VarP)[2, 1]

# Prediction accuracy 
Bv <- bv(G1_pop)
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


#### G2 #####
# Start of blind genomic selection  

# Generate progeny
# Each cross gets a random numner of progenies, with minimum 1
G_pops <- vector("list", nCross)
for (i in 1:nCross) {
  Prog_number <- max(round(rnorm(1, mean = nProg, sd = sdFam)), 1)
  G_pops[[i]] <- randCross(
    pop = Parents, 
    nCrosses = 1, 
    nProgeny = Prog_number, 
    ignoreSexes = TRUE
  )
}

# Combine all families into one pop-object
G2_pop <- do.call(c, G_pops)

# Update pedigree
Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G2_pop@id,])
Pedigree_All <- rbind(Pedigree_All, Pedigree_New)

A_Mat <- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                   dadid = Pedigree_All$father,
                   momid = Pedigree_All$mother)


####################### Calculate MAF Distributoin ###########################
G_new<-pullSnpGeno(G2_pop)
maf_values_new <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                             Gen = "G2")

maf_values_df <- rbind(maf_values_df,maf_values_new) 

M_G<-pullSnpGeno(G2_pop, snpChip = 1)-1
M <- rbind(M, M_G)
G = A.mat(M, impute.method = "mean", min.MAF = 0.01)

dim(G)

#### Blind selection of a subpart of the families
# No new phenotypes, only new marker data

#### Model ####
GLOBAL_Phenotypes <- GLOBAL_Phenotypes[match(rownames(G), GLOBAL_Phenotypes$ID), ]
GLOBAL_Phenotypes$ID <- row.names(G)

y<-as.matrix(GLOBAL_Phenotypes[,4:5])

new_model<-mkr(y,G)

# Save EBVs
EBVs<-as.data.frame(new_model$hat)
EBVs$ID <- rownames(G)

EBVs <- EBVs[EBVs$ID %in% G2_pop@id,]
EBVs <- EBVs[match(EBVs$ID, G2_pop@id), ]
G2_pop@ebv <- as.matrix(EBVs[,1:2])

Acc_blind <- cor(bv(G2_pop),G2_pop@ebv)

# The select function requires that phenotypes exist
G2_pop@pheno <- matrix(0, nrow = nInd(G2_pop), ncol = 2)


##### Select blindly ######
G2_pop <- selectWithinFam(G2_pop,                                           # If the number of progenies are less than nProg/2 
                          nInd = (nProg/2),                                # then that whole family will be selected  
                          use = "ebv",
                          selectTop = TRUE)                                 



################ Assign phenotypic values ##################################  
New_Pheno <- data.frame()

for (i in 1:Ramets) {
  pheno <- setPheno(G2_pop,
                    rep = 1,
                    corE = E,
                    simParam = NULL,
                    onlyPheno = T)
  New_Pheno <- rbind(New_Pheno,
                     data.frame(ID = G2_pop@id, Rep = i, Pheno = pheno))
}

New_Pheno$Rep<-NULL

New_Pheno_mean<- New_Pheno %>%
  group_by(ID) %>%
  summarise(
    Trait1 = mean(Pheno.Trait1, na.rm = TRUE),
    Trait2 = mean(Pheno.Trait2, na.rm = TRUE)
  )

New_Pheno_mean$ID <- as.character(New_Pheno_mean$ID)

# Reorder based on the rownames of K
New_Pheno_mean <- New_Pheno_mean[match((G2_pop@id), New_Pheno_mean$ID), ]

A_Mat_s <- A_Mat[G2_pop@id, G2_pop@id]

# Estimate inbreeding from A-matrix
Inb <- as.data.frame(diag(A_Mat_s)-1)

# Match the order of the pop object
Inb <- Inb[match(rownames(Inb), G2_pop@id), ]

New_Pheno_mean$Trait1_F <- New_Pheno_mean$Trait1 + New_Pheno_mean$Trait1 * (Inb * -0.35)
New_Pheno_mean$Trait2_F <- New_Pheno_mean$Trait2 + New_Pheno_mean$Trait2 * (Inb * -0.35)

# Find indices in df1 where Trait is NA
na_rows <- is.na(GLOBAL_Phenotypes$Trait1)

# Match IDs in df1 to df2 for those rows
matched_indices <- match(GLOBAL_Phenotypes$ID[na_rows], New_Pheno_mean$ID)

GLOBAL_Phenotypes$Trait1_F[na_rows] <- New_Pheno_mean$Trait1_F[matched_indices]
GLOBAL_Phenotypes$Trait2_F[na_rows] <- New_Pheno_mean$Trait2_F[matched_indices]
GLOBAL_Phenotypes$Trait1[na_rows] <- New_Pheno_mean$Trait1[matched_indices]
GLOBAL_Phenotypes$Trait2[na_rows] <- New_Pheno_mean$Trait2[matched_indices]

GLOBAL_Phenotypes <- GLOBAL_Phenotypes[match(rownames(G), GLOBAL_Phenotypes$ID), ]

#### Model 2 ####
y<-as.matrix(GLOBAL_Phenotypes[,4:5])

new_model<-mkr(y,G)

# Save EBVs
EBVs<-as.data.frame(new_model$hat)
EBVs$ID <- rownames(G)

EBVs <- EBVs[EBVs$ID %in% G2_pop@id,]
EBVs <- EBVs[match(EBVs$ID, G2_pop@id), ]
G2_pop@ebv <- as.matrix(EBVs[,1:2])


New_Pheno_mean <- New_Pheno_mean[match(New_Pheno_mean$ID, G2_pop@id), ]
G2_pop@pheno <- as.matrix(New_Pheno_mean[,4:5])

##### Select parents ######
Parents <- selectWithinFam(G2_pop,
                           nInd = 1,
                           use = "ebv",
                           selectTop = TRUE)



# Prediction accuracy 
Bv <- bv(G2_pop)
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



####                                G3                      #####

# Generate progeny
# Each cross gets a random numner of progenies, with minimum 1
G_pops <- vector("list", nCross)
for (i in 1:nCross) {
  Prog_number <- max(round(rnorm(1, mean = nProg, sd = sdFam)), 1)             # Need in intermediate G object?
  G_pops[[i]] <- randCross(
    pop = Parents, 
    nCrosses = 1, 
    nProgeny = Prog_number, 
    ignoreSexes = TRUE
  )
}

# Combine all families into one pop-object
G3_pop <- do.call(c, G_pops)

# Update pedigree
Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G3_pop@id,])
Pedigree_All <- rbind(Pedigree_All, Pedigree_New)

A_Mat <- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                   dadid = Pedigree_All$father,
                   momid = Pedigree_All$mother)


####################### Calculate MAF Distributoin ###########################
G_new<-pullSnpGeno(G3_pop)
maf_values_new <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                             Gen = "G3")

maf_values_df <- rbind(maf_values_df,maf_values_new) 

M_G<-pullSnpGeno(G3_pop, snpChip = 1)-1
M <- rbind(M, M_G)
G = A.mat(M, impute.method = "mean", min.MAF = 0.01)


#### Blind selection of a subpart of the families
# No new phenotypes, only new marker data

GLOBAL_Phenotypes <- GLOBAL_Phenotypes[match(rownames(G), GLOBAL_Phenotypes$ID), ]

GLOBAL_Phenotypes$ID <- colnames(G)

y<-as.matrix(GLOBAL_Phenotypes[,4:5])


#### Model ####
new_model<-mkr(y,G)

# Save EBVs
EBVs<-as.data.frame(new_model$hat)
EBVs$ID <- rownames(G)

EBVs <- EBVs[EBVs$ID %in% G3_pop@id,]
EBVs <- EBVs[match(EBVs$ID, G3_pop@id), ]
G3_pop@ebv <- as.matrix(EBVs[,1:2])

Acc_blind <- cor(bv(G3_pop),G3_pop@ebv)

# The select function requires that phenotypes exist
G3_pop@pheno <- matrix(0, nrow = nInd(G3_pop), ncol = 2)


##### Select blindly ######
G3_pop <- selectWithinFam(G3_pop,                                           # If the number of progenies are less than nProg/2 
                          nInd = (nProg/2),                                # then that whole family will be selected  
                          use = "ebv",
                          selectTop = TRUE)                                 


################ Assign phenotypic values ##################################  
New_Pheno <- data.frame()


for (i in 1:Ramets) {
  pheno <- setPheno(G3_pop,
                    rep = 1,
                    corE = E,
                    simParam = NULL,
                    onlyPheno = T)
  New_Pheno <- rbind(
    New_Pheno,
    data.frame(ID = G3_pop@id, Rep = i, Pheno = pheno)
  )
}

New_Pheno$Rep<-NULL

New_Pheno_mean<- New_Pheno %>%
  group_by(ID) %>%
  summarise(
    Trait1 = mean(Pheno.Trait1, na.rm = TRUE),
    Trait2 = mean(Pheno.Trait2, na.rm = TRUE)
  )

New_Pheno_mean$ID <- as.character(New_Pheno_mean$ID)

# Reorder based on the rownames of K
New_Pheno_mean <- New_Pheno_mean[match((G3_pop@id), New_Pheno_mean$ID), ]

A_Mat_s <- A_Mat[G3_pop@id, G3_pop@id]

# Estimate inbreeding from A-matrix
Inb <- as.data.frame(diag(A_Mat_s)-1)

# Match the order of the pop object
Inb <- Inb[match(rownames(Inb), G3_pop@id), ]

New_Pheno_mean$Trait1_F <- New_Pheno_mean$Trait1 + New_Pheno_mean$Trait1 * (Inb * -0.35)
New_Pheno_mean$Trait2_F <- New_Pheno_mean$Trait2 + New_Pheno_mean$Trait2 * (Inb * -0.35)

# Find indices in df1 where Trait is NA
na_rows <- is.na(GLOBAL_Phenotypes$Trait1)

# Match IDs in df1 to df2 for those rows
matched_indices <- match(GLOBAL_Phenotypes$ID[na_rows], New_Pheno_mean$ID)

GLOBAL_Phenotypes$Trait1_F[na_rows] <- New_Pheno_mean$Trait1_F[matched_indices]
GLOBAL_Phenotypes$Trait2_F[na_rows] <- New_Pheno_mean$Trait2_F[matched_indices]
GLOBAL_Phenotypes$Trait1[na_rows] <- New_Pheno_mean$Trait1[matched_indices]
GLOBAL_Phenotypes$Trait2[na_rows] <- New_Pheno_mean$Trait2[matched_indices]

GLOBAL_Phenotypes <- GLOBAL_Phenotypes[match(rownames(G), GLOBAL_Phenotypes$ID), ]


#GLOBAL_Phenotypes <- GLOBAL_Phenotypes[match(rownames(G), GLOBAL_Phenotypes$ID), ]

#### Model 2 ####
y<-as.matrix(GLOBAL_Phenotypes[,4:5])

G3_pheno <- GLOBAL_Phenotypes[GLOBAL_Phenotypes$ID %in% G3_pop@id,]

GLOBAL_Phenotypes[GLOBAL_Phenotypes$Trait2=="NA",]

new_model<-mkr(y,G)

# Save EBVs
EBVs<-as.data.frame(new_model$hat)
EBVs$ID <- rownames(G)

EBVs <- EBVs[EBVs$ID %in% G3_pop@id,]
EBVs <- EBVs[match(EBVs$ID, G3_pop@id), ]
G3_pop@ebv <- as.matrix(EBVs[,1:2])

New_Pheno_mean <- New_Pheno_mean[match(New_Pheno_mean$ID, G3_pop@id), ]
G3_pop@pheno <- as.matrix(New_Pheno_mean[,4:5])

##### Select parents ######
Parents <- selectWithinFam(G3_pop,
                           nInd = 1,
                           use = "ebv",
                           selectTop = TRUE)

plot(Parents@gv[,2],Parents@ebv[,2])

# Prediction accuracy 
Bv <- bv(G3_pop)
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

