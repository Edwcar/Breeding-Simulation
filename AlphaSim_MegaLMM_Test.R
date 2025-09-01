# Use MegaLMM with AlphaSim

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
library(MegaLMM)

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
sdFam = 5
# Mating schemes
MatePlan = "RandCross"

# Phenotyping efforts
ramets = 12 # Ramets per site

# Number of generations to simulate
# Define number of generations

# Genotyping Parameters
# SNP chip for GS
nSNP = 100     # Nr of SNPs per chromosome
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

# Create new pop object
G0_pop = newPop(founderPop)

Acc_df <- data.frame(Trait1 = 0,
                     Trait2 = 0)

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

