# Inbreeding test# AlphaSim Simulation 

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

# Rate of inbreeding deprssion
ID = -35 # Source:
meanDD = ID/(2*nQtl*PaChr) # 2 time the number of QTLs

initMeanG = c(200,800)  
CV = 0.3 # The standard deviation ~ 30 % of the mean 
initVarG = (CV * initMeanG)^2  # Initial trait genetic variance

initVarD = initVarG*0.25


initVarE = 5*(CV * initMeanG)^2

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
nParents = 20   # Number of parents to be selected
nCross = 30     # Number of families to create from Parents
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
              var          = initVarG,
              meanDD       = meanDD, 
              varDD        = initVarD, 
              gamma        = T,
              shape        = GShape,
              corA         = matrix(c(1, CorA,
                                      CorA, 1), nrow=2),
              corDD        = matrix(c(1, CorA,
                                      CorA, 1), nrow=2)) # Add AxA correlation between two height measurements

# Define heritabilities for both traits (this is REQUIRED before using varE)
SP$setVarE(h2 = c(0.2, 0.2))

# Set residual covariance matrix used by the setPheno function
SP$setVarE(varE = matrix(c(1, 0.8,
                           0.8, 1), nrow=2))

# Create new pop object
G0_pop = newPop(founderPop)

G0_pop <- setPheno(G0_pop,
                   h2 = 0.2)

F1_pop = selectCross(pop = G0_pop,
                     nInd = 5,
                     nCrosses = 25,
                     nProgeny = 20)

F1_pop<-setPheno(pop = F1_pop,
                 h2 = 0.2)

F2_pop = selectCross(pop = F1_pop,
                     nInd = 5,
                     nCrosses = 25,
                     nProgeny = 20)

F2_pop <- setPheno(F2_pop,
               h2 = 0.2)

F3_pop = selectCross(pop = F2_pop,
                     nInd = 5,
                     nCrosses = 25,
                     nProgeny = 20)

F3_pop <- setPheno(F3_pop,
                   h2 = 0.2)

F4_pop = selectCross(pop = F3_pop,
                     nInd = 5,
                     nCrosses = 35,
                     nProgeny = 10)

F4_pop <- setPheno(F4_pop,
                   h2 = 0.2)

F5_pop = selectCross(pop = F4_pop,
                     nInd = 10,
                     nCrosses = 100,
                     nProgeny = 5)

F5_pop <- setPheno(F5_pop,
                   h2 = 0.2)

F6_pop = selectCross(pop = F5_pop,
                     nInd = 100,
                     nCrosses = 100,
                     nProgeny = 20)

F6_pop <- setPheno(F6_pop,
                   h2 = 0.2)

Ped<-as.data.frame(SP$pedigree)

A<- 2*kinship(id = as.numeric(rownames(Ped)),
                           dadid = Ped$father,
                           momid = Ped$mother)

A_s <- A[F6_pop@id,F6_pop@id]

Inb<-(diag(A_s)-1)

Inb_df <- data.frame(I = Inb,
                     H = F6_pop@pheno[,1])

model <- lm(H ~ I, data = Inb_df)
slope <- round(coef(model)[2], 2)
slope/mean(Inb_df$H)

ggplot(Inb_df, aes(x = I, y = H)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Trait vs Inbreeding", x = "Inbreeding (F)", y = "Trait Value")

