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
library(MegaLMM)

# Set initial Simullation values #############################

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

run_megalm_analysis <- function(MegaLMM_state, priors, 
                                burn_iter = 1, sample_iter = 1, 
                                n_burn_runs = 1, n_sample_runs = 1,
                                drop_cor_threshold = 0.6,
                                verbose = TRUE) {
  
  # Setup missing data and priors
  if(verbose) cat("Setting up missing data map and priors...\n")
  maps <- make_Missing_data_map(MegaLMM_state)
  MegaLMM_state <- set_Missing_data_map(MegaLMM_state, maps$Missing_data_map)
  MegaLMM_state <- set_priors_MegaLMM(MegaLMM_state, priors)
  
  # Initialize variables and check memory
  if(verbose) cat("Initializing variables and checking memory...\n")
  MegaLMM_state <- initialize_variables_MegaLMM(MegaLMM_state)
  estimate_memory_initialization_MegaLMM(MegaLMM_state)
  MegaLMM_state <- initialize_MegaLMM(MegaLMM_state)
  MegaLMM_state <- clear_Posterior(MegaLMM_state)
  
  # Burn-in phase
  if(verbose) cat("Starting burn-in phase...\n")
  for(i in 1:n_burn_runs) {
    if(verbose) print(sprintf('Burn-in run %d', i))
    
    # Reorder factors to improve mixing
    MegaLMM_state <- reorder_factors(MegaLMM_state, drop_cor_threshold = drop_cor_threshold)
    MegaLMM_state <- clear_Posterior(MegaLMM_state)
    
    # Sample
    MegaLMM_state <- sample_MegaLMM(MegaLMM_state, burn_iter)
    
    # Optional: Check convergence with diagnostic plots
    # traceplot_array(MegaLMM_state$Posterior$Lambda, 
    #                 name = file.path(MegaLMM_state$run_ID, paste0('Lambda_burnin_', i, '.pdf')))
  }
  
  # Clear burn-in samples - keep current parameter values
  MegaLMM_state <- clear_Posterior(MegaLMM_state)
  
  # Sampling phase 
  if(verbose) cat("Starting sampling phase...\n")
  for(i in 1:n_sample_runs) {
    if(verbose) print(sprintf('Sampling run %d', i))
    MegaLMM_state <- sample_MegaLMM(MegaLMM_state, sample_iter)
    MegaLMM_state <- save_posterior_chunk(MegaLMM_state)
  }
  
  # Reload all posterior samples for analysis
  if(verbose) cat("Reloading posterior samples...\n")
  MegaLMM_state$Posterior <- reload_Posterior(MegaLMM_state)
  
  return(MegaLMM_state)
}

################################################################################
##### Set up Baysian Parameters ####
run_parameters = MegaLMM_control(
  max_NA_groups = 5,  
  scale_Y = T,         
  h2_divisions = 20,   
  burn = 0,     
  thin = 2,
  K = 5)

priors = MegaLMM_priors(
  tot_Y_var = list(V = 0.5 , nu = 3),  
  tot_F_var = list(V = 0.9, nu = 20), 
  
  Lambda_prior = list(                
    sampler = sample_Lambda_prec_ARD,  
    Lambda_df = 3,                    
    delta_1 = list(shape = 2, rate = 1), 
    delta_2 = list(shape = 3, rate = 1),
    
    Lambda_beta_var_shape = 3,
    Lambda_beta_var_rate = 1,
    delta_iterations_factor = 100
  ),
  
  h2_priors_resids_fun = function(h2s,n) 1, 
  h2_priors_factors_fun = function(h2s,n) 1)


#### Set up Baysian Sampling ####



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
SP$setVarE(varE = matrix(c(1, CorP,
                           CorP, 1), nrow=2))

# Create new pop object
G0_pop = newPop(founderPop)


################################## Burn-in #####################################

for (i in 1:nGenBurnIn) {
  G0_pop <- randCross(       # Perform random crossing
    pop = G0_pop,            # The population to cross
    nCrosses = nCrossBurnIn, # Total number of crosses
    nProgeny = nProgBurnIn,  # Number of progeny per cross
    ignoreSexes = TRUE)      # Ignore sexes in crossing
}



###################### G0 Founder generation ######################################

# Select clones for first real G0 generation
G0_pop <- selectInd(G0_pop, nG0, use = "rand")



# Save pedigree
G0_pop@mother <- rep("0",length(G0_pop@id))
G0_pop@father <- rep("0",length(G0_pop@id))

Pedigree_All <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G0_pop@id,])
Pedigree_All$mother<-0
Pedigree_All$father<-0

Kinship_matrix<- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                           dadid = Pedigree_All$father,
                           momid = Pedigree_All$mother)


#### Add SNP chips ####

# Randomly assigns eligible SNPs to SNP chip
# Use the latest population as reference
SP$addSnpChip(nSNP, minSnpFreq = minMAF,name = "GS", refPop = G0_pop)


# Add first phenotypes
G0_Pheno <- data.frame()

for (i in 1:1) {
  pheno <- setPheno(G0_pop,
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


Y <- as.matrix(GLOBAL_Phenotypes[, c("Pheno.Trait1", "Pheno.Trait2")])

MegaLMM_state = setup_model_MegaLMM(
  Y,  # your n x p data matrix
  ~1 + (1|ID),  # model formula for fixed and random effects
  data = G0_Pheno_mean,  # data frame with model variables
  relmat = list(ID = Kinship_matrix),  # relationship matrices for random effects
  run_parameters = run_parameters,
  run_ID = 'Sim_MegaLMM_output'  # name for output directory
)

MegaLMM_state <- run_megalm_analysis(
  MegaLMM_state = MegaLMM_state,
  priors = priors,
  burn_iter = 50,        # iterations per burn-in run
  sample_iter = 50,      # iterations per sampling run  
  n_burn_runs = 5,       # number of burn-in runs
  n_sample_runs = 5,     # number of sampling runs
  drop_cor_threshold = 0.6,
  verbose = TRUE)

GEBV_MegaLMM <- as.data.frame(get_posterior_mean(MegaLMM_state, U_R + U_F %*% Lambda))
rownames(GEBV_MegaLMM) <- sub("::ID", "", rownames(GEBV_MegaLMM))
EBVs <- GEBV_MegaLMM[match(G0_pop@id, rownames(GEBV_MegaLMM)), ]
G0_pop@ebv <- as.matrix(EBVs[,1:2])

Pheno <- G0_Pheno_mean[match(G0_pop@id, G0_Pheno_mean$ID), ]
G0_pop@pheno <- as.matrix(Pheno[,2:3])

##### Select parents ######
G1_pop <- selectCross(G0_pop,
                           nInd = 20,
                           use = "pheno",
                            nCrosses = 10,
                            nProgeny = 10,
                           selectTop = TRUE)

Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G1_pop@id,])
Pedigree_All <- rbind(Pedigree_All, Pedigree_New)

A_Mat <- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                   dadid = Pedigree_All$father,
                   momid = Pedigree_All$mother)

New_Pheno <- data.frame()

for (i in 1:10) {
  pheno <- setPheno(G1_pop,
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

New_Pheno_mean <- New_Pheno_mean[match(rownames(A_Mat), New_Pheno_mean$ID), ]

Y <- as.matrix(New_Pheno_mean[, c("Trait1", "Trait2")])

MegaLMM_state = setup_model_MegaLMM(
  Y,  # your n x p data matrix
  ~1 + (1|ID),  # model formula for fixed and random effects
  data = New_Pheno_mean,  # data frame with model variables
  relmat = list(ID = A_Mat),  # relationship matrices for random effects
  run_parameters = run_parameters,
  run_ID = 'Sim_MegaLMM_output'  # name for output directory
)

MegaLMM_state <- run_megalm_analysis(
  MegaLMM_state = MegaLMM_state,
  priors = priors,
  burn_iter = 50,        # iterations per burn-in run
  sample_iter = 50,      # iterations per sampling run  
  n_burn_runs = 5,       # number of burn-in runs
  n_sample_runs = 5,     # number of sampling runs
  drop_cor_threshold = 0.6,
  verbose = TRUE)

GEBV_MegaLMM <- as.data.frame(get_posterior_mean(MegaLMM_state, U_R + U_F %*% Lambda))
rownames(GEBV_MegaLMM) <- sub("::ID", "", rownames(GEBV_MegaLMM))
EBVs <- GEBV_MegaLMM[match(G1_pop@id, rownames(GEBV_MegaLMM)), ]
G1_pop@ebv <- as.matrix(EBVs[,1:2])

Pheno <- New_Pheno_mean[match(G1_pop@id, New_Pheno_mean$ID), ]
G1_pop@pheno <- as.matrix(Pheno[,2:3])

Parents <- selectWithinFam(G1_pop,
                           nInd = 1,
                           use = "ebv",
                           selectTop = TRUE)

G2_pop <- randCross(pop = Parents, 
                    nCrosses = 10, 
                    nProgeny = 10, 
                    ignoreSexes = TRUE
)

Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G2_pop@id,])
Pedigree_All <- rbind(Pedigree_All, Pedigree_New)

A_Mat <- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                   dadid = Pedigree_All$father,
                   momid = Pedigree_All$mother)

New_Pheno <- data.frame()

for (i in 1:10) {
  pheno <- setPheno(G2_pop,
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

New_Pheno_mean <- New_Pheno_mean[match(rownames(A_Mat), New_Pheno_mean$ID), ]

Y <- as.matrix(New_Pheno_mean[, c("Trait1", "Trait2")])

MegaLMM_state = setup_model_MegaLMM(
  Y,  # your n x p data matrix
  ~1 + (1|ID),  # model formula for fixed and random effects
  data = New_Pheno_mean,  # data frame with model variables
  relmat = list(ID = A_Mat),  # relationship matrices for random effects
  run_parameters = run_parameters,
  run_ID = 'Sim_MegaLMM_output'  # name for output directory
)

MegaLMM_state <- run_megalm_analysis(
  MegaLMM_state = MegaLMM_state,
  priors = priors,
  burn_iter = 50,        # iterations per burn-in run
  sample_iter = 50,      # iterations per sampling run  
  n_burn_runs = 5,       # number of burn-in runs
  n_sample_runs = 5,     # number of sampling runs
  drop_cor_threshold = 0.6,
  verbose = TRUE)

GEBV_MegaLMM <- as.data.frame(get_posterior_mean(MegaLMM_state, U_R + U_F %*% Lambda))
rownames(GEBV_MegaLMM) <- sub("::ID", "", rownames(GEBV_MegaLMM))
EBVs <- GEBV_MegaLMM[match(G2_pop@id, rownames(GEBV_MegaLMM)), ]
G2_pop@ebv <- as.matrix(EBVs[,1:2])

Pheno <- New_Pheno_mean[match(G2_pop@id, New_Pheno_mean$ID), ]
G2_pop@pheno <- as.matrix(Pheno[,2:3])

Parents <- selectWithinFam(G2_pop,
                           nInd = 1,
                           use = "ebv",
                           selectTop = TRUE)

