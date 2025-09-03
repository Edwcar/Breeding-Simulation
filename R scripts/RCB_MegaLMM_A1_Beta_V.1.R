library(AlphaSimR)
library(pedigree)
library(rrBLUP)
library(ASRgenomics)
library(dplyr)
library(pedigree)
library(ggplot2)
library(tidyr)
library(kinship2)
library(AGHmatrix)
library(MegaLMM)
library(gridExtra)

options(scipen = 999)  # large penalty against scientific notation
 
# Test change

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
SegSite = 100000

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
CV = 0.1 # The standard deviation ~ 30 % of the mean 
initVarG = (CV * initMeanG)^2  # Initial trait genetic variance

# Dominance and Inbreeding Depression
inb_penalty = -0.35 # Source:


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
nFounders = 200   # Number of founders in the population
neFounders = 100   # Effective population size of founders

nGenBurnIn <- 100   # Number of descrete burn-in generations
nCrossBurnIn <- 100 # Number of crosses within each burn-in generation
nProgBurnIn <- 10

# Breeding Parameters    
# Size of the G0 population
nG0 = 500

# Crosses
nParents = 10   # Number of parents to be selected
nCrosses = 20   # Number of families to create from Parents
nProgeny = 10   # Number of progeny per cross
nInd_Fam = 1    # Number of clones to select within each family

# Mating schemes
MatePlan = "RandCross"

# Phenotyping
Reps_G0 = 1
Reps = 10

# Genotyping Parameters
# SNP chip for GS
nSNP = 100     # Nr of SNPs per chromosome
minMAF = 0.01  # Minimum MAF for inclusion


##################### Number of generations to simulate ########################
nGen <- 5   # (after G0)

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


withinFamilyAccuracy <- function(pop, use = c("ebv","pheno")) {
  use <- match.arg(use)
  
  # Build dataframe
  df <- data.frame(
    ID    = pop@id,
    Mum   = pop@mother,
    Dad   = pop@father,
    TBV1  = AlphaSimR::gv(pop)[,1],
    TBV2  = AlphaSimR::gv(pop)[,2]
  )
  
  if(use == "ebv") {
    df$EBV1 <- pop@ebv[,1]
    df$EBV2 <- pop@ebv[,2]
  } else if(use == "pheno") {
    df$PHENO1 <- AlphaSimR::pheno(pop)[,1]
    df$PHENO2 <- AlphaSimR::pheno(pop)[,2]
  }
  
  # Family IDs
  df$Family <- as.integer(interaction(df$Mum, df$Dad, drop = TRUE))
  
  # Compute per-family correlations
  library(dplyr)
  summary_df <- df %>%
    group_by(Family) %>%
    summarize(
      famSize   = n(),
      accuracy1 = if(n() > 1) cor(TBV1, if(use=="ebv") EBV1 else PHENO1) else NA_real_,
      accuracy2 = if(n() > 1) cor(TBV2, if(use=="ebv") EBV2 else PHENO2) else NA_real_,
      .groups = "drop"
    )
  
  # Family-mean accuracy
  mean_acc1 <- mean(summary_df$accuracy1, na.rm = TRUE)
  mean_acc2 <- mean(summary_df$accuracy2, na.rm = TRUE)
  
  return(list(family_summary = summary_df,
              mean_accuracy_trait1 = mean_acc1,
              mean_accuracy_trait2 = mean_acc2))
}
################################################################################
##### Set up Baysian Parameters ####
run_parameters = MegaLMM_control(
  max_NA_groups = 5,  
  scale_Y = F,         
  h2_divisions = 20,   
  burn = 0,     
  thin = 2,
  K = 4)

# Sampling the the MCMC chain

burn_iter = 100
sample_iter = 50

# How many paralell samplings
n_burn_runs = 5
n_sample_runs = 5


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
SP$restrSegSites(overlap = T) # Ensure QTLs and SNPs do not overlap
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


############################ Burn-in Simulation ################################

for (i in 1:nGenBurnIn) {
  G0_pop <- randCross(       # Perform random crossing
    pop = G0_pop,            # The population to cross
    nCrosses = nCrossBurnIn, # Total number of crosses
    nProgeny = nProgBurnIn,  # Number of progeny per cross
    ignoreSexes = TRUE)      # Ignore sexes in crossing
}



# Select G0 founders
G0_pop <- selectInd(G0_pop, nG0, use = "rand")

# Start pedigree for founders
G0_pop@mother <- rep("0", length(G0_pop@id))
G0_pop@father <- rep("0", length(G0_pop@id))

Pedigree_All <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G0_pop@id, ])
Pedigree_All$mother <- 0
Pedigree_All$father <- 0

# Build SNP chip based on Founders
SP$addSnpChip(nSNP, minSnpFreq = minMAF, refPop = G0_pop)

# Container for populations
gen_list <- list()

############################ Generation 0 ############################
cat("=== Generation 0 ===\n")

# Create relationship matrix
A <- 2 * kinship(id = as.numeric(rownames(Pedigree_All)),
                 dadid = Pedigree_All$father,
                 momid = Pedigree_All$mother)

# Genotype Seedlings
M_All<-pullSnpGeno(G0_pop, snpChip = 1)-1 # 

# Set phenotypes for G0
G0_pop <- setPheno(G0_pop, varE = initVarE, corE =  matrix(c(1, CorP,
                                                             CorP, 1), nrow=2),rep = Reps_G0, simParam = NULL, onlyPheno = FALSE)
Pheno <- as.data.frame(G0_pop@pheno)
Pheno$ID <- G0_pop@id
Pheno_All <- Pheno
Pheno_All <- Pheno_All[match(rownames(A), Pheno_All$ID), ]


# Estimate EBVs with MegaLMM
Y <- as.matrix(Pheno_All[, c("Trait1", "Trait2")])
MegaLMM_state <- setup_model_MegaLMM(
  Y, ~1 + (1|ID), data = Pheno_All, relmat = list(ID = A),
  run_parameters = run_parameters, run_ID = "Gen0_MegaLMM")

MegaLMM_state <- run_megalm_analysis(
  MegaLMM_state = MegaLMM_state, priors = priors,
  burn_iter = burn_iter, sample_iter = sample_iter, n_burn_runs = n_burn_runs, n_sample_runs = n_sample_runs,
  drop_cor_threshold = 0.6, verbose = TRUE
)

Lambda_MCMC <- as.data.frame(MegaLMM_state$Posterior$Lambda)

Lambda_long <- Lambda_MCMC %>%
  mutate(iter = row_number()) %>%
  pivot_longer(-iter, names_to = "Factor_Trait", values_to = "Value") %>%
  separate(Factor_Trait, into = c("Factor", "Trait"), sep = "\\.") %>%
  mutate(Factor = as.integer(Factor))

ggplot(Lambda_long, aes(x = iter, y = Value, colour = Trait)) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ Factor, scales = "free_y") +
  labs(x = "Iteration", y = "Lambda Value", title = "MCMC Trace by Factor (G0)") +
  theme_minimal()

EBVs <- as.data.frame(get_posterior_mean(MegaLMM_state, U_R + U_F %*% Lambda))
rownames(EBVs) <- sub("::ID", "", rownames(EBVs))
EBVs <- EBVs[match(G0_pop@id, rownames(EBVs)), ]


G0_pop@ebv <- as.matrix(EBVs[, 1:2])
gen_list[[1]] <- G0_pop

# Select and cross for G1 (no families yet → use phenotypes or EBVs)
G1_pop <- selectCross(G0_pop,
                      nInd = 20, use = "pheno", 
                      nCrosses = nCrosses, nProgeny = nProgeny,
                      selectTop = TRUE)

# Update pedigree
Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G1_pop@id, ])
Pedigree_All <- rbind(Pedigree_All, Pedigree_New)

dim(A)
gen_list[[2]] <- G1_pop

############################ Loop G1 → Gn ############################
for (gen in 2:(nGen+1)) {
  cat("=== Generation", gen-1, "===\n")
  pop <- gen_list[[gen]]
  
  # Create relationship matrix
  A <- 2 * kinship(id = as.numeric(rownames(Pedigree_All)),
                   dadid = Pedigree_All$father,
                   momid = Pedigree_All$mother)
  
  # Genotype Seedlings
  # M<-pullSnpGeno(pop, snpChip = 1)-1
  #  M_All <- rbind(M_All, M)
  #  G = A.mat(M_All, impute.method = "mean", min.MAF = 0.01)
  
  # Set phenotypes
  pop <- setPheno(pop = pop,
                  rep = Reps,
                  varE = initVarE,
                  corE =  matrix(c(1, CorP, 
                                   CorP, 1), nrow=2),
                  simParam = NULL, onlyPheno = FALSE)
  
  Pheno <- as.data.frame(pop@pheno) # Save phenotypes
  Pheno$ID <- pop@id

  # Apply inbreeding penalty
  Inb <- data.frame(Inb = (diag(A)) - 1, ID = rownames(A))
  Inb <- Inb[match(pop@id, Inb$ID), ]
  Pheno$Trait1 <- Pheno$Trait1 + (Pheno$Trait1 * (Inb$Inb * inb_penalty))
  Pheno$Trait2 <- Pheno$Trait2 + (Pheno$Trait2 * (Inb$Inb * inb_penalty))
  pop@pheno <- as.matrix(Pheno[, 1:2])
  
  # Store phenotypes
  Pheno_All <- rbind(Pheno_All, Pheno)
  Pheno_All <- Pheno_All[match(rownames(A), Pheno_All$ID), ]
  
  # Create relationship matrix
  A <- 2 * kinship(id = as.numeric(rownames(Pedigree_All)),
                   dadid = Pedigree_All$father,
                   momid = Pedigree_All$mother)
  
  # Get GEBVs with MegaLMM
  Y <- as.matrix(Pheno_All[, c("Trait1", "Trait2")])
  MegaLMM_state <- setup_model_MegaLMM(
    Y, ~1 + (1|ID), data = Pheno_All, relmat = list(ID = A),
    run_parameters = run_parameters, run_ID = paste0("Gen", gen-1, "_MegaLMM"))
  
  MegaLMM_state <- run_megalm_analysis(
    MegaLMM_state = MegaLMM_state, priors = priors,
    burn_iter = burn_iter, sample_iter = sample_iter, n_burn_runs = n_burn_runs, n_sample_runs = n_sample_runs,
    drop_cor_threshold = 0.6, verbose = TRUE
  )
  
  # Track convergence
  Lambda_MCMC <- as.data.frame(MegaLMM_state$Posterior$Lambda)
  
  Lambda_long <- Lambda_MCMC %>%
    mutate(iter = row_number()) %>%
    pivot_longer(-iter, names_to = "Factor_Trait", values_to = "Value") %>%
    separate(Factor_Trait, into = c("Factor", "Trait"), sep = "\\.") %>%
    mutate(Factor = as.integer(Factor))
  
  Convergance<-ggplot(Lambda_long, aes(x = iter, y = Value, colour = Trait)) +
    geom_line(alpha = 0.6) +
    facet_wrap(~ Factor, scales = "free_y") +
    labs(x = "Iteration", y = "Lambda Value",    
    title = paste("MCMC Trace by Factor – Generation", gen-1) 
) +
    theme_minimal()
  
  print(Convergance)
  
  # Calculate Prediction Accuracy
  EBVs <- as.data.frame(get_posterior_mean(MegaLMM_state, U_R + U_F %*% Lambda))
  rownames(EBVs) <- sub("::ID", "", rownames(EBVs))
  EBVs <- EBVs[match(pop@id, rownames(EBVs)), ]
  pop@ebv <- as.matrix(EBVs[, 1:2])
  gen_list[[gen]] <- pop    
  
  # Select and cross for next generation 
  Parents <- selectWithinFam(pop, nInd = nInd_Fam,
                             use = "ebv",
                             selectTop = TRUE)
  next_pop <- randCross(Parents, nCrosses = nCrosses, nProgeny = nProgeny, ignoreSexes = TRUE)
  
  # Update pedigree
  Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% next_pop@id, ])
  Pedigree_All <- rbind(Pedigree_All, Pedigree_New)
  
  
  
  gen_list[[gen+1]] <- next_pop
}

############################ Output ############################
# Compute genetic variance per trait for all generations

gen_list[[gen]]<-NULL

results <- data.frame()
wf_results <- data.frame()   # per-family WF accuracies

for (gen in seq_along(gen_list)) {
  pop <- gen_list[[gen]]
  
  # Compute within-family accuracies (returns per-family + mean)
  wf_res <- withinFamilyAccuracy(pop, use = "ebv")
  
  for (t in 1:2) {   # loop over traits
    gv <- pop@gv[,t]
    pheno <- pop@pheno[,t]
    
    # Compute prediction accuracy (EBV vs TBV)
    if (!is.null(pop@ebv) && ncol(pop@ebv) >= t) {
      ebv <- pop@ebv[,t]
      valid <- complete.cases(gv, ebv)
      acc <- if (sum(valid) > 1) cor(ebv[valid], gv[valid]) else NA
    } else {
      acc <- NA
    }
    
    # Pull mean WF accuracy
    wf_acc <- if (t == 1) wf_res$mean_accuracy_trait1 else wf_res$mean_accuracy_trait2
    
    # Add to summary dataframe (for mean stats)
    results <- rbind(results, data.frame(
      Generation   = gen - 1,
      Trait        = paste0("Trait", t),
      Mean_GV      = mean(gv, na.rm = TRUE),
      Mean_Pheno   = mean(pheno, na.rm = TRUE),
      Accuracy     = acc,
      WF_Accuracy  = wf_acc
    ))
    
    # Add family-level accuracies
    fam_acc <- wf_res$family_summary
    fam_col <- if (t == 1) "TBV1_accuracy" else "TBV2_accuracy"
    
    if (!all(is.na(fam_acc[[fam_col]]))) {
      wf_results <- rbind(wf_results, data.frame(
        Generation = gen - 1,
        Trait      = paste0("Trait", t),
        Family     = fam_acc$Family,
        FamSize    = fam_acc$famSize,
        WF_Acc     = fam_acc[[fam_col]]
      ))
    }
  }  # closes trait loop
}    # closes generation loop

# ---------------- Collect Variances ---------------- #
gen_var_list <- lapply(seq_along(gen_list), function(gen) {
  pop <- gen_list[[gen]]
  data.frame(
    Generation   = gen - 1,
    Trait1_GV_var = var(pop@gv[,1], na.rm = TRUE),
    Trait2_GV_var = var(pop@gv[,2], na.rm = TRUE),
    Trait1_P_var  = var(pop@pheno[,1], na.rm = TRUE),
    Trait2_P_var  = var(pop@pheno[,2], na.rm = TRUE)
  )
})
gen_var_df <- do.call(rbind, gen_var_list)

# ---------------- Merge results per trait ---------------- #

# Trait1
results_trait1 <- results %>% 
  filter(Trait == "Trait1") %>%
  left_join(gen_var_df %>% select(Generation, Trait1_GV_var, Trait1_P_var), by = "Generation") %>%
  rename(Genetic_Var = Trait1_GV_var,
         Phenotypic_Var = Trait1_P_var)
wf_trait1 <- wf_results %>% filter(Trait == "Trait1")

# Trait2
results_trait2 <- results %>% 
  filter(Trait == "Trait2") %>%
  left_join(gen_var_df %>% select(Generation, Trait2_GV_var, Trait2_P_var), by = "Generation") %>%
  rename(Genetic_Var = Trait2_GV_var,
         Phenotypic_Var = Trait2_P_var)
wf_trait2 <- wf_results %>% filter(Trait == "Trait2")

# ---------------- Plotting ---------------- #


# ---------------- Trait1 Plots ---------------- #
# GV vs Phenotype
p1 <- ggplot(results_trait1, aes(x = Generation)) +
  geom_line(aes(y = Mean_GV, color = "GV"), size = 1.2) +
  geom_line(aes(y = Mean_Pheno, color = "Phenotype"), size = 1.2) +
  geom_point(aes(y = Mean_GV, color = "GV"), size = 2) +
  geom_point(aes(y = Mean_Pheno, color = "Phenotype"), size = 2) +
  labs(y = "Value", x = "Generation", title = "Trait1: GV vs Phenotype") +
  scale_color_manual(values = c("GV" = "turquoise", "Phenotype" = "darkblue")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Variances
p2 <- ggplot(results_trait1, aes(x = Generation)) +
  geom_line(aes(y = Genetic_Var, color = "Genetic"), size = 1.2) +
  geom_point(aes(y = Genetic_Var, color = "Genetic"), size = 2) +
  #geom_line(aes(y = Phenotypic_Var, color = "Phenotypic"), size = 1.2) +
  #geom_point(aes(y = Phenotypic_Var, color = "Phenotypic"), size = 2) +
  labs(y = "Variance", x = "Generation", title = "Trait1: Genetic vs Phenotypic Variance") +
  scale_color_manual(values = c("Genetic" = "darkgreen", "Phenotypic" = "limegreen")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Accuracy with boxplot (check if wf_trait1 has data)
p3 <- ggplot() +
  {if (nrow(wf_trait1) > 0) geom_boxplot(
    data = wf_trait1,
    aes(x = factor(Generation), y = WF_Acc, group = Generation),
    fill = "magenta", alpha = 0.3, outlier.shape = NA
  )} +
  geom_line(data = results_trait1, 
            aes(x = factor(Generation), y = WF_Accuracy, color = "WF", group = 1), 
            size = 1.2, linetype = "dashed") +
  geom_point(data = results_trait1, 
             aes(x = factor(Generation), y = WF_Accuracy, color = "WF"), 
             size = 2) +
  geom_line(data = results_trait1, 
            aes(x = factor(Generation), y = Accuracy, color = "Global", group = 1), 
            size = 1.2) +
  geom_point(data = results_trait1, 
             aes(x = factor(Generation), y = Accuracy, color = "Global"), 
             size = 2) +
  labs(y = "Accuracy", x = "Generation", title = "Trait1: Accuracy (EBV vs WF)") +
  ylim(0,1) +
  scale_color_manual(values = c("Global" = "purple", "WF" = "magenta")) +
  theme_minimal() +
  theme(legend.position = "bottom")
# ---------------- Trait2 Plots ---------------- #
# GV vs Phenotype
p4 <- ggplot(results_trait2, aes(x = Generation)) +
  geom_line(aes(y = Mean_GV, color = "GV"), size = 1.2) +
  geom_line(aes(y = Mean_Pheno, color = "Phenotype"), size = 1.2) +
  geom_point(aes(y = Mean_GV, color = "GV"), size = 2) +
  geom_point(aes(y = Mean_Pheno, color = "Phenotype"), size = 2) +
  labs(y = "Value", x = "Generation", title = "Trait2: GV vs Phenotype") +
  scale_color_manual(values = c("GV" = "turquoise", "Phenotype" = "darkblue")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Variances
p5 <- ggplot(results_trait2, aes(x = Generation)) +
  geom_line(aes(y = Genetic_Var, color = "Genetic"), size = 1.2) +
  geom_point(aes(y = Genetic_Var, color = "Genetic"), size = 2) +
  #geom_line(aes(y = Phenotypic_Var, color = "Phenotypic"), size = 1.2) +
  #geom_point(aes(y = Phenotypic_Var, color = "Phenotypic"), size = 2) +
  labs(y = "Variance", x = "Generation", title = "Trait2: Genetic vs Phenotypic Variance") +
  scale_color_manual(values = c("Genetic" = "darkgreen", "Phenotypic" = "limegreen")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Accuracy with boxplot (check if wf_trait2 has data)
p6 <- ggplot() +
  {if (nrow(wf_trait2) > 0) geom_boxplot(
    data = wf_trait2,
    aes(x = factor(Generation), y = WF_Acc, group = Generation),
    fill = "magenta", alpha = 0.3, outlier.shape = NA
  )} +
  geom_line(data = results_trait2, 
            aes(x = factor(Generation), y = WF_Accuracy, color = "WF", group = 1), 
            size = 1.2, linetype = "dashed") +
  geom_point(data = results_trait2, 
             aes(x = factor(Generation), y = WF_Accuracy, color = "WF"), 
             size = 2) +
  geom_line(data = results_trait2, 
            aes(x = factor(Generation), y = Accuracy, color = "Global", group = 1), 
            size = 1.2) +
  geom_point(data = results_trait2, 
             aes(x = factor(Generation), y = Accuracy, color = "Global"), 
             size = 2) +
  labs(y = "Accuracy", x = "Generation", title = "Trait2: Accuracy (EBV vs WF)") +
  ylim(0,1) +
  scale_color_manual(values = c("Global" = "purple", "WF" = "magenta")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# ---------------- Arrange all plots ---------------- #
grid.arrange(p1, p2, p3,
             p4, p5, p6,
             ncol = 3)
