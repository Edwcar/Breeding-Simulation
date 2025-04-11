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
    # Height at 10 years
    initMeanG = c(200,800)     # Initial trait mean genetic value
    CV = 0.3 # The standard deviation ~ 30 % of the mean 
    initVarG = (CV * initMeanG)^2  # Initial trait genetic variance
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
    nFounders = 100   # Number of founders in the population
    neFounders = 10   # Effective population size of founders

    nGenBurnIn <- 100  # Number of descrete burn-in generations
    nCrossBurnIn <- 100 # Number of crosses within each burn-in generation
    nProgBurnIn <- 10

# Breeding Parameters    
  # Size of the G0 population
    nG0 = 100

  # Crosses
    nParents = 10   # Number of parents to be selected
    nCross = 10     # Number of families to create from Parents
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
    num_generations <- 10  
    
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
SP$addTraitA(nQtlPerChr   = nQtl,
             mean         = initMeanG,
             var          = initVarG,
             gamma        = T,
             shape        = GShape,
             corA         = matrix(c(1, CorA,
                                     CorA, 1), nrow=2)) # Add AxA correlation between two height measurements

# Define heritabilities for both traits (this is REQUIRED before using varE)
SP$setVarE(h2 = c(0.2, 0.2))

# Set residual covariance matrix used by the setPheno function
#SP$setVarE(varE = matrix(c(1, 0.8,
#                           0.8, 1), nrow=2))
# Create new pop object
G0_pop = newPop(founderPop)


Acc_df <- data.frame(Trait1 = 0,
                     Trait2 = 0)

################################## Burn-in #####################################
for (i in 1:nGenBurnIn) {
  # Perform random crossing
  G0_pop <- randCross(
    pop = G0_pop,            # The population to cross
    nCrosses = nCrossBurnIn, # Total number of crosses
    nProgeny = nProgBurnIn,  # Number of progeny per cross
    ignoreSexes = TRUE)      # Ignore sexes in crossing
}

###################### Founder generation ######################################
                    
# Select clones for first real G0 generation
G0_pop <- selectInd(G0_pop, nG0, use = "rand")

#### Add SNP chips ####

# Randomly assigns eligible SNPs to SNP chip
# Use the latest population as reference
SP$addSnpChip(nSNP, minSnpFreq = minMAF,name = "GS", refPop = G0_pop)
SP$addSnpChip(nSNP2, minSnpFreq = minMAF,name = "Inb", refPop = G0_pop)

# Genotype with imaginary high-density SNP chip to estimate an inbreeding baseline
G_G0<-pullSnpGeno(G0_pop, snpChip = "Inb")-1
  G<-A.mat(G_G0, min.MAF = 0.01)
  F_0<-mean(diag(G))

  F_new <-(diag(G))-F_0

  Inbreeding<-data.frame(F = F_new,
                       ID = as.numeric(G0_pop@id))

    # Save mean self-kinship 
   F_df <- data.frame(Inbreeding = mean(diag(G)))


  ###### Calculate MAF distribution ####
    G_G0<-pullSnpGeno(G0_pop, snpChip = "Inb")
    maf_values_df <- data.frame(Frequency = apply(G_G0, 2, calculate_maf),
                                Gen = 0)
  #maf_values0.05 <-maf_values[maf_values$Frequency > 0.01, ]
  
# Add first phenotypes
    G0_Pheno <- data.frame(genotype_id = (G0_pop@id))  # Start with IDs

    for (i in 1:n_reps) {
      pheno <- setPheno(G0_pop,
                      varE = initVarE,
                      corE = matrix(c(1, CorP,
                                      CorP, 1), nrow = 2),
                      simParam = NULL)
  
      G0_Pheno[[paste0("trait1_rep", i)]] <- pheno@pheno[, 1]
      G0_Pheno[[paste0("trait2_rep", i)]] <- pheno@pheno[, 2]
    }

  # Add missing data
  num_values <- prod(dim(G0_Pheno[, -1]))  # Exclude genotype_id
  num_na <- round(0.1 * num_values)

  rows <- sample(nrow(G0_Pheno), num_na, replace = TRUE)
  cols <- sample(2:ncol(G0_Pheno), num_na, replace = TRUE)
  
  for (i in seq_along(rows)) {
    G0_Pheno[rows[i], cols[i]] <- NA
  }

  # Add missing data
  num_values <- prod(dim(G0_Pheno[, -1]))  # Exclude genotype_id
  num_na <- round(0.1 * num_values)

  rows <- sample(nrow(G0_Pheno), num_na, replace = TRUE)
  cols <- sample(2:ncol(G0_Pheno), num_na, replace = TRUE)

  for (i in seq_along(rows)) {
    G0_Pheno[rows[i], cols[i]] <- NA
  }

  # Trait 1 columns (those starting with "trait1_")
  trait1_cols <- grep("^trait1_", names(G0_Pheno), value = TRUE)
  trait1_df <- G0_Pheno[, c("genotype_id", trait1_cols)]

  # Trait 2 columns (those starting with "trait2_")
  trait2_cols <- grep("^trait2_", names(G0_Pheno), value = TRUE)
  trait2_df <- G0_Pheno[, c("genotype_id", trait2_cols)]

  trait1_df$Mean <- rowMeans(trait1_df[ , -1], na.rm = T)
  trait2_df$Mean <- rowMeans(trait2_df[ , -1], na.rm = T)

  New_Pheno_long_t1<-trait1_df %>%
    pivot_longer(cols = c(-genotype_id, -Mean), names_to = "Rep", values_to = "Pheno")
  New_Pheno_long_t1$Rep <- gsub("^V[0-9]+\\.|_rep[0-9]+$", "", New_Pheno_long_t1$Rep)
  New_Pheno_long_t1$ID <- as.numeric(New_Pheno_long_t1$genotype_id) 
  New_Pheno_long_t1 <- na.omit(New_Pheno_long_t1)
  New_Pheno_long_t1$`G0_pop@id`<-NULL
  New_Pheno_long_t1$genotype_id<-NULL

  New_Pheno_long_t2<-trait2_df %>%
    pivot_longer(cols = c(-genotype_id, -Mean), names_to = "Rep", values_to = "Pheno")
  New_Pheno_long_t2$Rep <- gsub("^V[0-9]+\\.|_rep[0-9]+$", "", New_Pheno_long_t2$Rep)
  New_Pheno_long_t2$ID <- as.numeric(New_Pheno_long_t2$genotype_id)
  New_Pheno_long_t2 <- na.omit(New_Pheno_long_t2)
  New_Pheno_long_t2$`G0_pop@id`<-NULL
  New_Pheno_long_t2$genotype_id<-NULL
  

  
  
  # Add linear penalty for inbreeding depression only if F > 0
  New_Pheno_long_t1 <- left_join(New_Pheno_long_t1, Inbreeding, by = "ID")
  New_Pheno_long_t2 <- left_join(New_Pheno_long_t2, Inbreeding, by = "ID")
  
  # Apply the penalty only if F > 0
  New_Pheno_long_t1$Pheno_Inb <- New_Pheno_long_t1$Pheno + 
    abs(New_Pheno_long_t1$Pheno) * ifelse(New_Pheno_long_t1$F > 0, (InbDepr * New_Pheno_long_t1$F), 0)
  
  # Apply the penalty only if F > 0
  New_Pheno_long_t2$Pheno_Inb <- New_Pheno_long_t2$Pheno + 
    abs(New_Pheno_long_t2$Pheno) * ifelse(New_Pheno_long_t2$F > 0, (InbDepr * New_Pheno_long_t2$F), 0)
  

  All_Pheno_long <- rbind(New_Pheno_long_t1, New_Pheno_long_t2)
  

  

  # Summerize to mean values to get same number of values per trait
  Trait1_mean <- New_Pheno_long_t1 %>%
    group_by(ID) %>%
    summarise(mean_pheno = mean(Pheno_Inb, na.rm = TRUE))
  
  Trait2_mean <- New_Pheno_long_t2 %>%
    group_by(ID) %>%
    summarise(mean_pheno = mean(Pheno_Inb, na.rm = TRUE))
  
  
    G0_Pheno <- data.frame(Trait1 = Trait1_mean$mean_pheno,
                           Trait2 = Trait2_mean$mean_pheno)

  G0_pop@pheno <- as.matrix(G0_Pheno)

  # Save G0 parameters
  MeanG <- meanG(G0_pop)  
  VarA <- varG(G0_pop)
  Cor_A = cov2cor(VarA)[2, 1]
  Additive_df <- data.frame(Trait1_Mean = MeanG[1],
                            Trait1_Var = VarA[1,1],
                            Trait2_Mean = MeanG[2],
                            Trait2_Var = VarA[2,2],
                            Covariance = VarA[1,2],
                            Correlation = Cor_A)
  
  
  
  Phen.Mean <- meanP(G0_pop)  
  VarP <- varP(G0_pop)
  Cor_P = cov2cor(VarP)[2, 1]
  Pheno_df <- data.frame(Trait1_Mean = Phen.Mean[1],
                         Trait1_Var = VarP[1,1],
                         Trait2_Mean = Phen.Mean[2],
                         Trait2_Var = VarP[2,2],
                         Covariance = VarP[1,2],
                         Correlation = Cor_P) 
  
  Acc_df <- data.frame(Trait1 = cor(G0_pop@gv[,1], G0_pop@pheno[,1]),
                       Trait2 = cor(G0_pop@gv[,2], G0_pop@pheno[,2]))
  
  Bias_df <- data.frame(Trait1 = 1,
                        Trait2 = 1)
  

  
####################### Large Loop Start #######################################

# Save Parameters during looping
# Initialize populations
  Pheno_mean_All <- data.frame(ID = G0_Pheno$`G0_pop@id`, Trait = G0_Pheno$Mean)
  Pedigree_All <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G0_pop@id,])
  Pedigree_All$mother <- 0 # Reset pedigree for Plus trees
  Pedigree_All$father <- 0 # Reset pedigree for Plus trees
  Pedigree<-as.data.frame(SP$pedigree)
  
  A <- data.frame(ID = as.numeric(rownames(Pedigree_All)),
                  Sire = Pedigree_All$father,
                  Dam = Pedigree_All$mother)
  
  model <- remlf90(
    fixed = Pheno_Inb ~ Rep,  # trait-specific means
    genetic = list(
      model = 'add_animal',
      pedigree = A,
      id = All_Pheno_long$ID
    ),
    data = All_Pheno_long,
  )
  
  BV.full <- ranef(model)$genetic
  BV <- model.matrix(model)$genetic %*% BV.full
  All_Pheno_long$EBV <- BV@x
  EBV_df <- All_Pheno_long %>% distinct(ID, EBV)
  All_Pheno_long$EBV <- NULL
  
  EBV_new <- EBV_df[EBV_df$ID %in% G0_pop@id,]
  G0_pop@ebv <- as.matrix(EBV_new$EBV)
  
  
  current_pop <- G0_pop
  # Store populations for downstream analysis
  pop_list <- list(G0 = G0_pop)


######################### Loop through generations  ############################
for (gen in 1:num_generations) {
  cat("Generation:", gen, "\n")
  
  # Select best individuals based on phenotype or EBV
  if (gen == 1) {
    
    # Calculate the Smith Hazel Index
    #b = smithHazel(econWt = c(1,1), 
    #               varG = varG(current_pop),
    #               varP = varP(current_pop))
    
       # Select the first parent based on phenotype
    Parents <- selectInd(current_pop,
                         nInd = nParents,
                         use = "ebv",
                         selectTop = TRUE)
    
    # Calculate selection intensity
    mean_EBV_new <- mean(Parents@ebv)
    mean_EBV_all <- mean(G0_pop@ebv)
    sd_EBV_all <- sd(G0_pop@ebv)
    intensity_EBV_df <- data.frame(Intensity = ((mean_EBV_new - mean_EBV_all) / sd_EBV_all),
                                   Gen = gen)
    
    mean_TBV_new <- mean(Parents@gv)
    mean_TBV_all <- mean(G0_pop@gv)
    sd_TBV_all <- sd(G0_pop@gv)
    intensity_TBV_df <- data.frame(Intensity = ((mean_TBV_new - mean_TBV_all) / sd_TBV_all),
                                   Gen = gen)
    
  } else {

    
    Parents <- selectWithinFam(current_pop,
                         nInd = 1,
                         use = "ebv",
                         selectTop = TRUE
                         )
  }
  
  # Generate progeny
  # Each cross gets a random numner of progenies, with minimum 1
  G_pops <- vector("list", nParents)
  for (i in 1:nParents) {
    Prog_number <- max(round(rnorm(1, mean = nProg, sd = 10)), 1)
    G_pops[[i]] <- randCross(
      pop = Parents, 
      nCrosses = 1, 
      nProgeny = Prog_number, 
      ignoreSexes = TRUE
    )
  }
  
  # Combine all families into one pop-object
  new_pop <- do.call(c, G_pops)
  
  # Save population object
  pop_list[[paste0("G", gen)]] <- new_pop
  
  # Update pedigree
  Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% new_pop@id,])
  Pedigree_All <- rbind(Pedigree_All, Pedigree_New)
  
  if (gen > 5) {  # Ensure at least 5 previous generations exist
    gen_ids_to_keep <- unlist(lapply((gen-5):gen, function(g) pop_list[[paste0("G", g)]]@id))
  }

  if (gen < 6) {
    A <- data.frame(ID = rownames(Pedigree_All),
                    Sire = Pedigree_All$father,
                    Dam = Pedigree_All$mother)
  
  } else {
    A <- Pedigree_All[rownames(Pedigree_All) %in% gen_ids_to_keep, ]
  
    A <- data.frame(ID = rownames(A),
                    sire = A$father,
                    dam = A$mother)

  }
  
  A$ID <-as.numeric(A$ID)
  
  
  ################ To do ######################
  # Add mechanism for inbreeding depression
  # Calculate the A matrix first
  
  #A_matrix<-2*kinship(id = A$ID, dadid = A$sire, momid = A$dam)
  
  #Inbreeding<-data.frame(F = (diag(A_matrix)-1),
   #            ID = A$ID)
  

  G_new<-pullSnpGeno(new_pop, snpChip = "Inb")-1
  G<-A.mat(G_new, min.MAF = 0.01)
  F_new <-(diag(G))-F_0
  
  F_df <- rbind(F_df,mean(diag(G)))
  
  Inbreeding<-data.frame(F = F_new,
              ID = as.numeric(new_pop@id))
  
  
  ####################### Calculate MAF Distributoin ###########################
  G_new<-pullSnpGeno(new_pop, snpChip = "Inb")
  maf_values_new <- data.frame(Frequency = apply(G_new, 2, calculate_maf),
                               Gen = gen)

  maf_values_df <- rbind(maf_values_df,maf_values_new) 
  
  ################ Assign phenotypic values ####################################
  new_Pheno <- as.data.frame(new_pop@id)
  for (i in 1:n_reps) {
    new_Pheno[[paste0("V", i)]] <- setPheno(new_pop,
                                            varE = initVarE,
                                            corE = matrix(c(1, CorP,
                                                            CorP, 1), nrow = 2),
                                            simParam = NULL,
                                            onlyPheno = T)
  } 
 
  # Add missing data
  num_values <- prod(dim(new_Pheno[, -1]))  # Exclude genotype_id
  num_na <- round(0.1 * num_values)
  
  rows <- sample(nrow(new_Pheno), num_na, replace = TRUE)
  cols <- sample(2:ncol(new_Pheno), num_na, replace = TRUE)
  
  for (i in seq_along(rows)) {
    new_Pheno[rows[i], cols[i]] <- NA
  }
  
  # Add missing data
  num_values <- prod(dim(new_Pheno[, -1]))  # Exclude genotype_id
  num_na <- round(0.1 * num_values)
  
  rows <- sample(nrow(new_Pheno), num_na, replace = TRUE)
  cols <- sample(2:ncol(new_Pheno), num_na, replace = TRUE)
  
  for (i in seq_along(rows)) {
    new_Pheno[rows[i], cols[i]] <- NA
  }

  # Flatten the data.frame
  new_Pheno <- do.call(cbind, lapply(new_Pheno, as.data.frame))

  # Subset the phenotypes based on trait
      trait1_cols <- grep("\\.Trait1$", names(new_Pheno), value = TRUE)
      trait1_df <- new_Pheno[, c("X[[i]]", trait1_cols)]

      trait2_cols <- grep("\\.Trait2$", names(new_Pheno), value = TRUE)
      trait2_df <- new_Pheno[, c("X[[i]]", trait2_cols)]


  # Calculate mean values per clone per trait
     trait1_df$Mean <- rowMeans(trait1_df[ , -1], na.rm = T)
     trait2_df$Mean <- rowMeans(trait2_df[ , -1], na.rm = T)
  
  # Save individual measurements for both traits for downstream analysis   
    All_Pheno_long_t1<-trait1_df %>%
      pivot_longer(cols = c(-`X[[i]]`, -Mean), names_to = "Rep", values_to = "Pheno")
    All_Pheno_long_t1$Rep <- gsub("^V[0-9]+\\.|_rep[0-9]+$", "", All_Pheno_long_t1$Rep)
    All_Pheno_long_t1$ID <- as.numeric(All_Pheno_long_t1$`X[[i]]`)
    All_Pheno_long_t1 <- na.omit(All_Pheno_long_t1)
    All_Pheno_long_t1$`new_Pheno@id`<-NULL
    All_Pheno_long_t1$`X[[i]]`<-NULL
  
    All_Pheno_long_t2<-trait2_df %>%
      pivot_longer(cols = c(-`X[[i]]`, -Mean), names_to = "Rep", values_to = "Pheno")
    All_Pheno_long_t2$Rep <- gsub("^V[0-9]+\\.|_rep[0-9]+$", "", All_Pheno_long_t2$Rep)
    All_Pheno_long_t2$ID <- as.numeric(All_Pheno_long_t2$`X[[i]]`)
    All_Pheno_long_t2 <- na.omit(All_Pheno_long_t2)
    All_Pheno_long_t2$`new_Pheno@id`<-NULL
    All_Pheno_long_t2$`X[[i]]`<-NULL
    
    # Merge phenotype data with their inbreeding coefficient
    All_Pheno_long_t1 <- left_join(All_Pheno_long_t1, Inbreeding, by = "ID")
    All_Pheno_long_t2 <- left_join(All_Pheno_long_t2, Inbreeding, by = "ID")
    
    # Apply inbreeding penalty only if F > 0
    
      # Trait 1
      All_Pheno_long_t1$Pheno_Inb <- All_Pheno_long_t1$Pheno + 
      abs(All_Pheno_long_t1$Pheno) * ifelse(All_Pheno_long_t1$F > 0, (InbDepr * All_Pheno_long_t1$F), 0)
    
      # Trait 2    
      All_Pheno_long_t2$Pheno_Inb <- All_Pheno_long_t2$Pheno + 
      abs(All_Pheno_long_t2$Pheno) * ifelse(All_Pheno_long_t2$F > 0, (InbDepr * All_Pheno_long_t2$F), 0)
    
    # Add to preexisting data
      All_Pheno_long <- rbind(All_Pheno_long,All_Pheno_long_t1,All_Pheno_long_t2)
  
      # Summerize to mean values to get same number of values per trait
      Trait1_mean <- All_Pheno_long_t1 %>%
        group_by(ID) %>%
        summarise(mean_pheno = mean(Pheno_Inb, na.rm = TRUE))
      
      Trait2_mean <- All_Pheno_long_t2 %>%
        group_by(ID) %>%
        summarise(mean_pheno = mean(Pheno_Inb, na.rm = TRUE))
      
      
      New_Pheno <- data.frame(Trait1 = Trait1_mean$mean_pheno,
                              Trait2 = Trait2_mean$mean_pheno)
      
      new_pop@pheno <- as.matrix(New_Pheno)
      
      
  # Keep only the last 5 generations of data
  if (gen < 6) {
    
 ############################## Estimate BVs ###################################
    model <- remlf90(
      fixed = Pheno_Inb ~ Rep,  # trait-specific means
      genetic = list(
        model = 'add_animal',
        pedigree = A,
        id = All_Pheno_long$ID
      ),
      data = All_Pheno_long,
    )
    
    BV.full <- ranef(model)$genetic
    BV <- model.matrix(model)$genetic %*% BV.full
    All_Pheno_long$EBV <- BV@x
    EBV_df <- All_Pheno_long %>% distinct(ID, EBV)
    All_Pheno_long$EBV <- NULL
  } else {
    # Ensure at least 5 previous generations exist

    G5_Pheno_long <- All_Pheno_long[All_Pheno_long$ID %in% gen_ids_to_keep,]
    

    # Estimate BVs using a multi-trait model
    model <- remlf90(
      fixed = Pheno_Inb ~ Rep,  # trait-specific means
      genetic = list(
        model = 'add_animal',
        pedigree = A,
        id = G5_Pheno_long$ID
      ),
      data = G5_Pheno_long,
    )
    
    BV.full <- ranef(model)$genetic
    BV <- model.matrix(model)$genetic %*% BV.full
    G5_Pheno_long$EBV <- BV@x
    EBV_df <- G5_Pheno_long %>% distinct(ID, EBV)
    
  }
  
  
  #EBV <- data.frame(ebv = model$g, ID = rownames(model$g))
  EBV_new <- EBV_df[EBV_df$ID %in% new_pop@id,]
  new_pop@ebv <- as.matrix(EBV_new$EBV)
  
  # Save population object
  pop_list[[paste0("G", gen)]] <- new_pop
  
  # Save generation data
    Acc <- as.data.frame(cor(new_pop@gv,new_pop@ebv)) # Prediction Accuracy
   
    # Estimated Selection intensity
    all_ebvs <- lapply(pop_list, function(pop) {
      ebv(pop)  # This will extract the EBVs from each Pop object
    })
    
    ebv_df <- data.frame(
      Generation = rep(seq_along(pop_list), times = sapply(all_ebvs, length)),  # Repeat generation number
      EBV = unlist(all_ebvs)  # Flatten the list of EBVs into a single vector
    )
    
    mean_EBV_new <- mean(Parents@ebv)
    mean_EBV_all <- mean(ebv_df$EBV)
    sd_EBV_all <- sd(ebv_df$EBV)
    intensity_EBV <- data.frame(Intensity = ((mean_EBV_new - mean_EBV_all) / sd_EBV_all),
                               Gen = gen)
    
    intensity_EBV_df <- rbind(intensity_EBV_df,intensity_EBV)
    
    # True Selection intensity
    all_tbvs <- lapply(pop_list, function(pop) {
      gv(pop)  # This will extract the EBVs from each Pop object
    })
    
    tbv_df <- data.frame(
      Generation = rep(seq_along(pop_list), times = sapply(all_tbvs, length)),  # Repeat generation number
      TBV = unlist(all_tbvs)  # Flatten the list of EBVs into a single vector
    )
    
    mean_TBV_new <- mean(Parents@gv)
    mean_TBV_all <- mean(tbv_df$TBV)
    sd_TBV_all <- sd(tbv_df$TBV)
    intensity_TBV <- data.frame(Intensity = ((mean_TBV_new - mean_TBV_all) / sd_TBV_all),
                                Gen = gen)
    
    intensity_TBV_df <- rbind(intensity_TBV_df,intensity_TBV)
    
    # Additive Genetic value (Breeding value)
    MeanG <- meanG(new_pop)  
    VarA <- varG(new_pop)
    Cor_A = cov2cor(VarA)[2, 1]
    Additive_new <- data.frame(Trait1_Mean = MeanG[1],
                              Trait1_Var = VarA[1,1],
                              Trait2_Mean = MeanG[2],
                              Trait2_Var = VarA[2,2],
                              Covariance = VarA[1,2],
                              Correlation = Cor_A)  
    Additive_df <- rbind(Additive_df, Additive_new)
    
    # Phenotype
    Phen.Mean <- meanP(new_pop)  
    VarP <- varP(new_pop)
    Cor_P = cov2cor(VarP)[2, 1]
    Pheno_new <- data.frame(Trait1_Mean = Phen.Mean[1],
                           Trait1_Var = VarP[1,1],
                           Trait2_Mean = Phen.Mean[2],
                           Trait2_Var = VarP[2,2],
                           Covariance = VarP[1,2],
                           Correlation = Cor_P) 
    
    Pheno_df <- rbind(Pheno_df, Pheno_new)        

    Acc <- as.data.frame(cor(new_pop@gv,new_pop@ebv)) # Prediction Accuracy
    Acc_new <- data.frame(Trait1 = Acc[1,1],
                          Trait2 = Acc[2,1])
    Acc_df <- rbind(Acc_df, Acc_new)
    
  #  Bias_new <- data.frame(Bias = mean((new_pop@gv[,1]*0.5+new_pop@gv[,2]*0.5)-new_pop@ebv))
    
    #                                   -new_pop@ebv),
    #                       Triat1 = mean(new_pop@gv[,2]-new_pop@ebv))
   # Bias_df <- rbind()

  # Update current population for next iteration
  current_pop <- new_pop
}
####### Loop End #########



  
  
  
  
  
  
######################### Visualize results ####################################
  plot(x = 1:(1+num_generations), y = F_df$Inbreeding, type = "l", col = "black", lwd = 3,
       xlab = "Generation", ylab = "Mean Self-kinship")
  
    plot(x = 1:(1+num_generations), y = Pheno_df$Correlation, type = "l", col = "black", lwd = 3,
       xlab = "Generation", ylab = "Phenotypic correlation")
  
  plot(x = 1:(1+num_generations), y = Additive_df$Correlation, type = "l", col = "black", lwd = 3,
       xlab = "Generation", ylab = "Genetic correlation")
  
  
  
  varRanges = range(c(Additive_df$Trait1_Var, Additive_df$Trait2_Var))
  plot(x = 1:(1+num_generations), y = Additive_df$Trait1_Var, type = "l", col = "black", lwd = 3,
       xlab = "Generation", ylab = "Variance of genetic values", ylim = varRanges) 
  lines(x = 1:(1+num_generations), y = Additive_df$Trait2_Var, type = "l", col = "purple", lty = 2, lwd = 3) 
  legend(x = "topleft", legend = c("1", "2"), title = "Trait",
         lwd = 3, lty = c(1, 2), col = c("black", "purple"))
  
  varRanges = range(c(Additive_df$Trait1_Mean, Additive_df$Trait2_Mean))
  plot(x = 1:(1+num_generations), y = Additive_df$Trait1_Mean, type = "l", col = "black", lwd = 3,
       xlab = "Generation", ylab = "Mean breeding values", ylim = varRanges) 
  lines(x = 1:(1+num_generations), y = Additive_df$Trait2_Mean, type = "l", col = "purple", lty = 2, lwd = 3) 
  legend(x = "topleft", legend = c("1", "2"), title = "Trait",
         lwd = 3, lty = c(1, 2), col = c("black", "purple"))
  
  
  varRanges = range(c(Pheno_df$Trait1_Var, Pheno_df$Trait2_Var))
  plot(x = 1:(1+num_generations), y = Pheno_df$Trait1_Var, type = "l", col = "black", lwd = 3,
       xlab = "Generation", ylab = "Variance of phenotypiv values", ylim = varRanges) 
  lines(x = 1:(1+num_generations), y = Pheno_df$Trait2_Var, type = "l", col = "purple", lty = 2, lwd = 3) 
  legend(x = "topleft", legend = c("1", "2"), title = "Trait",
         lwd = 3, lty = c(1, 2), col = c("black", "purple"))
  
  varRanges = range(c(Pheno_df$Trait1_Mean, Pheno_df$Trait2_Mean))
  plot(x = 1:(1+num_generations), y = Pheno_df$Trait1_Mean, type = "l", col = "black", lwd = 3,
       xlab = "Generation", ylab = "Mean phenotypiv values", ylim = varRanges) 
  lines(x = 1:(1+num_generations), y = Pheno_df$Trait2_Mean, type = "l", col = "purple", lty = 2, lwd = 3) 
  legend(x = "topleft", legend = c("1", "2"), title = "Trait",
         lwd = 3, lty = c(1, 2), col = c("black", "purple"))
  
  varRanges = range(c(Acc_df$Trait1, Acc_df$Trait2))
  plot(x = 1:(1+num_generations), y = Acc_df$Trait1, type = "l", col = "black", lwd = 3,
       xlab = "Generation", ylab = "Accuracy", ylim = varRanges) 
  lines(x = 1:(1+num_generations), y = Acc_df$Trait2, type = "l", col = "purple", lty = 2, lwd = 3) 
  legend(x = "topleft", legend = c("1", "2"), title = "Trait",
         lwd = 3, lty = c(1, 2), col = c("black", "purple"))
  
  # Selection intensities
    varRanges = range(c(intensity_EBV_df$Intensity, intensity_TBV_df$Intensity))
  plot(x = 1:(1+num_generations), y = intensity_EBV_df$Intensity, type = "l", col = "black", lwd = 3,
       xlab = "Generation", ylab = "Selection intensity", ylim = varRanges) 
    lines(x = 1:(1+num_generations), y = intensity_TBV_df$Intensity, type = "l", col = "purple", lty = 2, lwd = 3) 
  legend(x = "topright", legend = c("Estimated selection intensity", "True selection intensity"), title = "Trait",
         lwd = 3, lty = c(1, 2), col = c("black", "purple"))
  
  
  
  

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
                  Chr.Size.bp = ChrLen,
                  Chr.Size.mo = PG_Chr_M,
                  Mut.Rate = ConMR,
                  FoundersNr = nFounders,
                  Founders.ne = neFounders,
                  BurnInGen = nGenBurnIn,
                  BurnInCross = nCrossBurnIn,
                  BurnInProg = nProgBurnIn,
                  ParentsNr = nParents,
                  FamiliesNr = nCross,
                  ProgenyNr = nProg,
                  Trait = Trait,
                  Trait.mean = initMeanG,
                  Trait.varA = initVarG,
                  h2_G0 = h2_G0,
                  h2_G1 = h2_G1,
                  MatePlan = MatePlan,
                  Gamma = GAMMA,
                  GammaShape = GShape,
                  SegSite = SegSite)

# Simulation Parameters
filename <- paste0("Sim_Parameters_Ne1000_MultiGen", today_date, ".txt")
write.table(Sim, filename, quote = F, col.names = T, row.names = F, sep = "\t")