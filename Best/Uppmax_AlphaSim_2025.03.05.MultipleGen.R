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
    initMeanG = 200     # Initial trait mean genetic value
    CV = 0.3 # The standard deviation ~ 30 % of the mean 
    initVarG = (CV * initMeanG)^2  # Initial trait genetic variance
    
    # Heritability
    h2_G0 = 0.2
    h2_G1 = 0.2
    

# Founder Parameters    
  # Founder Population
    nFounders = 100   # Number of founders in the population
    neFounders = 10    # Effective population size of founders

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
  # SNP chip
    nSNP = 100     # Nr of SNPs per chromosome
    minMAF = 0.01  # Minimum MAF for inclusion

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
             shape        = GShape)

# Create new pop object
G0_pop = newPop(founderPop)
rm(founderPop)

# Keep track of the fluctuations in genetic mean and variance
VarA_list<-list()
VarA_G0 <- varG(G0_pop)
VarA_list <- c(VarA_list, VarA_G0)

MeanA_list<-list()
MeanA_G0 <- meanG(G0_pop)
MeanA_list <- c(MeanA_list, MeanA_G0)

################################## Burn-in #####################################
for (i in 1:nGenBurnIn) {
  # Perform random crossing
  G0_pop <- randCross(
    pop = G0_pop,            # The population to cross
    nCrosses = nCrossBurnIn, # Total number of crosses
    nProgeny = nProgBurnIn,  # Number of progeny per cross
    ignoreSexes = TRUE)      # Ignore sexes in crossing
  VarA_G0 <- varG(G0_pop)
  VarA_list <- c(VarA_list, VarA_G0)
  MeanA_G0 <- meanG(G0_pop)
  MeanA_list <- c(MeanA_list, MeanA_G0)
}

###################### Founder generation ######################################
                    
G0_pop <- selectInd(G0_pop, nG0, use = "rand")

G0_Pheno <- as.data.frame(G0_pop@id)

# Loop through and assign phenotypic values
for (i in 1:n_reps) {
  G0_Pheno[[paste0("V", i)]] <- setPheno(G0_pop,
                                         h2 = h2_G0, 
                                         fixEff = 1L,
                                         onlyPheno = TRUE,
                                         traits = 1,
                                         reps = 1,
                                         simParam = NULL)
}

# Get total number of values excluding genotype_id
num_values <- prod(dim(G0_Pheno[, -1]))
num_na <- round(0.1 * num_values)
rows <- sample(nrow(G0_Pheno), num_na, replace = TRUE)
cols <- sample(2:ncol(G0_Pheno), num_na, replace = TRUE) 
# Assign NA to selected positions
for (i in seq_along(rows)) {
  G0_Pheno[rows[i], cols[i]] <- NA
}

for (i in seq_along(cols)) {
  G0_Pheno[rows[i], cols[i]] <- NA
}



G0_Pheno$Mean <- rowMeans(G0_Pheno[ , -1], na.rm = T)

All_Pheno_long<-G0_Pheno %>%
  pivot_longer(cols = c(-`G0_pop@id`, -Mean), names_to = "Rep", values_to = "Pheno")
All_Pheno_long$ID <- All_Pheno_long$`G0_pop@id`
All_Pheno_long <- na.omit(All_Pheno_long)
All_Pheno_long$`G0_pop@id`<-NULL

G0_pop@pheno <- as.matrix(G0_Pheno$Mean)

# Randomly assigns eligible SNPs to SNP chip
# Use the latest population as reference
SP$addSnpChip(nSNP, minSnpFreq = minMAF, refPop = G0_pop)

####################### Large Loop Start #######################################

# Save Parameters during looping

Acc_List <- list()
Add.Var_List <- list()
Phen.Mean_List <- list()
Phen.Var_List <- list()
BV.Mean_List <- list()
BV.Var_List <- list()


# Initialize populations
current_pop <- G0_pop
Pheno_mean_All <- data.frame(ID = G0_Pheno$`G0_pop@id`, Trait = G0_Pheno$Mean)
Pedigree_All <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G0_pop@id,])
Pedigree_All$mother <- 0 # Reset pedigree for Plus trees
Pedigree_All$father <- 0 # Reset pedigree for Plus trees

# Store populations for downstream analysis
pop_list <- list(G0 = G0_pop)

# Loop through generations
for (gen in 1:num_generations) {
  cat("Generation:", gen, "\n")
  
  # Select best individuals based on phenotype or EBV
  if (gen == 1) {
    Parents <- selectInd(current_pop,
                         nInd = nParents,
                         trait = 1,
                         use = "pheno",
                         selectTop = TRUE)
  } else {
    Parents <- selectWithinFam(pop = current_pop,
                               nInd = 1, # Select the top within each family
                               use = "ebv")
  }
  
  # Generate progeny
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
  
    # Calculate relationship matrix
    #A_matrix <- 2 * kinship(id = A$ID,
    #                        dadid = A$Sire,
    #                        momid = A$Dam)
  } else {
    A <- Pedigree_All[rownames(Pedigree_All) %in% gen_ids_to_keep, ]
  
    A <- data.frame(ID = rownames(A),
                    sire = A$father,
                    dam = A$mother)
  
    # Calculate relationship matrix
    #A_matrix <- 2 * kinship(id = A$id,
    #                        dadid = A$sire,
    #                        momid = A$dam)
  }
  
  A$ID <-as.numeric(A$ID)
  
  #### Assign phenotypic values ####
  new_Pheno <- as.data.frame(new_pop@id)
  for (i in 1:n_reps) {
    new_Pheno[[paste0("V", i)]] <- setPheno(new_pop,
                                            h2 = h2_G1,
                                            fixEff = 1L,
                                            onlyPheno = TRUE,
                                            traits = 1,
                                            reps = 1)
  }
  
  # Introduce missing values
  num_values <- prod(dim(new_Pheno[, -1]))
  num_na <- round(0.1 * num_values)
  rows <- sample(nrow(new_Pheno), num_na, replace = TRUE)
  cols <- sample(2:ncol(new_Pheno), num_na, replace = TRUE)
  for (i in seq_along(rows)) new_Pheno[rows[i], cols[i]] <- NA
  new_Pheno$Mean <- rowMeans(new_Pheno[, -1], na.rm = TRUE)
  new_pop@pheno <- as.matrix(new_Pheno$Mean)
  
  # Store phenotypic values
  new_Pheno_long<-new_Pheno %>%
    pivot_longer(cols = c(-`new_pop@id`, -Mean), names_to = "Rep", values_to = "Pheno")
    new_Pheno_long$ID <- new_Pheno_long$`new_pop@id`
    new_Pheno_long <- na.omit(new_Pheno_long)
    new_Pheno_long$`new_pop@id`<-NULL
    
    All_Pheno_long <- rbind(All_Pheno_long,new_Pheno_long)
    
    All_Pheno_long$ID <- as.numeric(All_Pheno_long$ID)
  
  # Update phenotype dataset
  Pheno_mean_New <- data.frame(ID = new_Pheno$`new_pop@id`, Trait = new_Pheno$Mean)
  Pheno_mean_All <- rbind(Pheno_mean_All, Pheno_mean_New)

  # Keep only the last 5 generations of data
  if (gen < 6) {
    #model <- kin.blup(data = Pheno_mean_All,
    ##                  geno = "ID",
    #                  pheno = "Trait",
    #                  K = A_matrix)
    
    model<-remlf90(fixed = Pheno ~ 1,
                        genetic = list(model = 'add_animal',
                                       pedigree = A,
                                       id = All_Pheno_long$ID),
                        data = All_Pheno_long)
    
    BV.full <- ranef(model)$genetic
    BV <- model.matrix(model)$genetic %*% BV.full
    All_Pheno_long$EBV <- BV@x
    EBV_df <- All_Pheno_long %>% distinct(ID, EBV)
    All_Pheno_long$EBV <- NULL
  } else {
    # Ensure at least 5 previous generations exist

    G5_Pheno <- All_Pheno_long[All_Pheno_long$ID %in% gen_ids_to_keep,]
    
   # model <- kin.blup(data = Pheno_mean_5G,
  #                geno = "ID",
  #                  pheno = "Trait",
  #                  K = A_matrix)
    
    model<-remlf90(fixed = Pheno ~ 1,
                        genetic = list(model = 'add_animal',
                                       pedigree = A,
                                       id = G5_Pheno$ID),
                        data = G5_Pheno)
    
    BV.full <- ranef(model)$genetic
    BV <- model.matrix(model)$genetic %*% BV.full
    G5_Pheno$EBV <- BV@x
    EBV_df <- G5_Pheno %>% distinct(ID, EBV)
    
  }
  
  
  #EBV <- data.frame(ebv = model$g, ID = rownames(model$g))
  EBV_new <- EBV_df[EBV_df$ID %in% new_pop@id,]
  new_pop@ebv <- as.matrix(EBV_new$EBV)
  
  # Save population object
  pop_list[[paste0("G", gen)]] <- new_pop
  
  Acc <- cor(new_pop@gv,new_pop@ebv)
  Add.Var<-varA(new_pop)
  Phen.Mean <- mean(new_Pheno_long$Pheno)
  Phen.Var <- var(new_Pheno_long$Pheno)
  BV.Mean <- mean(new_pop@gv)
  BV.Var <- var(new_pop@gv)
  
  Acc_List <- list(Acc_List, Acc)
  Add.Var_List <- list(Add.Var_List, Add.Var)
  Phen.Mean_List <- list(Phen.Mean_List, Phen.Mean)
  Phen.Var_List <- list(Phen.Var_List,Phen.Var)
  BV.Mean_List <- list(BV.Mean_List, BV.Mean)
  BV.Var_List <- list(BV.Var_List, BV.Var)
  

  # Update current population for next iteration
  current_pop <- new_pop
}

Generations <- data.frame(Generation = 1:num_generations,
                          Accuracy = unlist(Acc_List),
                          AdditiveVariance = unlist(Add.Var_List),
                          PhenotypeMean = unlist(Phen.Mean_List),
                          PhenotypeVar = unlist(Phen.Var_List),
                          BV_Mean = unlist(BV.Mean_List),
                          BV_Var = unlist(BV.Var_List))

1-330/440

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