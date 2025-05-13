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
    nCross = 50     # Number of families to create from Parents
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

# Save pedigree
  Pedigree_All <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G0_pop@id,])
  Pedigree_All$mother <- 0 # Reset pedigree for Plus trees
  Pedigree_All$father <- 0 # Reset pedigree for Plus trees
  Pedigree<-as.data.frame(SP$pedigree)
  
  # Save as data frame
  A <- data.frame(ID = as.numeric(rownames(Pedigree_All)),
                  Sire = Pedigree_All$father,
                  Dam = Pedigree_All$mother)

#### Add SNP chips ####

# Randomly assigns eligible SNPs to SNP chip
# Use the latest population as reference
  SP$addSnpChip(nSNP, minSnpFreq = minMAF,name = "GS", refPop = G0_pop)
  SP$addSnpChip(nSNP2, minSnpFreq = minMAF,name = "Inb", refPop = G0_pop)

# Calculate initial inbreeding
                    # Alternative 1: A-matrix
  Kinship_matrix<- 2*kinship(id = as.numeric(rownames(Pedigree_All)),
                        dadid = Pedigree_All$father,
                        momid = Pedigree_All$mother)


                    # Alternative 2: G-matrix
# Genotype with imaginary high-density SNP chip to estimate an inbreeding baseline
# G_G0<-pullSnpGeno(G0_pop, snpChip = "Inb")-1
# Kinship_matrix<-A.mat(G_G0, min.MAF = 0.01)


  F_0<-mean(diag(Kinship_matrix))-1

  # Save first estimate of inbreeding
    Inbreeding<-data.frame(F = F_0,
                          ID = as.numeric(G0_pop@id))

# Save mean self-kinship 
  F_df <- data.frame(Inbreeding = mean(diag(Kinship_matrix)))

# Initialize coancestry 
  Ns_df <- data.frame(Status_number = 0,
                         Coancestry = 0,
                                Gen = 0)
# Add first phenotypes
    G0_Pheno <- data.frame(genotype_id = (G0_pop@id))  # Start with IDs

    
      pheno <- setPheno(G0_pop,
                      varE = initVarE,
                      corE = matrix(c(1, CorP,
                                      CorP, 1), nrow = 2),
                      simParam = NULL)
  
      G0_Pheno[[paste0("trait1_rep", i)]] <- pheno@pheno[, 1]
      G0_Pheno[[paste0("trait2_rep", i)]] <- pheno@pheno[, 2]

  
  # Trait 1 columns (those starting with "trait1_")
  trait1_cols <- grep("^trait1_", names(G0_Pheno), value = TRUE)
  trait1_df <- G0_Pheno[, c("genotype_id", trait1_cols)]

  # Trait 2 columns (those starting with "trait2_")
  trait2_cols <- grep("^trait2_", names(G0_Pheno), value = TRUE)
  trait2_df <- G0_Pheno[, c("genotype_id", trait2_cols)]


  New_Pheno_long_t1<-trait1_df %>%
    pivot_longer(cols = c(-genotype_id), names_to = "Rep", values_to = "Pheno")
  New_Pheno_long_t1$Rep <- gsub("^V[0-9]+\\.|_rep[0-9]+$", "", New_Pheno_long_t1$Rep)
  New_Pheno_long_t1$ID <- as.numeric(New_Pheno_long_t1$genotype_id) 
  New_Pheno_long_t1 <- na.omit(New_Pheno_long_t1)
  New_Pheno_long_t1$`G0_pop@id`<-NULL
  New_Pheno_long_t1$genotype_id<-NULL

  New_Pheno_long_t2<-trait2_df %>%
    pivot_longer(cols = c(-genotype_id), names_to = "Rep", values_to = "Pheno")
  New_Pheno_long_t2$Rep <- gsub("^V[0-9]+\\.|_rep[0-9]+$", "", New_Pheno_long_t2$Rep)
  New_Pheno_long_t2$ID <- as.numeric(New_Pheno_long_t2$genotype_id)
  New_Pheno_long_t2 <- na.omit(New_Pheno_long_t2)
  New_Pheno_long_t2$`G0_pop@id`<-NULL
  New_Pheno_long_t2$genotype_id<-NULL
  

  # Add linear penalty for inbreeding depression only if F > 0
  New_Pheno_long_t1 <- left_join(New_Pheno_long_t1, Inbreeding, by = "ID")
  New_Pheno_long_t2 <- left_join(New_Pheno_long_t2, Inbreeding, by = "ID")
  

    New_Pheno_long_t1$Pheno_Inb <- New_Pheno_long_t1$Pheno + 
      abs(New_Pheno_long_t1$Pheno) * ifelse(New_Pheno_long_t1$F > 0, (InbDepr * New_Pheno_long_t1$F), 0)
  

    New_Pheno_long_t2$Pheno_Inb <- New_Pheno_long_t2$Pheno + 
      abs(New_Pheno_long_t2$Pheno) * ifelse(New_Pheno_long_t2$F > 0, (InbDepr * New_Pheno_long_t2$F), 0)
  

  All_Pheno_long <- rbind(New_Pheno_long_t1, New_Pheno_long_t2)
  All_Pheno_long$Mean <- All_Pheno_long$Pheno_Inb
  
All_Pheno_long <- All_Pheno_long %>%                                            # Adding standardizing the phenotypes, as to not give too high importance of later measurements                                                                                                                                                                                                               
   mutate(Pheno_std = scale(Pheno_Inb)) %>%                                     # even if there could be an argument for keeping it that way, as the older tree is closer              
   ungroup()     
  
  # Summerize to mean values to get same number of values per trait
  Trait1_mean <- New_Pheno_long_t1 %>%
    group_by(ID) %>%
    summarise(mean_pheno = mean(Pheno_Inb, na.rm = TRUE))
  
  Trait2_mean <- New_Pheno_long_t2 %>%
    group_by(ID) %>%
    summarise(mean_pheno = mean(Pheno_Inb, na.rm = TRUE))
  
  
    G0_Pheno <- data.frame(Trait1 = Trait1_mean$mean_pheno,
                           Trait2 = Trait2_mean$mean_pheno)

    # Add the phenotypes to the pop-object
      G0_pop@pheno <- as.matrix(G0_Pheno)
  
  # Save object to use in modeling
     Pheno_mean_All <- data.frame(ID = G0_Pheno$`G0_pop@id`, Trait = G0_Pheno$Mean)


  
  # Save data to estimate breeding values
  Pheno_mean_All <- data.frame(ID = G0_Pheno$`G0_pop@id`,
                               Trait = G0_Pheno$Mean)

############################# Mixed model 1 ####################################
 model_reml <- remlf90(
    fixed = Pheno_std ~ Rep,  # trait-specific means
    genetic = list(
      model = 'add_animal',
      pedigree = A,
      id = All_Pheno_long$ID
    ),
    data = All_Pheno_long
  )


  
 #####################
  # It is not a true multi-trait model. breedR cannot model covaraince structures 
  # between traits. Now they are assumed to be independendent, which is very not true
  
  # Can try sommer again
    # Need to format my data to handle sommer
  

  ################################
   
  library(sommer)
  
  #All_Pheno_long$Rep<-as.factor(All_Pheno_long$Rep)
  #All_Pheno_long$ID<-as.factor(All_Pheno_long$ID)
  #All_Pheno_long$Pheno_Inb<-as.numeric(All_Pheno_long$Pheno_Inb)
  #All_Pheno_long$Trait <- All_Pheno_long$Rep
  #All_Pheno_long$Trait<-as.factor(All_Pheno_long$Trait)
  
#Kinship_matrix <- as.matrix(Kinship_matrix)
  
#### Sommer model ####
#Ai <- solve(Kinship_matrix + diag(1e-4,ncol(Kinship_matrix),ncol(Kinship_matrix)))
#Ai <- as(as(as( Ai, "dMatrix"), "generalMatrix"), "CsparseMatrix")

#All_Pheno_long <- All_Pheno_long %>%
#  arrange(by=Trait)


#All_Pheno_long <- All_Pheno_long %>%                                            # Adding standardizing the phenotypes, as to not give too high importance of later measurements                                                                                                                                                                                                               
#   mutate(Pheno_std = scale(Pheno_Inb)) %>%                                     # even if there could be an argument for keeping it that way, as the older tree is closer              
#  ungroup()                                                                     # to clear cut age 

#   model_sommer <- mmes(
#    Pheno_std ~ Rep-1,                                                          # Trait-specific intercepts (can add trait:block etc. too)
#    random = ~ vsm(usm(Rep),ism(ID), Gu=Ai) ,                                   # Random genetic effects by trait using A matrix
#    rcov   = ~ units,
#    data = All_Pheno_long
#  )
  
#    Gmat <-model_sommer$theta[[1]]
#    Rmat<-model_sommer$theta[[2]]
    
    # Genetic correlations
#    Gcor <- cov2cor(Gmat)
    
    # Residual correlations
#    Rcor <- cov2cor(Rmat)



#EBVs<-randef(model_sommer)
#EBVsTrait1<-data.frame(EBV=EBVs[1:100],
#                       ID = Ai@Dimnames[[1]])

#EBVsTrait2<-data.frame(EBV=EBVs[101:200],
#                       ID = Ai@Dimnames[[1]])

#EBVsTrait1 <- EBVsTrait1[match(G0_pop@iid, EBVsTrait1$ID), ]
#EBVsTrait2 <- EBVsTrait2[match(G0_pop@iid, EBVsTrait2$ID), ]


#GV_1<-data.frame(GV = G0_pop@gv[,1],
#                 ID = G0_pop@id)

#GV_2<-data.frame(GV = G0_pop@gv[,2],
#                 ID = G0_pop@id)


#cor(GV_1$GV,EBVsTrait1$EBV)
#cor(GV_2$GV,EBVsTrait2$EBV)

#EBVs <- as.matrix(EBVsTrait1$EBV,EBVsTrait2$EBV)

#G0_pop@ebv<-as.matrix(EBVsTrait1$EBV)
#G0_pop@ebv<-as.matrix(EBVsTrait2$EBV)

    
###########################
  
  # Estimate heritability
  var_comp<-summary(model_reml)
  Va <- var_comp$var[1]
  Ve <- var_comp$var[2]
  h2_df<-data.frame(h2 = Va/(Va+Ve))
  
  
  BV.full <- ranef(model_reml)$genetic
  BV <- model.matrix(model_reml)$genetic %*% BV.full
  All_Pheno_long$EBV <- BV@x
  EBV_df <- All_Pheno_long %>% distinct(ID, EBV)
  All_Pheno_long$EBV <- NULL
  
  EBV_new <- EBV_df[EBV_df$ID %in% G0_pop@id,]
  EBV_new <- EBV_new[match(G0_pop@iid, EBV_new$ID), ]
  
  cor(GV_2$GV,EBV_new$EBV)
  
  # The remlf90 model gives higher accuracy...
  # How to go forward...
  # Something may not be completely correctly set up in the sommer model
  
  G0_pop@ebv <- as.matrix(EBV_new$EBV)
  

  ###### Calculate MAF distribution ####
  G_G0<-pullSnpGeno(G0_pop, snpChip = "Inb")
  maf_values_df <- data.frame(Frequency = apply(G_G0, 2, calculate_maf),
                              Gen = 0)
  #maf_values0.05 <-maf_values[maf_values$Frequency > 0.01, ]
  
  
  current_pop <- G0_pop
  # Store populations for downstream analysis
  pop_list <- list(G0 = G0_pop)
  
######################### Save initial parameters ##############################
  # Additive genetic value   
  MeanG <- meanG(G0_pop)  
  VarA <- varG(G0_pop)
  Cor_A = cov2cor(VarA)[2, 1]
  Additive_df <- data.frame(Trait1_Mean = MeanG[1],
                            Trait1_Var = VarA[1,1],
                            Trait2_Mean = MeanG[2],
                            Trait2_Var = VarA[2,2],
                            Covariance = VarA[1,2],
                            Correlation = Cor_A)
  
  
  # Phenotypes
  Phen.Mean <- meanP(G0_pop)  
  VarP <- varP(G0_pop)
  Cor_P = cov2cor(VarP)[2, 1]
  Pheno_df <- data.frame(Trait1_Mean = Phen.Mean[1],
                         Trait1_Var = VarP[1,1],
                         Trait2_Mean = Phen.Mean[2],
                         Trait2_Var = VarP[2,2],
                         Covariance = VarP[1,2],
                         Correlation = Cor_P) 
  
  # Prediction accuracy 
  Acc <- as.data.frame(cor(G0_pop@gv,G0_pop@ebv)) # Prediction Accuracy
  Acc_df <- data.frame(Trait1 = Acc[1,1],
                       Trait2 = Acc[2,1])
  # Estimation bias
  Bias_df <- data.frame(Bias = mean((G0_pop@ebv) - ((G0_pop@gv[,1]+G0_pop@gv[,2])
                                                      /2)))
  # Estimate disversion 
  gv1 <- G0_pop@gv[,1]
  gv2 <- G0_pop@gv[,2]
  ebv <- G0_pop@ebv[,1]
  
  # Fit regression model
  bias_model1 <- lm(gv1 ~ ebv)
  bias_model2 <- lm(gv2 ~ ebv)
  
  # Optional: view slope too
  slope1 <- coef(bias_model1)[2]
  slope2 <- coef(bias_model2)[2]
  
  
  Dispersion_df <- data.frame(Trait1 = slope1,
                           Trait2 = slope2)

  gen= 1
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
  new_pop <- do.call(c, G_pops)
  
  # Save population object
  pop_list[[paste0("G", gen)]] <- new_pop
  
                             # # Record distributoin of family sizes
                            #        df <- data.frame(
                            #  ID = new_pop@id,        # Individual ID (as string)
                            #  father = new_pop@father, # Father's ID
                            #  mother = new_pop@mother,
                            #  pheno = new_pop@pheno) # Mother's ID # Phenotype for the first trait
   
                           # df$Family <- paste(df$father, df$mother, sep = "_")
                            
                          ##  ggplot(df, aes(x = Family, y = pheno.Trait1, fill = Family)) +
                          ##    geom_boxplot() +
                          #    theme_minimal() +
                          #    labs(title = "Phenotype Distribution by Family",
                          #         x = "Family",
                          #         y = "Phenotype") +
                          #    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                            
  # Update pedigree
  Pedigree_New <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% new_pop@id,])
  Pedigree_All <- rbind(Pedigree_All, Pedigree_New)
  
  if (gen > 5) {  # Ensure at least 5 previous generations exist
    gen_ids_to_keep <- unlist(lapply((gen-5):gen, function(g) pop_list[[paste0("G", g)]]@id))
  }

  if (gen < 6) {
    A <- data.frame(ID = as.numeric(rownames(Pedigree_All)),
                    Sire = Pedigree_All$father,
                    Dam = Pedigree_All$mother)
  
  } else {
    A <- Pedigree_All[rownames(Pedigree_All) %in% gen_ids_to_keep, ]
  
    A <- data.frame(ID = as.numeric(rownames(A)),
                    Sire = A$father,
                    Dam = A$mother)

  }
  
  
  A_Mat <- kinship(id = as.numeric(rownames(Pedigree_All)),
                  dadid = Pedigree_All$father,
                  momid = Pedigree_All$mother)
                          
                    indices<-which((rownames(A_Mat) %in% new_pop@id))
                    A_Mat <-A_Mat[indices,indices]
                    n <- nrow(A_Mat)
                          
                    group_coancestry <- sum(A_Mat) / (n^2)
                    
                    status_number <- 1 / (2*group_coancestry)
                    
                    Ns_new <- data.frame(Status_number = status_number,
                                                    Coancestry = group_coancestry,
                                                    Gen = gen)
                    Ns_df <- rbind(Ns_df,Ns_new)
  
  
  
  ########################## Estimate inbreeding ###############################
                    # Add mechanism for inbreeding depression
                          
                      # Using the A matrix 
        Kinship_matrix<-2*kinship(id = A$ID, dadid = A$Sire, momid = A$Dam)
        
        Inbreeding<-data.frame(F = (diag(Kinship_matrix)-1),
                     ID = A$ID)
  

                      # Using the G-matrix
        #G_new<-pullSnpGeno(new_pop, snpChip = "Inb")-1
        #Kinship_matrix<-A.mat(G_new, min.MAF = 0.01)
        #F_new <-(diag(Kinship_matrix))-F_0
        
        
        #Inbreeding<-data.frame(F = F_new,
        #            ID = as.numeric(new_pop@id))
        
        F_df <- rbind(F_df,mean(diag(Kinship_matrix)))
  
  
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
    
    ######### Apply inbreeding penalty #################
    
      # Trait 1
      All_Pheno_long_t1$Pheno_Inb <- All_Pheno_long_t1$Pheno + 
      abs(All_Pheno_long_t1$Pheno) * ifelse(All_Pheno_long_t1$F > 0, (InbDepr * All_Pheno_long_t1$F), 0)
    
      # Trait 2    
      All_Pheno_long_t2$Pheno_Inb <- All_Pheno_long_t2$Pheno + 
      abs(All_Pheno_long_t2$Pheno) * ifelse(All_Pheno_long_t2$F > 0, (InbDepr * All_Pheno_long_t2$F), 0)
      
      # Standardize
      All_Pheno_long_t1 <- All_Pheno_long_t1 %>%                                                                                                                                                                                                                                                          
        mutate(Pheno_std = scale(Pheno_Inb)) %>%                                                   
        ungroup()
        
      All_Pheno_long_t2 <- All_Pheno_long_t2 %>%                                                                                                                                                                                                                                                 
        mutate(Pheno_std = scale(Pheno_Inb)) %>%                                                  
        ungroup()
    
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
    
 ############################ Mixed model 2 ####################################
    
    model <- remlf90(
      fixed = Pheno_std ~ Rep,  # trait-specific means
      genetic = list(
        model = 'add_animal',
        pedigree = A,
        id = All_Pheno_long$ID
      ),
      data = All_Pheno_long,
    )
    
    # Estimate heritability
    var_comp<-summary(model)
    Va <- var_comp$var[1]
    Ve <- var_comp$var[2]
    h2_new<-data.frame(h2 = Va/(Va+Ve))
    
    h2_df <- rbind(h2_df, h2_new)
    
    BV.full <- ranef(model)$genetic
    BV <- model.matrix(model)$genetic %*% BV.full
    All_Pheno_long$EBV <- BV@x
    EBV_df <- All_Pheno_long %>% distinct(ID, EBV)
    All_Pheno_long$EBV <- NULL
  } else {
    # Ensure at least 5 previous generations exist

    G5_Pheno_long <- All_Pheno_long[All_Pheno_long$ID %in% gen_ids_to_keep,]
    

####################### Mixed model 3 ##########################################    
    model <- remlf90(
      fixed = Pheno_std ~ Rep,  # trait-specific means
      genetic = list(
        model = 'add_animal',
        pedigree = A,
        id = G5_Pheno_long$ID
      ),
      data = G5_Pheno_long,
    )
    
    var_comp<-summary(model)
    Va <- var_comp$var[1]
    Ve <- var_comp$var[2]
    h2_new<-data.frame(h2 = Va/(Va+Ve))
    
    h2_df <- rbind(h2_df, h2_new)
    
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
    
    Bias_new <- data.frame(Bias = mean((new_pop@ebv) - ((new_pop@gv[,1]+new_pop@gv[,2])
                                                        /2)))
    
    Bias_df <- rbind(Bias_df,Bias_new)
    
    # Estimate disversion 
    gv1 <- new_pop@gv[,1]
    gv2 <- new_pop@gv[,2]
    ebv <- new_pop@ebv[,1]
    
    # Fit regression model
    bias_model1 <- lm(gv1 ~ ebv)
    bias_model2 <- lm(gv2 ~ ebv)
    
    # Optional: view slope too
    slope1 <- coef(bias_model1)[2]
    slope2 <- coef(bias_model2)[2]
    
    
    Dispersion_new <- data.frame(Trait1 = slope1,
                                 Trait2 = slope2)
    
    Dispersion_df <- rbind(Dispersion_df,Dispersion_new)
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
  
  plot(x = 1:(1+num_generations), y = h2_df$h2, type = "l", col = "black", lwd = 3,
       xlab = "Generation", ylab = "Estimated narrow sense heritability")
  
  plot(x = 1:(1+num_generations), y = Ns_df$Status_number, type = "l", col = "black", lwd = 3,
       xlab = "Generation", ylab = "Status Number")
  
  plot(x = 1:(1+num_generations), y = Ns_df$Coancestry, type = "l", col = "black", lwd = 3,
       xlab = "Generation", ylab = "Group Coancestry")
  
  plot(x = 1:(1+num_generations), y = Bias_df$Bias, type = "l", col = "black", lwd = 3,
       xlab = "Generation", ylab = "Bias")
  
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
  
  varRanges = range(c(Dispersion_df$Trait1, Dispersion_df$Trait2))
  plot(x = 1:(1+num_generations), y = Dispersion_df$Trait1, type = "l", col = "black", lwd = 3,
       xlab = "Generation", ylab = "Dispersion Bias", ylim = varRanges) 
  lines(x = 1:(1+num_generations), y = Dispersion_df$Trait2, type = "l", col = "purple", lty = 2, lwd = 3) 
  legend(x = "topright", legend = c("1", "2"), title = "Trait",
         lwd = 3, lty = c(1, 2), col = c("black", "purple"))
  
  # Selection intensities
    varRanges = range(c(intensity_EBV_df$Intensity, intensity_TBV_df$Intensity))
  plot(x = 1:(1+num_generations), y = intensity_EBV_df$Intensity, type = "l", col = "black", lwd = 3,
       xlab = "Generation", ylab = "Selection intensity", ylim = varRanges) 
    lines(x = 1:(1+num_generations), y = intensity_TBV_df$Intensity, type = "l", col = "darkred", lty = 2, lwd = 3) 
  legend(x = "topright", legend = c("Estimated SI", "True SI"), title = "Trait",
         lwd = 3, lty = c(1, 2), col = c("black", "darkred"))
  
  
  
  

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