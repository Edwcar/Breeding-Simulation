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

######## Set initial values ########
# Nr of chromosomes 
PaChr = 12 # Source: Nystedt et al., 2013 https://doi.org/10.1038/nature12211

# Chromosome length in bp 
ChrLen = ((19.6*10^9)/PaChr) # Source: Source: Nystedt et al., 2013 https://doi.org/10.1038/nature12211

# Chromosome length in morgans 
PG_Chr_M <- ((3556/100)/PaChr) # Source: Bernhardsson et. al. 2019. https://doi.org/10.1534/g3.118.200840

MatAge = 45   # Average age age at reproduction during burn-in      

# Mutation rate between generations
ConMR = (2.2*10^-9) * MatAge # Source: Nystedt et al., 2013 https://doi.org/10.1038/nature12211

# Founder Population
nFounders = 100   # Number of founders in the population
neFounders = 10    # Effective population size of founders

nGenBurnIn <- 100  # Number of descrete burn-in generations
nCrossBurnIn <- 100 # Number of crosses within each burn-in generation
nProgBurnIn <- 10

# Size of the G0 population
nG0 = 100

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
SD = 0.3 # Size of the standard deviation ~ 30 % of the mean 
initVarG = (SD * initMeanG)^2  # Initial trait genetic variance

# Clonal Heritability*
h2_G0 = 0.2
h2_G1 = 0.2

# Crosses
nParents = 10       # Number of parents to be crossed
nCross = 10     # Number of families to create from Parents
nProg = 10      # Number of progeny per cross
# Comment: Create ~ 1000 individuals

# Segregating sites
SegSite = 1000

# SNP chip
nSNP = 100     # Nr of SNPs per chromosome
# Comment: Replicate the number of real SNPs after filtering (~ 39600) 

# Mating schemes
MatePlan = "RandCross"

RandDe = 0.1   # Percentage random drift

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

####### Burn-in #########
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

G0_pop <- selectInd(G0_pop, nG0, use = "rand")

G0_Pheno <- as.data.frame(G0_pop@id)

# Number of replicates
n_reps <- 20

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
G0_pop@pheno <- as.matrix(G0_Pheno$Mean)

# Randomly assigns eligible SNPs to SNP chip
# Use the latest population as reference
SP$addSnpChip(nSNP, minSnpFreq = 0.01, refPop = G0_pop)

# Generate Phenotype Data

# Define number of generations
num_generations <- 3  # Change as needed

# Initialize populations
current_pop <- G0_pop
G_Pheno_mean <- data.frame(ID = G0_Pheno$`G0_pop@id`, Trait = G0_Pheno$Mean)
Ped_G <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% G0_pop@id,])
Ped_G$mother <- 0
Ped_G$father <- 0

# Store populations for downstream analysis
pop_list <- list(G0 = G0_pop)

# Loop through generations
for (gen in 1:num_generations) {
  cat("Generation:", gen, "\n")
  
  # Select best individuals based on phenotype or EBV
  if (gen == 1) {
    Parents <- selectInd(current_pop, nParents, trait = 1, use = "pheno", selectTop = TRUE)
  } else {
    Parents <- selectWithinFam(pop = current_pop, nInd = 1, use = "ebv")
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
  Ped_new <- as.data.frame(SP$pedigree[rownames(SP$pedigree) %in% new_pop@id,])
  Ped_G <- rbind(Ped_G, Ped_new)
  Pedigree <- data.frame(ID = rownames(Ped_G), Sire = Ped_G$father, Dam = Ped_G$mother)
  
  # Calculate relationship matrix
  A_matrix <- 2 * kinship(id = Pedigree$ID, dadid = Pedigree$Sire, momid = Pedigree$Dam)
  
  # Assign phenotypic values
  new_Pheno <- as.data.frame(new_pop@id)
  for (i in 1:10) {
    new_Pheno[[paste0("V", i)]] <- setPheno(new_pop, h2 = h2_G1, fixEff = 1L, onlyPheno = TRUE, traits = 1, reps = 1)
  }
  
  # Introduce missing values
  num_values <- prod(dim(new_Pheno[, -1]))
  num_na <- round(0.1 * num_values)
  rows <- sample(nrow(new_Pheno), num_na, replace = TRUE)
  cols <- sample(2:ncol(new_Pheno), num_na, replace = TRUE)
  for (i in seq_along(rows)) new_Pheno[rows[i], cols[i]] <- NA
  new_Pheno$Mean <- rowMeans(new_Pheno[, -1], na.rm = TRUE)
  new_pop@pheno <- as.matrix(new_Pheno$Mean)
  
  

  # Update phenotype dataset
  new_Pheno_mean <- data.frame(ID = new_Pheno$`new_pop@id`, Trait = new_Pheno$Mean)
  G_Pheno_mean <- rbind(G_Pheno_mean, new_Pheno_mean)
  
  # Calculate EBVs
  model <- kin.blup(data = G_Pheno_mean, geno = "ID", pheno = "Trait", K = A_matrix)
  EBV <- data.frame(ebv = model$g, ID = G_Pheno_mean$ID)
  EBV_new <- EBV[EBV$ID %in% new_pop@id,]
  new_pop@ebv <- as.matrix(EBV_new$ebv)
  
  # Save population object
  pop_list[[paste0("G", gen)]] <- new_pop
  
  # Update current population for next iteration
  current_pop <- new_pop
}














####################################### G3 #####################################
G_pops <- vector("list", nParents)

for (i in 1:nParents) {
  # Generate random number of progeny for this cross (with mean = 30 and sd = 10)
  Prog_number <- max(round(rnorm(1, mean = nProg, sd = 10)), 1)
  
  # Create a population for this iteration
  G_pops[[i]] <- randCross(
    pop = New_candidates,   # The population to cross
    nCrosses = 1,    # Total number of crosses
    nProgeny = Prog_number,  # Number of progeny per cross (randomly generated)
    balance = NULL,   # If using sexes, this option will balance the number of progeny per parent
    parents = NULL,
    ignoreSexes = TRUE  # Should sexes be ignored?
  )
}

G3_pop <- do.call(c, G_pops)

Ped<-SP$pedigree

Ped_G3 <- as.data.frame(Ped[rownames(Ped) %in% G3_pop@id,])

Ped_G<-rbind(Ped_G0,Ped_G1,Ped_G2,Ped_G3)

Pedigree <- data.frame(ID = rownames(Ped_G),
                       Sire = Ped_G$father,
                       Dam = Ped_G$mother)

# Calcualte the A matrix
A_matrix<-2*kinship(id = Pedigree$ID, dadid = Pedigree$Sire, momid = Pedigree$Dam, )

G3_Pheno <- as.data.frame(G3_pop@id)

# Number of replicates
n_reps <- 10

# Loop through and assign phenotypic values
for (i in 1:n_reps) {
  G3_Pheno[[paste0("V", i)]] <- setPheno(G3_pop,
                                         h2 = h2_G1, 
                                         fixEff = 1L,
                                         onlyPheno = TRUE,
                                         traits = 1,
                                         reps = 1,
                                         simParam = NULL)
}

num_values <- prod(dim(G3_Pheno[, -1]))
num_na <- round(0.1 * num_values)
rows <- sample(nrow(G3_Pheno), num_na, replace = TRUE)
cols <- sample(2:ncol(G3_Pheno), num_na, replace = TRUE) 
# Assign NA to selected positions
for (i in seq_along(rows)) {
  G3_Pheno[rows[i], cols[i]] <- NA
}

for (i in seq_along(cols)) {
  G3_Pheno[rows[i], cols[i]] <- NA
}

G3_Pheno$Mean <- rowMeans(G3_Pheno[ , -1], na.rm = T)

G3_pop@pheno <- as.matrix(G3_Pheno$Mean)

G3_Pheno_mean <- data.frame(
  ID = G3_Pheno$`G3_pop@id`,
  Trait = G3_Pheno$Mean)

G_Pheno_mean <- rbind(G0_Pheno_mean, G1_Pheno_mean, G2_Pheno_mean, G3_Pheno_mean)


# Calcualte EBVs
model<-kin.blup(data = G_Pheno_mean,
                geno = "ID",
                pheno = "Trait",
                K = A_matrix)

EBV <- data.frame(ebv = model$g,
                  ID = G_Pheno_mean$ID)

EBV_G3 <- EBV[EBV$ID %in% G3_pop@id,]

G3_pop@ebv <- as.matrix(EBV_G3$ebv)

New_candidates<-selectWithinFam(pop = G3_pop,
                                nInd = 1, # Number Selected within each family
                                use = "ebv")

##################################### G4 #######################################
G_pops <- vector("list", nParents)

for (i in 1:nParents) {
  # Generate random number of progeny for this cross (with mean = 30 and sd = 10)
  Prog_number <- max(round(rnorm(1, mean = nProg, sd = 10)), 1)
  
  # Create a population for this iteration
  G_pops[[i]] <- randCross(
    pop = New_candidates,   # The population to cross
    nCrosses = 1,    # Total number of crosses
    nProgeny = Prog_number,  # Number of progeny per cross (randomly generated)
    balance = NULL,   # If using sexes, this option will balance the number of progeny per parent
    parents = NULL,
    ignoreSexes = TRUE  # Should sexes be ignored?
  )
}

G4_pop <- do.call(c, G_pops)

Ped<-SP$pedigree

Ped_G4 <- as.data.frame(Ped[rownames(Ped) %in% G4_pop@id,])

Ped_G<-rbind(Ped_G0,Ped_G1,Ped_G2,Ped_G3,Ped_G4)

Pedigree <- data.frame(ID = rownames(Ped_G),
                       Sire = Ped_G$father,
                       Dam = Ped_G$mother)

# Calcualte the A matrix
A_matrix<-2*kinship(id = Pedigree$ID, dadid = Pedigree$Sire, momid = Pedigree$Dam, )

G4_Pheno <- as.data.frame(G4_pop@id)

# Number of replicates
n_reps <- 10

# Loop through and assign phenotypic values
for (i in 1:n_reps) {
  G4_Pheno[[paste0("V", i)]] <- setPheno(G4_pop,
                                         h2 = h2_G1, 
                                         fixEff = 1L,
                                         onlyPheno = TRUE,
                                         traits = 1,
                                         reps = 1,
                                         simParam = NULL)
}

num_values <- prod(dim(G4_Pheno[, -1]))
num_na <- round(0.1 * num_values)
rows <- sample(nrow(G4_Pheno), num_na, replace = TRUE)
cols <- sample(2:ncol(G4_Pheno), num_na, replace = TRUE) 
# Assign NA to selected positions
for (i in seq_along(rows)) {
  G4_Pheno[rows[i], cols[i]] <- NA
}

for (i in seq_along(cols)) {
  G4_Pheno[rows[i], cols[i]] <- NA
}


G4_Pheno$Mean <- rowMeans(G4_Pheno[ , -1], na.rm = T)

G4_pop@pheno <- as.matrix(G4_Pheno$Mean)

G4_Pheno_mean <- data.frame(
  ID = G4_Pheno$`G4_pop@id`,
  Trait = G4_Pheno$Mean)

G_Pheno_mean <- rbind(G0_Pheno_mean, G1_Pheno_mean, G2_Pheno_mean, G3_Pheno_mean, G4_Pheno_mean)


# Calcualte EBVs
model<-kin.blup(data = G_Pheno_mean,
                geno = "ID",
                pheno = "Trait",
                K = A_matrix)

EBV <- data.frame(ebv = model$g,
                  ID = G_Pheno_mean$ID)

EBV_G4 <- EBV[EBV$ID %in% G4_pop@id,]

G4_pop@ebv <- as.matrix(EBV_G4$ebv)

New_candidates<-selectWithinFam(pop = G4_pop,
                                nInd = 1, # Number Selected within each family
                                use = "ebv")

################################## G5 ##########################################
G_pops <- vector("list", nParents)

for (i in 1:nParents) {
  # Generate random number of progeny for this cross (with mean = 30 and sd = 10)
  Prog_number <- max(round(rnorm(1, mean = nProg, sd = 10)), 1)
  
  # Create a population for this iteration
  G_pops[[i]] <- randCross(
    pop = New_candidates,   # The population to cross
    nCrosses = 1,    # Total number of crosses
    nProgeny = Prog_number,  # Number of progeny per cross (randomly generated)
    balance = NULL,   # If using sexes, this option will balance the number of progeny per parent
    parents = NULL,
    ignoreSexes = TRUE  # Should sexes be ignored?
  )
}

G5_pop <- do.call(c, G_pops)

Ped<-SP$pedigree

Ped_G5 <- as.data.frame(Ped[rownames(Ped) %in% G5_pop@id,])

Ped_G<-rbind(Ped_G0,Ped_G1,Ped_G2,Ped_G3,Ped_G4,Ped_G5)

Pedigree <- data.frame(ID = rownames(Ped_G),
                       Sire = Ped_G$father,
                       Dam = Ped_G$mother)

# Calcualte the A matrix
A_matrix<-2*kinship(id = Pedigree$ID, dadid = Pedigree$Sire, momid = Pedigree$Dam, )

G5_Pheno <- as.data.frame(G5_pop@id)

# Number of replicates
n_reps <- 10

# Loop through and assign phenotypic values
for (i in 1:n_reps) {
  G5_Pheno[[paste0("V", i)]] <- setPheno(G5_pop,
                                         h2 = h2_G1, 
                                         fixEff = 1L,
                                         onlyPheno = TRUE,
                                         traits = 1,
                                         reps = 1,
                                         simParam = NULL)
}

num_values <- prod(dim(G5_Pheno[, -1]))
num_na <- round(0.1 * num_values)
rows <- sample(nrow(G5_Pheno), num_na, replace = TRUE)
cols <- sample(2:ncol(G5_Pheno), num_na, replace = TRUE) 
# Assign NA to selected positions
for (i in seq_along(rows)) {
  G5_Pheno[rows[i], cols[i]] <- NA
}

for (i in seq_along(cols)) {
  G5_Pheno[rows[i], cols[i]] <- NA
}

G5_Pheno$Mean <- rowMeans(G5_Pheno[ , -1], na.rm = T)

G5_pop@pheno <- as.matrix(G5_Pheno$Mean)

G5_Pheno_mean <- data.frame(
  ID = G5_Pheno$`G5_pop@id`,
  Trait = G5_Pheno$Mean)

G_Pheno_mean <- rbind(G0_Pheno_mean, G1_Pheno_mean, G2_Pheno_mean, G3_Pheno_mean, G4_Pheno_mean, G5_Pheno_mean)

# Calcualte EBVs
model<-kin.blup(data = G_Pheno_mean,
                geno = "ID",
                pheno = "Trait",
                K = A_matrix)

EBV <- data.frame(ebv = model$g,
                  ID = G_Pheno_mean$ID)

EBV_G5 <- EBV[EBV$ID %in% G5_pop@id,]

G5_pop@ebv <- as.matrix(EBV_G5$ebv)

##### Images #####
# Plot the data
df <- data.frame(
  VarA = c(varA(G0_pop),varA(G1_pop),varA(G2_pop),varA(G3_pop),varA(G4_pop),varA(G5_pop)),
  Gen = c(0,1,2,3,4,5))

VarA_Image<-ggplot(df, aes(x = Gen, y = VarA)) +
  geom_line() +
  labs(title = "Additive variance",
       x = "Generation",
       y = "Value",
  ) +
  theme_minimal()

filename <- paste0("VarA_gen.", today_date, ".png")
ggsave(filename, VarA_Image, width = 8, height = 6, dpi = 300)


df <- data.frame(
  Bv = c(meanP(G0_pop),meanP(G1_pop),meanP(G2_pop),meanP(G3_pop),meanP(G4_pop),meanP(G5_pop)),
  Gen = c(0,1,2,3,4,5))

meanP_Image<-ggplot(df, aes(x = Gen, y = Bv)) +
  geom_line() +
  labs(title = "Mean Phenotype",
       x = "Generation",
       y = "Value",
  ) +
  theme_minimal()

filename <- paste0("meanPheno_gen.", today_date, ".png")
ggsave(filename, meanP_Image, width = 8, height = 6, dpi = 300)

df <- data.frame(
  Bv = c(meanG(G0_pop),meanG(G1_pop),meanG(G2_pop),meanG(G3_pop),meanG(G4_pop),meanG(G5_pop)),
  Gen = c(0,1,2,3,4,5))

MeanG_Image<-ggplot(df, aes(x = Gen, y = Bv)) +
  geom_line() +
  labs(title = "Mean BV",
       x = "Generation",
       y = "Value",
  ) +
  theme_minimal()

filename <- paste0("meanG_gen.", today_date, ".png")
ggsave(filename, MeanG_Image, width = 8, height = 6, dpi = 300)


############################## Merge Generations ###############################


Data_G0 <- data.frame(Genotype_ID = as.character(G0_pop@iid),
                      Mum = G0_pop@mother,
                      Dad = G0_pop@father,
                      Tbv = G0_pop@gv[,"Trait1"],
                      #Ebv = G0_pop@ebv,
                      Pheno.Mean = G0_pop@pheno,
                      Gen = "G0")

Data_G1 <- data.frame(Genotype_ID = as.character(G1_pop@iid),
                      Mum = G1_pop@mother,
                      Dad = G1_pop@father,
                      Tbv = G1_pop@gv[,"Trait1"],
                      #Ebv = G1_pop@ebv,
                      Pheno.Mean = G1_pop@pheno,
                      Gen = "G1")

Data_G2 <- data.frame(Genotype_ID = as.character(G2_pop@iid),
                      Mum = G2_pop@mother,
                      Dad = G2_pop@father,
                      Tbv = G2_pop@gv[,"Trait1"],
                      #Ebv = G2_pop@ebv,
                      Pheno.Mean = G2_pop@pheno,
                      Gen = "G2")

Data_G3 <- data.frame(Genotype_ID = as.character(G3_pop@iid),
                      Mum = G3_pop@mother,
                      Dad = G3_pop@father,
                      Tbv = G3_pop@gv[,"Trait1"],
                      #Ebv = G3_pop@ebv,
                      Pheno.Mean = G3_pop@pheno,
                      Gen = "G3")

Data_G4 <- data.frame(Genotype_ID = as.character(G4_pop@iid),
                      Mum = G4_pop@mother,
                      Dad = G4_pop@father,
                      Tbv = G4_pop@gv[,"Trait1"],
                      #Ebv = G4_pop@ebv,
                      Pheno.Mean = G4_pop@pheno,
                      Gen = "G4")

Data_G5 <- data.frame(Genotype_ID = as.character(G5_pop@iid),
                      Mum = G5_pop@mother,
                      Dad = G5_pop@father,
                      Tbv = G5_pop@gv[,"Trait1"],
                      #Ebv = G5_pop@ebv,
                      Pheno.Mean = G5_pop@pheno,
                      Gen = "G5")

Data <- rbind(Data_G1, Data_G2, Data_G3, Data_G4, Data_G5)

# Set genotype_id as factor
Data$Genotype_ID <- as.factor(Data$Genotype_ID)

# Create Families
Data$Family <- with(Data, paste(Mum, Dad, sep = "_"))
Data$Family <- as.numeric(as.factor(Data$Family))

Data_G0$Family <- 0
Data <- rbind(Data_G0,Data)

# Add all phenotype data
G0_long <- G0_Pheno %>%
  pivot_longer(cols = c(-`G0_pop@id`, -Mean), names_to = "Replicate", values_to = "Phenotype")
G0_long$Genotype_ID<-G0_long$`G0_pop@id`
G0_long$`G0_pop@id` <- NULL
G0_long$Replicate <- NULL

G1_long <- G1_Pheno %>%
  pivot_longer(cols = c(-`G1_pop@id`, -Mean), names_to = "Replicate", values_to = "Phenotype")
G1_long$Genotype_ID<-G1_long$`G1_pop@id`
G1_long$`G1_pop@id` <- NULL
G1_long$Replicate <- NULL

G2_long <- G2_Pheno %>%
  pivot_longer(cols = c(-`G2_pop@id`, -Mean), names_to = "Replicate", values_to = "Phenotype")
G2_long$Genotype_ID<-G2_long$`G2_pop@id`
G2_long$`G2_pop@id` <- NULL
G2_long$Replicate <- NULL

G3_long <- G3_Pheno %>%
  pivot_longer(cols = c(-`G3_pop@id`, -Mean), names_to = "Replicate", values_to = "Phenotype")
G3_long$Genotype_ID<-G3_long$`G3_pop@id`
G3_long$`G3_pop@id` <- NULL
G3_long$Replicate <- NULL

G4_long <- G4_Pheno %>%
  pivot_longer(cols = c(-`G4_pop@id`, -Mean), names_to = "Replicate", values_to = "Phenotype")
G4_long$Genotype_ID<-G4_long$`G4_pop@id`
G4_long$`G4_pop@id` <- NULL
G4_long$Replicate <- NULL

G5_long <- G5_Pheno %>%
  pivot_longer(cols = c(-`G5_pop@id`, -Mean), names_to = "Replicate", values_to = "Phenotype")
G5_long$Genotype_ID<-G5_long$`G5_pop@id`
G5_long$`G5_pop@id` <- NULL
G5_long$Replicate <- NULL

G_long <- rbind(G0_long, G1_long, G2_long, G3_long, G4_long, G5_long)

Data_merged <- left_join(Data, G_long, by = "Genotype_ID")

# Sort the Data after Genotyoe_ID
Data_merged <- Data_merged %>%
  arrange(by = Genotype_ID)

# Fix the column names
Data_merged$Mean <- NULL

Data_merged<-na.omit(Data_merged)

#### Create G matricies

G_G0<-pullSnpGeno(G0_pop, snpChip = 1)-1
G_G1<-pullSnpGeno(G1_pop, snpChip = 1)-1
G_G2<-pullSnpGeno(G2_pop, snpChip = 1)-1
G_G3<-pullSnpGeno(G3_pop, snpChip = 1)-1
G_G4<-pullSnpGeno(G4_pop, snpChip = 1)-1
G_G5<-pullSnpGeno(G5_pop, snpChip = 1)-1

cleaned_rownames <- gsub('"', '', rownames(G_G0))
G_G0<-G_G0[cleaned_rownames %in% Data$Genotype_ID, ]
cleaned_rownames <- gsub('"', '', rownames(G_G1))
G_G1<-G_G1[cleaned_rownames %in% Data$Genotype_ID, ]
cleaned_rownames <- gsub('"', '', rownames(G_G2))
G_G2<-G_G2[cleaned_rownames %in% Data$Genotype_ID, ]
cleaned_rownames <- gsub('"', '', rownames(G_G3))
G_G3<-G_G3[cleaned_rownames %in% Data$Genotype_ID, ]
cleaned_rownames <- gsub('"', '', rownames(G_G4))
G_G4<-G_G4[cleaned_rownames %in% Data$Genotype_ID, ]
cleaned_rownames <- gsub('"', '', rownames(G_G5))
G_G5<-G_G5[cleaned_rownames %in% Data$Genotype_ID, ]

# Create several G_matricies, as they are too large to compute in R
G01_G_matrix <- A.mat(rbind(G_G0,G_G1), min.MAF = 0.01)
G12_G_matrix <- A.mat(rbind(G_G1,G_G2), min.MAF = 0.01)
G23_G_matrix <- A.mat(rbind(G_G2,G_G3), min.MAF = 0.01)
G34_G_matrix <- A.mat(rbind(G_G3,G_G4), min.MAF = 0.01)
G45_G_matrix <- A.mat(rbind(G_G4,G_G5), min.MAF = 0.01)

G012_G_matrix <- A.mat(rbind(G_G0,G_G1,G_G2), min.MAF = 0.01)
G0123_G_matrix <- A.mat(rbind(G_G0,G_G1,G_G2,G_G3), min.MAF = 0.01)
G01234_G_matrix <- A.mat(rbind(G_G0,G_G1,G_G2,G_G3,G_G4), min.MAF = 0.01)
G012345_G_matrix <- A.mat(rbind(G_G0,G_G1,G_G2,G_G3,G_G4,G_G5), min.MAF = 0.01)

# Marker matrix
filename <- paste0("Sim_G01_G_Matrix_", today_date, ".txt")
write.table(G01_G_matrix, filename, quote = F, col.names = T, row.names = T, sep = "\t")

filename <- paste0("Sim_G12_G_Matrix_", today_date, ".txt")
write.table(G12_G_matrix, filename, quote = F, col.names = T, row.names = T, sep = "\t")

filename <- paste0("Sim_G23_G_Matrix_", today_date, ".txt")
write.table(G23_G_matrix, filename, quote = F, col.names = T, row.names = T, sep = "\t")

filename <- paste0("Sim_G34_G_Matrix_", today_date, ".txt")
write.table(G34_G_matrix, filename, quote = F, col.names = T, row.names = T, sep = "\t")

filename <- paste0("Sim_G45_G_Matrix_", today_date, ".txt")
write.table(G45_G_matrix, filename, quote = F, col.names = T, row.names = T, sep = "\t")


filename <- paste0("Sim_G012_G_Matrix_", today_date, ".txt")
write.table(G012_G_matrix, filename, quote = F, col.names = T, row.names = T, sep = "\t")

filename <- paste0("Sim_G0123_G_Matrix_", today_date, ".txt")
write.table(G0123_G_matrix, filename, quote = F, col.names = T, row.names = T, sep = "\t")

filename <- paste0("Sim_G01234_G_Matrix_", today_date, ".txt")
write.table(G01234_G_matrix, filename, quote = F, col.names = T, row.names = T, sep = "\t")

filename <- paste0("Sim_G012345_G_Matrix_", today_date, ".txt")
write.table(G012345_G_matrix, filename, quote = F, col.names = T, row.names = T, sep = "\t")


# Raw data produced
filename <- paste0("Sim_Data_", today_date, ".txt")
write.table(Data_merged, filename, quote = F, col.names = T, row.names = F, sep = "\t")


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