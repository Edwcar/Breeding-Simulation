# AlphaSim Simulation 

# Try shorter burn in after 

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
nFounders = 1000    # Number of founders in the population
neFounders = 1000    # Effective population size of founders

nGenBurnIn <- 100  # Number of descrete burn-in generations
nCrossBurnIn <- 1000 # Number of crosses within each burn-in generation
nProgBurnIn <- 2

# Size of the G0 population
nG0 = 1000

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
initVarG = 0.2 * (SD * initMeanG)^2  # Initial trait genetic variance

# Set to mimik real trail design
# Field trials
nSite_G0 = 4
nSite_F1 = 3


# Clonal Heritability*
h2_G0 = 0.7
h2_F1 = 0.7

# Crosses
nParents = 30       # Number of parents to be crossed
nCross = 35     # Number of families to create from Parents
nProg = 30      # Number of progeny per cross
# Comment: Create ~ 1000 individuals


# SNP chip
nSNP = round(40470/12)     # Nr of SNPs per chromosome
# Comment: Replicate the number of real SNPs after filtering (~ 39600) 

# Mating schemes
MatePlan = "RandCross"

RandDe = 0.05   # Percentage random drift

###### Initiate Founder Population ######
founderPop = runMacs2(nInd = nFounders,         
                      nChr = PaChr,
                      segSites = NULL,
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

VarA_G0 <- varG(G0_pop)
VarA_list <- c(VarA_list, VarA_G0)
MeanA_G0 <- meanG(G0_pop)
MeanA_list <- c(MeanA_list, MeanA_G0)

time <- seq_along(VarA_list)  # Time series

# Create a data frame
df <- data.frame(
  Time = rep(time, 2),  # Repeat time for both parameters
  Value = c(unlist(VarA_list), unlist(MeanA_list)),  # Combine the values of both parameters
  Parameter = rep(c("VarG", "MeanG"), each = length(time))  # Label each parameter
)

# Plot the data
Burn_Image<-ggplot(df, aes(x = Time, y = Value, color = Parameter)) +
  geom_line(linewidth = 0.5) +         # Lines for both parameters
  geom_point(size = 1) + 
  scale_color_manual( values = c("VarG" = "darkgreen", "MeanG" = "darkred")  # Custom colors
  ) +# Points for both parameters
  labs(title = "",
       x = "Generation",
       y = "Value",
       color = "Parameter"
  ) +
  theme_minimal()       

filename <- paste0("Burn_in_effect.", today_date, ".png")
ggsave(filename, Burn_Image, width = 8, height = 6, dpi = 300)

# Randomly assigns eligible SNPs to SNP chip
# Use the latest population as reference
SP$addSnpChip(nSNP, minSnpFreq = 0.01, refPop = G0_pop)

# Can remove after
geno_matrix <- pullSnpGeno(G0_pop)

# Function for facculating MAF
calculate_maf <- function(column) {
  # Total number of alleles
  total_alleles <- 2 * length(column)
  # Count the number of minor alleles (i.e., the count of '1's and '2's)
  minor_allele_count <- sum(column)
  # Calculate allele frequency
  minor_allele_frequency <- minor_allele_count / total_alleles
  # Ensure MAF is the smaller of the two allele frequencies
  maf <- min(minor_allele_frequency, 1 - minor_allele_frequency)
  return(maf)
}

# Get number of segregating sites
denom <- G0_pop@nInd * G0_pop@ploidy

# Efficient calculation of F1 segregating sites
SegSite <- sum(sapply(G0_pop@geno[1:12], length)) / denom

# Apply the function. The 2 indicates that it should be applied column-wise
maf_values <- as.data.frame(apply(geno_matrix, 2, calculate_maf))
maf_values$Frequency <- maf_values$`apply(geno_matrix, 2, calculate_maf)`
rm(geno_matrix)

# Show in density plot
MAF_after<-ggplot(maf_values, aes(x = Frequency)) +
  geom_density(fill = "brown", alpha = 0.5) +
  labs(title = "G0",
       x = "Minor Allele Frequency",
       y = "Density") +
  theme_minimal()

G0_pop<-setPheno(G0_pop,
                 h2 = h2_G0, 
                 fixEff = 1L,
                 onlyPheno = FALSE,
                 traits = 1,
                 reps = nSite_G0,
                 simParam = NULL)

df <- G0_pop@pheno

Pheno_G0<-ggplot(df, aes(x = Trait1)) +
  geom_histogram(binwidth = 5, fill = "brown", color = "black", alpha = 0.7) +
  labs(title = "G0",
       x = "H_10 dm",
       y = "Frequency") +
  theme_minimal()


# Select best plus trees based on phenotype
Parents<-selectInd(G0_pop, nParents, trait = 1, use = "pheno", selectTop = T)

#### Create the F1 population #####

# Cross Parents to make the next generation
F1_pop <- randCross(pop = Parents,   # The population to cross
                    nCrosses = nCross, # Total number of crosses
                    nProgeny = nProg , # Number of progeny per cross
                    balance = NULL,    # if using sexes, this option will balance the number of progeny per parent
                    parents = NULL,
                    ignoreSexes = T)  # should sexes be ignored

F1_pop<-setPheno(F1_pop,
                 h2 = h2_F1,
                 fixEff = 1L,       # Fixed effect to assign the population
                 onlyPheno = FALSE, # Should only phenotypes be returned?
                 traits = 1,
                 reps = nSite_F1,
                 simParam = NULL)

df <- F1_pop@pheno

Pheno_F1<-ggplot(df, aes(x = Trait1)) +
  geom_histogram(binwidth = 5, fill = "forestgreen", color = "black", alpha = 0.7) +
  labs(title = "G1",
       x = "H_10 dm",
       y = "Frequency") +
  theme_minimal()


pheno_image <- grid.arrange(Pheno_G0, Pheno_F1, nrow = 1, ncol = 2)
filename <- paste0("Phenotype_distributions.", today_date, ".png")
ggsave(filename, pheno_image, width = 8, height = 6, dpi = 300)

geno_matrix <- pullSnpGeno(F1_pop)

# Apply the function. The 2 indicates that it should be applied column-wise
maf_values <- as.data.frame(apply(geno_matrix, 2, calculate_maf))
maf_values$Frequency <- maf_values$`apply(geno_matrix, 2, calculate_maf)`
rm(geno_matrix)

# Show in density plot
MAF_F1<-ggplot(maf_values, aes(x = Frequency)) +
  geom_density(fill = "forestgreen", alpha = 0.5) +
  labs(title = "G1",
       x = "Minor Allele Frequency",
       y = "Density") +
  theme_minimal()


maf_image <- grid.arrange(MAF_after, MAF_F1, nrow = 1, ncol = 2)
filename <- paste0("MAF_distributions.", today_date, ".png")
ggsave(filename, maf_image, width = 8, height = 6, dpi = 300)

#### Record genetic data ####

# Record additive genetic variance
VarA_G0<-as.numeric(genicVarA(pop = G0_pop))
VarA_F1<-as.numeric(genicVarA(pop = F1_pop))

MeanA_G0<-meanG(G0_pop)
MeanA_F1<-meanG(F1_pop)

Phen_G0 <- meanP(G0_pop)
Phen_Pa <- meanP(Parents)
Phen_F1 <- meanP(F1_pop)

#### Merge populations ####


Data_G0 <- data.frame(Genotype_ID = G0_pop@iid,
                      Mum = G0_pop@mother,
                      Dad = G0_pop@father,
                      Tbv = G0_pop@gv[, "Trait1"],
                      Pheno = G0_pop@pheno[, "Trait1"],
                      Gen = "G0")

Data_F1 <- data.frame(Genotype_ID = F1_pop@iid,
                      Mum = F1_pop@mother,
                      Dad = F1_pop@father,
                      Tbv = F1_pop@gv[, "Trait1"],
                      Pheno = F1_pop@pheno[, "Trait1"],
                      Gen = "F1")

Data <- rbind(Data_F1, Data_G0)

# Set genotype_id as factor
Data$Genotype_ID <- as.factor(Data$Genotype_ID)

# Create Families
Data$Family <- with(Data, paste(Mum, Dad, sep = "_"))
Data$Family <- as.numeric(as.factor(Data$Family))

# Sort the Data after Genotyoe_ID
Data <- Data %>%
  arrange(by = Genotype_ID)

# Identify genotypes that are parents (mum or dad) in F1
parent_ids <- Parents@id

# Add a new column indicating if the genotype is a parent
Data$IsParent <- ifelse(Data$Genotype_ID %in% parent_ids, 1, 0)

# Calculate the number of rows to remove
n_to_remove <- round(0.05 * nrow(Data))

# Randomly sample rows to remove
rows_to_remove <- sample(1:nrow(Data), n_to_remove)

# Remove the selected rows
Data <- Data[-rows_to_remove, ]

Parents_removed <- length(Parents@id)-sum(Data$IsParent)

G_F1<-pullSnpGeno(F1_pop, snpChip = 1)-1
G_G0<-pullSnpGeno(G0_pop, snpChip = 1)-1

G_all<-rbind(G_G0,G_F1)

cleaned_rownames <- gsub('"', '', rownames(G_all))

G_all<-G_all[cleaned_rownames %in% Data$Genotype_ID, ]

# Marker matrix
filename <- paste0("Sim_M_Matrix_", today_date, ".txt")
write.table(G_all, filename, quote = F, col.names = T, row.names = F, sep = "\t")

# Create the G-matrix
G = A.mat(G_all)
rm(G_all)

# G matrix
filename <- paste0("Sim_G_Matrix_", today_date, ".txt")
write.table(G, filename, quote = F, col.names = T, row.names = F, sep = "\t")
rm(G)

# Raw data produced
filename <- paste0("Sim_Data_", today_date, ".txt")
write.table(Data, filename, quote = F, col.names = T, row.names = F, sep = "\t")

###### Check LD ###### 
# G0
writePlink(G0_pop, baseName="G0_plink", use = "gv")
rm(G0_pop)
system("plink --file G0_plink --r2 --ld-window-r2 0 --ld-window 1500 --ld-window-kb 100000 --out G0_Sim_r2")  

r2<-read.table("G0_Sim_r2.ld", header = T)

r2 <- r2 %>%
  mutate(markerDistance = (abs(BP_A - BP_B)/(1000)))

r2$markerDistance <- r2$markerDistance # kb

# Sample a random subset of rows  to make the picture more clear


# Calculate average LD and save for the simulation parameters
r2_s <- r2[sample(nrow(r2), 100000), ]

distance<-r2_s$markerDistance 
LD.data<-r2_s$R2
n<-1000
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance))) 

LD_Image<-ggplot(r2_s, aes(x = distance)) +
  geom_point(aes(y = R2, color = "Observed LD"), color = "gray", size = 1, alpha = 0.6) +  # Observed data points (size adjusted)
  geom_line(aes(x = distance, y = fpoints), color = "blue") +  
  labs(title = "G0", x = "Distance between SNPs (kb)", y = expression(paste("LD (", r^2, ")"))) +
  theme(
    legend.title = element_blank(),  # Remove legend title
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Black border around the plot
    panel.background = element_rect(fill = "white")  # White background inside the plot
  )

# Define the equation to solve
find_distance <- function(d) {
  # Equation for fpoints
  fpoint <- ((10 + new.rho * d) / ((2 + new.rho * d) * (11 + new.rho * d))) *
    (1 + ((3 + new.rho * d) * (12 + 12 * new.rho * d + (new.rho * d)^2)) /
       (n * (2 + new.rho * d) * (11 + new.rho * d)))
  return(fpoint - 0.2)  # Find where fpoints is approximately 0.2
}

# Use uniroot to find the distance
result <- uniroot(find_distance, lower = min(distance), upper = max(distance))

# Distance at which fpoints is ~0.2
G0_distance_at_0.2 <- result$root

filename <- paste0("G0_Sim_LD.", today_date, ".png")
ggsave(filename, LD_Image, width = 8, height = 6, dpi = 300)

LD_G0<-mean(r2_s$R2)

rm(r2)
rm(r2_s)
rm(distance)
rm(LD.data)
rm(tt)
rm(LD_Image)
gc()

# Check LD for the F1
writePlink(F1_pop, baseName="F1_plink", use = "gv")
rm(F1_pop)
system("plink --file F1_plink --r2 --ld-window-r2 0 --ld-window 1500 --ld-window-kb 100000 --out F1_Sim_r2")
r2<-read.table("F1_Sim_r2.ld", header = T)
r2 <- r2 %>%
  mutate(markerDistance = (abs(BP_A - BP_B)/1000))

r2_s <- r2[sample(nrow(r2), 100000), ]

distance<-r2_s$markerDistance 
LD.data<-r2_s$R2
n<-(nCross*nProg)
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance))) 

LD_Image<-ggplot(r2_s, aes(x = distance)) +
  geom_point(aes(y = R2, color = "Observed LD"), color = "gray", size = 1, alpha = 0.6) +  # Observed data points (size adjusted)
  geom_line(aes(x = distance, y = fpoints), color = "blue") +  
  labs(title = "G1", x = "Distance between SNPs (kb)", y = expression(paste("LD (", r^2, ")"))) +
  theme(
    legend.title = element_blank(),  # Remove legend title
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Black border around the plot
    panel.background = element_rect(fill = "white")  # White background inside the plot
  )

filename <- paste0("F1_Sim_LD.", today_date, ".png")
ggsave(filename, LD_Image, width = 8, height = 6, dpi = 300)

# Define the equation to solve
find_distance <- function(d) {
  # Equation for fpoints
  fpoint <- ((10 + new.rho * d) / ((2 + new.rho * d) * (11 + new.rho * d))) *
    (1 + ((3 + new.rho * d) * (12 + 12 * new.rho * d + (new.rho * d)^2)) /
       (n * (2 + new.rho * d) * (11 + new.rho * d)))
  return(fpoint - 0.2)  # Find where fpoints is approximately 0.2
}

# Use uniroot to find the distance
result <- uniroot(find_distance, lower = min(distance), upper = max(distance))

# Distance at which fpoints is ~0.2
F1_distance_at_0.2 <- result$root

LD_F1<-mean(r2_s$R2)

rm(r2)
rm(r2_s)
rm(distance)
rm(LD.data)
rm(tt)
rm(LD_Image)


###############################################################################
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
                  Parents_removed = Parents_removed,
                  Trait = Trait,
                  Trait.mean = initMeanG,
                  Trait.varA = initVarG,
                  Trials_G1 = nSite_F1,
                  Trials_G0 = nSite_G0,
                  h2_G0 = h2_G0,
                  h2_G1 = h2_F1,
                  MatePlan = MatePlan,
                  Gamma = GAMMA,
                  GammaShape = GShape,
                  SegSites = SegSite,
                  LD_G0 = LD_G0,
                  LD_G1 = LD_F1,
                  LD_G0_dist_0.2 = G0_distance_at_0.2,
                  LD_G1_dist_0.2 = F1_distance_at_0.2,
                  VarA_G0 = VarA_G0,
                  VarA_G1 = VarA_F1,
                  MeanG_G0 = MeanA_G0,
                  MeanG_G1 = MeanA_F1,
                  Phen_G0 = Phen_G0,
                  Phen_Pa = Phen_Pa,
                  Phen_G1 = Phen_F1)

# Simulation Parameters
filename <- paste0("Sim_Parameters_Ne1000", today_date, ".txt")
write.table(Sim, filename, quote = F, col.names = T, row.names = F, sep = "\t")