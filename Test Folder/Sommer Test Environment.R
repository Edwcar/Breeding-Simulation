library(sommer)

# How to set up a multi-trait model

  # Two traits (juvenile and old measurement)
  # Trait-specific intercepts
  # Trait-specific additive variance
  # Trait specific residual variance
# Too much to fit into one model?


##### Generate Data #####
library(AlphaSimR)
library(dplyr)


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
nFounders = 1000   # Number of founders in the population
neFounders = 10   # Effective population size of founders

# Breeding Parameters    
# Size of the G0 population
nG0 = 1000

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

# Genotyping Parameters
# SNP chip for GS
nSNP = 100     # Nr of SNPs per chromosome
minMAF = 0.01  # Minimum MAF for inclusion
# "Imaginary" SNp chip for estimating true inbreeding coefficients
nSNP2 = 300



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

G0_Pheno <- data.frame()

for (i in 1:1) {
  pheno <- setPheno(G0_pop,
                    varE = initVarE,
                    corE = matrix(c(1, CorP,
                                    CorP, 1), nrow = 2),
                    rep = 1,
                    simParam = NULL,
                    onlyPheno = T)
  G0_Pheno <- rbind(
                    G0_Pheno,
                    data.frame(ID = G0_pop@id, Rep = i, Pheno = pheno)
              )
}

G1_pop<-randCross(pop = G0_pop,
          nCrosses = 100,
          nProgeny = 10)

G1_Pheno <- data.frame()

for (i in 1:10) {
  pheno <- setPheno(G1_pop,
                    varE = initVarE,
                    corE = matrix(c(1, CorP,
                                    CorP, 1), nrow = 2),
                    rep = 1,
                    simParam = NULL,
                    onlyPheno = T)
  G1_Pheno <- rbind(
    G1_Pheno,
    data.frame(ID = G1_pop@id, Rep = i, Pheno = pheno)
  )
}

Data_X <- rbind(G0_Pheno,G1_Pheno)

# Convert to long format

Pedigree<-as.data.frame(SP$pedigree)

# Save as data frame
A <- data.frame(ID = as.numeric(rownames(Pedigree)),
                Sire = Pedigree$father,
                Dam = Pedigree$mother)

library(kinship2)
Kinship_matrix<- 2*kinship(id = as.numeric(rownames(Pedigree)),
                           dadid = Pedigree$father,
                           momid = Pedigree$mother)


##### Model 1 ######
Data_X$ID <- as.factor(Data_X$ID)

DataX2 <- stackTrait(data = Data_X, traits = c("Pheno.Trait1","Pheno.Trait2"))


Ai <- solve(Kinship_matrix + diag(1e-4,ncol(Kinship_matrix),ncol(Kinship_matrix)))
Ai <- as(as(as( Ai, "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Ai, 'inverse')=TRUE

DataX2$long=DataX2$long[with(DataX2$long, order(trait)), ]

model<-mmes(fixed = value ~ trait,
            random = ~vsm(ism(ID), Gu = Ai),
            rcov = ~ units,
            data = DataX2$long,
            nIters = 25,
            henderson = T,
            dateWarning = T)

# 162 sec 20000 obs

model2<-mmes(fixed = value ~ trait,
     random = ~vsm(dsm(trait), ism(ID), Gu = Ai),
     rcov = ~vsm(dsm(trait),ism(units)),
     data = DataX2$long,
     nIters = 25,
     henderson = T,
     dateWarning = T)

# 2912 sec 20000 obs (not converged)

model3<-mmes(fixed = value ~ trait,
             random = ~vsm(dsm(trait), ism(ID), Gu = Ai),
             rcov = ~units,
             data = DataX2$long,
             nIters = 25,
             henderson = T,
             dateWarning = T)

# 95 sec for one iteration

model4<-mmes(fixed = value ~ trait,
             random = ~vsm(usm(trait), ism(ID), Gu = Ai),
             rcov = ~units,
             data = DataX2$long,
             nIters = 25,
             henderson = T,
             dateWarning = T)

# 224 sec for one iteration



model$AIC
model2$AIC
model3$AIC
model4$AIC

#### calculate the genetic correlation
 cov2cor(model4$theta[[1]])

vpredict(model, h2 ~ V1 / ( V1 + V2 ) )
vpredict(model2, h2 ~ V1 / ( V1 + V2 ) )
vpredict(model3, h2 ~ V1 / ( V1 + V2 ) )
vpredict(model4, h2 ~ V1 / ( V1 + V2 ) )


BLUP1 <- as.data.frame(model4$u)
BLUP1_t1 <- as.data.frame(BLUP1[0:((nFounders+(length(G1_pop@id)))),])
BLUP1_t2 <- as.data.frame(BLUP1[(1+((nFounders+(length(G1_pop@id))))):(2*(nFounders+(length(G1_pop@id)))),])

BLUP1$Id <- as.numeric(rownames(BLUP1))


BLUP1_t1$Id <- as.numeric(rownames(BLUP1_t1))
BLUP1_t2$Id <- as.numeric(rownames(BLUP1_t2))

BLUP1_t1 <- BLUP1_t1 %>%
  arrange(Id)
BLUP1_t2 <- BLUP1_t2 %>%
  arrange(Id)

BLUP1 <- BLUP1 %>%
  arrange(Id)


G0_t1<-BLUP1_t1[0:100,]
G0_t2<-BLUP1_t1[101:200,]

cor(G0_pop@gv[,1],G0_t1$`BLUP1[0:200, ]`)
plot(G0_pop@gv[,2],G0_t2$`BLUP1[0:200, ]`)


BLUP1_G0 <- BLUP1[0:100,]

cor(G0_pop@gv[,1],BLUP1_G0$V1)
cor(G0_pop@gv[,2],BLUP1_G0$V1)

plot(G0_pop@gv[,2],BLUP1_G0$V1)

BLUP1_G1 <- BLUP1[101:200,]

plot(G1_pop@gv[,1],BLUP1_G1$V1)
cor(G1_pop@gv[,2],BLUP1_G1$V1)

plot(G1_pop@gv[,1],BLUP1_G1$V1)

library(bWGR)

Data_0 <-Data_X
Data_0$Mum<-NULL
Data_0$Dad<-NULL
Data_0$id<-NULL
Data_0$Rep<-NULL
Data_0$Pheno.Trait1 <- as.numeric(Data_0$Pheno.Trait1)
Data_0$Pheno.Trait2 <- as.numeric(Data_0$Pheno.Trait2)

df_means <- Data_0 %>%
  group_by(ID) %>%
  summarise(
    Trait1 = mean(Pheno.Trait1, na.rm = TRUE),
    Trait2 = mean(Pheno.Trait2, na.rm = TRUE)
  )

df_means$ID <- as.character(df_means$ID)

# Reorder df_means based on the rownames of K
df_means_ordered <- df_means[match(rownames(Kinship_matrix), df_means$ID), ]
df_means_ordered$ID<-NULL


y<-as.matrix(df_means_ordered)

library(tictoc)
tic()
new_model<-mkr(y,Kinship_matrix)
toc()

# Genetic correlations
new_model$h2
new_model$GC
new_model$vb
new_model$mu

EBVs<-new_muEBVs<-new_model$hat
G0<-EBVs[1:1000,]
G1<-EBVs[1001:2000,]

cor(G0_pop@gv[,1],G0[,1])
cor(G0_pop@gv[,2],G0[,2])

plot(G0_pop@gv[,1],G0[,1])
plot(G0_pop@gv[,2],G0[,2])

cor(G1_pop@gv[,1],G1[,1])
cor(G1_pop@gv[,2],G1[,2])

plot(G1_pop@gv[,1],G1[,1])
plot(G1_pop@gv[,2],G1[,2])



# Genetic correlations
new_model$h2

new_model$GC


#### Model 2 ####
# Add trait specific residual variance

model<-mmes(fixed = valueS ~ trait-1,
            random = ~vsm(ism(id), Gu = Ai),
            rcov = ~vsm(ism(units, trait)),
            data = Data2$long,
            nIters = 25,
            henderson = T,
            dateWarning = T)

### Example ####
data("DT_wheat")
 rownames(GT_wheat) <- rownames(DT_wheat)
 G <- A.mat(GT_wheat)
 Y <- data.frame(DT_wheat)
 
# # make the decomposition
 UD<-eigen(G) # get the decomposition: G = UDU'
 U<-UD$vectors
 D<-diag(UD$values)# This will be our new 'relationship-matrix'
 rownames(D) <- colnames(D) <- rownames(G)
 X<-model.matrix(~1, data=Y) # here: only one fixed effect (intercept)
 UX<-t(U)%*%X # premultiply X and y by U' 
 UY <- t(U) %*% as.matrix(Y) # multivariate


 
# # dataset for decomposed model
 DTd<-data.frame(id = rownames(G) ,UY, UX =UX[,1])
 DTd$id<-as.character(DTd$id)
 
 names <- c("X1","X2")
 Stacked_DTd<-stackTrait(DTd, names)
 Stacked_DTd$long
# 
 modeld <- mmes(fixed = trait ~ UX - 1, 
               random = ~vsm(ism(id),Gu=D), 
               rcov = ~units,
               data=Stacked_DTd, verbose = T)

 
 head(DTd$id)
 length(unique(DTd$id))  # Ensure it has 599 unique values
 
 min(UD$values)  # Check for the smallest eigenvalue
 max(UD$values)  # Check for the largest eigenvalue
 summary(UD$values)  # Get a quick summary of the eigenvalues
 
 
 
# # dataset for normal model
 DTn<-data.frame(id = rownames(G) , DT_wheat)
 DTn$id<-as.character(DTn$id)
 
 modeln <- mmes(cbind(X1,X2) ~ 1, 
               random = ~vsm(ism(id),Gu=G), 
               rcov = ~vsm(ism(units)),
               data=DTn, verbose = FALSE)
 
# ## compare regular and transformed blups
 plot(x=(solve(t(U)))%*%modeld$U$`u:id`$X2[colnames(D)], 
      y=modeln$U$`u:id`$X2[colnames(D)], xlab="UDU blup",
      ylab="blup")
 
 data(DT_example)
 DT <- DT_example
 head(DT)
 DT2 <- stackTrait(DT, traits = c("Yield","Weight"))
 head(DT2)
 
 
  data(DT_btdata)
  DT <- DT_btdata
  head(DT)
  mix4 <- mmes(tarsus ~ sex,
  random = ~ dam + fosternest,
  rcov=~units,
  data = DT)
  summary(mix4)$varcomp
 # MULTI-TRAIT EXAMPLE
  head(DT)
  DT2 <- stackTrait(DT, traits = c("tarsus","back"))
  head(DT2$long)
  DT2$long=DT2$long[with(DT2$long, order(trait)), ]
 
  mix3 <- mmes(valueS ~ trait:sex - 1, henderson=TRUE,
  random = ~ vsm(usm(trait),ism(dam)) +
  vsm(usm(trait), ism(fosternest)),
  rcov= ~ vsm(dsm(trait),ism(units)),
  data = DT2$long)
 #
 #
  summary(mix3)$varcomp
 # #### calculate the genetic correlation
  cov2cor(mix3$theta[[1]])
  cov2cor(mix3$theta[[2]])
 


