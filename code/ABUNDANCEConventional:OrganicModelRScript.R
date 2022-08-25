########################################################################################################
#############ABUNDANCE MODEL SCRIPT: CONVENTIONAL:ORGANIC###############################################
########################################################################################################

#Loading required packages and importing data
library(tidyverse) # for data manipulation and quick data reading and writing
library(metafor)
library(dplyr)
library(esc)

#Need to set working directory to code directory before running this
data <- read.csv("../data/RevisedDataLong.csv", stringsAsFactors = FALSE)

#Cleaning data and preparing for analysis
data <- subset(data, Abundance != "n/a") #removing n/as
data <- subset(data, Variance != "n/a")
data <- subset(data, n1i != "n/a") 

#These columns were originally of the class character, cant pass these to binary operators
data <- transform(data, Abundance = as.numeric(Abundance))
data <- transform(data, Variance = as.numeric(Variance))

#Defining a function to produce all possible combinations of pairs in a dataframe
calc.pairs <- function(data) {
  pairs = t(combn(nrow(data), m = 2))
  pairs <- as.data.frame(pairs)
  colnames(pairs) <- c("SiteOneRowNumber", "SiteTwoRowNumber")
  pairs$LRR <- 0
  pairs$LRR_var <- 0
  return(pairs)
}

#Defining a function to calculate the effect size and variance
calc.effect <- function(siteone, sitetwo) {
  
  if (siteone$Management == "Conventional" && sitetwo$Management == "Organic"){
    n2i <- siteone$n1i
    m2i <- siteone$Abundance
    sd2i <- siteone$Variance
    n1i <- sitetwo$n1i
    m1i <- sitetwo$Abundance
    sd1i <- sitetwo$Variance
    # print("Condition met.")
    
    #Calculating effect sizes with assigned variables
    effectsize <- escalc(measure="ROM", 
                         n1i=n1i,
                         n2i=n2i,
                         m1i=m1i,
                         m2i=m2i,
                         sd1i=sd1i,
                         sd2i=sd2i,
                         var.names=c("LRR", "LRR_var"))
    return(effectsize)
  } else {
    # print("Condition not met.")
    out <- data.frame(LRR=c(0), LRR_var=c(0))
    return(out)
  }
}

#Initialising an empty dataframe to store results in
pairs_final <- data.frame(
  ID=as.character(c()),
  SiteOneRowNumber=as.integer(c()),
  SiteTwoRowNumber=as.integer(c()),
  LRR=as.numeric(c()),
  LRR_var=as.numeric(c())
)

#Calculating all effect sizes
for (x in unique(data$ID)){
  # if (x == "24.1a"){
  #   print(x)
  #   break
  # }
  print(x)
  df <- subset(data, ID == x)
  #Creating a matrix of all possible combinations of pairs in order to do pairwise comparisons on all of the sites
  pairs = t(combn(nrow(df), m = 2))
  #Some more data wrangling
  pairs <- as.data.frame(pairs)
  colnames(pairs) <- c("SiteOneRowNumber", "SiteTwoRowNumber")
  pairs$LRR <- 0
  pairs$LRR_var <- 0
  pairs$ID <- x
  pairs$SoilType[x] <- df$Soil.Type
  pairs$Length.of.Experiment[x] <- df$Length.of.Experiment
  pairs$Climate[x] <- df$Climate
  pairs$Crop.Type[x] <- df$Crop
  #pairs$SoilDepth[x] <- df$Soil.Level
  
  for (i in 1:nrow(pairs)) {
    # print(i)
    #Assigning Paper IDs to variables
    a <- pairs[i,1]
    b <- pairs[i,2]
    # print(a)
    # print(b)
    siteone <- df[a,]
    sitetwo <- df[b,]
    #print(paperone)
    #print(papertwo)
    
    #Inputting variables into calc.effect function and saving the output
    effect.size <- calc.effect(siteone, sitetwo)
    # print(effect.size)
    pairs$LRR[i] <- effect.size$LRR
    pairs$LRR_var[i] <- effect.size$LRR_var
    
    pairs_final <- rbind(pairs_final, pairs)
    
  }
}

#Remove replicates
pairs_final <- pairs_final[!duplicated(pairs_final), ]

#Remove rows that have not satisfied the conditions
pairs_final[pairs_final==0] <- NA
pairs_final <- pairs_final[complete.cases(pairs_final),]

#########################################################################################################################
#################################################TOTALCOMMUNITYMODEL#####################################################
#########################################################################################################################

#Constructing a model for the total community abundance between the two conditions
community.model <- rma.mv(yi=pairs_final$LRR,
                          V=pairs_final$LRR_var,
                          random=list(~ 1 | ID),
                          slab=ID,
                          data=pairs_final)
summary(community.model)

#Making a forest plot of the community model
forest(community.model)

#Pooling estimates and variances from the model to see overall effect
estimates <- c(coef(community.model))
variances <- c(vcov(community.model))

#Making a forest plot containing estimates and variances for the bacteria model
forest(estimates, variances, slab="Pooled effect")

#Running the model again with moderators
community.model2 <- rma.mv(yi=pairs_final$LRR,
                           V=pairs_final$LRR_var,
                           random=list(~ 1 | ID),
                           mods = ~ SoilType,
                           slab=ID,
                           data=pairs_final)
summary(community.model2)

#Forest plot of model with moderators
forest(community.model2, addpred = TRUE, header = TRUE)


#########################################################################################################################
#################################################PROKARYOTEMODEL#########################################################
#########################################################################################################################

#Subsetting by species 
prokaryotes <- c("26", "27", "27.2", "29","32a", "32b", "37.1", "37a", "37b", "37c", "37d", "39.1", "39a", "39b", "39c", "43")

#Deciding what data to include in the model
prokaryotes_subset <- filter(pairs_final, ID %in% prokaryotes)

#Constructing the model for the prokaryote subset
prokaryotes.model <- rma.mv(yi=prokaryotes_subset$LRR,
                            V=prokaryotes_subset$LRR_var,
                            random=list(~ 1 | ID),
                            slab=ID,
                            data=prokaryotes_subset)
summary(prokaryotes.model)

#Making a forest plot of the prokaryotes model
forest(prokaryotes.model)

#Pooling estimates and variances from the model to see overall effect
estimates <- c(coef(prokaryotes.model))
variances <- c(vcov(prokaryotes.model))

#Making a forest plot containing estimates and variances for the prokaryotes model
forest(estimates, variances, slab="Pooled effect")

#Running the model again with moderators
prokaryotes.model2 <- rma.mv(yi=prokaryotes_subset$LRR,
                             V=prokaryotes_subset$LRR_var,
                             random=list(~ 1 | ID),
                             mods = ~ SoilType,
                             slab=ID,
                             data=prokaryotes_subset)
summary(prokaryotes.model2)

#Forest plot of model with moderators
forest(prokaryotes.model2, addpred = TRUE, header = TRUE)


#########################################################################################################################
#################################################EUKARYOTEMODEL##########################################################
#########################################################################################################################

#Subsetting by species 
eukaryotes <- c("26.1", "27.1", "29.1", "32.1a", "32.1b", "39.2", "43.1")

#Deciding what data to include in the model
eukaryotes_subset <- filter(pairs_final, ID %in% eukaryotes)

#Constructing the model for the prokaryote subset
eukaryotes.model <- rma.mv(yi=eukaryotes_subset$LRR,
                           V=eukaryotes_subset$LRR_var,
                           random=list(~ 1 | ID),
                           slab=ID,
                           data=eukaryotes_subset)
summary(eukaryotes.model)

#Making a forest plot of the eukaryotes model
forest(eukaryotes.model)

#Pooling estimates and variances from the model to see overall effect
estimates <- c(coef(eukaryotes.model))
variances <- c(vcov(eukaryotes.model))

#Making a forest plot containing estimates and variances for the eukaryotes model
forest(estimates, variances, slab="Pooled effect")

#Running the model again with moderators
eukaryotes.model2 <- rma.mv(yi=eukaryotes_subset$LRR,
                            V=eukaryotes_subset$LRR_var,
                            random=list(~ 1 | ID),
                            mods = ~ SoilType,
                            slab=ID,
                            data=eukaryotes_subset)
summary(eukaryotes.model2)

#Forest plot of model with moderators
forest(eukaryotes.model2, addpred = TRUE, header = TRUE)

#########################################################################################################################
#################################################BACTERIAMODEL###########################################################
#########################################################################################################################

#Subsetting by species 
bacteria <- c("26", "27", "29", "32a", "32b", "37a", "37b", "37c", "37d", "39.1", "39a", "39b", "39c", "43")

#Deciding what data to include in the model
bacteria_subset <- filter(pairs_final, ID %in% bacteria)

#Constructing the model for the bacterial subset
bacteria.model <- rma.mv(yi=bacteria_subset$LRR,
                         V=bacteria_subset$LRR_var,
                         random=list(~ 1 | ID),
                         slab=ID,
                         data=bacteria_subset)
summary(bacteria.model)

#Making a forest plot of the bacteria model
forest(bacteria.model)

#Pooling estimates and variances from the model to see overall effect
estimates <- c(coef(bacteria.model))
variances <- c(vcov(bacteria.model))

#Making a forest plot containing estimates and variances for the bacteria model
forest(estimates, variances, slab="Pooled effect")

#Running the model again with moderators
bacteria.model2 <- rma.mv(yi=bacteria_subset$LRR,
                          V=bacteria_subset$LRR_var,
                          random=list(~ 1 | ID),
                          mods = ~ SoilType,
                          slab=ID,
                          data=bacteria_subset)
summary(bacteria.model2)

#Forest plot of model with moderators
forest(bacteria.model2, addpred = TRUE, header = TRUE)

#########################################################################################################################
#################################################ARCHAEAMODEL############################################################
#########################################################################################################################

#Subsetting by species 
archaea <- c("27.2", "37.1")

#Deciding what data to include in the model
archaea_subset <- filter(pairs_final, ID %in% archaea)

#Constructing the model for the archaea subset
archaea.model <- rma.mv(yi=archaea_subset$LRR,
                        V=archaea_subset$LRR_var,
                        random=list(~ 1 | ID),
                        slab=ID,
                        data=archaea_subset)
summary(archaea.model)

#Making a forest plot of the archaea model
forest(archaea.model)

#Pooling estimates and variances from the model to see overall effect
estimates <- c(coef(archaea.model))
variances <- c(vcov(archaea.model))

#Making a forest plot containing estimates and variances for the archaea model
forest(estimates, variances, slab="Pooled effect")

#Running the model again with moderators
archaea.model2 <- rma.mv(yi=archaea_subset$LRR,
                         V=archaea_subset$LRR_var,
                         random=list(~ 1 | ID),
                         mods = ~ SoilType,
                         slab=ID,
                         data=archaea_subset)
summary(archaea.model2)

#Forest plot of model with moderators
forest(archaea.model2, addpred = TRUE, header = TRUE)
