########################################################################################################
#############DIVERSITY MODEL SCRIPT: CONVENTIONAL:CONSERVATION##########################################
########################################################################################################

#Loading required packages and importing data
library(tidyverse) # for data manipulation and quick data reading and writing
library(metafor)
library(dplyr)
library(esc)

#Need to set working directory to code directory before running this
data <- read.csv("../data/RevisedDataLong.csv", stringsAsFactors = FALSE)

#Cleaning data and preparing for analysis
data <- subset(data, sd1i != "n/a") #removing n/as
data <- subset(data, m1i != "n/a")
data <- subset(data, n1i != "n/a")  

#These columns were originally of the class character, cant pass these to binary operators
data <- transform(data, m1i = as.numeric(m1i))
data <- transform(data, sd1i = as.numeric(sd1i))

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
  
  if (siteone$Management == "Conventional" && sitetwo$Management == "Conservation"){
    n2i <- siteone$n1i
    m2i <- siteone$m1i
    sd2i <- siteone$sd1i
    n1i <- sitetwo$n1i
    m1i <- sitetwo$m1i
    sd1i <- sitetwo$sd1i
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
  #pairs$SoilType[x] <- df$Soil.Type
  #pairs$Length.of.Experiment <- df$Length.of.Experiment
  #pairs$Crop.Type <- df$Crop
  #pairs$Climate <- df$Climate
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

#Constructing a model for the total community diversity between the two conditions
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

#Making a forest plot containing estimates and variances for the community model
forest(estimates, variances, slab="Pooled effect")
