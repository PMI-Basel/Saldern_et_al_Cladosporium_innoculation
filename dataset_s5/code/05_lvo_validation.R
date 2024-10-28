#!/usr/bin/env Rscript

# Lodding of libraries
library(glmulti)
library(dplyr)
library(phyloseq)
library(randomForest)

set.seed(1968)

# Data importation and preparation
load("./mi13_norm.rds")

# loading soil_data
soil_data <- data.frame(sample_data(mi13.norm)) %>% dplyr::select(!c("Silt", "Nitrate", "SOM"))

# loading otus_data
otus <- data.frame(otu_table(mi13.norm))
otus_data <- merge(otus, soil_data %>% dplyr::select(BRI), by = 0) # adding BRI data to otus
rownames(otus_data) = otus_data$Row.names # putting the row names back
otus_data <- otus_data[,-1] # removing column row.names

# combined data
combined_data <- merge(soil_data, otus, by = 0)
rownames(combined_data) = combined_data$Row.names
combined_data = combined_data[,-1]

# Defining the dataframe collecting the predictions
predictions <- data.frame(matrix(ncol = 3))
colnames(predictions) = c("Obs_bri", "Pred_bri", "field")

# Defining the dataframe collecting the predictions
best.predictors <- data.frame(matrix(ncol = 2))
colnames(best.predictors) = c("field", "predictors")

# Function that return the predictors selected by the best model
findBestSoilModel <- function(trainData) {
  # glmulti
  soil.glm <-do.call("glmulti", list(BRI ~., data=trainData, level = 1, plotty = F, report = F, crit = aicc, fitfunction='lm'))
  # compute the model
  print(summary(soil.glm@objects[[1]]))
  predictors <- gsub("BRI ~ 1 \\+ ", "", (summary(soil.glm)$bestmodel)) # getting the selected predictors from the model
  soil.predictors <- strsplit(predictors, " \\+ ")[[1]] # get them into a vector for subsequent analysis
  print(soil.predictors)
  return(soil.predictors)
}

findBestOTUModel <- function(trainData) {
  # Selection of the top 20 OTUs with RF
  rf.asv <- randomForest(BRI ~., 
                         corr.bias = T,
                         importance = T,
                         ntree = 1000,
                         data = trainData)
  
  # get OTU importance
  imp.asv <- importance(rf.asv)
  imp.asv <- data.frame(imp.asv) %>% filter(IncNodePurity > 0) %>% dplyr::arrange(-IncNodePurity) # sort by DESC importance
  top20 <- rownames(imp.asv[1:20,])
  print(top20)
 
  # compute the model
  otus.formula <- formula(paste0("BRI ~", paste0(top20, collapse = "+")))
  otus.lm <- lm(otus.formula, data = trainData)
  
  # exhaustive screening 
  otus.glm <-do.call("glmulti", list(otus.formula, data = trainData, level = 1, plotty = F, report = F, crit = aicc, fitfunction='lm'))
  print(summary(otus.glm@objects[[1]]))
  predictors <- gsub("BRI ~ 1 \\+ ", "", (summary(otus.glm)$bestmodel)) # getting the selected predictors from the model
  otus.predictors <- strsplit(predictors, " \\+ ")[[1]] # get them into a vector for subsequent analysis
  return(otus.predictors)
}

findBestCombined <- function(trainData, soil_predictors, otus_predictors) {
  # build formula from previous selection
  combined.formula <- formula(paste0("BRI ~", paste0(soil_predictors, collapse = "+"), "+", paste0(otus_predictors, collapse ="+")))
  
  # compute the model
  combined.lm <- lm(combined.formula, data = trainData)
  
  # Exhaustive screening 
  combined.glm <-do.call("glmulti", list(combined.formula, data = trainData, level = 1, plotty = F, report = F, crit = aicc, fitfunction='lm'))
 print("glm combined done")
  print(summary(combined.glm@objects[[1]]))
  
  predictors <- gsub("BRI ~ 1 \\+ ", "", (summary(combined.glm)$bestmodel)) # getting the selected predictors from the model
  combined.predictors <- strsplit(predictors, " \\+ ")[[1]] # get them into a vector for subsequent analysis
 print("before return of combine function")  
return(combined.predictors)
}

# For loop for the LVO

for(field in rownames(otus_data)) {
  print(paste0("Starting ", field)) 
  testData <- combined_data[field,]
  
  soil.preds <- findBestSoilModel(soil_data[rownames(soil_data) != field, ]) # getting best soil predictors
  otus.preds <- findBestOTUModel(otus_data[rownames(otus_data) != field, ]) # getting best otus predictors
  combined.preds <- findBestCombined(combined_data[rownames(combined_data) != field, ], soil.preds, otus.preds)
  print("3 functions succesfully applied in the for loop") 
  best.formula <- formula(paste0("BRI ~", paste0(combined.preds, collapse = "+")))
  best.lm <- lm(best.formula, data = combined_data[rownames(combined_data) != field, ])
  summary(best.lm) # to have a look on the model 
 
 print("before saving predictors") 
  # save predictors used for this field
  best.predictors <- rbind(best.predictors, c(field, paste0(combined.preds, collapse = "+")))
  print(best.predictors) 
  # make predictions
  pred <-as.data.frame(predict(best.lm, newdata = testData))
  predictions <- rbind(predictions, c(combined_data[field, 1], pred[1,1])) # save predictions in a file
  print(predictions)
}

save(best.predictors, file = "./best_predictors.rda")
save(predictions, file = "./predictions.rda")




