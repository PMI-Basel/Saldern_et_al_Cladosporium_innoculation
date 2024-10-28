#!/usr/bin/env Rscript



library(glmulti)

load(file = "./BRI_soil.rda")

print("Loaded and ready")
print("Starting with OTUs model")

soil.lm <- lm(BRI ~ ., data = BRI.soil)
soil.glm <- glmulti(soil.lm, level = 1, plotty = F, trace = F, crit = aicc)


summary(soil.glm@objects[[1]])
save(soil.glm, file = "./soil_glm.rda")














