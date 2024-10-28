#!/usr/bin/env Rscript



library(glmulti)

load(file = "./combined_data.rda")

print("Loaded and ready")
print("Starting with OTUs model")


 combined.lm <- lm(BRI ~ OTU25 + OTU576 + OTU2203 + OTU97 + OTU110 + Ca_H2O + P_AAE + WHC + Corg, data = combined_data)
combined.glm <- glmulti(combined.lm, level = 1, plotty = F, trace = F, crit = aicc)


summary(combined.glm@objects[[1]])
save(combined.glm, file = "./combined_glm.rda")














