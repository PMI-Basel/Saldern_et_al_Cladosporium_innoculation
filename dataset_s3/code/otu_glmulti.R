#!/usr/bin/env Rscript



library(glmulti)

load(file = "./BRI_otus.rda")

print("Loaded and ready")
print("Starting with OTUs model")

otus.lm <- lm(BRI ~ OTU25+OTU576+OTU543+OTU241+OTU14+OTU5+OTU2+OTU3+OTU2203+OTU287+OTU17+OTU305+OTU4+OTU97+OTU15+OTU163+OTU110+OTU47+OTU11+OTU10, data = BRI_otus1)
otus.glm <- glmulti(otus.lm, level = 1, plotty = F, trace = F, crit = aicc)


summary(otus.glm@objects[[1]])
save(otus.glm, file = "./otus_glm.rda")














