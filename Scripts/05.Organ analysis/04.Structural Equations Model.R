library(lavaan)
library(semPlot)
library(psych)

## 1. Import organ matrix
organ <- read.csv(file.choose(), row.names = 1, check.names = FALSE)

# Import group information
info <- read.delim(file.choose())  # Sample | Group (AR / HC)

common_samples <- intersect(rownames(organ), info$ID)
organ <- organ[common_samples, ]
info  <- info[match(common_samples, info$ID), ]

## 2. Z-score normalization
organ_z <- scale(organ)
organ_z <- as.data.frame(organ_z)

# 3. Group：NC = 0, AR = 1
info$Group_bin <- ifelse(info$Group == "2AR", 1, 0)

dat <- cbind(
  organ_z,
  Group = info$Group_bin
)

##4. Construct SEM 
model2 <- '
Lymphoidtissue ~ Bonemarrow
Lymphoidtissue ~ Gastrointestinal
Hepatobiliary ~ Lymphoidtissue
Hepatobiliary ~ Brain
Hepatobiliary ~ Endocrine
Brain ~ Gastrointestinal
Brain ~ Pancreas
Brain ~ Hepatobiliary
Brain ~ Endocrine
Gastrointestinal ~ Hepatobiliary
Gastrointestinal ~ Pancreas
Gastrointestinal ~ Endocrine
Bonemarrow ~ Brain
Bonemarrow ~ Gastrointestinal
Pancreas ~ Gastrointestinal
Endocrine ~ Lymphoidtissue
'

##5. Model assessment
fit <- sem(model2, data = dat, estimator = "MLR", missing = "fiml", meanstructure = TRUE)
summary(fit, fit.measures = TRUE, standardized = TRUE)