library(glmnet)
library(Metrics)
library(ModelMetrics)
library(caret)

#Import multi-omics data (121 features with group information)
data <- read.delim(file.choose(), header = T, row.names = 1)

#Training and Test dataset
set.seed(123)
d_index <- createDataPartition(data$Group,p = 0.7)
train_d <- data[d_index$Resample1,]
test_d <- data[-d_index$Resample1,]

#Training LASSO Model
lambdas <- seq(0, 2, length.out = 100)
X <- as.matrix(train_d[,1:13])
Y <- train_d[,14]
set.seed(123)
lasso_model <- cv.glmnet(X,Y, alpha = 1, lambda = lambdas, nfolds =5)

#Plotting 
plot(lasso_model)

plot(lasso_model$glmnet.fit, "lambda", label = T)

#Select the optimal regularization parameters, train the optimal model, and output the Lasso coefficients of each variable.
lasso_min <- lasso_model$lambda.min
lasso_best <- glmnet(X, Y, alpha = 1,lambda = lasso_min)
write.csv(as.matrix(coef(lasso_best)), "Lasso_feature.csv")