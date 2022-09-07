# insert data and check for Missing values 

file.choose()
data <- read.csv("D:\\Downloads D\\Ergasia2-2020 (1)\\airpoll.csv")

sum(is.na(data))

#############################
### Best Subset Selection ###
#############################

library(ISLR)
library(leaps)

# Perform best subset selection for 15 variables
BSS <- regsubsets(mort~., data = data, nvmax = 15)
summary(BSS)

# Examine R2, RSS, adjusted R2, Cp, and BIC for the models
Bss.summary <- summary(BSS)
names(Bss.summary)
Bss.summary$rsq

# Plotting RSS, adjusted R2, Cp, and BIC for all of the models at once
par(mfrow=c(2,2))
plot(Bss.summary$rss, xlab = "Number of Variables", ylab = "RSS", type = "l")
plot(Bss.summary$adjr2, xlab = "Number of Variables", ylab = "Adjusted RSq", type = "l")
which.max(Bss.summary$adjr2)
points(10, Bss.summary$adjr2[10], col="Maroon", cex=3, pch=20)
plot(Bss.summary$cp, xlab = "Number of Variables", ylab = "Cp", type = "l")
which.min(Bss.summary$cp)
points(6, Bss.summary$cp[6],col="Maroon", cex=3, pch = 20)
plot(Bss.summary$bic, xlab = "Number of Variables", ylab= "BIC", type = "l")
which.min(Bss.summary$bic)
points(4, Bss.summary$bic[4], col = "Maroon", cex = 3, pch = 20)

# coefficient estimates for the four-variable model (lowest BIC)
coef(BSS, 4)


######## validation set approach 

x <- model.matrix (mort~.,data )[,-16]
y <- data$mort


#split the dataset
set.seed(1)
train.rows <- sample(1:nrow(x), .80 * nrow(x))
x.train <- x[train.rows,]
x.test <- x[-train.rows,]
y.train <- y[train.rows]
y.test <- y[ -train.rows]
xy.train.d <- data[train.rows,]
xy.test.d <- data[-train.rows,]

#perform best subset selection to the training set 
BSS.vld <- regsubsets(mort~., data = xy.train.d , nvmax = 15)


#compute the validation set error for the best model of each model size
#model matrix from the test data
test.mat <- model.matrix(mort~., data = xy.test.d)

val.errors <- rep(NA, 15)
for (i in 1 : 15){
  coefi <- coef(BSS.vld, id = i)
  pred <- test.mat[, names(coefi)]%*%coefi
  val.errors[i] <- mean((xy.test.d$mort-pred)^2)
}

val.errors
which.min(val.errors)
val.errors[4]
coef(BSS.vld,4)

BSS.vld <- regsubsets(mort~., data = data, nvmax = 15)
coef(BSS.vld, 4)

# predict method function
predict.regsubsets <- function(object, newdata, id){
  form <- as.formula(object$call[[2]])
  mat <- model.matrix(form, newdata)
  coefi <- coef(object, id = id)
  xvars <- names(coefi)
  mat[, xvars]%*%coefi
}


########## k-fold cross-validation

k = 10
folds <-sample(1:k, nrow(xy.train.d), replace = TRUE)
cv.errors <- matrix(NA, k, 15, dimnames=list(NULL, paste(1:15)))


for (j in 1:k){
  best.fit <- regsubsets(mort~., data = xy.train.d[folds!=j,], nvmax=15)
  for (i in 1:15){
    pred <- predict.regsubsets(best.fit, xy.test.d[folds==j,], id = i)
    cv.errors[j, i] <- mean((xy.test.d$mort[folds==j]-pred)^2)
  }
}


mean.cv.errors  <- apply(cv.errors,2, mean)
mean.cv.errors
plot(mean.cv.errors, type = "b")
min(mean.cv.errors)
BSS.vld <- regsubsets(mort~., data = data, nvmax = 15)
coef(BSS.vld,6)



#################################
### Ridge Regression & Lasso  ###
#################################

 
# produce a matrix corresponding to the 19 predictors
x <- model.matrix (mort~., data)[,-16]
y <- data$mort

### Ridge Regression
library(glmnet)
grid  <- 10 ^ seq(10,-2, length =100)
Ridge <- glmnet(x, y, alpha =0, lambda =grid)

ridge.pred <- predict(Ridge, s=6, newx=x.test)
mean((ridge.pred - y.test)^2)

dim(coef(Ridge))

# split the data 
set.seed(1)
train.rows <- sample(1:nrow(x), .80 * nrow(x))
x.train <- x[train.rows,]
x.test <- x[-train.rows,]
y.train <- y[train.rows]
y.test <- y[ -train.rows]



# Use cross-validation to choose the tuning parameter 
set.seed(1)
cv.out <- cv.glmnet(x.train, y.train, alpha=0)
plot(cv.out)
bestlam <- cv.out$lambda.min
bestlam

# Test MSE associated with bestlam
ridge.pred <- predict(Ridge, s=bestlam, newx=x.test)
mean((ridge.pred-y.test)^2)


# refit our ridge regression model on the full data set and examine the coefficient estimates
out <- glmnet (x,y,alpha =0)
predict(out,type="coefficients",s=6 )[1:16 ,]



### Lasso 
lasso.mod <- glmnet(x.train, y.train, alpha = 1, lambda = grid)
plot(lasso.mod)

# perform cross-validation and compute the associated test error

set.seed(1)
cv.out <- cv.glmnet(x.train,y.train,alpha =1)
plot(cv.out)
bestlam <- cv.out$lambda.min
bestlam
lasso.pred <- predict(lasso.mod, s=bestlam, newx=x.test)
mean((lasso.pred - y.test)^2)

# The lasso model with lambda chosen by cross-validation contains 
#only eight variables:

out <- glmnet(x, y, alpha=1, lambda = grid)
lasso.coef <- predict(out, type = "coefficients", s = bestlam)[1:16,]
lasso.coef

#######################################
### Principal Components Regression ###
#######################################

library(pls)
set.seed(2)
pcr.fit <- pcr(mort~., data=data ,scale=TRUE ,validation ="CV")

summary(pcr.fit)

# Plot the cross-validation scores
validationplot(pcr.fit, val.type="MSEP")

# Perform PCR on the training data and evaluate its test set performance
set.seed(1)
pcr.fit <- pcr(mort~., data=data, subset =x.train, scale =TRUE, validation ="CV")
validationplot(pcr.fit ,val.type="MSEP")
which.min(MSEP(pcr.fit)$val[1,1,])

# The lowest cross-validation error occurs when 15 component are used
# Test MSE

pcr.pred <- predict(pcr.fit,x.test, ncomp =15)
mean((pcr.pred -y.test)^2)

# Fit PCR on the full data set, using M = 15 the number of components identified by cross-validation.
pcr.fit <- pcr(mort~., data=data, scale = TRUE, ncomp = 15)
summary(pcr.fit)

### Partial least squares

set.seed(1)
pls.fit<- plsr(mort~., data=data ,subset =x.train ,scale=TRUE, validation ="CV")
summary(pls.fit )
validationplot(pls.fit, val.type="MSEP")
which.min(MSEP(pls.fit)$val[1,1,])

pls.pred <- predict(pls.fit, x.test, ncomp = 14)
mean((pls.pred - y.test)^2)




#######################################################################
#######################################################################


datab <- data[,c(1,2,9,14,16)]
data.train <- xy.train.d[ , c(1,2,9,14,16)]
data.test <- xy.test.d[ , c(1,2,9,14,16)]
mean.mort.train <- data.train$mort > mean(data$mort)
mean.mort.test <- data.test$mort > mean(data$mort)
### Linear Discriminant Analysis (LDA)

library(MASS)
lda.fit <- lda(mean.mort.train~prec+temp1+noncauc+so2, data = data.train, family= "gaussian")
lda.fit
lda.pred <- predict(lda.fit, data.test)
names(lda.pred)


lda.class <- lda.pred$class
table(lda.class, mean.mort.test)
mean(lda.class==mean.mort.test)
##The LDA predictions are accurate 83% of the time


## Quadratic Discriminant Analysis (QDA)
qda.fit <- qda(mean.mort.train~prec+temp1+noncauc+so2, data = data.train)
qda.fit
qda.class <- predict(qda.fit,  data.test)$class
table(qda.class, mean.mort.test)
mean(qda.class== mean.mort.test)
##The QDA predictions are accurate almost 84% 


## K-Nearest Neighbors (KNN)
library(class)
train.X <- as.matrix(data.train[,-5])
test.X <- as.matrix(data.test[,-5])

# kNN fit for K = 1
set.seed(123) 
knn.pred <- knn(train.X, test.X, mean.mort.train, k = 1)
table(knn.pred, mean.mort.test)
mean(knn.pred==mean.mort.test)
knn.pred
#Try K = 3
knn.pred <- knn(train.X, test.X, mean.mort.train, k = 3)
table(knn.pred, mean.mort.test)
mean(knn.pred==mean.mort.test)
     

glm.fit <- glm(mean.mort.train ~ prec+temp1+noncauc+so2, data = data.train, family="binomial")
summary(glm.fit)






































































