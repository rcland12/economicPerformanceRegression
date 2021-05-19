setwd("~/Documents/school/Semester 10/Applied Regression Analysis 7134 - A/Final")

# Load data
data <- read.csv("EconomicPerformance.csv", header = FALSE)[,-1]
colnames(data) <- c("ROR", "GM", "GOS", "OGMR", "CFR", "EBITM", "TDT")

# Simple linear regression analysis
x <- as.matrix(data[,6])
y <- data[,1]
plot(x, y, main = "EBIT Margin vs. Economic RoR")

# Calculate parameters using definitions
n <- length(x)
Sxx <- sum(x^2)-1/n*sum(x)^2
Sxy <- sum(x*y)-1/n*sum(x)*sum(y)
beta1 <- Sxy/Sxx
beta0 <- mean(y)-beta1*mean(x)
beta0+beta1*0.1
abline(beta0,beta1)

# ANOVA table test
residuals <- y-(beta0+beta1*x)
SSE <- sum(residuals^2)
SST <- sum((y-mean(y))^2)
SSR <- SST-SSE
MSR <- SSR
MSE <- SSE/(n-2)
Fstat <- MSR/MSE
Pvalue <- 1-pf(Fstat, 1, n-2)

# Confidence interval for significance of regression
s_eBeta1 <- sqrt(MSE/Sxx)
critBeta1 <- qt(0.995, n-2)
CIBeta1 <- c(beta1-critBeta1*s_eBeta1, beta1+critBeta1*s_eBeta1)

# Making predictions and confidence interval bounds
beta0+beta1*(0.1)
val <- seq(0, 0.2, 0.01)
s_eY1 <- sqrt((1/n+(val-mean(x))^2/Sxx)*MSE)
s_eY2 <- sqrt((1+1/n+(val-mean(x))^2/Sxx)*MSE)
upper1 <- (beta0+beta1*val)+qt(0.995,n-2)*s_eY1
lower1 <- (beta0+beta1*val)-qt(0.995,n-2)*s_eY1
upper2 <- (beta0+beta1*val)+qt(0.995,n-2)*s_eY2
lower2 <- (beta0+beta1*val)-qt(0.995,n-2)*s_eY2
lines(val, upper1, type = "l", lty = 2)
lines(val, lower1, type = "l", lty = 2)
lines(val, upper2, type = "l", lty = 3)
lines(val, lower2, type = "l", lty = 3)
legend(0, 0.285, legend=c("ERoR", "Future ERoR"), lty=2:3, cex = 0.8)

# R-squared value
R2 <- SSR/SST



# Multiple Linear Regression
X <- as.matrix(data[,-c(1)])
Xs <- as.matrix(cbind(1, scale(data[,-c(1)])))
Y <- data[,1]
betaStandardized <- solve(t(Xs) %*% Xs) %*% t(Xs) %*% Y
beta <- rep(0,7)
beta[1] <- betaStandardized[1] - sum(betaStandardized[2:(7)]*colMeans(X)/apply(X, 2, sd))
beta[2:(7)] <- betaStandardized[2:(7)]/apply(X, 2, sd)
ERROR <- solve(t(X) %*% X)           ### SINGULARITY ERROR

# Multicollinearity
sink('cor.txt')
diag(solve(cor(data[,-1])))
sink()

# Rebuild the model
X <- as.matrix(data[,-c(1,4,7)])
Xs <- as.matrix(cbind(1, scale(data[,-c(1,4,7)])))
Y <- data[,1]
betaStandardized <- solve(t(Xs) %*% Xs) %*% t(Xs) %*% Y
beta <- rep(0,5)
beta[1] <- betaStandardized[1] - sum(betaStandardized[2:(5)]*colMeans(X)/apply(X, 2, sd))
beta[2:(5)] <- betaStandardized[2:(5)]/apply(X, 2, sd)

# ANOVA table
n <- length(Y)
SSR <- t(beta) %*% t(cbind(1, X)) %*% Y - sum(Y)^2/n
SSE <- t(Y) %*% Y-beta %*% t(cbind(1, X)) %*% Y
SST <- SSR + SSE
MSR <- SSR/4
MSE <- SSE/(n-5)
F0 <- MSR/MSE
pValue <- 1-pf(F0, 4, n-5)

# R-squared 
R2 <- SSR/SST
R2adj <- 1-MSE/SST*(n-1)

# Automatic model
MLRM <- lm(Y ~ X)

sink('test.txt')
summary(MLRM)$coefficients
sink()

# Extra sum of squares
EEX1 <- as.matrix(data[,-1])
EEmodel1 <- lm(Y ~ EEX1)
SSR1 <- anova(EEmodel1)[1,2]
EEX2 <- as.matrix(data[,-c(1,2,3,6)])
EEmodel2 <- lm(Y ~ EEX2)
SSR2 <- anova(EEmodel2)[1,2]
SSR3 <- SSR1-SSR2
MSR1 <- SSR3/3
EEF <- MSR1/anova(EEmodel1)[2,3]
result <- 1-pf(EEF, 3, 21)

# Simultaneous confidence interval
X <- as.matrix(data[,-c(1,4,5,7)])
newModel <- lm(Y ~ X)
critVal <- qt(1-(0.05/2/4), n-4)
SCIlower <- newModel$coefficients-critVal*summary(newModel)$coefficients[,2]
SCIupper <- newModel$coefficients+critVal*summary(newModel)$coefficients[,2]
SCI <- cbind(SCIlower, SCIupper)

sink('SCI.txt')
SCI
sink()

# Checking model assumptions
sink('newModel.txt')
summary(newModel)
sink()

plot(log(data[,2]),asin(sqrt(data[,1])),xlab = "Good Money",ylab = "Rate of Return")
plot(data[,3],asin(sqrt(data[,1])),xlab = "Gross Operating Surplus",ylab = "Rate of Return")
plot(data[,6]^2,asin(sqrt(data[,1])),xlab = "EBIT Margin",ylab = "Rate of Return")

plot(data[,2], rstudent(newModel),xlab = "Good Money", ylab = "Studentized Residuals")
plot(data[,3], rstudent(newModel),xlab = "Gross Operating Surplus", ylab = "Studentized Residuals")
plot(data[,4], rstudent(newModel),xlab = "EBIT Margin", ylab = "Studentized Residuals")

qqnorm(rstudent(newModel))
qqline(rstudent(newModel))


# Applying a transformation
nX <- cbind(log(data[,2]), data[,6])
nY <- asin(sqrt(data[,1]))
nmodel <- lm(nY ~ nX)
summary(nmodel)

sink('nmodel.txt')
summary(nmodel)
sink()

plot(log(data[,2]),asin(sqrt(data[,1])),xlab = "Good Money",ylab = "Rate of Return")
plot(data[,3],asin(sqrt(data[,1])),xlab = "Gross Operating Surplus",ylab = "Rate of Return")
plot(data[,6]^2,asin(sqrt(data[,1])),xlab = "EBIT Margin",ylab = "Rate of Return")

plot(data[,2], rstudent(nmodel),xlab = "Good Money", ylab = "Studentized Residuals")
plot(data[,3], rstudent(nmodel),xlab = "Gross Operating Surplus", ylab = "Studentized Residuals")
plot(data[,4], rstudent(nmodel),xlab = "EBIT Margin", ylab = "Studentized Residuals")

qqnorm(rstudent(nmodel))
qqline(rstudent(nmodel))

resd <- cbind(nmodel$residuals, rstandard(nmodel), rstudent(nmodel))
colnames(resd) <- c("Residuals", "Standardized", "Studentized")

# Residuals of transformed model
sink('resid.txt')
resd
sink()



