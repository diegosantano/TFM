library(dplyr)
library(fda)
library("factoextra")
library(refund)
library(fda.usc)
library(qlcMatrix)
library(ggplot2)

# We read the tables resulting from the preprocessing
PersonVarsfil <- readRDS("PersonVarsfil.rds")
AllAct <- readRDS("AllAct.rds")
# We compute again the unique identifiers we have
SEQN <- unique(AllAct$SEQN)

# Functional PCA

# Case taking all the observations (not used in the final work)

Actcols <- grep("^MIN", names(AllAct), value = TRUE)
Actdays <- AllAct[, Actcols]
Actdays <- as.matrix(Actdays)
Actdays <- t(Actdays)


# We create the basis used
bspline_basis <- create.bspline.basis(
  rangeval = c(1 / 60, 1440 / 60),
  nbasis = 41, norder = 4
)
smooth_func <- smooth.basis(
  argvals = (1:1440) / 60, y = Actdays,
  fdParobj = bspline_basis
)

# We perform the pca
pca <- pca.fd(smooth_func$fd,
  nharm = 10,
  harmfdPar = fdPar(smooth_func$fd)
)

# Table with eigenvalues, variability, and cumulative variability
eig <- cbind(pca$values[1:10], pca$varprop, cumsum(pca$varprop))
eig

# We plot the cumulative variability
plot(1:10, eig[, 3],
  pch = 19, col = "blue", type = "b", lwd = 2,
  main = "Accumulated proportion of variability (daily case)",
  xlab = "Number of eigenfunction", ylab = "Value"
)


# We plot the functions of the first 4 components
plot(pca$harmonics[1],
  xlab = "Time", ylab = "Value",
  main = "First FPC (daily case)",
  lty = 1, lwd = 2, col = "blue"
)

plot(pca$harmonics[2],
  xlab = "Time", ylab = "Value",
  main = "Second FPC (daily case)",
  lty = 1, lwd = 2, col = "blue"
)

plot(pca$harmonics[3],
  xlab = "Time", ylab = "Value",
  main = "Third FPC (daily case)",
  lty = 1, lwd = 2, col = "blue"
)

plot(pca$harmonics[4],
  xlab = "Time", ylab = "Value",
  main = "Fourth FPC (daily case)",
  lty = 1, lwd = 2, col = "blue"
)

# We show the observations according to their first 2 scores
plot(pca$scores[, c(1, 2)],
  pch = 20, col = "blue", main = "FPCs scores (daily case)",
  xlab = "First score", ylab = "Second score"
)



# We try with the mean of the day for each person (method 3 in the TFM document)
ActPerson <- readRDS("ActPerson.rds")

# We create the basis
bspline_basis <- create.bspline.basis(
  rangeval = c(1 / 60, 1440 / 60),
  nbasis = 41, norder = 4
)
# Penalization addition
pen_basis <- fdPar(bspline_basis, Lfdobj = 2, lambda = 100)


smooth_func3 <- smooth.basis(
  argvals = (1:1440) / 60, y = ActPerson,
  fdParobj = pen_basis
)

# You can test with any other individual curve, in order to see how this penalization affects the curve
par(mfrow = c(1, 1))
plot(1:1440 / 60, ActPerson[, 100],
  col = "red",
  main = "Physical activity of an individual along with its smoothed curve",
  xlab = "Hour of the day", ylab = "Measure of physical activity"
)
lines(smooth_func3$fd[100])

plot(smooth_func3,
  main = "Smoothed curves",
  xlab = "Hour of the day", ylab = "Measure of physical activity"
)

# We perform the pca
pca3 <- pca.fd(smooth_func3$fd,
  nharm = 10,
  harmfdPar = fdPar(smooth_func3$fd)
)

# Table with eigenvalues, variability, and cumulative variability
eig3 <- cbind(pca3$values[1:10], pca3$varprop, cumsum(pca3$varprop))
eig3



# Plot of the cumulative variability
plot(1:10, eig3[, 3],
  pch = 19, col = "blue", type = "b", lwd = 2,
  main = "Accumulated proportion of variability (mean case)",
  xlab = "Number of eigenfunction", ylab = "Value"
)

# Functions of the first 4 PCs
par(mfrow = c(2, 2))
plot(pca3$harmonics[1],
  xlab = "Time", ylab = "Value",
  main = "First FPC (mean case)",
  lty = 1, lwd = 2, col = "blue"
)

plot(pca3$harmonics[2],
  xlab = "Time", ylab = "Value",
  main = "Second FPC (mean case)",
  lty = 1, lwd = 2, col = "blue"
)

plot(pca3$harmonics[3],
  xlab = "Time", ylab = "Value",
  main = "Third FPC (mean case)",
  lty = 1, lwd = 2, col = "blue"
)

plot(pca3$harmonics[4],
  xlab = "Time", ylab = "Value",
  main = "Fourth FPC (mean case)",
  lty = 1, lwd = 2, col = "blue"
)


# Observations according to their 2 first scores
par(mfrow = c(1, 1))
plot(pca3$scores[, c(1, 2)],
  pch = 20, col = "blue", main = "FPCs scores (mean case)",
  xlab = "First score", ylab = "Second score"
)



# We look at the curves with minimum and maximum scores for the first 3 components

smooth_funcmin1 <- smooth_func3$fd[which.min(pca3$scores[, 1])]


smooth_funcmax1 <- smooth_func3$fd[which.max(pca3$scores[, 1])]

plot(smooth_funcmax1,
  main = "Physical activity of the maximum first score and an intermediate case (mean of the 7 days)",
  xlab = "Hour of the day", ylab = "Measure of physical activity", col = "black"
)
lines(smooth_funcmin1, lty = 1, lwd = 2, col = "red")
legend("topleft",
  legend = c("Minimum", "Maximum"), col = c("red", "black"),
  lty = c(1, 1)
)



smooth_funcmin2 <- smooth_func3$fd[which.min(pca3$scores[, 2])]
smooth_funcmax2 <- smooth_func3$fd[which.max(pca3$scores[, 2])]

plot(smooth_funcmin2,
  main = "Physical activity of the maximum and minimum second score (mean of the 7 days)",
  xlab = "Hour of the day", ylab = "Measure of physical activity", col = "red"
)
lines(smooth_funcmax2, lty = 1, lwd = 2, col = "black")
legend("topleft",
  legend = c("Minimum", "Maximum"), col = c("red", "black"),
  lty = c(1, 1)
)


smooth_funcmin3 <- smooth_func3$fd[which.min(pca3$scores[, 3])]
smooth_funcmax3 <- smooth_func3$fd[which.max(pca3$scores[, 3])]

plot(smooth_funcmin3,
  main = "Physical activity of the maximum and minimum third score (mean of the 7 days)",
  xlab = "Hour of the day", ylab = "Measure of physical activity", col = "red"
)
lines(smooth_funcmax3, lty = 1, lwd = 2, col = "black")
legend("topleft",
  legend = c("Minimum", "Maximum"), col = c("red", "black"),
  lty = c(1, 1)
)



# We now repeat but separating between weekends and not weekends

ActPersonEnd <- readRDS("ActPersonEnd.rds")
ActPersonnotEnd <- readRDS("ActPersonnotEnd.rds")
SEQNEnd <- readRDS("SEQNEnd.rds")

# We create the new basis (same one as previous case)
bspline_basis <- create.bspline.basis(
  rangeval = c(1 / 60, 1440 / 60),
  nbasis = 41, norder = 4
)

pen_basis <- fdPar(bspline_basis, Lfdobj = 2, lambda = 100)


smooth_funcEnd <- smooth.basis(
  argvals = (1:1440) / 60, y = ActPersonEnd,
  fdParobj = pen_basis
)

smooth_funcnotEnd <- smooth.basis(
  argvals = (1:1440) / 60, y = ActPersonnotEnd,
  fdParobj = pen_basis
)

# We compute the PCA of both cases
pcaEnd <- pca.fd(smooth_funcEnd$fd,
  nharm = 10,
  harmfdPar = fdPar(smooth_funcEnd$fd)
)

pcanotEnd <- pca.fd(smooth_funcnotEnd$fd,
  nharm = 10,
  harmfdPar = fdPar(smooth_funcnotEnd$fd)
)


# Tables with the same columns as before (eigenvalues, variability and cumulative variability)
eigEnd <- cbind(pcaEnd$values[1:10], pcaEnd$varprop, cumsum(pcaEnd$varprop))
eigEnd

eignotEnd <- cbind(pcanotEnd$values[1:10], pcanotEnd$varprop, cumsum(pcanotEnd$varprop))
eignotEnd


# Cumulative variability of both cases
plot(1:10, eigEnd[, 3],
  pch = 19, col = "blue", type = "b", lwd = 2,
  main = "Accumulated proportion of variability (Weekend and weekday case)",
  xlab = "Number of eigenfunction", ylab = "Value"
)
lines(1:10, eignotEnd[, 3], pch = 19, col = "black", type = "b", lwd = 2)
legend("topleft",
  legend = c("Weekend", "Weekday"), col = c("blue", "black"),
  lty = c(1, 1)
)

# Functions of the first 4 components for both cases
par(mfrow = c(2, 2))
plot(pcaEnd$harmonics[1],
  xlab = "Time", ylab = "Value",
  main = "First FPC (Weekend and weekday case)",
  lty = 1, lwd = 2, col = "blue"
)
lines(pcanotEnd$harmonics[1], lty = 1, lwd = 2, col = "black")
legend("topleft",
  legend = c("Weekend", "Weekday"), col = c("blue", "black"),
  lty = c(1, 1)
)

plot(pcaEnd$harmonics[2],
  xlab = "Time", ylab = "Value",
  main = "Second FPC (Weekend and weekday case)",
  lty = 1, lwd = 2, col = "blue", ylim = c(-0.5, 0.5)
)
lines(pcanotEnd$harmonics[2], lty = 1, lwd = 2, col = "black")
legend("topleft",
  legend = c("Weekend", "Weekday"), col = c("blue", "black"),
  lty = c(1, 1)
)

plot(pcaEnd$harmonics[3],
  xlab = "Time", ylab = "Value",
  main = "Third FPC (Weekend and weekday case)",
  lty = 1, lwd = 2, col = "blue", ylim = c(-0.4, 0.4)
)
lines(pcanotEnd$harmonics[3], lty = 1, lwd = 2, col = "black")
legend("bottomleft",
  legend = c("Weekend", "Weekday"), col = c("blue", "black"),
  lty = c(1, 1)
)

plot(pcaEnd$harmonics[4],
  xlab = "Time", ylab = "Value",
  main = "Fourth FPC (Weekend and weekday case)",
  lty = 1, lwd = 2, col = "blue", ylim = c(-0.5, 0.5)
)
lines(pcanotEnd$harmonics[4], lty = 1, lwd = 2, col = "black")
legend("topleft",
  legend = c("Weekend", "Weekday"), col = c("blue", "black"),
  lty = c(1, 1)
)

# Didtribution of scores
par(mfrow = c(2, 1))
plot(pcaEnd$scores[, c(1, 2)],
  pch = 20, col = "blue", main = "FPCs scores (Weekend case)",
  xlab = "First score", ylab = "Second score"
)
plot(pcanotEnd$scores[, c(1, 2)],
  pch = 20, col = "blue", main = "FPCs scores (Not weekend case)",
  xlab = "First score", ylab = "Second score"
)


# ANOVA
# We select only the categorical variables
Categories <- select_if(PersonVarsfil, is.factor)

# We create an array where the p-values of each category will be stored
pvalue <- array(1, ncol(Categories))
# We perform the ANOVA for all the categorical variables and store their p-values
# It takes a while with 500 bootstrap samples, reduce if necessary
# Seed for replicability
set.seed(100)
for (i in 1:ncol(Categories)){
  pvalue[i] = fanova.onefactor(fdata(smooth_func3$fd),Categories[,i], nboot = 500)$pvalue
}
print(pvalue)


# Model (means for each person)

summary(PersonVarsfil)
regressors <- PersonVarsfil
# We remove the identifying variables and the mean variable (it makes no sense to include it)
regressors$SEQN <- NULL
regressors$SDDSRVYR <- NULL
regressors$MeanActDay <- NULL


# We only need the numeric representation of BMI
regressors$BMI_cat <- NULL

# We will use the categoric form of the number of drinks
regressors$DrinksPerWeek <- NULL

summary(regressors)


# We need regressors to be numeric
numreg <- select_if(regressors, is.numeric)
catreg <- select_if(regressors, is.factor)
summary(catreg)


# We need to convert these factors into binary variables through dummy indices
# We will need to take one of the categories as the base category for the rest.
# It will be explained which category is selected

# We include all these in a new data frame
dummyvars <- data.frame(
  # Race (White as base)
  RaceMexAm = catreg$Race == "Mexican American",
  RaceOtHis = catreg$Race == "Other Hispanic",
  RaceBlack = catreg$Race == "Black",
  RaceOther = catreg$Race == "Other",

  # Gender (female as base)
  GenderMale = catreg$Gender == "Male",

  # All these will have No as base
  # Diabetes
  DiabetesYes = catreg$Diabetes == "Yes",
  DiabetesBorderline = catreg$Diabetes == "Borderline",

  # CHF
  CHFYes = catreg$CHF == "Yes",

  # CHD
  CHDYes = catreg$CHD == "Yes",

  # Cancer
  CancerYes = catreg$Cancer == "Yes",

  # Stroke
  StrokeYes = catreg$Stroke == "Yes",

  # Education (We take the least possible education as base, 'Less than 9th grade')
  Edu911 = catreg$EducationAdult == "9-11th grade",
  EduHighscGrad = catreg$EducationAdult == "High school grad/GED or equivalent",
  EduCol = catreg$EducationAdult == "Some College or AA degree",
  EduColGrad = catreg$EducationAdult == "College graduate or above",

  # MobilityProblem (No Difficulty as base)
  MobProb = catreg$MobilityProblem == "Any Difficulty",

  # DrinkStatus (Non-Drinker as base)
  DrinkMod = catreg$DrinkStatus == "Moderate Drinker",
  DrinkHeavy = catreg$DrinkStatus == "Heavy Drinker",

  # SmokeCigs (Never as base)
  SmokeFormer = catreg$SmokeCigs == "Former",
  SmokeCurrent = catreg$SmokeCigs == "Current",

  # mortstat (Deceased (1) as base)
  Alive = catreg$mortstat == "0",

  # ucod_leading (Alive left out, categories not being 001, 002 and 010 left out) (check)
  ucod_Cardiovascular = catreg$ucod_leading == "001",
  ucod_Neoplasms = catreg$ucod_leading == "002",
  ucod_Respiratory = catreg$ucod_leading == "003",
  ucod_Cerebrovascular = catreg$ucod_leading == "005",


  # diabetes_mcod (0 as base)
  diab_mcod = catreg$diabetes_mcod == "1"
)


# We group again both tables
regressorsnew <- cbind(numreg, dummyvars)



regressorList <- vector("list", ncol(regressorsnew) + 1)
# Intercept term
regressorList[[1]] <- rep(1, nrow(regressorsnew))

# Rest of the terms
for (i in 2:(ncol(regressorsnew) + 1)) {
  regressorList[[i]] <- as.numeric(regressorsnew[, i - 1])
}

# We take the fd object of the smoothing
meanactfd <- smooth_func3$fd

# We create an fd object for the beta coefficients
betabasis <- create.bspline.basis(
  rangeval = c(1 / 60, 1440 / 60),
  nbasis = 11, norder = 4
)
betafdPar <- fdPar(betabasis)

# We create a list, where all the fdobjects of the beta coefficients will be stored
betaList <- vector("list", ncol(regressorsnew) + 1)
for (j in 1:(ncol(regressorsnew) + 1)) betaList[[j]] <- betafdPar

# We compute the regression through the fRegress function
regresslist <- fRegress(meanactfd, regressorList, betaList)

# We obtain the estimated beta functions
betaListEst <- regresslist$betaestlist

# We plot all of them for visualization
par(mfrow = c(1, 1))
plot(betaListEst[[1]]$fd,
  lwd = 2, col = "red",
  xlab = "Hour of the day",
  ylab = "Physical activity", main = "Intercept"
)

names <- colnames(regressorsnew)
for (i in 1:ncol(regressorsnew)) {
  plot(betaListEst[[i + 1]]$fd,
    lwd = 2, col = "red",
    xlab = "Hour of the day",
    ylab = "Physical activity", main = names[i]
  )
}

# Example with the two diabetes functions
plot(betaListEst[[13]]$fd,
  lwd = 2, col = "black",
  xlab = "Hour of the day",
  ylab = "Physical activity", main = "Diabetes", ylim = c(-50, 10)
)
lines(betaListEst[[14]]$fd, lwd = 2, col = "blue")
legend("top",
  legend = c("Yes", "Borderline"),
  col = c("black", "blue"), lty = c(1, 1)
)




predy <- regresslist$yhatfdobj


predvalues <- eval.fd(seq(1, 1440, 10) / 60, predy)
obsvalues <- eval.fd(seq(1, 1440, 10) / 60, smooth_func3$fd)

res <- obsvalues - predvalues

matplot(seq(1, 1440, 10) / 60, res, type = "l", lty = 1, col = "blue")


MSE <- apply(res^2, 1, mean)
RMSE <- sqrt(MSE)

plot(seq(1, 1440, 10) / 60, RMSE)

SSE <- apply(res^2, 1, sum)

func_mean <- mean.fd(smooth_func3$fd)
plot(func_mean)
meanvalues <- eval.fd(seq(1, 1440, 10) / 60, func_mean)
meanvaluesmat <- meanvalues %*% matrix(1, 1, nrow(PersonVarsfil))
SSY <- apply((obsvalues - meanvaluesmat)^2, 1, sum)

Rsq <- (SSY - SSE) / SSY

plot(seq(1, 1440, 10) / 60, Rsq)


# Models (means for weekend and workdays)

summary(PersonVarsfil)
regressors <- PersonVarsfil
# # We remove SDDSRVYR and the mean variable (it makes no sense to include it)
regressors$SDDSRVYR <- NULL
regressors$MeanActDay <- NULL

# We only need the numeric representation of BMI
regressors$BMI_cat <- NULL

# We will use the categoric form of the number of drinks
regressors$DrinksPerWeek <- NULL

summary(regressors)


# We need regressors to be numeric
numreg <- select_if(regressors, is.numeric)
catreg <- select_if(regressors, is.factor)
summary(catreg)


# We need to convert these factors into binary variables through dummy indices
# The dummy indices are just the same as in the previous model

# We group again both tables
regressorsnew2 <- cbind(numreg, dummyvars)
# Note why we didn't remove SEQN here, as we need it in order to filter the weekend case
regressorsnewEnd <- regressorsnew2[regressorsnew2$SEQN %in% SEQNEnd, ]
# We now just remove the SEQN case so that it isn't used in the model
regressorsnewEnd$SEQN <- NULL

regressorsnewnotEnd <- regressorsnew2
regressorsnewnotEnd$SEQN <- NULL

# We now just repeat the whole procedure we did in the previous model in the two cases we have
regressorListEnd <- vector("list", ncol(regressorsnewEnd) + 1)
regressorListnotEnd <- vector("list", ncol(regressorsnewnotEnd) + 1)


# Intercept term
regressorListEnd[[1]] <- rep(1, nrow(regressorsnewEnd))
regressorListnotEnd[[1]] <- rep(1, nrow(regressorsnewnotEnd))

# Rest of the terms
for (i in 2:(ncol(regressorsnewEnd) + 1)) {
  regressorListEnd[[i]] <- as.numeric(regressorsnewEnd[, i - 1])
  regressorListnotEnd[[i]] <- as.numeric(regressorsnewnotEnd[, i - 1])
}

# We take the fd object of the smoothing
meanactfdEnd <- smooth_funcEnd$fd
meanactfdnotEnd <- smooth_funcnotEnd$fd


betabasis <- create.bspline.basis(
  rangeval = c(1 / 60, 1440 / 60),
  nbasis = 11, norder = 4
)
betafdPar <- fdPar(betabasis)
betaListEnd <- vector("list", ncol(regressorsnewEnd) + 1)
betaListnotEnd <- vector("list", ncol(regressorsnewnotEnd) + 1)
for (j in 1:(ncol(regressorsnewEnd) + 1)) {
  betaListEnd[[j]] <- betafdPar
  betaListnotEnd[[j]] <- betafdPar
}


regresslistEnd <- fRegress(meanactfdEnd, regressorListEnd, betaListEnd)
regresslistnotEnd <- fRegress(meanactfdnotEnd, regressorListnotEnd, betaListnotEnd)

betaListEstEnd <- regresslistEnd$betaestlist
betaListEstnotEnd <- regresslistnotEnd$betaestlist

names <- c("Intercept", colnames(regressorsnewEnd))
plots <- vector("list", ncol(regressorsnewEnd) + 1)

par(mfrow = c(1, 1))
for (i in 1:(ncol(regressorsnewEnd) + 1)) {
  table <- data.frame(
    time = seq(1, 1440, 1) / 60,
    Mean = eval.fd(seq(1, 1440, 1) / 60, betaListEst[[i]]$fd),
    Weekend = eval.fd(seq(1, 1440, 1) / 60, betaListEstEnd[[i]]$fd),
    Weekday = eval.fd(seq(1, 1440, 1) / 60, betaListEstnotEnd[[i]]$fd)
  )

  plots[[i]] <- ggplot(table, aes(x = time)) +
    geom_line(aes(y = Mean), color = "black") +
    geom_line(aes(y = Weekend), color = "red") +
    geom_line(aes(y = Weekday), color = "blue") +
    geom_hline(yintercept = 0, lty = "dashed") +
    theme(
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    ) +
    theme(panel.background = element_rect(fill = NA)) +
    labs(title = names[i], x = "Time of the day (hour)", y = "Physical activity")
  plot(plots[[i]])
}




predyEnd <- regresslistEnd$yhatfdobj[-(length(SEQNEnd) + 1)]
predynotEnd <- regresslistnotEnd$yhatfdobj[-(nrow(PersonVarsfil) + 1)]


predvaluesEnd <- eval.fd(seq(1, 1440, 10) / 60, predyEnd)
obsvaluesEnd <- eval.fd(seq(1, 1440, 10) / 60, smooth_funcEnd$fd)

predvaluesnotEnd <- eval.fd(seq(1, 1440, 10) / 60, predynotEnd)
obsvaluesnotEnd <- eval.fd(seq(1, 1440, 10) / 60, smooth_funcnotEnd$fd)

resEnd <- obsvaluesEnd - predvaluesEnd
resnotEnd <- obsvaluesnotEnd - predvaluesnotEnd

matplot(seq(1, 1440, 10) / 60, resEnd, type = "l", lty = 1, col = "blue")

matplot(seq(1, 1440, 10) / 60, resnotEnd, type = "l", lty = 1, col = "blue")




SSEEnd <- apply(resEnd^2, 1, sum)
SSEnotEnd <- apply(resnotEnd^2, 1, sum)

func_meanEnd <- mean.fd(smooth_funcEnd$fd)
func_meannotEnd <- mean.fd(smooth_funcnotEnd$fd)


meanvaluesEnd <- eval.fd(seq(1, 1440, 10) / 60, func_meanEnd)
meanvaluesmatEnd <- meanvaluesEnd %*% matrix(1, 1, length(SEQNEnd))
SSYEnd <- apply((obsvaluesEnd - meanvaluesmatEnd)^2, 1, sum)
RsqEnd <- (SSYEnd - SSEEnd) / SSYEnd

meanvaluesnotEnd <- eval.fd(seq(1, 1440, 10) / 60, func_meannotEnd)
meanvaluesmatnotEnd <- meanvaluesnotEnd %*% matrix(1, 1, nrow(PersonVarsfil))
SSYnotEnd <- apply((obsvaluesnotEnd - meanvaluesmatnotEnd)^2, 1, sum)
RsqnotEnd <- (SSYnotEnd - SSEnotEnd) / SSYnotEnd

plot(seq(1, 1440, 10) / 60, RsqnotEnd,
  type = "l", col = "blue", ylim = c(0, 0.3),
  xlab = "Time of the day (Hour)", ylab = "Square multiple correlation"
)
lines(seq(1, 1440, 10) / 60, RsqEnd, col = "red")
lines(seq(1, 1440, 10) / 60, Rsq, col = "black")
legend("topleft",
  legend = c("Mean case", "Weekend", "Not weekend"),
  col = c("black", "red", "blue"), lty = c(1, 1)
)


# Clustering (functional)

X <- t(smooth_func3$fd$coefs)

fviz_nbclust(X, kmeans, method = "silhouette", k.max = 5)

kmeans_X <- kmeans(X, centers = 2, iter.max = 1000, nstart = 100)
colors_kmeans_X <- c("skyblue", "red")[kmeans_X$cluster]

par(mfrow = c(1, 1))
plot(smooth_func3,
  xlab = "Hour of the day", ylab = "Measurement of physical activity",
  lty = 1, lwd = 2, col = colors_kmeans_X
)
title("Division in clusterings of the smoothed data")
legend("topleft",
  legend = c("First cluster", "Second cluster"),
  col = c("skyblue", "red"), lty = c(1, 1)
)


plot(smooth_func3$fd[which(kmeans_X$cluster == 1)],
  xlab = "Hour of the day", ylab = "Measurement of physical activity",
  lty = 1, lwd = 2, col = "lightblue"
)
title("First clustering of the smoothed data")


plot(smooth_func3$fd[which(kmeans_X$cluster == 2)],
  xlab = "Hour of the day", ylab = "Measurement of physical activity",
  lty = 1, lwd = 2, col = "red"
)
title("Second clustering of the smoothed data")


plot(pca3$scores[, c(1, 2)],
  pch = 20, col = colors_kmeans_X, main = "FPCs scores (mean case)",
  xlab = "First score", ylab = "Second score"
)
legend("topleft",
  legend = c("First cluster", "Second cluster"),
  col = c("skyblue", "red"), lty = c(1, 1)
)


# Clustering (covariates and scores)


Clusttable <- cbind(regressorsnew, pca3$scores[, c(1, 2, 3, 4)])
fviz_nbclust(Clusttable, kmeans, method = "silhouette", k.max = 5)
kmeans_cov <- kmeans(Clusttable, centers = 2, iter.max = 1000, nstart = 100)

sum(kmeans_cov$cluster == kmeans_X$cluster) / length(kmeans_X$cluster)

colors_kmeans_cov <- c("skyblue", "red")[kmeans_cov$cluster]

plot(smooth_func3,
  xlab = "Hour of the day", ylab = "Measurement of physical activity",
  lty = 1, lwd = 2, col = colors_kmeans_cov
)
title("Division in clusterings of the smoothed data (covariates case)")


plot(smooth_func3$fd[which(kmeans_cov$cluster == 1)],
  xlab = "Hour of the day", ylab = "Measurement of physical activity",
  lty = 1, lwd = 2, col = "lightblue"
)
title("First clustering of the smoothed data (covariates case)")


plot(smooth_func3$fd[which(kmeans_cov$cluster == 2)],
  xlab = "Hour of the day", ylab = "Measurement of physical activity",
  lty = 1, lwd = 2, col = "red"
)
title("Second clustering of the smoothed data (covariates case)")





# Correlation

# For this, we need to use only people who have observations at both
# weekends and work days

# Function computing the correlation for a given ordering
funccor <- function(X, Y) {
  conc <- 0
  for (i in 2:length(X)) {
    for (j in 1:(i - 1)) {
      if ((X[i] > X[j] & Y[i] > Y[j]) | (X[i] < X[j] & Y[i] < Y[j])) {
        conc <- conc + 1
      } else {
        conc <- conc - 1
      }
    }
  }
  corrs <- 2 * conc / (length(X) * (length(X) - 1))
  return(corrs)
}


X1 <- pcaEnd$scores[, 1]
Y1 <- pcanotEnd$scores[, 1]
Y1 <- Y1[PersonVarsfil$SEQN %in% SEQNEnd]


corrs1 <- funccor(X1, Y1)
corrs1

X2 <- pcaEnd$scores[, 2]
Y2 <- pcanotEnd$scores[, 2]
Y2 <- Y2[PersonVarsfil$SEQN %in% SEQNEnd]

corrs2 <- funccor(X2, Y2)
corrs2

X3 <- pcaEnd$scores[, 3]
Y3 <- pcanotEnd$scores[, 3]
Y3 <- Y3[PersonVarsfil$SEQN %in% SEQNEnd]

corrs3 <- funccor(X3, Y3)
corrs3

X4 <- pcaEnd$scores[, 4]
Y4 <- pcanotEnd$scores[, 4]
Y4 <- Y4[PersonVarsfil$SEQN %in% SEQNEnd]

corrs4 <- funccor(X4, Y4)
corrs4


Endpoints <- eval.fd(seq(1, 1440, 1) / 60, smooth_funcEnd$fd)
notEndpoints <- eval.fd(seq(1, 1440, 1) / 60, smooth_funcnotEnd$fd)
notEndpoints <- notEndpoints[, PersonVarsfil$SEQN %in% SEQNEnd]

Xmax <- colMax(Endpoints)
Xmax <- as.vector(Xmax)
Ymax <- colMax(notEndpoints)
Ymax <- as.vector(Ymax)

corrsmax <- funccor(Xmax, Ymax)
corrsmax
