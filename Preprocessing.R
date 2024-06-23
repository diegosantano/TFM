library(rnhanesdata)
library(dplyr)
library(tibble)
library(ggplot2)
library(GGally)


# Create a (local) temporary directory
# where lab measurement (cholesterol, blood presure) data will be
# downloaded from the CDC website and then loaded into R. These files need to be
# downloaded separately as the raw files associated with these lab measurements
# are not included in the rnhanesdata package.
dir_tmp <- tempfile()
dir.create(dir_tmp)

if (!dir.exists(dir_tmp)) {
  dir.create(dir_tmp, showWarnings = FALSE)
}
dl_file <- function(url) {
  bn <- basename(url)
  destfile <- file.path(dir_tmp, bn)
  if (!file.exists(destfile)) {
    out <- download.file(url, destfile = destfile, mode = "wb")
  }
  stopifnot(file.exists(destfile))
}
## download the lab measurement data for the cohort 2003-2004
# Cholesterol - Total & HDL: LBXTC and LBXHDD
dl_file("https://wwwn.cdc.gov/Nchs/Nhanes/2003-2004/L13_C.XPT")
# Blood Pressure: BPXSY1 , BPXSY2, BPXSY3 and BPXSY4
dl_file("https://wwwn.cdc.gov/Nchs/Nhanes/2003-2004/BPX_C.XPT")

## download the lab measurement data for the cohort 2005-2006
# Total Cholesterol: LBXTC
dl_file("https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/TCHOL_D.XPT")
# HDL Cholesterol: LBDHDD
dl_file("https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/HDL_D.XPT")
# Blood Pressure, up to 4 measurements per person: BPXSY1 , BPXSY2, BPXSY3 and BPXSY4
dl_file("https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/BPX_D.XPT")


varnames <- c(
  "LBXTC", "LBXHDD", "LBDHDD", ## 1. cholesterol. Note LBXHDD and LBDHDD are the same variable,
  ##    but different names for 2003-2004 and 2005-2006 cohorts
  "BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4" ## 2. blood pressure measurements
)

## load and merge the lab data
lab_data <- process_covar(varnames = varnames, localpath = dir_tmp)

## change column name for cholesterol variable that changed names
colnames(lab_data$Covariate_C)[colnames(lab_data$Covariate_C) == "LBXHDD"] <- "LBDHDD"

## combine waves
CVMarkers <- bind_rows(lab_data$Covariate_C, lab_data$Covariate_D)

rm(list = c("lab_data", "dir_tmp", "varnames"))


## load the data
data("PAXINTEN_C")
data("PAXINTEN_D") # processed physical activity data 2003-2004 and 2005-2006, respectively
data("Flags_C")
data("Flags_D")
data("Mortality_2015_C")
data("Mortality_2015_D")
data("Covariate_C")
data("Covariate_D")

## re-code activity counts which are considered "non-wear" to be 0
## this doesn't impact many data points, most estimated non-wear times correspond to 0 counts
PAXINTEN_C[, paste0("MIN", 1:1440)] <- PAXINTEN_C[, paste0("MIN", 1:1440)] * Flags_C[, paste0("MIN", 1:1440)]
PAXINTEN_D[, paste0("MIN", 1:1440)] <- PAXINTEN_D[, paste0("MIN", 1:1440)] * Flags_D[, paste0("MIN", 1:1440)]

# Change name for better interpretability
AllAct_C <- PAXINTEN_C
AllAct_D <- PAXINTEN_D

AllFlags_C <- Flags_C
AllFlags_D <- Flags_D

# Merge covariate and mortality data
PersonVarsC <- left_join(Covariate_C, Mortality_2015_C, by = "SEQN")
PersonVarsD <- left_join(Covariate_D, Mortality_2015_D, by = "SEQN")

## clean up the workspace for memory purposes
rm(list = c(paste0(
  c("PAXINTEN_", "Covariate_", "Mortality_2015_", "Flags_"),
  rep(LETTERS[3:4], each = 4)
)))

## combine data for the two waves
AllAct <- bind_rows(AllAct_C, AllAct_D)
AllFlags <- bind_rows(AllFlags_C, AllFlags_D)
PersonVars <- bind_rows(PersonVarsC, PersonVarsD)

# We also need to merge with cardiovascular markers
PersonVars <- left_join(PersonVars, CVMarkers, by = "SEQN")

## clean up the workspace again
rm(list = c(
  "AllAct_C", "AllAct_D", "AllFlags_C", "AllFlags_D", "CVMarkers",
  "PersonVarsC", "PersonVarsD"
))


# It may be interesting to create a variable stating the mean activity per
# day of a person
Actcols <- grep("^MIN", names(AllAct), value = TRUE)
Actdays <- AllAct[, Actcols]

# We take NA values to be 0
Actdays[is.na(Actdays) == TRUE] <- 0

AllAct[, Actcols] <- Actdays
Minsum <- rowSums(AllAct[, Actcols])

for (i in 1:length(PersonVars$SEQN)) {
  PersonVars$MeanActDay[i] <- mean(Minsum[AllAct$SEQN == PersonVars$SEQN[i]])
}

length(unique(AllAct$SEQN))
length(PersonVars$SEQN)

# There are people where no activity registers are found, but its covariates
# are found in the other charts, we remove said observations as it makes no
# sense to keep observations with no response

# Care: we're removing a lot of rows
length(PersonVars$SEQN[PersonVars$MeanActDay == "NaN"])

PersonVarsfil <- PersonVars[!PersonVars$MeanActDay == "NaN", ]

# Let's check the variables we have
summary(PersonVarsfil)

# There are many individuals with all values of physical activity having
# repeated sequences of 32767 values at some point, spoiling the results.
# This must be an error in the measurement instrument used, so they are removed
# As a first option, we remove all variables that have either PAXCAL or PAXSTAT
# not equal to 1
length(AllAct$PAXCAL[AllAct$PAXCAL == 1]) / length(AllAct$PAXCAL)
length(AllAct$PAXSTAT[AllAct$PAXSTAT == 1]) / length(AllAct$PAXSTAT)
AllAct <- AllAct[AllAct$PAXCAL == 1, ]
AllAct <- AllAct[AllAct$PAXSTAT == 1, ]

# We filter the values in PersonVarsfil
SEQN <- unique(AllAct$SEQN)
AllFlags <- AllFlags[AllFlags$SEQN %in% SEQN, ]
PersonVarsfil <- PersonVars[PersonVars$SEQN %in% SEQN, ]

summary(PersonVarsfil)



hist(log(PersonVarsfil$MeanActDay + 1))




# Remove variables

# SEQN not removed as it is needed to compare the covariate with the observations
# in the functional data

# SDDSRVYR states the wave. We keep it by now as it expresses which of the two
# waves we are using

# Variables SDMVPSU, SDMVSTRA, WTINT2YR, WTMEC2YR don't seem relevant
PersonVarsfil$SDMVPSU <- NULL
PersonVarsfil$SDMVSTRA <- NULL
PersonVarsfil$WTINT2YR <- NULL
PersonVarsfil$WTMEC2YR <- NULL

# permth_exm and permth_int in the mortality data also seem irrelevant
PersonVarsfil$permth_exm <- NULL
PersonVarsfil$permth_int <- NULL

# NA treatment, possible transformations and possible further removal of variables

# RIDAGEMN NA: More than 85 years are coded as NA, we substitute the NAs by 85 years (in months)
PersonVarsfil$RIDAGEMN[is.na(PersonVarsfil$RIDAGEMN) == TRUE] <- 85 * 12

# RIDAGEEX NA: Same as the previous
PersonVarsfil$RIDAGEEX[is.na(PersonVarsfil$RIDAGEEX) == TRUE] <- 85 * 12

# These variables provide very similar information, probably one should be removed

# BMI and BMI_cat: we substitute BMI NAs by the median, and BMI_cat by overweight,
# which is the category corresponding to both the median and the mean
PersonVarsfil$BMI[is.na(PersonVarsfil$BMI) == TRUE] <- median(PersonVarsfil$BMI, na.rm = TRUE)
PersonVarsfil$BMI_cat[is.na(PersonVarsfil$BMI_cat) == TRUE] <- "Overweight"
# In the model, we should only use one of these two variables

# From here onwards, variables start to have tons of NAs

test <- PersonVarsfil[is.na(PersonVarsfil$CHF) == FALSE, ]
summary(test)

# We see that the observations with NAs are the same for CHF, CHD, Cancer, Stroke, EducationAdult
# and MobilityProblem, and the reason for that seem to be the age, as all the observations taken
# for minor people have been removed when removing the NA observations. Most of the NAs of the other covariate
# variables (DrinkStatus, DrinksPerWeek and SmokeCigs) have been removed as well, so they seem
# to have a similar restriction. Let's try with DrinkStatus
test <- PersonVarsfil[is.na(PersonVarsfil$DrinkStatus) == FALSE, ]
summary(test)

# Checking in the documentation, we see that in fact, only data of individuals over 20 years 
# old has been published, and refused or don't know information have been taken as NA as well. Therefore, 
# we just remove these observations as we don't have info to consider people smaller than 20 years old
PersonVarsfil <- PersonVarsfil[is.na(PersonVarsfil$CHF) == FALSE, ]
PersonVarsfil <- PersonVarsfil[is.na(PersonVarsfil$DrinkStatus) == FALSE, ]
PersonVarsfil <- PersonVarsfil[is.na(PersonVarsfil$SmokeCigs) == FALSE, ]



# We can also remove the few NAs in mortstat now (only 7), that are due to being in category 3
# in eligstat (non eligible)
PersonVarsfil <- PersonVarsfil[is.na(PersonVarsfil$mortstat) == FALSE, ]

# eligstat is useless now as all the observations left are eligible (all observations have value 1)
sum(PersonVarsfil$eligstat == 1) / nrow(PersonVarsfil)
PersonVarsfil$eligstat <- NULL


# However, we have NA problems with the mortality causes measurements

# This seems to be the case because these patients remain alive while these variables are implied for dead patients
# Let's convert all these into factors, as well as change the values for alive people to alive

PersonVarsfil$mortstat <- factor(PersonVarsfil$mortstat)

PersonVarsfil$ucod_leading[PersonVarsfil$mortstat == "0"] <- "Alive"
PersonVarsfil$ucod_leading <- factor(PersonVarsfil$ucod_leading)

PersonVarsfil$diabetes_mcod[PersonVarsfil$mortstat == "0"] <- "Alive"
PersonVarsfil$diabetes_mcod <- factor(PersonVarsfil$diabetes_mcod)

PersonVarsfil$hyperten_mcod[PersonVarsfil$mortstat == "0"] <- "Alive"
PersonVarsfil$hyperten_mcod <- factor(PersonVarsfil$diabetes_mcod)

summary(PersonVarsfil)

# We see that now only one NA is found in each of these, indeed explaining that this is the source of the NAs
# We can just remove this observation
PersonVarsfil <- PersonVarsfil[is.na(PersonVarsfil$ucod_leading) == FALSE, ]

# We still have a problem with the laboratory measurements, which has a lot of NAs
# The problem here is that 4 attempts were made at recording the blood pressure, with some of them failing,
# hence the NAs. We can just take the mean of the 4 observations for each patient, excluding NAs

BPXSY <- rowMeans(PersonVarsfil[, c("BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4")],
  na.rm = TRUE
)

PersonVarsfil <- PersonVarsfil %>% add_column(BPXSY, .after = "BPXSY4")

PersonVarsfil$BPXSY1 <- NULL
PersonVarsfil$BPXSY2 <- NULL
PersonVarsfil$BPXSY3 <- NULL
PersonVarsfil$BPXSY4 <- NULL

# The few NAs left can now be removed
PersonVarsfil <- PersonVarsfil[is.na(PersonVarsfil$LBDHDD) == FALSE, ]
PersonVarsfil <- PersonVarsfil[is.na(PersonVarsfil$BPXSY) == FALSE, ]

# It is also nice to remove the people with Don't know or refuse answers, to avoid future
# problems modelling. There remains few variables with this problem so it shouldn't be too problematic
PersonVarsfil <- PersonVarsfil[PersonVarsfil$CHD != "Don't know", ]
PersonVarsfil <- PersonVarsfil[PersonVarsfil$CHF != "Don't know", ]
PersonVarsfil <- PersonVarsfil[PersonVarsfil$Cancer != "Don't know", ]
PersonVarsfil <- PersonVarsfil[PersonVarsfil$Stroke != "Don't know", ]
PersonVarsfil <- PersonVarsfil[PersonVarsfil$EducationAdult != "Don't know", ]
PersonVarsfil <- PersonVarsfil[PersonVarsfil$EducationAdult != "Refused", ]
PersonVarsfil <- PersonVarsfil[PersonVarsfil$Diabetes != "Don't know", ]

# We drop all the unused levels
PersonVarsfil <- droplevels(PersonVarsfil)

summary(PersonVarsfil)


# We filter the values in AllAct and AllFlags so that they have the same obs. as in PersonVarsfil
SEQN <- unique(PersonVarsfil$SEQN)
AllAct <- AllAct[AllAct$SEQN %in% SEQN, ]
AllFlags <- AllFlags[AllFlags$SEQN %in% SEQN, ]


# We make a last filter, selecting only the days with more than 10 hours
Flagcols <- grep("^MIN", names(AllFlags), value = TRUE)
Flagdays <- AllFlags[, Flagcols]
Flagdays[is.na(Flagdays) == TRUE] <- 0
Flagtime <- rowSums(Flagdays) / 60
mean(Flagtime)
length(Flagtime[Flagtime > 10]) / length(Flagtime)

AllFlags <- AllFlags[Flagtime > 10, ]
AllAct <- AllAct[Flagtime > 10, ]

# Again, we repeat the filter in the scalar variables
SEQN <- unique(AllAct$SEQN)
PersonVarsfil <- PersonVarsfil[PersonVarsfil$SEQN %in% SEQN, ]

# To have enough significant data of each person, we select only people with more
# than 4 valid days
PersonVarsfil <- PersonVarsfil[table(AllAct$SEQN) >= 4, ]
# Repeat for the other two tables
SEQN <- unique(PersonVarsfil$SEQN)
AllAct <- AllAct[AllAct$SEQN %in% SEQN, ]
AllFlags <- AllFlags[AllFlags$SEQN %in% SEQN, ]

# Let's rename some of the variables that may be hard to understand
names(PersonVarsfil)[names(PersonVarsfil) == "BPXSY"] <- "sysBP"
names(PersonVarsfil)[names(PersonVarsfil) == "LBXTC"] <- "TotChol"
names(PersonVarsfil)[names(PersonVarsfil) == "LBDHDD"] <- "HDLChol"

# Let's make a last check, this time looking for possible outliers in the continuous variables
summary(PersonVarsfil)

# There are two variables with clear outliers, BMI and DrinksPerWeek.
# DrinksPerWeek have several very big amounts. However, they will be kept as there seems
# to be enough to be considered just a small group of people who consume a lot of drinks
sort(PersonVarsfil$DrinksPerWeek, decreasing = TRUE)[1:10]

# BMI, on the contrary, is just one observation which is clearly above the rest
sort(PersonVarsfil$BMI, decreasing = TRUE)[1:10]

# The best thing seems to just remove this observation, as even if it's not wrong, it is
# so big of an outlier that it may alter a lot future models
PersonVarsfil <- PersonVarsfil[-which.max(PersonVarsfil$BMI), ]

# Once again, apply this filter to the two other tables
SEQN <- unique(PersonVarsfil$SEQN)
AllAct <- AllAct[AllAct$SEQN %in% SEQN, ]
AllFlags <- AllFlags[AllFlags$SEQN %in% SEQN, ]


# We look for multicollinearity in the continuous variables
numeric <- select_if(PersonVarsfil, is.numeric)
numeric <- numeric[, c(-1, -2)]

ggcorr(numeric, label = T)

# Not high correlations as expected, except for the variables that we suspected in the previous analysis
# that should be very similar, that indeed have correlation of 1, therefore, we remove the variables
# that provide repeated information

PersonVarsfil$RIDAGEEX <- NULL
PersonVarsfil$RIDAGEYR <- NULL

# We also rename the now unique age variable
names(PersonVarsfil)[names(PersonVarsfil) == "RIDAGEMN"] <- "AgeMonths"


# Let's make also a dataframe with categorical variables (+ MeanActDay)

factorvars <- select_if(PersonVarsfil, is.factor)
factorvars$MeanActDay <- PersonVarsfil$MeanActDay

# Let's compare some variables that we can suspect to be very similar

# CHF and CHD
sum(factorvars$CHF == factorvars$CHD) / nrow(factorvars)
ggplot(factorvars, aes(x = CHF, fill = CHD)) +
  geom_bar(position = "dodge")
# The extreme similarity is due to the extreme disbalance towards the No categories,
# probably should be kept because the Yes cases in each seem to differ

# CHD and stroke
sum(factorvars$CHD == factorvars$Stroke) / nrow(factorvars)
ggplot(factorvars, aes(x = Stroke, fill = CHD)) +
  geom_bar(position = "dodge")
# Similar to the previous case, the Yes cases do seem to differ for both cases

# diabetes_mcod and hyperten_mcod
sum(factorvars$diabetes_mcod == factorvars$hyperten_mcod) / nrow(factorvars)
### Exactly the same, remove one
PersonVarsfil$hyperten_mcod <- NULL

# Let's make some boxplots and histograms looking at the differences by categories

# Gender
ggplot(factorvars, aes(x = Gender, y = log(MeanActDay + 1), fill = Gender)) +
  geom_boxplot() +
  theme(panel.background = element_rect(fill = NA)) +
  labs(x = "Gender", y = "Activity per day")

ggplot(factorvars, aes(log(MeanActDay + 1))) +
  geom_density(aes(group = Gender, colour = Gender, fill = Gender), alpha = 0.1) +
  theme(panel.background = element_rect(fill = NA)) +
  labs(x = "Activity per day", y = "Density")
# Slightly higher activity for men

# BMI_cat
ggplot(factorvars, aes(x = BMI_cat, y = log(MeanActDay + 1), fill = BMI_cat)) +
  geom_boxplot() +
  theme(panel.background = element_rect(fill = NA)) +
  labs(x = "BMI_cat", y = "Activity per day")

ggplot(factorvars, aes(log(MeanActDay + 1))) +
  geom_density(aes(group = BMI_cat, colour = BMI_cat, fill = BMI_cat), alpha = 0.1) +
  theme(panel.background = element_rect(fill = NA)) +
  labs(x = "Activity per day", y = "Density")
# Overweight and obese seem to have less average activity


# DrinkStatus
ggplot(factorvars, aes(x = DrinkStatus, y = log(MeanActDay + 1), fill = DrinkStatus)) +
  geom_boxplot() +
  theme(panel.background = element_rect(fill = NA)) +
  labs(x = "DrinkStatus", y = "Activity per day")

ggplot(factorvars, aes(log(MeanActDay + 1))) +
  geom_density(aes(group = DrinkStatus, colour = DrinkStatus, fill = DrinkStatus), alpha = 0.1) +
  theme(panel.background = element_rect(fill = NA)) +
  labs(x = "Activity per day", y = "Density")
# No clear conclusions, but heavy drinkers surprisingly seem to have higher physical activity.
# The relation seems non-linear, which may suggest that using this variable in the model is better
# than using the numeric version.

# Cancer
ggplot(factorvars, aes(x = Cancer, y = log(MeanActDay + 1), fill = Cancer)) +
  geom_boxplot() +
  theme(panel.background = element_rect(fill = NA)) +
  labs(x = "Cancer", y = "Activity per day")

ggplot(factorvars, aes(log(MeanActDay + 1))) +
  geom_density(aes(group = Cancer, colour = Cancer, fill = Cancer), alpha = 0.1) +
  theme(panel.background = element_rect(fill = NA)) +
  labs(x = "Activity per day", y = "Density")
# More daily activity for people without cancer as expected

# CHD
ggplot(factorvars, aes(x = CHD, y = log(MeanActDay + 1), fill = CHD)) +
  geom_boxplot() +
  theme(panel.background = element_rect(fill = NA)) +
  labs(x = "CHD", y = "Activity per day")

ggplot(factorvars, aes(log(MeanActDay + 1))) +
  geom_density(aes(group = CHD, colour = CHD, fill = CHD), alpha = 0.1) +
  theme(panel.background = element_rect(fill = NA)) +
  labs(x = "Activity per day", y = "Density")
# Again more daily activity for people without CHD


# We save the results of the preprocesing
saveRDS(PersonVarsfil, "PersonVarsfil.rds")
saveRDS(AllAct, "AllAct.rds")
