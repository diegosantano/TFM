# We read our two main tables
PersonVarsfil <- readRDS("PersonVarsfil.rds")
AllAct <- readRDS("AllAct.rds")
# We calculate the identifiers we have
SEQN <- unique(AllAct$SEQN)

# Mean case
# We extract only the physical activity columns
Actcols <- grep("^MIN", names(AllAct), value = TRUE)
holder <- AllAct[, Actcols]

# We execute the for loop which executes the mean for each identifier (it takes a while)
ActPerson <- matrix(0, length(SEQN), 1440)
for (i in 1:length(SEQN)) {
  ActPerson[i, ] <- colMeans(holder[AllAct$SEQN == SEQN[i], ])
}
# We just traspose the resulting matrix to make it have the required shape
ActPerson <- t(ActPerson)
# We finally save it as a file
saveRDS(ActPerson, "ActPerson.rds")


# Weekend and week day case
# We separate each case
AllActEnd <- AllAct[AllAct$WEEKDAY == 1 | AllAct$WEEKDAY == 7, ]
ActEnd <- AllActEnd[, Actcols]
# We take the identifiers in the weekend case (should be less than the general case)
SEQNEnd <- unique(AllActEnd$SEQN)
PersonVarsfilEnd <- PersonVarsfil[PersonVarsfil$SEQN %in% SEQNEnd, ]

# Repeat with the weekday case
AllActnotEnd <- AllAct[AllAct$WEEKDAY != 1 & AllAct$WEEKDAY != 7, ]
ActnotEnd <- AllActnotEnd[, Actcols]

# We now perform the two means through two loops (both take a while)
ActPersonEnd <- matrix(0, length(SEQNEnd), 1440)
ActPersonnotEnd <- matrix(0, length(SEQN), 1440)

for (i in 1:length(SEQN)) {
  ActPersonnotEnd[i, ] <- colMeans(ActnotEnd[AllActnotEnd$SEQN == SEQN[i], ])
}

for (i in 1:length(SEQNEnd)) {
  ActPersonEnd[i, ] <- colMeans(ActEnd[AllActEnd$SEQN == SEQNEnd[i], ])
}

# We transpose both matrices
ActPersonEnd <- t(ActPersonEnd)
ActPersonnotEnd <- t(ActPersonnotEnd)

# We finally save the result
saveRDS(ActPersonEnd, "ActPersonEnd.rds")
saveRDS(ActPersonnotEnd, "ActPersonnotEnd.rds")
saveRDS(SEQNEnd, "SEQNEnd.rds")
