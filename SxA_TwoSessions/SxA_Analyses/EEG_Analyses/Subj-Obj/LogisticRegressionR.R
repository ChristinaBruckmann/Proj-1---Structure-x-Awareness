# Logistic Regression Alpha-Contrast-Perception for one subject


# Import Data from Matlab
library(R.matlab)

alldata <- readMat('~/PhD Projects/Proj1 - StructurexAwareness/SxA_TwoSessions/Data Analysis/SxA_EEG_Analyses_Current/Results/SubjObj/DataforLogReg_Subj20.mat')

powerdata_rhy <- alldata[["data.power"]][[1]][[1]] # frequencies x timepoints x trials
powerdata_int <- alldata[["data.power"]][[2]][[1]]
behavdata_rhy <- alldata[["all.behav"]][[1]][[1]]
behavdata_int <- alldata[["all.behav"]][[2]][[1]]
rm(alldata)



# Convert behavioural data into data frames and 
behavdata_rhy=as.data.frame(behavdata_rhy)
colnames(behavdata_rhy) <- c('con', 'obj', 'subj')

behavdata_int=as.data.frame(behavdata_int)
colnames(behavdata_int) <- c('con', 'obj', 'subj')


# Initialize Output Matrices
contrast_coeff_objrhy <- matrix(, nrow = dim(powerdata_rhy)[1], ncol = dim(powerdata_rhy)[2])
power_coeff_objrhy <- matrix(, nrow = dim(powerdata_rhy)[1], ncol = dim(powerdata_rhy)[2])

contrast_coeff_objint <- matrix(, nrow = dim(powerdata_int)[1], ncol = dim(powerdata_int)[2])
power_coeff_objint <- matrix(, nrow = dim(powerdata_int)[1], ncol = dim(powerdata_int)[2])

contrast_coeff_subjrhy <- matrix(, nrow = dim(powerdata_rhy)[1], ncol = dim(powerdata_rhy)[2])
power_coeff_subjrhy <- matrix(, nrow = dim(powerdata_rhy)[1], ncol = dim(powerdata_rhy)[2])

contrast_coeff_subjint <- matrix(, nrow = dim(powerdata_int)[1], ncol = dim(powerdata_int)[2])
power_coeff_subjint <- matrix(, nrow = dim(powerdata_int)[1], ncol = dim(powerdata_int)[2])

# Iterate across frequencies
for (freq in (1:dim(powerdata_int)[1])) {  # currently takes interval trial n, change!
  
  # Iterate across time points
  for (tp in (1:dim(powerdata_int)[2])) {
    
    # Select correct power data
    currpow_rhy = powerdata_rhy[freq,tp,]
    currpow_int = powerdata_int[freq,tp,]
    
    # Run Logistic Regression Models
    
    objmodel_rhy <- glm(behavdata_rhy$obj ~ behavdata_rhy$con + currpow_rhy, family="binomial")
    contrast_coeff_objrhy[[freq,tp]] <- objmodel_rhy[["coefficients"]][["behavdata_rhy$con"]]
    power_coeff_objrhy[[freq,tp]] <- objmodel_rhy[["coefficients"]][["currpow_rhy"]]
    
    objmodel_int <- glm(behavdata_int$obj ~ behavdata_int$con + currpow_int, family="binomial")
    contrast_coeff_objint[[freq,tp]] <- objmodel_int[["coefficients"]][["behavdata_int$con"]]
    power_coeff_objint[[freq,tp]] <- objmodel_int[["coefficients"]][["currpow_int"]]
    
    subjmodel_rhy <- glm(behavdata_rhy$subj ~ behavdata_rhy$con + currpow_rhy, family="binomial")
    contrast_coeff_subjrhy[[freq,tp]] <- subjmodel_rhy[["coefficients"]][["behavdata_rhy$con"]]
    power_coeff_subjrhy[[freq,tp]] <- subjmodel_rhy[["coefficients"]][["currpow_rhy"]]
    
    subjmodel_int <- glm(behavdata_int$subj ~ behavdata_int$con + currpow_int, family="binomial")
    contrast_coeff_subjint[[freq,tp]] <- subjmodel_int[["coefficients"]][["behavdata_int$con"]]
    power_coeff_subjint[[freq,tp]] <- subjmodel_int[["coefficients"]][["currpow_int"]]
  }
}

# Plot
library(plot3D)
require(plot3D)
persp3D(z = contrast_coeff_objrhy, theta = 120)
persp3D(z = power_coeff_objrhy, theta = 120)

save(contrast_coeff_objrhy, power_coeff_objrhy, file = "ObjectiveRhyhmAlpha_Subj22.Rdata")
save(contrast_coeff_objint, power_coeff_objint, file = "ObjectiveIntervalRhyhmAlpha_Subj22.Rdata")
save(contrast_coeff_subjrhy, power_coeff_subjrhy, file = "SubjectiveRhyhmAlpha_Subj22.Rdata")
save(contrast_coeff_subjint, power_coeff_subjint, file = "SubjectiveIntervalAlpha_Subj22.Rdata")
