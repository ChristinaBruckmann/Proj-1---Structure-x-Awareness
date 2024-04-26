# Logistic Regression Alpha-Contrast-Perception


# Import Data from Matlab

install.packages("R.matlab")
library("R.matlab")

#subj<-c(106, 107, 108, 111, 113, 114)
subj<-c(111)

for (s in subj) {
  readname=sprintf('~/PhD Projects/Proj1 - StructurexAwareness/SxA_TwoSessions/SxA_Data/EEG Results/SubjObj/RInput_LogReg_S%i.mat',s)
  
  alldata <- readMat(readname)
  
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
  # library(plot3D)
  # require(plot3D)
  # persp3D(z = contrast_coeff_objrhy, theta = 120)
  # persp3D(z = power_coeff_objrhy, theta = 120)
  
  filename1=sprintf("Subj%i_ObjectiveRhyhmAlpha.Rdata",s)
  filename2=sprintf("Subj%i_ObjectiveIntervalAlpha.Rdata",s)
  filename3=sprintf("Subj%i_SubjectiveRhyhmAlpha.Rdata",s)
  filename4=sprintf("Subj%i_SubjectiveIntervalAlpha.Rdata",s)
  # save(contrast_coeff_objrhy, power_coeff_objrhy, file = filename1)
  # save(contrast_coeff_objint, power_coeff_objint, file = filename2)
  # save(contrast_coeff_subjrhy, power_coeff_subjrhy, file = filename3)
  # save(contrast_coeff_subjint, power_coeff_subjint, file = filename4)
  
  #Save Objective Rhythm
  savefilename=sprintf('Subj%i_ObjectiveRhyhmAlpha_contrast.txt',s)
  write.table(contrast_coeff_objrhy, file = savefilename, row.names=FALSE)
  
  savefilename=sprintf('Subj%i_ObjectiveRhyhmAlpha_power.txt',s)
  write.table(power_coeff_objrhy, file = savefilename, row.names=FALSE)
  
  #Save Objective Interval
  
  savefilename=sprintf('Subj%i_ObjectiveIntervalAlpha_contrast.txt',s)
  write.table(contrast_coeff_objint, file = savefilename, row.names=FALSE)
  
  savefilename=sprintf('Subj%i_ObjectiveIntervalAlpha_power.txt',s)
  write.table(power_coeff_objint, file = savefilename, row.names=FALSE)
  
  # Save Subjective Rhythm
  
  savefilename=sprintf('Subj%i_SubjectiveRhyhmAlpha_contrast.txt',s)
  write.table(contrast_coeff_subjrhy, file = savefilename, row.names=FALSE)
  
  savefilename=sprintf('Subj%i_SubjectiveRhyhmAlpha_power.txt',s)
  write.table(power_coeff_subjrhy, file = savefilename, row.names=FALSE)
  
  # Save Subjective Interval
  
  savefilename=sprintf('Subj%i_SubjectiveIntervalAlpha_contrast.txt',s)
  write.table(contrast_coeff_subjint, file = savefilename, row.names=FALSE)
  
  savefilename=sprintf('Subj%i_SubjectiveIntervalAlpha_power.txt',s)
  write.table(power_coeff_subjint, file = savefilename, row.names=FALSE)
}
