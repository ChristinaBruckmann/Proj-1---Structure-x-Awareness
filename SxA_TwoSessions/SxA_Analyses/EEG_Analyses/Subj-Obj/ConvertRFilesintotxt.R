# Convert R file into matlab-readable file

subj<-c(105:108, 111, 113, 114)

setwd("C:/Users/cbruckmann/Documents/PhD Projects/Proj1 - StructurexAwareness/SxA_TwoSessions/SxA_Data/EEG Results/SubjObj")

for (s in subj){
  loadfilename=sprintf("~/PhD Projects/Proj1 - StructurexAwareness/SxA_TwoSessions/SxA_Data/EEG Results/SubjObj/Subj%i_ObjectiveRhyhmAlpha.Rdata",s)
  load(loadfilename)
  
  savefilename=sprintf('Subj%i_ObjectiveRhyhmAlpha_contrast.txt',s)
  write.table(contrast_coeff_objrhy, file = savefilename, row.names=FALSE)
  
  savefilename=sprintf('Subj%i_ObjectiveRhyhmAlpha_power.txt',s)
  write.table(power_coeff_objrhy, file = savefilename, row.names=FALSE)

}

for (s in subj){
  loadfilename=sprintf("~/PhD Projects/Proj1 - StructurexAwareness/SxA_TwoSessions/SxA_Data/EEG Results/SubjObj/Subj%i_SubjectiveRhyhmAlpha.Rdata",s)
  load(loadfilename)
  
  savefilename=sprintf('Subj%i_SubjectiveRhyhmAlpha_contrast.txt',s)
  write.table(contrast_coeff_subjrhy, file = savefilename, row.names=FALSE)
  
  savefilename=sprintf('Subj%i_SubjectiveRhyhmAlpha_power.txt',s)
  write.table(power_coeff_subjrhy, file = savefilename, row.names=FALSE)
  
}

for (s in subj){
  loadfilename=sprintf("~/PhD Projects/Proj1 - StructurexAwareness/SxA_TwoSessions/SxA_Data/EEG Results/SubjObj/Subj%i_ObjectiveIntervalAlpha.Rdata",s)
  load(loadfilename)
  
  savefilename=sprintf('Subj%i_ObjectiveIntervalAlpha_contrast.txt',s)
  write.table(contrast_coeff_objint, file = savefilename, row.names=FALSE)
  
  savefilename=sprintf('Subj%i_ObjectiveIntervalAlpha_power.txt',s)
  write.table(power_coeff_objint, file = savefilename, row.names=FALSE)
  
}

for (s in subj){
  loadfilename=sprintf("~/PhD Projects/Proj1 - StructurexAwareness/SxA_TwoSessions/SxA_Data/EEG Results/SubjObj/Subj%i_SubjectiveIntervalAlpha.Rdata",s)
  load(loadfilename)
  
  savefilename=sprintf('Subj%i_SubjectiveIntervalAlpha_contrast.txt',s)
  write.table(contrast_coeff_subjint, file = savefilename, row.names=FALSE)
  
  savefilename=sprintf('Subj%i_SubjectiveIntervalAlpha_power.txt',s)
  write.table(power_coeff_subjint, file = savefilename, row.names=FALSE)
  
}
