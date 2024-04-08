## Simple Logistic Regression
# Integrates to the Matlab Code Comparing different fitting methods.

# Import Data from Matlab
setwd("C:/Program Files/R/R-4.1.3/library")
install.packages('R.matlab')
install.packages('faraway')
library(R.matlab)
library(faraway)

setwd("~/PhD Projects/Proj1 - StructurexAwareness/Other Files/Simulations")

alldata <- readMat("C:/Users/cbruckmann/Documents/PhD Projects/Proj1 - StructurexAwareness/Other Files/Simulations/SimDataforR.mat")

Rdata=alldata[["rdata"]]
rm(alldata)


# Convert behavioural data into data frames and 
Rdata=as.data.frame(Rdata)
colnames(Rdata) <- c('con','resp')

    
# Run Logistic Regression Models
    
R_model <- glm(Rdata$resp ~ Rdata$con, family="binomial")
R_res<- c(R_model[["coefficients"]][["(Intercept)"]],R_model[["coefficients"]][["Rdata$con"]])

summary(R_model)

write.table(R_res, file = 'R_output.txt', row.names=FALSE)

