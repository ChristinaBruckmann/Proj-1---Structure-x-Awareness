% Merge 107
clear
clc

% Load 
load('SxA1.1_s107a-session1.mat')
subresults1=subresults;
load('SxA1.1_s107b-session1.mat')
subresults2=subresults;
clear subresults

% Create Merged Version (First block from 107a, the rest from 107b)
subresults.trialmatrix3=[subresults1.trialmatrix3(1:50,:);subresults2.trialmatrix3(51:end,:)];
subresults.data=[subresults1.data;subresults2.data];

subresults.irrcuetimes=[subresults1.irrcuetimes;subresults2.irrcuetimes];

% Add the rest
subresults.trialmatrix1=subresults2.trialmatrix1;
subresults.trialmatrix2=subresults2.trialmatrix2;
subresults.thresholdresults=subresults1.thresholdresults;
subresults.firstreversal=subresults1.firstreversal;
subresults.calculatedthreshold=subresults1.calculatedthreshold;
subresults.conditions=subresults1.conditions;
subresults.visStimifo=subresults1.visStimifo;
subresults.mixinfo=subresults1.mixinfo;
subresults.timeinfo=subresults1.timeinfo;


% Save
save('SxA1.1_s107-session1.mat',"subresults")
