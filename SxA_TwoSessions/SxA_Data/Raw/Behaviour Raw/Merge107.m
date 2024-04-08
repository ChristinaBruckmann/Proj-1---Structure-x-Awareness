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

% Add the rest
subresults.trialmatrix1=subresults2.trialmatrix1;
subresults.trialmatrix2=subresults2.trialmatrix2;

% Save
save('SxA1.1_s107-session1.mat',"subresults")
