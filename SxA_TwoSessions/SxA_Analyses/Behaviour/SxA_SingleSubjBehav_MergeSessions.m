%% Merge Session 1 and 2
% This code takes the raw behavioural data for a single subject from two sessions of the SxA
% experiment and merges the two files into 1 big one.
%% Load and clean up each session
clear
clc

cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data\Raw'
% Session 1
subj_n=input('Subject Number? ');
filename1=sprintf('SxA1.1_s%i-session1.mat',subj_n);
filename2=sprintf('SxA1.1_s%i-session2.mat',subj_n);
savename=sprintf('SxA1.1_s%i-merged.mat',subj_n);

subresults1=load(filename1,'subresults');
subresults2=load(filename2,'subresults');
subresults.trialmatrix1=[subresults1.subresults.trialmatrix1;subresults2.subresults.trialmatrix1];
subresults.trialmatrix2=[subresults1.subresults.trialmatrix2;subresults2.subresults.trialmatrix2];
subresults.trialmatrix3=[subresults1.subresults.trialmatrix3;subresults2.subresults.trialmatrix3];
subresults.conditions=[subresults1.subresults.conditions;subresults2.subresults.conditions];
subresults.data=[subresults1.subresults.data;subresults2.subresults.data];

cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data\Raw\Behaviour Raw'

save(savename,'subresults')

clear
clc