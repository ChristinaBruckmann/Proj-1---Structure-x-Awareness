%% Single Subject Pipeline SxA
% Every new subject for which we have behavioural data from both session, as well as preprocessed EEG-Data from the second one
% Should be run through this pipeline to create the result files necessary for all GL analyses
clear
clc
subj=input("Subject Number: ");

%% Behavioural Analysis
% Merge
SxA_SingleSubjBehav_MergeSessions(subj) % Merge behavioural data from session 1 and 2

% Analyse
disp('Starting Behavioural Analyses.')
SxA_SingleSubjBehav_Analysis(subj,1,0) % Pre-Processs and Analyse session 1
SxA_SingleSubjBehav_Analysis(subj,2,0) % Session 2
SxA_SingleSubjBehav_Analysis(subj,3,1) % Both sessions together
disp('Behavioural Analyses Done.')
%% EEG Analyses

% Time-Frequency
disp('Starting TF Analysis.')
sxa_TF(subj,0) % Time-Frequency Analysis (All Freqs, Saving Single Trials for later Behav-EEG LogReg)
sxa_TF(subj,1) % Baseline (pre-trial)
disp('TF Analysis Done.')

% Alpha Amplitude
disp('Starting Alpha Analysis.')
sxa_TF_alpaamp(subj,1,0) % Alpha Amplitude All Trials
sxa_TF_alpaamp(subj,1,1) % Alpha Amplitude Catch Trials Only
disp('Alpha Analysis Done.')

% Delta ITPC
SxA_DeltaPhaseExtraction(subj)
disp('Starting Delta Sementation.')
SxA_DeltaSegmentation(subj,1,0) % Occipital Cluster, All Trials
SxA_DeltaSegmentation(subj,2,0) % Central Cluster, All Trials
SxA_DeltaSegmentation(subj,1,1) % Occipital Cluster, Catch Only
SxA_DeltaSegmentation(subj,2,1) % Central Cluster, Catch Only

%% Behaviour-EEG Analyses

disp('All analyses done and saved.')