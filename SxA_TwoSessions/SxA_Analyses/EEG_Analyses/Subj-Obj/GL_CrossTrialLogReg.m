%% Cross_Trial_Logregression_GL
% (misspelled 'Rhythm' as 'rhyhm' when saving, keeping the change here)
clear
clc
subj=[101:103, 105:108, 111, 113, 114];

cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\SubjObj\Routput'
%% Load Data
% Data Structure (Subject, Obj(1)-Sub(2),Rhy(1)-Int(2), ContrastCoeff(1)-PowerCoeff(2))
for s=1:length(subj)

    % Objective Rhythm
    loadfilename1=sprintf('Subj%i_ObjectiveRhyhmAlpha_contrast.txt',subj(s));
    loadfilename2=sprintf('Subj%i_ObjectiveRhyhmAlpha_power.txt',subj(s));
    logregcoeff(s,1,1,1,:,:)=readmatrix(loadfilename1); % read and save contrast coefficients
    logregcoeff(s,1,1,2,:,:)=readmatrix(loadfilename2); % read and save power coefficients

    % Objective Interval
    loadfilename1=sprintf('Subj%i_ObjectiveIntervalAlpha_contrast.txt',subj(s));
    loadfilename2=sprintf('Subj%i_ObjectiveIntervalAlpha_power.txt',subj(s));
    logregcoeff(s,1,2,1,:,:)=readmatrix(loadfilename1); % read and save contrast coefficients
    logregcoeff(s,1,2,2,:,:)=readmatrix(loadfilename2); % read and save power coefficients

    % Subjective Rhythm
    loadfilename1=sprintf('Subj%i_SubjectiveRhyhmAlpha_contrast.txt',subj(s));
    loadfilename2=sprintf('Subj%i_SubjectiveRhyhmAlpha_power.txt',subj(s));
    logregcoeff(s,2,1,1,:,:)=readmatrix(loadfilename1); % read and save contrast coefficients
    logregcoeff(s,2,1,2,:,:)=readmatrix(loadfilename2); % read and save power coefficients

    % Subjective Interval
    loadfilename1=sprintf('Subj%i_SubjectiveIntervalAlpha_contrast.txt',subj(s));
    loadfilename2=sprintf('Subj%i_SubjectiveIntervalAlpha_power.txt',subj(s));
    logregcoeff(s,2,2,1,:,:)=readmatrix(loadfilename1); % read and save contrast coefficients
    logregcoeff(s,2,2,2,:,:)=readmatrix(loadfilename2); % read and save power coefficients
end


% Average across subjects
logregcoeff_mean=mean(logregcoeff,1);
size(logregcoeff_mean)

% Average across conditions (for now, as both are predictive)
logregcoeff_mean=mean(logregcoeff_mean,3);
size(logregcoeff_mean)

% Plot
for objsubj=[1:2]
    figure;
    subplot(1,2,1)
    heatmap(squeeze(logregcoeff_mean(1,objsubj,1,1,:,:))) % Plot Contrast Coefficients
    subplot(1,2,1)
    heatmap(squeeze(logregcoeff_mean(1,objsubj,1,2,:,:))) % Plot Power Coefficients
    title('Objective Power Coefficients - Predictive Conditions')
end
