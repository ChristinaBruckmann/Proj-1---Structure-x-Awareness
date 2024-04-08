%% Correlation Alpha Baseline and Trough (SxA)
% Checks if the intensity of the suppression depends on the baseline
clc
clear

subj=[17:22 101:103 105:106 108];
catchonly=0;
alpharange=8:12;

% For each subject load and store baseline and through for each condition
for s=1:length(subj)
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data Analysis\SxA_EEG_Analyses_Current\Results\AlphaRes'
    if catchonly
        loadfilename=sprintf('EEG_SxA_Subj%i_AlphaResults_Catch.mat',subj(s));
        baserange=[-900 -800]; % baseline time points in relation to Target (Target at 0)
        troughrange=[0 800];
    else
        loadfilename=sprintf('EEG_SxA_Subj%i_AlphaResults.mat',subj(s));
        baserange=[-100 0]; % baseline time points in relation to WS (WS at time point 0)
        troughrange=[-800 0];
    end

    load(loadfilename,'alpha_Results'); % timepoints x freq x electrodes
    load(loadfilename,'alpha_timeVecTotal'); %timevector
    bigtimevec=alpha_timeVecTotal;

    for c=1:3 % For each condition
    
    % Select data
    timevec=bigtimevec{c};
    data=alpha_Results{c};
    data=squeeze(mean(data(:,alpharange,:),3)); % Select alpha range and average across electrodes
    data=squeeze(mean(data,2)); % Average across alpha range

    % Extract baseline for each condition
    baseidx=(timevec>baserange(1))&(baserange(2)>timevec);
    baseline(s,c)=mean(data(baseidx)); % Calculate average alpha activity across alpha range at baseline

    % Search for lowest trough
    troughidx=(timevec>troughrange(1))&(troughrange(2)>timevec);
    troughval(s,c)=min(data(troughidx)); % find lowest point in data
    troughdepth(s,c)=baseline(s,c)-troughval(s,c);
    end
end

%% Correlate Baseline and Trough

% General
corrcoef(troughdepth,baseline)
