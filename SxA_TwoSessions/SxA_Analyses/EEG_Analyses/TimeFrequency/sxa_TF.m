% TF Analysis
function []=sxa_TF(subj,baseline)
%% TF Analysis
% Input: Subject Number, generate plots (1/0)
disp('Starting Time Frequency Analysis')
singletrials=0; %(save single trials?)
artrej=1; % already remove artifact trials?

% Parameters (real analysis)
if baseline % baseline means, that the whole code will extract the TF at the baseline (beginning of trial), note that it doesnt mean that the data itself is already baseline corrected.
    % Pre-Target Baseline
    triggercodes={31;32;33}; % Start of trial
    timerange=[0 1000]; % The shortest ITI is 1 second (jittered in experiment between 1s and 1.6s)
else
    % Around WS
    triggercodes={71;72;73}; % Warning Signals per condition
    timerange=[-200 1500];
end

paddingLength=500; % in ms
TF_wavFreqs=2.^[0:1/6:5]; % Log Range
%alpharange=8:12;

% Load Data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Preprocessed'
loadfilename=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj);

if singletrials
    if baseline
        savefilename=sprintf('EEG_SxA_Subj%i_TF_SingleTrials_BL.mat',subj);
    else
        savefilename=sprintf('EEG_SxA_Subj%i_TF_SingleTrials.mat',subj);
    end
else
    if baseline
        savefilename=sprintf('EEG_SxA_Subj%i_TF_Results_BL.mat',subj);
    else
        savefilename=sprintf('EEG_SxA_Subj%i_TF_Results.mat',subj);
    end
end
load(loadfilename)

% Extract Info and Data from File
data=SDATA.data;
artifacts=SDATA.metadata.artifacts;
triggers=SDATA.events.triggerChannel;
srate=SDATA.info.sampling_rate;

% Segmentation
for c=1:size(triggercodes,1) % For conditions
    [segmentedData, isNotArtifact, timeVec]=segmentContEEGdata(triggercodes{c}, timerange+[-paddingLength paddingLength], data, triggers, artifacts, srate);

    % Check how many trials are artifact-free
    sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))

    % Remove trials with artifacts
    if artrej==1
        segmentedData=segmentedData(:,:,isNotArtifact==1); 
    end

    % TF Analysis
    tic
    if ~singletrials
        parfor el=1:width(data)
            [wvlt_amp, ~] = morletwave(TF_wavFreqs, 12, squeeze(segmentedData(:,el,:))', srate, 0, 'waitbar', 'on'); % time points x frequnencies x trials
            inducedMat=squeeze(mean(wvlt_amp,3)); % average over trials
            condResults(:,:,el)=inducedMat'; % reformat to: time points x frequencies x electrodes
        end
    else
        parfor el=1:width(data)
            [wvlt_amp, ~] = morletwave(TF_wavFreqs, 12, squeeze(segmentedData(:,el,:))', srate, 0, 'waitbar', 'on'); % time points x frequnencies x trials
            inducedMat=wvlt_amp; % do not average over trials
            condResults(:,:,:,el)=permute(inducedMat,[2,1,3]); % reformat to: time points x frequencies x trials x electrodes
        end
    end
    toc

    % Remove Padding
    if ~singletrials
        condResults=condResults(timeVec>=timerange(1) & timeVec<=timerange(2) ,:,:);
    else
        condResults=condResults(timeVec>=timerange(1) & timeVec<=timerange(2) ,:,:,:);
    end

    timeVec=timeVec(timeVec>=timerange(1) & timeVec<=timerange(2));
    TF_Results{c}=condResults; % Time Points, Frequencies, Trials, Electrodes
    TF_timeVecTotal{c}=timeVec;
    TF_NotArtifact{c}=isNotArtifact;
    clear condResults segmentedData timeVec subj
end
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\TimeFreq'
if ~singletrials
    save(savefilename, 'TF_Results', 'TF_timeVecTotal', 'TF_wavFreqs');
else
    TF_Results_Trial=TF_Results;
    TF_trial_timeVec=TF_timeVecTotal;
    save(savefilename, 'TF_Results_Trial','TF_trial_timeVec','TF_NotArtifact','-v7.3');
end
disp('TF Results Saved')
end