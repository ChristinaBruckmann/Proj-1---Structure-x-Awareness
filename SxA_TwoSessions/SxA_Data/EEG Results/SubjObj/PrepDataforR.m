%%  Prepare Data for R 
% To run logistic regression in R
clear
clc

% Parameters
subj=[22];
freqrange=[1:40];
electrodes=[25:30 62:64]; % occipital
timep=[250 1000]; %(data aligned to warning signal, target at 800)
cond=[1 2]; % which conditions (rhythm, interval, irregular)
%% Load Data

% Load TF Data
for s=subj
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data Analysis\SxA_EEG_Analyses_Current\Results'
    loadfilename=sprintf('EEG_SxA_Subj%i_Results_New.mat',s);
    load(loadfilename,'TF_Results_Trial','TF_trial_timeVec','TF_NotArtifact')
    savefilename=sprintf('EEG_SxA_Subj%i_Results_SubjObj.mat',s);
    TF_res=TF_Results_Trial; % Size: frequencies timepoints trials electrodes
    TF_time=TF_trial_timeVec;
    TF_ArtVectors=TF_NotArtifact;

    % Convert to power and average across electrodes
    for c=cond % for each condition (only predictable for now)
        % Select Data
        data=TF_res{1, c};
        timevec=TF_time{1, c};
        idx=timevec>=timep(1)&timevec<=timep(2);
        data=mean(data(idx,freqrange,:,electrodes),4); % average across electrodes (time, freq, trials)
        % Convert to Power
        data_power{c}=squeeze(abs(data).^2);
    end

    % Load Behavioural Data
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data Analysis\Behaviour'
    loadfilename=sprintf('SxA_ResultsSubject%i_Session2',s);
    load(loadfilename,'alldataclean')
    behaviour=alldataclean;
    clear alldataclean

    % For each subject, create a table with contrast level, subj, obj for each
    % trial (artifacts removed)
    for c=cond
        logidx=logical(behaviour{:,'Condition'}==c);
        obj_resp=behaviour{logidx,'Correct/Incorrect'};
        subj_resp=behaviour{logidx,'Binary Visibility'};
        contrast=behaviour{logidx,'Contrast Level'};
        all_behav{c}=[contrast(logical(TF_NotArtifact{c})) obj_resp(logical(TF_NotArtifact{c}))  subj_resp(logical(TF_NotArtifact{c}))]; % save trial list with artifacts removed

        % Remove Catch Trials
        catchidx=all_behav{c}(:,1)~=0; % select non-catch trials only
        all_behav{c}=all_behav{c}(catchidx,:); % remove catch trials
        data_power{c}=data_power{c}(:,:,catchidx);

        % For irreular, remove trials with target before 800
        if c==3
            targettime=behaviour{logidx,'Irregular Target Time'};
            targettime=targettime(logical(TF_NotArtifact{c})); % remove artifacts
            targettime=targettime(catchidx,:); % remove catch trials
            targetidx=ismember(targettime,[3,4,5]); % find trials with target at 800ms or later
            all_behav{c}=all_behav{c}(targetidx,:);
            data_power{c}=data_power{c}(:,:,targetidx); % remove early target trials from data
        end
    end
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data Analysis\SxA_EEG_Analyses_Current\Results\SubjObj'
    save(sprintf('DataforLogReg_Subj%i',s),"data_power",'all_behav')
end

% After Running code in R, come back here.
Routput= readmatrix('R_output.txt');


