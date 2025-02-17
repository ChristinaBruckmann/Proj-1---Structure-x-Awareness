%% SxA Multivariate Log Reg Code
clear
clc

%% Parameters
freqs=8:12;
elecs=[25:30 62:64]; % occipital
subj=[101:103 105:108 110 112:114 117:119 121 122 124 126 127 129 130]; 
downsample=1; % downsample time to 30ms windows for faster processing?
irrtartimes=[3 4 5]; % when irregular targets can appear (1 and 2 before 900ms, 3 is 900ms, 4 and 5 after 900ms)

%% Run for each subject
for s=1:length(subj)
    %% Load Data
    % EEG
    % cd 'Y:\el-Christina\SxA\SxA_Results\New TF Results'
    cd '/home/cbruckmann/LogRegStuff/New TF Results'
    loadfilename=sprintf('EEG_SxA_Subj%i_TF_SingleTrials.mat',subj(s));
    savefilename=sprintf('EEG_SxA_LogRes_GL_subj%i.mat',subj(s));
    load(loadfilename); % three conditions: time points x frequencies x trials x electrodes

    % Behaviour
    %cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\Behavioural Preprocessed'
    cd '/home/cbruckmann/LogRegStuff/LogRegBehaviour'
    loadfilename=sprintf('SxA_ResultsSubject%i_Session2.mat',subj(s));
    load(loadfilename,'alldataclean'); % three conditions: time points x frequencies x trials x electrodes

    %% Run for each condition
    for c=1:3
        EEG_data=TF_Results_Trial{1, c}; % time points x frequencies x trials x electrodes
        timeVec=TF_trial_timeVec{1, c}; % choose time vector
        notart=logical(TF_NotArtifact{1, c}); % choose artifact vector
        %% Select Electrodes and average
        EEG_data=squeeze(mean(EEG_data(:,:,:,elecs),4)); % output: tp x freq x trials
        %% Select Frequencies
        EEG_data=EEG_data(:,freqs,:); % output: tp x freqs x trials
        %% Downsample
        if downsample
            % Define window size in milliseconds
            window_size_ms = 30;

            % Calculate the number of points per window
            fs = length(timeVec) / (timeVec(end) - timeVec(1)) * 1000; % Sampling frequency in Hz (points per ms)
            points_per_window = round(window_size_ms * fs / 1000);

            % Reshape the EEG_data to group time points into 30ms windows
            nWindows = floor(size(EEG_data, 1) / points_per_window); % Number of full windows
            reshaped_data = reshape(EEG_data(1:nWindows*points_per_window, :), points_per_window, nWindows, size(EEG_data, 2),size(EEG_data, 3));

            % Average across the first dimension (time points in each window)
            EEG_data = squeeze(mean(reshaped_data, 1));

            % Downsampled time vector:
            timeVec = timeVec(1:points_per_window:end);
            timeVec = timeVec(1:nWindows);  % Truncate to the number of windows
        end
        %% Format and Clean (Behavioural) Data
        % Select
        behav_data_obj=alldataclean(alldataclean{:,"Condition"}==c,{'Correct/Incorrect'}); % objective performance for current condition trials
        behav_data_subj=alldataclean(alldataclean{:,"Condition"}==c,{'Binary Visibility'}); % subjective visibility for current condition trials
        contrast_level=alldataclean(alldataclean{:,"Condition"}==c,{'Contrast Level'}); % get contrast level for each trial

        % Remove artifacts from behaviour (already removed from EEG data)
        behav_data_obj=table2array(behav_data_obj(notart,:));
        behav_data_subj=table2array(behav_data_subj(notart,:));
        contrast_level=table2array(contrast_level(notart,:));
        EEG_data=EEG_data(:,:,notart,:); % Dimensions: tp x freq x trials x elec

         % Remove trials with unwanted irregular target times
        clean_behavdata=alldataclean{(alldataclean{:,'Condition'}==c)&notart,:}; % first remove artifacts from behavioural data to match indeces of trials for EEG
       
        if c==3
            idx_notar=ismember(clean_behavdata{(clean_behavdata{:,'Condition'}==c),'Irregular Target Time'},irrtartimes);
        else
            idx_notar=ones(size(behav_data_obj,1));
        end

        % Remove Catch Trials
        idx=~isnan(behav_data_obj); % Find catch trials and mark as 0
        behav_data_obj=behav_data_obj(idx);
        behav_data_subj=behav_data_subj(idx);
        contrast_level=contrast_level(idx);
        EEG_data=EEG_data(:,:,idx);
        %% Run Log Reg for each  frequency and time point
        for f=1:size(EEG_data,2) % freqs
            for tp=1:size(EEG_data,1) %time points
                curr_power=double(squeeze(EEG_data(tp,f,:))); % Select power for each trial at current time point and frequency

                % Run for objective and save
                [bestFitParams,~, ~, ~]=fitPsychometric_binaryData([contrast_level,curr_power], behav_data_obj, [0.4, 0.6, 0.9, 1]);
                res_obj(s, c, f, tp , :) = bestFitParams;

                % Run for subjective and save
                [bestFitParams,~, ~,~]=fitPsychometric_binaryData([contrast_level,curr_power], behav_data_subj, [0, 0.3, 0.9, 1]);
                res_subj(s, c, f, tp , :)= bestFitParams; % subj x cond x freq x time point x params
            end
        end
    end
    %% Save
    %cd 'Y:\el-Christina\SxA\SxA_Results\LogRegResults'
    fprintf("Completed Subject %i/%i.\n",s,length(subj)) 
    save(savefilename,"res_subj","res_obj","timeVec","freqs","elecs")
end
save('GL_Res_LogReg',"res_subj","res_obj","timeVec","freqs","elecs")
disp("LogReg Results Completed")