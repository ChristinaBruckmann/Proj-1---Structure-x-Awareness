%% CNV Analysis - Each Condition Separately
% Adding of option to only chose visible trials to also see resulution is
% planned

clear
clc
subj=[15 17:22];

% Single Subject
for s=1:length(subj)
    disp('Starting CNV Analysis')
    % Load data
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\EEG Data'
    loadfilename=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj(s));
    savefilename=sprintf('EEG_SxA_Subj%i_Results_New.mat',subj(s));
    load(loadfilename)

    % Parameters
    triggercodes={71;72;73}; % Warning Signals per condition
    timerange=[-200 1500];

    % Extract Info and Data from File
    data=SDATA.data;
    triggers=SDATA.events.triggerChannel;
    CNV_channelnames=SDATA.info.channel_labels;
    srate=SDATA.info.sampling_rate;
    artifacts=SDATA.metadata.artifacts;
    %timepoints=(0:length(data)-1)/srate; % time points in seconds
    figure;

    for c=1:length(triggercodes)
        % Segment Data
        [CNV_SingleTrials, isNotArtifact, CNV_timeVec]=segmentContEEGdata(triggercodes{c}, timerange, data, triggers, artifacts, srate);

        % Select Specific Trials Only
        onlyvis=1;
        vislevels=[6:10]; % Which trials are classified as visible?
        if onlyvis
            load(sprintf('SxA_ResultsSubject%i_Session2.mat',subj(s)),'alldataclean')
            idx_vis=ismember(alldataclean{(alldataclean{:,'Condition'}==c),'Contrast Level'},vislevels);
        else
            idx_vis=ones(size(CNV_SingleTrials,3));
        end

        % Remove trials with artifacts
        %artrej=input("Remove artifact trials (1=yes)? ");
        artrej=1;
        if artrej==1
            CNV_SingleTrials=CNV_SingleTrials(:,:,(isNotArtifact & idx_vis));
        end

        % Select region of interest (some occipital electrodes)
        CNV_electodes=["FC1 (11)","FCz (47)","FC2 (46)", "C1 (12)", "Cz (48)", "C2 (49)", "CP1 (19)", "CPz (32)", "CP2 (56)"];

        roiidx=contains(CNV_channelnames,CNV_electodes,'IgnoreCase',true);

        % Average across trials and across ROI electrodes
        CNV_allroiidx=find(roiidx==1);

        CNV_perElec=mean(CNV_SingleTrials(:,CNV_allroiidx,:),3);

        CNV_Mean(c,:)=mean(CNV_perElec,2);

        cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\SxA_EEG_Analyses_Current\Results'
        save(savefilename,'CNV_Mean','CNV_electodes','CNV_timeVec',"-append")

        % disp('CNV results saved')
        % Plot

        %         figure;
        %         for elecidx=1:length(CNV_electodes)
        %             subplot(2,round(length(CNV_electodes)/2),elecidx); plot(CNV_timeVec,CNV_perElec(:,elecidx))
        %             title(CNV_channelnames{CNV_allroiidx(elecidx),1})
        %             x1=xline(0,'r--','Warning Signal');
        %             x2=xline(800,'b--','Predicted Target');
        %             x1.FontSize = 8;
        %             x2.FontSize = 8;
        %         end

        subplot(1,3,c); plot(CNV_timeVec,CNV_Mean(c,:))
        title(sprintf('Average Activity from ROI channels %s', strcat(CNV_channelnames{CNV_allroiidx,1})))
        xline(0,'r--','Warning Signal');
        xline(800,'b--','Predicted Target');

    end
end

%% Group Level 

% Load Data
for s=1:length(subj)
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\SxA_EEG_Analyses_Current\Results'
    loadfilename=sprintf('EEG_SxA_Subj%i_Results_New.mat',subj(s));
    load(loadfilename,'CNV_Mean')
    for c=1:3 % conditions
    CNV_GL(s,c,:)=CNV_Mean(c,:);
    end
end

% Average
figure;
for c=1:3 % conditions
    CNV_GL_Mean(c,:)=mean(squeeze(CNV_GL(:,c,:)),1);

    % Plot
    subplot(1,3,c); plot(CNV_timeVec,CNV_GL_Mean(c,:))
    title(sprintf('Average Activity from ROI channels.'))
    xline(0,'r--','Warning Signal');
    xline(800,'b--','Predicted Target');
end

% %% Jackknifing
% for s=1:length(subj)
%     figure;
%     currentsubjs=1:length(subj);
%     currentsubjs(s)=[]; % leave out one subject
%     for c=1:3 % conditions
%         CNV_GL_Mean_JK(s,c,:)=mean(squeeze(CNV_GL(currentsubjs,c,:)),1); % calculate mean across other subjects
% 
%         % Plot
%         subplot(1,3,c); plot(CNV_timeVec,squeeze(CNV_GL_Mean_JK(s,c,:)))
%         title(sprintf('CNVs without Subject %i.',subj(s)))
%         xline(0,'r--','Warning Signal');
%         xline(800,'b--','Predicted Target');
%     end
% end