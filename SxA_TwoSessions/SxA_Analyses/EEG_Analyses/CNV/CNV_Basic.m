%% CNV Analysis - Each Condition Separately
% Adding of option to only chose visible trials to also see resulution is
% planned

clear
clc
subj=[101:103 105:108 110:114 116:119 121:122 124 126 127 129 131 132];
%subj=[124 126 127 129 131 132];
lpf_cutoff=20;
irrtartimes=[3]; % when irregular targets can appear (1 and 2 before 800ms, 3 is 800ms, 4 and 5 after 800ms)

onlyvis=0;
vislevels=[6:10]; % Which trials are classified as visible?

tarfree_before800=1; % only choose irregular trials which have a target at or after 800
tar_at800=1; % only choose irregular trials which have a target at 800
%% Single Subject
for s=1:length(subj)
    disp('Starting CNV Analysis')
    % Load data
    cd 'Y:\el-Christina\SxA\SxA_Data\EEG Preprocessed'
    loadfilename=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj(s));
    savefilename=sprintf('EEG_SxA_Subj%i_CNV.mat',subj(s));
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

        % Select Only Visible Trials
        if onlyvis
            load(sprintf('SxA_ResultsSubject%i_Session2.mat',subj(s)),'alldataclean')
            idx_vis=ismember(alldataclean{(alldataclean{:,'Condition'}==c),'Contrast Level'},vislevels);
        else
            idx_vis=ones(size(CNV_SingleTrials,3),1);
        end

        % Select only irregular trials with targets at moments specified above (tendentially only at or after 800 to reduce smearing)
        if c==3
            cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\Behavioural Preprocessed'
            load(sprintf('SxA_ResultsSubject%i_Session2.mat',subj(s)),'alldataclean','subresults')
            idx_notar=ismember(alldataclean{(alldataclean{:,'Condition'}==c),'Irregular Target Time'},irrtartimes);
        else
            idx_notar=ones(size(CNV_SingleTrials,3),1);
        end

        % Remove trials with artifacts
        %artrej=input("Remove artifact trials (1=yes)? ");
        artrej=1;
        if artrej==1
            CNV_SingleTrials=CNV_SingleTrials(:,:,(isNotArtifact & idx_vis& idx_notar));
        end

        % Select region of interest (some occipital electrodes)
        CNV_electodes=["FC1 (11)","FCz (47)","FC2 (46)", "C1 (12)", "Cz (48)", "C2 (49)", "CP1 (19)", "CPz (32)", "CP2 (56)"];

        roiidx=contains(CNV_channelnames,CNV_electodes,'IgnoreCase',true);

        % Average across trials and across ROI electrodes
        CNV_allroiidx=find(roiidx==1);

        CNV_perElec=mean(CNV_SingleTrials(:,CNV_allroiidx,:),3);

        CNV_Mean(c,:)=mean(CNV_perElec,2);
        CNV_Mean(c,:)=LPF( CNV_Mean(c,:), 1024, lpf_cutoff); % low pass filter the mean
        cd 'Y:\el-Christina\SxA\SxA_Results\CNV Results'
        save(savefilename,'CNV_Mean','CNV_electodes','CNV_timeVec')

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

bc=1; % baseline correct
bc_win=[-200 0]; % 0 is warning signal

% Load Data
for s=1:length(subj)
    cd 'Y:\el-Christina\SxA\SxA_Results\CNV Results'
    loadfilename=sprintf('EEG_SxA_Subj%i_CNV.mat',subj(s));
    load(loadfilename,'CNV_Mean','CNV_timeVec')
    CNV_GL(s,:,:)=CNV_Mean;
    timeVec_GL(s,:)=CNV_timeVec;
end

% Average 
timeVec=timeVec_GL(1,:);
CNV_GL_Mean=squeeze(mean(CNV_GL,1));

% Baseline correct
if bc
    bl_activity=CNV_GL_Mean(:,timeVec>=bc_win(1) &timeVec<=bc_win(2));
    bl=mean(bl_activity,2);
    CNV_GL_Mean=CNV_GL_Mean-bl;
end

% Plot

figure;
for c=1:3 % conditions
    plot(CNV_timeVec,CNV_GL_Mean(c,:),"LineWidth",3)
    hold on
end
    title(sprintf('Average Activity from ROI channels.'))
    xline(0,'r--','Warning Signal');
    xline(900,'b--','Predicted Target');
    legend(["Rhythm","Interval","Irregular"])

figure;
tiledlayout(1,3);
for c=1:3 % conditions
    %CNV_GL_Mean(c,:)=mean(squeeze(CNV_GL(:,c,:)),1);
    %CNV_GL_Mean(c,:)=LPF( CNV_GL_Mean(c,:), 1024, lpf_cutoff);
    % Plot
    nexttile
    %subplot(1,3,c); 
    plot(CNV_timeVec,CNV_GL_Mean(c,:))
    title(sprintf('Average Activity from ROI channels.'))
    xline(0,'r--','Warning Signal');
    xline(900,'b--','Predicted Target');
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