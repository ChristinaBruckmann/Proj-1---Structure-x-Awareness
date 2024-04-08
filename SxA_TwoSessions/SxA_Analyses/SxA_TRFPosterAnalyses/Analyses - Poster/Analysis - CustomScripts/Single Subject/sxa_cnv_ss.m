function []=sxa_cnv_ss(subj,plots)
% ERP: CNV ANalysis
% clear
% clc

disp('Starting CNV Analysis')
% Load data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG'
%subj=input("Subject Number? ");
loadfilename=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj);
savefilename=sprintf('EEG_SxA_Subj%i_Results.mat',subj);
load(loadfilename)

% Parameters
triggercodes=[71 72 73]; % Warning Signals per condition
timerange=[-200 1500];

% Extract Info and Data from File
data=SDATA.data;
triggers=SDATA.events.triggerChannel;
CNV_channelnames=SDATA.info.channel_labels;
srate=SDATA.info.sampling_rate;
artifacts=SDATA.metadata.artifacts;
%timepoints=(0:length(data)-1)/srate; % time points in seconds

% Segment Data
[CNV_SingleTrials, isNotArtifact, CNV_timeVec]=segmentContEEGdata(triggercodes, timerange, data, triggers, artifacts, srate);

% Remove trials with artifacts
artrej=input("Remove artifact trials (1=yes)? ");
if artrej==1
    CNV_SingleTrials=CNV_SingleTrials(:,:,isNotArtifact==1);
end

% Select region of interest (some occipital electrodes)
CNV_electodes=["FC1 (11)","FCz (47)","FC2 (46)", "C1 (12)", "Cz (48)", "C2 (49)", "CP1 (19)", "CPz (32)", "CP2 (56)"];

roiidx=contains(CNV_channelnames,CNV_electodes,'IgnoreCase',true);

% Average across trials and across ROI electrodes
CNV_allroiidx=find(roiidx==1);

CNV_perElec=mean(CNV_SingleTrials(:,CNV_allroiidx,:),3);

CNV_Mean=mean(CNV_perElec,2);

cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
save(savefilename,'CNV_Mean', 'CNV_perElec','CNV_electodes','CNV_timeVec','CNV_allroiidx','CNV_channelnames')

disp('CNV results saved')
% Plot
if plots
    figure;
    for elecidx=1:length(CNV_electodes)
        subplot(2,round(length(CNV_electodes)/2),elecidx); plot(CNV_timeVec,CNV_perElec(:,elecidx))
        title(CNV_channelnames{CNV_allroiidx(elecidx),1})
        x1=xline(0,'r--','Warning Signal');
        x2=xline(800,'b--','Predicted Target');
        x1.FontSize = 8;
        x2.FontSize = 8;
    end

    figure; plot(CNV_timeVec,CNV_Mean)
    title(sprintf('Average Activity from ROI channels %s', strcat(CNV_channelnames{CNV_allroiidx,1})))
    xline(0,'r--','Warning Signal');
    xline(800,'b--','Predicted Target');
end
end