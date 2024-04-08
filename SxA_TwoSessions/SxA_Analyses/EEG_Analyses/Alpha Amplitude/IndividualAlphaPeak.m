%% Extract every persons alpha peak frequency (Welch across whole data), take a range around that (+-2) and run alpha analysis with that.
% Run this after TF analysis
clear
clc

% Individual Alpha Peak Frequency vs. General Average Across 8-12Hz
wavFreqs=1:40; % Range taken from Elmira, Assaf uses 1:30 in plos bio
alpharange=6:16; % Range to check for peak
subj=[14 15 17:22];

% Welch parameters
window_seconds=1; % window size in seconds
overlap_propotion=0.5; % propotion of window size that overlaps 

remart=input('Remove Artifacts? Yes-1: ');

%figure; 
for s=1:length(subj)
% Load Data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\EEG Data'
% subj=input("Subject Number? ");
occelonly=1;
loadfilename=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj(s));
%savefilename=sprintf('EEG_SxA_Subj%i_Results.mat',subj(s));
load(loadfilename)

% Extract Info and Data from File
data=SDATA.data;
artifacts=SDATA.metadata.artifacts;
triggers=SDATA.events.triggerChannel;
channelnames=SDATA.info.channel_labels;
srate=SDATA.info.sampling_rate;

if remart
data=data(~artifacts,:);
end

if occelonly
alpha_tfElectrodes=[25:30 62:64]; % Occipital only
else
alpha_tfElectrodes=1:64; % All electrodes
end

% Get Spectral information
data=data(:,alpha_tfElectrodes,:);
win=window_seconds*srate;
overlap=overlap_propotion*win;
[frequencyspec,Hz]=pwelch(data,win, overlap,[],srate);
plot(Hz,frequencyspec)

% Select Peak
averagespec=mean(frequencyspec,2);
alphaidx=sum(Hz==alpharange,2);
[maxvalue,~]=max(averagespec(logical(alphaidx)));
peakfreq=Hz(find(averagespec==maxvalue)); % only works if the value doesn't show up twice, find a better way to do this

alpha_peakrange=[peakfreq-2:peakfreq+2];
peakidx=sum(Hz==alpha_peakrange,2);

% subplot(2,4,s); plot(Hz,frequencyspec)
% gcf;xlim([0 30])

% % Load TF results
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\SxA_EEG_Analyses_Current\Results'
occelonly=1;
loadfilename=sprintf('EEG_SxA_Subj%i_Results_New.mat',subj(s));
savefilename=sprintf('EEG_SxA_Subj%i_Results_New.mat',subj(s));
load(loadfilename)

% Extract Power
for c=1:3 % for each condition
amp_data(c,:,:,:)=TF_Results{c}; % data(conditionxtimepointsxfrequenciesxelectrode)
end
timeVec=TF_timeVecTotal{1};

if occelonly
    alphaelec=[25:30 62:64]; % Occipital only
else
    alphaelec=1:64; % all electrodes
end

for c=1:3 % separately for each condition
alpha_amp_peak_ss(c,:)=mean(squeeze(mean(amp_data(c,:,alpha_peakrange,alphaelec),3)),2); % individual peak frequency, average across frequencies and electrodes
alpha_amp_gen_ss(c,:)=mean(squeeze(mean(amp_data(c,:,8:12,alphaelec),3)),2); % same alpha for everyone

alpha_amp_peak_gl(s,c,:)=mean(squeeze(mean(amp_data(c,:,alpha_peakrange,alphaelec),3)),2); % individual peak frequency, average across frequencies and electrodes
alpha_amp_gen_gl(s,c,:)=mean(squeeze(mean(amp_data(c,:,8:12,alphaelec),3)),2); % same alpha for everyone
end

% Save
save(savefilename, 'alpha_amp_peak_ss', 'alpha_amp_gen_ss', 'alpha_peakrange', 'frequencyspec','Hz',"-append");
end

% Plot Individual Alpha Amps
for s=1:8
    figure;
    for c=1:3
        subplot(1,2,1)
        plot(timeVec,squeeze(alpha_amp_peak_gl(s,c,:)))
        hold on
        subplot(1,2,2)
        plot(timeVec,squeeze(alpha_amp_gen_gl(s,c,:)))
        hold on
    end
end

% Plot Average
basec=1;
baselinetp=[0 100];

figure; 
for c=1:3
    if basec
        datatoplot1=baselineCorrectSegmentedData(squeeze(mean(alpha_amp_peak_gl(:,c,:),1)), timeVec, baselinetp);
        datatoplot2=baselineCorrectSegmentedData(squeeze(mean(alpha_amp_gen_gl(:,c,:),1)), timeVec, baselinetp);
    else
        datatoplot1=squeeze(mean(alpha_amp_peak_gl(:,c,:),1));
        datatoplot2=squeeze(mean(alpha_amp_gen_gl(:,c,:),1));
    end

    subplot(1,2,1)
    plot(timeVec,datatoplot1) % average across subjects
    title('Individual Alpha Peak Frequency')
    hold on

    subplot(1,2,2)
    plot(timeVec,datatoplot2)
    hold on
    title('General Alpha Peak Frequency')
end

% Get SD
for c=1:3
    standdev_peak(c)=mean(std(squeeze(alpha_amp_peak_gl(:,c,:))));
    standdev_gen(c)=mean(std(squeeze(alpha_amp_gen_gl(:,c,:))));
end

subplot(1,2,1);bar(standdev_peak)
subplot(1,2,2);bar(standdev_gen)

% Plot with SD
figure; 
for c=1:3
    subplot(1,2,1)
    varplot(squeeze(alpha_amp_peak_gl(:,c,:))') % average across subjects
    title('Individual Alpha Peak Frequency')
    hold on

    subplot(1,2,2)
    varplot(squeeze(alpha_amp_gen_gl(:,c,:))') 
    hold on
    title('General Alpha Peak Frequency')
end
