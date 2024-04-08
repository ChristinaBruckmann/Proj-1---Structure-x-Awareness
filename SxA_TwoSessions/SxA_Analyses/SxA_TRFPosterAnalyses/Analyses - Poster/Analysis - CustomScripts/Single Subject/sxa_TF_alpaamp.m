function []=sxa_TF_alpaamp(subj,plots)
%% TF Analysis - Alpha Amplitude - Compare Conditions
% Still missing statistical analysis
% Input: Subject Number, generate plots (1/0)

% clc
% clear 
disp('Starting Alpha Amplitude Analysis')
% Parameters
triggercodes={71;72;73}; % Warning Signals per condition
%triggercodes=[71 72 73]; % All Warning Signals (rhythm, interval, irregular)
timerange=[-200 1500];
paddingLength=500; % in ms
alpha_wavFreqs=1:40; % Range taken from Elmira, Assaf uses 1:30 in plos bio
alpharange=8:12;

% Load Data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG'
% subj=input("Subject Number? ");
occelonly=input("Display occipital electrodes only (1-yes, 0-all)? ");
loadfilename=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj);
savefilename=sprintf('EEG_SxA_Subj%i_Results.mat',subj);
load(loadfilename)

% Extract Info and Data from File
data=SDATA.data;
artifacts=SDATA.metadata.artifacts;
triggers=SDATA.events.triggerChannel;
channelnames=SDATA.info.channel_labels;
srate=SDATA.info.sampling_rate;
if occelonly
alpha_tfElectrodes=[25:30 62:64]; % Occipital only
else
alpha_tfElectrodes=1:64; % All electrodes
end

% Segmentation
for c=1:size(triggercodes,1) % For conditions
    [segmentedData, isNotArtifact, timeVec]=segmentContEEGdata(triggercodes{c}, timerange+[-paddingLength paddingLength], data, triggers, artifacts, srate);

    % Check how many trials are artifact-free
    sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))

    % Remove trials with artifacts
    artrej=input("Remove artifact trials (1=yes)? ");
    if artrej==1
        segmentedData=segmentedData(:,:,isNotArtifact==1); 
    end

    % TF Analysis
    tic
    parfor el=1:length(alpha_tfElectrodes)
        [wvlt_amp, ~] = morletwave(alpha_wavFreqs, 12, squeeze(segmentedData(:,alpha_tfElectrodes(el),:))', srate, 0, 'waitbar', 'on'); % time points x frequnencies x trials
        inducedMat=squeeze(mean(wvlt_amp,3)); % average over trials
        condResults(:,:,el)=inducedMat'; % time points x frequencies x electrodes
    end
    toc

    % Remove Padding
    condResults=condResults(timeVec>=timerange(1) & timeVec<=timerange(2) ,:,:);
    timeVec=timeVec(timeVec>=timerange(1) & timeVec<=timerange(2));
    alpha_Results{c}=condResults; % Time Points, Frequencies, Electrodes
    alpha_timeVecTotal{c}=timeVec;

    clear condResults segmentedData timeVec
end
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
save(savefilename, 'alpha_Results', 'alpha_timeVecTotal', 'alpha_wavFreqs', 'alpha_tfElectrodes',"-append");
disp('Alpha Amplitude Results Saved')
%% Plot TF
if plots
% Baseline Correct?
%basec=input("Baseline Correction (1-yes)? ");
basec=1;
baselinetp=[0 100]; %baseline time points in relation to WS (WS at time point 0)

for c=1:size(triggercodes,1) % For conditions

    if basec
        %     [~,startbaseidx]=(min(abs(timeVec-baselinetp(1))));
        %     [~,endbaseidx]=(min(abs(timeVec-baselinetp(2))));
        %     baseline=mean(tfResults(startbaseidx:endbaseidx,:,:), 1);
        %     datatoplot=tfResults-baseline;
        datatoplot=baselineCorrectSegmentedData(alpha_Results{c}, alpha_timeVecTotal{c}, baselinetp);
    else
        datatoplot=alpha_Results{c};
    end

    %minmaxscale= [min(mean(datatoplot,3),[],"all") max(mean(datatoplot,3),[],"all")];
    minmaxscale=[-3 3];

    figure;
    imagesc(alpha_timeVecTotal{c}, alpha_wavFreqs, squeeze(mean(datatoplot,3))',minmaxscale); axis xy; hold all
    xline(0,'r--','Warning Signal');
    xline(800,'b--','Predicted Target');
    colorbar();
    xlabel('Time (ms)', 'FontWeight','bold', 'FontSize', 10);
    ylabel('Frequency (Hz)', 'FontWeight','bold', 'FontSize', 10);
    if c==1
        title('Rhythm - Pre-Target Alpha Amplitude');
    elseif c==2
        title('Interval - Pre-Target Alpha Amplitude');
    elseif c==3
        title('Irregular - Pre-Target Alpha Amplitude');
    end
end

figure;
for c=1:size(triggercodes,1) % For conditions
    if basec
        %     [~,startbaseidx]=(min(abs(timeVec-baselinetp(1))));
        %     [~,endbaseidx]=(min(abs(timeVec-baselinetp(2))));
        %     baseline=mean(tfResults(startbaseidx:endbaseidx,:,:), 1);
        %     datatoplot=tfResults-baseline;
        datatoplot=baselineCorrectSegmentedData(alpha_Results{c}, alpha_timeVecTotal{c}, baselinetp);
    else
        datatoplot=alpha_Results{c};
    end
    % Plot Alpha Only (8-12 hz)
    plot(alpha_timeVecTotal{c},mean(mean(datatoplot(:,alpharange,:),3),2), 'LineWidth', 2)
    % alphaamp=tfResults(:,alpharange,:);
    % meanelec=mean(alphaamp,3);
    % varplot(mean(alphaamp,3))
    legend('Rhythm','Interval', 'Irregular')
    title("Alpha Power between WS and Target")
    hold on;
end
xline(0,'r--','Warning Signal');
xline(800,'b--','Predicted Target');
end
end