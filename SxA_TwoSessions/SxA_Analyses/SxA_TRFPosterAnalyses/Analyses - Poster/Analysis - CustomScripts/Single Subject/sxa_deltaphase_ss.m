function []=sxa_deltaphase_ss(subj, plots)
% ITPC Pre-Target
disp('Starting Delta Phase Analysis')
% Parameters
triggercodes={71;72;73}; % Warning Signals per condition
timerange=[-600 1600];
paddingLength=500; % in ms
% Freq range should be around 1/1.7= 0.59 till 1/0.7= 1.43
% Geometric mean of 0.59 and 1.43 (1/f logarithmic relation) = sqrt(1.43*0.59)=0.9185 
% calclating frquencies in octave: taking base 2 logarithm of the frequency in Hz ----> 2^(freq in octave) = freq in Hz
%ITPCdelta_wavFreqs= sqrt(1.43*0.59)*[2^-1  2^1];
ITPCdelta_wavFreqs=[1.25*2^-0.5 1.25*2^0.5];

% Load Data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG'
%subj=input("Subject Number? ");
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
ITPCdelta_elec = 1:71;

% band-pass filter 0.5-2 Hz butterworth, 24dB/octave
bpFilteredData = bandPassFilter(min(ITPCdelta_wavFreqs),max(ITPCdelta_wavFreqs),data(:,ITPCdelta_elec),srate);
% Hilbert Transform
hil = hilbert(bpFilteredData);

% Extract Phase
delta_phase = angle(hil); % phase

for c=1:size(triggercodes,1) % For conditions
% Segmentation
    [Delta_phaseSingleTrials, isNotArtifact, timeVecITPC]=segmentContEEGdata(triggercodes{c} , timerange,...
        delta_phase , triggers, artifacts, srate);

    % Check how many trials are artifact-free
    sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))

    % Remove trials with artifacts
    artrej=input("Remove artifact trials (1=yes)? ");
    if artrej==1
       Delta_phaseSingleTrials=Delta_phaseSingleTrials(:,:,isNotArtifact==1);
    end

    ITPCdelta_nTrialsCond(c)=size(Delta_phaseSingleTrials,3);
    ITPCdelta_res{c}=circ_r(Delta_phaseSingleTrials,[], [], 3);
    ITPCdelta_timevec{c}=timeVecITPC;
    clear timeVecITPC Delta_phaseSingleTrials isNotArtifact artrej
end
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
save(savefilename, 'ITPCdelta_res', 'ITPCdelta_timevec', 'ITPCdelta_wavFreqs', 'ITPCdelta_elec', 'ITPCdelta_nTrialsCond',"-append");

disp('Delta Phase Results Saved')


%% check ITPC time course and topography
if plots
    timeROI = [700 750];
    baseRange = [0 100];

    if occelonly
        dispElectrodes = [25:30 62:64]; % Occipital only
    else
        dispElectrodes = ITPCdelta_elec;
    end

    % mean time course across occipital electrodes
    figure;
    for c=1:size(triggercodes,1) % For conditions
        subplot(1,3,c)
        plot(ITPCdelta_timevec{c}, squeeze(mean(ITPCdelta_res{c}(:,dispElectrodes), 2)));
        ylim([0 0.5])
        line([0 0],ylim,'LineStyle', '-', 'Color', [0 0 0])
        line([700 700],ylim,'LineStyle', '- -', 'Color', [0 0 0])
        if c==1
            title('ITPC Rhythm')
        elseif c==2
            title('ITPC Interval')
        elseif c==3
            title('ITPC Irregular')
        end
    end
    
    figure;
    for c=1:size(triggercodes,1) % For conditions
        % topography (in relation to baseline)
        subplot(1,3,c)
        bp_phase= baselineCorrectSegmentedData(ITPCdelta_res{c}, ITPCdelta_timevec{c}, baseRange);
        topoplot(squeeze(mean(bp_phase(ITPCdelta_timevec{c}>=timeROI(1) & ITPCdelta_timevec{c}<=timeROI(2), :))), 'head71.locs');
        colorbar
        if c==1
            title('Rhythm')
        elseif c==2
            title('Interval')
        elseif c==3
            title('Irregular')
        end
    end
end
end