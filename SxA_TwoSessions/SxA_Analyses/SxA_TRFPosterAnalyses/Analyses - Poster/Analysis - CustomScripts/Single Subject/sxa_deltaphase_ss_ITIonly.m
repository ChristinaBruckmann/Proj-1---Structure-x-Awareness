% Inter trial interval
function []=sxa_deltaphase_ss_ITIonly(subj, plots)
% ITPC Pre-Target
disp('Starting Delta Phase Analysis')
% Parameters
triggercodes={31;32;33}; % Warning Signals per condition
timerange=[-900 1400]; % ITI is one second, but just to be sure there is no smearing I take a bit less. from mask onset to S1 is 1 second
ITPCdelta_wavFreqs_ITI=[1.25*2^-0.5 1.25*2^0.5];

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
ITPCdelta_elec_ITI = 1:71;

% band-pass filter 0.5-2 Hz butterworth, 24dB/octave
bpFilteredData = bandPassFilter(min(ITPCdelta_wavFreqs_ITI),max(ITPCdelta_wavFreqs_ITI),data(:,ITPCdelta_elec_ITI),srate);
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

    ITPCdelta_nTrialsCond_ITI(c)=size(Delta_phaseSingleTrials,3);
    ITPCdelta_res_ITI{c}=circ_r(Delta_phaseSingleTrials,[], [], 3);
    ITPCdelta_timevec_ITI{c}=timeVecITPC;
    clear timeVecITPC Delta_phaseSingleTrials isNotArtifact artrej
end
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
save(savefilename, 'ITPCdelta_res_ITI', 'ITPCdelta_timevec_ITI', 'ITPCdelta_wavFreqs_ITI', 'ITPCdelta_elec_ITI', 'ITPCdelta_nTrialsCond_ITI',"-append");

disp('Delta Phase Results Saved')


%% check ITPC time course and topography
if plots
    if occelonly
        dispElectrodes = [25:30 62:64]; % Occipital only
    else
        dispElectrodes = ITPCdelta_elec_ITI;
    end

    % mean time course across occipital electrodes
    figure;
    for c=1:size(triggercodes,1) % For conditions
        subplot(1,3,c)
        plot(ITPCdelta_timevec_ITI{c}, squeeze(mean(ITPCdelta_res_ITI{c}(:,dispElectrodes), 2)));
        ylim([0 0.5])
        if c==1
            title('ITPC Rhythm')
        elseif c==2
            title('ITPC Interval')
        elseif c==3
            title('ITPC Irregular')
        end
    end
end
end