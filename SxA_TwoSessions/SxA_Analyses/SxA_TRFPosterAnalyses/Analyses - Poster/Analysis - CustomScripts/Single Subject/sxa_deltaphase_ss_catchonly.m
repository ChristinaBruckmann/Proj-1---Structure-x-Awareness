% Catch trials only (triggers 1,2,3 for rhythm, interval and irregular)
function []=sxa_deltaphase_ss_catchonly(subj, plots)
% ITPC Pre-Target
disp('Starting Delta Phase Analysis')
% Parameters
triggercodes={1;2;3}; % Warning Signals per condition
timerange=[-1300 500];
ITPCdelta_wavFreqs_catch=[1.25*2^-0.75 1.25*2^0.75];

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
ITPCdelta_elec_catch = 1:71;

% band-pass filter 0.5-2 Hz butterworth, 24dB/octave
bpFilteredData = bandPassFilter(min(ITPCdelta_wavFreqs_catch),max(ITPCdelta_wavFreqs_catch),data(:,ITPCdelta_elec_catch),srate);
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

    ITPCdelta_nTrialsCond_catch(c)=size(Delta_phaseSingleTrials,3);
    ITPCdelta_res_catch{c}=circ_r(Delta_phaseSingleTrials,[], [], 3);
    ITPCdelta_timevec_catch{c}=timeVecITPC;
    clear timeVecITPC Delta_phaseSingleTrials isNotArtifact artrej
end
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
save(savefilename, 'ITPCdelta_res_catch', 'ITPCdelta_timevec_catch', 'ITPCdelta_wavFreqs_catch', 'ITPCdelta_elec_catch', 'ITPCdelta_nTrialsCond_catch',"-append");

disp('Delta Phase Results Saved')


%% check ITPC time course and topography
if plots
    timeROI = [-100 -50];
    baseRange = [-800 -700];

    if occelonly
        dispElectrodes = [25:30 62:64]; % Occipital only
    else
        dispElectrodes = ITPCdelta_elec_catch;
    end

    % mean time course across occipital electrodes
    figure;
    for c=1:size(triggercodes,1) % For conditions
        subplot(1,3,c)
        plot(ITPCdelta_timevec_catch{c}, squeeze(mean(ITPCdelta_res_catch{c}(:,dispElectrodes), 2)));
        ylim([0 0.5])
        line([-800 -800],ylim,'LineStyle', '-', 'Color', [0 0 0])
        line([0 0],ylim,'LineStyle', '- -', 'Color', [0 0 0])
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
        bp_phase= baselineCorrectSegmentedData(ITPCdelta_res_catch{c}, ITPCdelta_timevec_catch{c}, baseRange);
        topoplot(squeeze(mean(bp_phase(ITPCdelta_timevec_catch{c}>=timeROI(1) & ITPCdelta_timevec_catch{c}<=timeROI(2), :))), 'head71.locs');
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