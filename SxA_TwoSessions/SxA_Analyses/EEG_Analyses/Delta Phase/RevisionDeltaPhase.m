%% Re-Analysis of Delta Phase from scratch to make sure there are no bugs in the first script. Only occipital for now.
% Addition: Individual Trial data gets saved to that they can all be
% plotted and compared.
clear 
clc
%subj=[17:22];
subj=[113];
%% Extract ITPC from whole data (with different reference depending on cluster)

for s=1:length(subj)
disp('Starting Delta Phase Analysis')

% Load Data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data\EEG Preprocessed'

loadfilename=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj(s));
savefilename=sprintf('EEG_SxA_Subj%i_Results_NewDelta.mat',subj(s));
load(loadfilename)

% Parameters
ITPCdelta_wavFreqs=[1.25*2^-0.75 1.25*2^0.75];
srate=SDATA.info.sampling_rate;
ITPCdelta_elec= 1:64;

% Re-Reference for different clusters
if subj(s)==101 % has a noisy nose, always refernce to mastoids
    data_occ=referenceContEEGdata(SDATA.data,[69 71]);
else
    data_occ=referenceContEEGdata(SDATA.data,71); %nose
end

if subj==19 % has a noisy mastoid, always reference to nose
    data_cen=data_occ;
else
    data_cen=referenceContEEGdata(SDATA.data,[69 71]); % mastoids
end

% band-pass filter 0.5-2 Hz butterworth, 24dB/octave
bpFilteredData_occ = bandPassFilter(min(ITPCdelta_wavFreqs),max(ITPCdelta_wavFreqs),data_occ(:,ITPCdelta_elec),srate);
bpFilteredData_cen = bandPassFilter(min(ITPCdelta_wavFreqs),max(ITPCdelta_wavFreqs),data_cen(:,ITPCdelta_elec),srate);

% Hilbert Transform
hil_occ = hilbert(bpFilteredData_occ);
hil_cen = hilbert(bpFilteredData_cen);


% Extract Phase
delta_ITPC_whole_occ = angle(hil_occ); % phase
delta_ITPC_whole_cen = angle(hil_cen); % phase

cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data Analysis\SxA_EEG_Analyses_Current\Results\DeltaRes'
save(savefilename, 'delta_ITPC_whole_occ','delta_ITPC_whole_cen','-v7.3')
end

%% Segmentation into Trials
clear
clc
subj=[107 108 113];
cluster=2; %(1 occipital, 2 central)
catchonly=0; %(1-yes)

if catchonly
    triggercodes={1;2;3}; % Catch Target Signals per condition (cut around target) 
    timerange=[-1500 400];
else
    irrtartimes=[3 4 5]; % only choose trials where the target does not appear before 800ms
    triggercodes={71;72;73}; % From Warning Signal
    timerange=[-700 1200];
end

if cluster==1
    elec=[25:30 62:64]; % occipital electrodes
else
    elec=[11:12 46:49]; % central electrodes
end

for s=1:length(subj)
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data Analysis\SxA_EEG_Analyses_Current\Results\DeltaRes'
    loadfilename1=sprintf('EEG_SxA_Subj%i_Results_NewDelta.mat',subj(s));
    loadfilename2=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj(s));
    %savefilename=sprintf('EEG_SxA_Subj%i_Results_NewDelta.mat',subj(s));
    if cluster==1
        load(loadfilename1,'delta_ITPC_whole_occ')
        delta_data=delta_ITPC_whole_occ;
    else
        load(loadfilename1,'delta_ITPC_whole_cen')
        delta_data=delta_ITPC_whole_cen;
    end
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data\EEG Preprocessed'
    load(loadfilename2,'SDATA')

    artifacts=SDATA.metadata.artifacts;
    triggers=SDATA.events.triggerChannel;
    channelnames=SDATA.info.channel_labels;
    srate=SDATA.info.sampling_rate;

    for c=1:size(triggercodes,1) % For conditions

        % Segmentation
        [Delta_phaseSingleTrials, isNotArtifact, timeVecITPC]=segmentContEEGdata(triggercodes{c} , timerange,...
            delta_data , triggers, artifacts, srate);

        % Get index for irregular trials without target in time window if not catch trials
        if ~catchonly
            if c==3
                cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data\Behavioural Preprocessed'
                load(sprintf('SxA_ResultsSubject%i_Session2.mat',subj(s)),'alldataclean')
                idx_notar=ismember(alldataclean{(alldataclean{:,'Condition'}==c),'Irregular Target Time'},irrtartimes);
            else
                idx_notar=ones(size(Delta_phaseSingleTrials,3));
            end
        end

        % Check how many trials are artifact-free
        sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))

        % Remove trials with artifacts
        artrej=1;
        if artrej==1 && c==3 && ~catchonly
            Delta_SingleTrials{c}=Delta_phaseSingleTrials(:,:,logical(isNotArtifact & idx_notar));
        else
            Delta_SingleTrials{c}=Delta_phaseSingleTrials(:,:,logical(isNotArtifact));
        end

        ITPCdelta_nTrials(c)=size(Delta_SingleTrials{c},3);
        ITPCdelta_timevec{c}=timeVecITPC;

        % Jackknifing to get single trial ITPC
        for tr=1:ITPCdelta_nTrials(c)
            tridx=ones(1,ITPCdelta_nTrials(c));
            tridx(tr)=0;
            currtrials=Delta_SingleTrials{c}(:,elec,logical(tridx)); % leave out one trial (time points x electrodes x trials)
            ITPCdelta_jack(s,c,tr,:,:)=circ_r(currtrials,[], [], 3); % (subject, condition, trial, ITPC, electrode)
            % Plot (each trial for each condition for each subject)
            subplot(1,3,c); plot(ITPCdelta_timevec{c},squeeze(mean(ITPCdelta_jack(s,c,tr,:,:),5)))
            hold on
        end

        if c==1
            title(sprintf('ITPC Rhythm', subj(s)))
        elseif c==2
            title(sprintf('ITPC Interval', subj(s)))
        else
            title(sprintf('ITPC Irregular', subj(s)))
        end

        ylim([0 0.4]);

        if catchonly %aligned to target
        xline(-800)
        xline(0)
        else %aligned to WS
        xline(800)
        xline(0)
        end
    end

    % Save
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data Analysis\SxA_EEG_Analyses_Current\Results\DeltaRes'
    if catchonly
        if cluster==1
            savefilename=sprintf('OccipitalCatch_Subj%i',subj(s));
            Delta_SingleTrials_occ_catch=Delta_SingleTrials;
            ITPCdelta_nTrials_occ_catch=ITPCdelta_nTrials;
            ITPCdelta_timevec_occ_catch=ITPCdelta_timevec;
            ITPCdelta_jack_occ_catch=ITPCdelta_jack;
            save(savefilename, "Delta_SingleTrials_occ_catch","ITPCdelta_nTrials_occ_catch","ITPCdelta_timevec_occ_catch","ITPCdelta_jack_occ_catch")
        else
            savefilename=sprintf('CentralCatch_Subj%i',subj(s));
            Delta_SingleTrials_cen_catch=Delta_SingleTrials;
            ITPCdelta_nTrials_cen_catch=ITPCdelta_nTrials;
            ITPCdelta_timevec_cen_catch=ITPCdelta_timevec;
            ITPCdelta_jack_cen_catch=ITPCdelta_jack;
            save(savefilename, "Delta_SingleTrials_cen_catch","ITPCdelta_nTrials_cen_catch","ITPCdelta_timevec_cen_catch","ITPCdelta_jack_cen_catch")
        end
    else
        if cluster==1
            savefilename=sprintf('OccipitalAll_Subj%i',subj(s));
            Delta_SingleTrials_occ=Delta_SingleTrials;
            ITPCdelta_nTrials_occ=ITPCdelta_nTrials;
            ITPCdelta_timevec_occ=ITPCdelta_timevec;
            ITPCdelta_jack_occ=ITPCdelta_jack;
            save(savefilename, "Delta_SingleTrials_occ","ITPCdelta_nTrials_occ","ITPCdelta_timevec_occ","ITPCdelta_jack_occ")
        else
            savefilename=sprintf('CentralAll_Subj%i',subj(s));
            Delta_SingleTrials_cen=Delta_SingleTrials;
            ITPCdelta_nTrials_cen=ITPCdelta_nTrials;
            ITPCdelta_timevec_cen=ITPCdelta_timevec;
            ITPCdelta_jack_cen=ITPCdelta_jack;
            save(savefilename, "Delta_SingleTrials_cen","ITPCdelta_nTrials_cen","ITPCdelta_timevec_cen","ITPCdelta_jack_cen")
        end
    end
end

%% Delta Average
clear
clc
%subj=[17:22 101:103 105:106];
subj=[101:103 105:108 113];
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data Analysis\SxA_EEG_Analyses_Current\Results\DeltaRes'
cluster=1; %(1 occipital, 2 central)
catchonly=1; %(1-yes)

if cluster==1
    elec=[25:30 62:64]; % occipital electrodes
else
    elec=[11:12 46:49]; % central electrodes
end

% Load
for s=1:length(subj)
    if cluster==1 && ~catchonly
        loadfilename=sprintf('OccipitalAll_Subj%i',subj(s));
        load(loadfilename,"Delta_SingleTrials_occ","ITPCdelta_timevec_occ","ITPCdelta_nTrials_occ") %(time points x electrodes x trials)
        for c=1:3
        Delta_SingleTrials(s,c,:,:)=circ_r(Delta_SingleTrials_occ{c},[], [], 3);
        Delta_TimeVec(s,c,:)=ITPCdelta_timevec_occ{c};
        Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_occ(c);
        end
    elseif cluster==1 && catchonly
        loadfilename=sprintf('OccipitalCatch_Subj%i',subj(s));
        load(loadfilename,"Delta_SingleTrials_occ_catch","ITPCdelta_timevec_occ_catch","ITPCdelta_nTrials_occ_catch") %(time points x electrodes x trials)
        for c=1:3
        Delta_SingleTrials(s,c,:,:)=circ_r(Delta_SingleTrials_occ_catch{c},[], [], 3);
        Delta_TimeVec(s,c,:)=ITPCdelta_timevec_occ_catch{c};
        Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_occ_catch(c);
        end
    elseif cluster==2 && ~catchonly
        loadfilename=sprintf('CentralAll_Subj%i',subj(s));
        load(loadfilename,"Delta_SingleTrials_cen","ITPCdelta_timevec_cen","ITPCdelta_nTrials_cen") %(time points x electrodes x trials)
        for c=1:3
        Delta_SingleTrials(s,c,:,:)=circ_r(Delta_SingleTrials_cen{c},[], [], 3);
        Delta_TimeVec(s,c,:)=ITPCdelta_timevec_cen{c};
        Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_cen(c);
        end
    elseif cluster==2 && catchonly
        loadfilename=sprintf('CentralCatch_Subj%i',subj(s));
        load(loadfilename,"Delta_SingleTrials_cen_catch","ITPCdelta_timevec_cen_catch","ITPCdelta_nTrials_cen_catch") %(time points x electrodes x trials)
        for c=1:3
        Delta_SingleTrials(s,c,:,:)=circ_r(Delta_SingleTrials_cen_catch{c},[], [], 3);
        Delta_TimeVec(s,c,:)=ITPCdelta_timevec_cen_catch{c};
        Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_cen_catch(c);
        end
    end
end

%Average
timeVec=squeeze(Delta_TimeVec(1,1,:)); % All TV are identical
nTrials=sum(Delta_Ntrials,1); % nTrials for chance ITPC calculation
for c=1:3
    SubjMean=squeeze(mean(Delta_SingleTrials(:,c,:,elec),1));
    Delta(c,:)=squeeze(mean(SubjMean,2)); % mean across subjects and electrodes
end

% Plot
figure;
for c=1:3
    plot(timeVec,Delta(c,:),"LineWidth",2)
    ylim([0 0.6])
    hold on
    %     if c==1
    %         title('Rhythm')
    %     elseif c==2
    %         title('Interval')
    %     else
    %         title('Irregular')
    %     end
end
if catchonly
    xline(-800)
    title('Catch Trials')
else
    xline(800)
    title('All Trials')
end
xline(0)

% Save
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data Analysis\SxA_EEG_Analyses_Current\Results\DeltaRes'
if cluster==1 && catchonly
save("DeltaGroupLevel_Occ_Catch",'Delta','subj','timeVec','nTrials')
elseif cluster==1 && ~catchonly
save("DeltaGroupLevel_Occ_All",'Delta','subj','timeVec','nTrials')
elseif cluster==2 && catchonly
save("DeltaGroupLevel_Cen_Catch",'Delta','subj','timeVec','nTrials')
elseif cluster==2 && ~catchonly
save("DeltaGroupLevel_Cen_All",'Delta','subj','timeVec','nTrials')
end
