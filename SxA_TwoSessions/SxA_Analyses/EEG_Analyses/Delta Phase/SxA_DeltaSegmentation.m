%% Delta Segmentation into Trials
function SxA_DeltaSegmentation(subj,cluster,catchonly)
% cluster=2; %(1 occipital, 2 central)
% catchonly=0; %(1-yes)

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
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\DeltaRes'
    loadfilename1=sprintf('EEG_SxA_Subj%i_DeltaPhaseSingleTrials.mat',subj(s));
    loadfilename2=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj(s));
    %savefilename=sprintf('EEG_SxA_Subj%i_Results_NewDelta.mat',subj(s));
    if cluster==1
        load(loadfilename1,'delta_phase_whole_occ')
        delta_data=delta_phase_whole_occ;
    else
        load(loadfilename1,'delta_phase_whole_cen')
        delta_data=delta_phase_whole_cen;
    end
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Preprocessed'
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
                cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\Behavioural Preprocessed'
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
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\DeltaRes'
    if catchonly
        if cluster==1
            savefilename=sprintf('OccipitalCatch_Subj%i',subj(s));
            Delta_SingleTrials_occ_catch=Delta_SingleTrials;
            ITPCdelta_nTrials_occ_catch=ITPCdelta_nTrials;
            ITPCdelta_timevec_occ_catch=ITPCdelta_timevec;
            ITPCdelta_jack_occ_catch=ITPCdelta_jack;
            save(savefilename, "Delta_SingleTrials_occ_catch","ITPCdelta_nTrials_occ_catch","ITPCdelta_timevec_occ_catch","ITPCdelta_jack_occ_catch",'-v7.3')
        else
            savefilename=sprintf('CentralCatch_Subj%i',subj(s));
            Delta_SingleTrials_cen_catch=Delta_SingleTrials;
            ITPCdelta_nTrials_cen_catch=ITPCdelta_nTrials;
            ITPCdelta_timevec_cen_catch=ITPCdelta_timevec;
            ITPCdelta_jack_cen_catch=ITPCdelta_jack;
            save(savefilename, "Delta_SingleTrials_cen_catch","ITPCdelta_nTrials_cen_catch","ITPCdelta_timevec_cen_catch","ITPCdelta_jack_cen_catch",'-v7.3')
        end
    else
        if cluster==1
            savefilename=sprintf('OccipitalAll_Subj%i',subj(s));
            Delta_SingleTrials_occ=Delta_SingleTrials;
            ITPCdelta_nTrials_occ=ITPCdelta_nTrials;
            ITPCdelta_timevec_occ=ITPCdelta_timevec;
            ITPCdelta_jack_occ=ITPCdelta_jack;
            save(savefilename, "Delta_SingleTrials_occ","ITPCdelta_nTrials_occ","ITPCdelta_timevec_occ","ITPCdelta_jack_occ",'-v7.3')
        else
            savefilename=sprintf('CentralAll_Subj%i',subj(s));
            Delta_SingleTrials_cen=Delta_SingleTrials;
            ITPCdelta_nTrials_cen=ITPCdelta_nTrials;
            ITPCdelta_timevec_cen=ITPCdelta_timevec;
            ITPCdelta_jack_cen=ITPCdelta_jack;
            save(savefilename, "Delta_SingleTrials_cen","ITPCdelta_nTrials_cen","ITPCdelta_timevec_cen","ITPCdelta_jack_cen",'-v7.3')
        end
    end
end
end