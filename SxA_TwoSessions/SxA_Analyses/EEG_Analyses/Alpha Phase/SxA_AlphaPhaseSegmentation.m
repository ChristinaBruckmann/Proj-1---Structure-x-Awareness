%% Alpha Segmentation into Trials
function SxA_AlphaPhaseSegmentation(subj,cluster,catchonly)
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
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\AlphaPhaseRes'
    loadfilename1=sprintf('EEG_SxA_Subj%i_AlphaPhaseSingleTrials.mat',subj(s));
    loadfilename2=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj(s));
    %savefilename=sprintf('EEG_SxA_Subj%i_Results_NewAlpha.mat',subj(s));
    if cluster==1
        load(loadfilename1,'alpha_phase_whole_occ')
        alpha_data=alpha_phase_whole_occ;
    else
        load(loadfilename1,'alpha_phase_whole_cen')
        alpha_data=alpha_phase_whole_cen;
    end
   
    %cd 'Z:\el-Christina\SxA\SxA_Data\EEG Preprocessed'
    cd 'D:\'
    load(loadfilename2,'SDATA')

    artifacts=SDATA.metadata.artifacts;
    triggers=SDATA.events.triggerChannel;
    channelnames=SDATA.info.channel_labels;
    srate=SDATA.info.sampling_rate;

    for c=1:size(triggercodes,1) % For conditions

        % Segmentation
        [Alpha_phaseSingleTrials, isNotArtifact, timeVecITPC]=segmentContEEGdata(triggercodes{c} , timerange,...
            alpha_data , triggers, artifacts, srate);

        % Get index for irregular trials without target in time window if not catch trials
        if ~catchonly
            if c==3
                cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\Behavioural Preprocessed'
                load(sprintf('SxA_ResultsSubject%i_Session2.mat',subj(s)),'alldataclean')
                idx_notar=ismember(alldataclean{(alldataclean{:,'Condition'}==c),'Irregular Target Time'},irrtartimes);
            else
                idx_notar=ones(size(Alpha_phaseSingleTrials,3));
            end
        end

        % Check how many trials are artifact-free
        sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))

        % Remove trials with artifacts
        artrej=1;
        if artrej==1 && c==3 && ~catchonly
            Alpha_SingleTrials{c}=Alpha_phaseSingleTrials(:,:,logical(isNotArtifact & idx_notar));
        else
            Alpha_SingleTrials{c}=Alpha_phaseSingleTrials(:,:,logical(isNotArtifact));
        end

        ITPCalpha_nTrials(c)=size(Alpha_SingleTrials{c},3);
        ITPCalpha_timevec{c}=timeVecITPC;

        % Jackknifing to get single trial ITPC
        for tr=1:ITPCalpha_nTrials(c)
            tridx=ones(1,ITPCalpha_nTrials(c));
            tridx(tr)=0;
            currtrials=Alpha_SingleTrials{c}(:,elec,logical(tridx)); % leave out one trial (time points x electrodes x trials)
            ITPCalpha_jack(s,c,tr,:,:)=circ_r(currtrials,[], [], 3); % (subject, condition, trial, ITPC, electrode)
            % Plot (each trial for each condition for each subject)
            subplot(1,3,c); plot(ITPCalpha_timevec{c},squeeze(mean(ITPCalpha_jack(s,c,tr,:,:),5)))
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
        xline(-900)
        xline(0)
        else %aligned to WS
        xline(900)
        xline(0)
        end
    end

    % Save
     cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\AlphaPhaseRes'
   % cd 'Z:\el-Christina\SxA\SxA_Results\AlphaPhaseRes'
    if catchonly
        if cluster==1
            savefilename=sprintf('Alpha_OccipitalCatch_Subj%i',subj(s));
            Alpha_SingleTrials_occ_catch=Alpha_SingleTrials;
            ITPCalpha_nTrials_occ_catch=ITPCalpha_nTrials;
            ITPCalpha_timevec_occ_catch=ITPCalpha_timevec;
            ITPCalpha_jack_occ_catch=ITPCalpha_jack;
            save(savefilename, "Alpha_SingleTrials_occ_catch","ITPCalpha_nTrials_occ_catch","ITPCalpha_timevec_occ_catch","ITPCalpha_jack_occ_catch",'-v7.3')
        else
            savefilename=sprintf('Alpha_CentralCatch_Subj%i',subj(s));
            Alpha_SingleTrials_cen_catch=Alpha_SingleTrials;
            ITPCalpha_nTrials_cen_catch=ITPCalpha_nTrials;
            ITPCalpha_timevec_cen_catch=ITPCalpha_timevec;
            ITPCalpha_jack_cen_catch=ITPCalpha_jack;
            save(savefilename, "Alpha_SingleTrials_cen_catch","ITPCalpha_nTrials_cen_catch","ITPCalpha_timevec_cen_catch","ITPCalpha_jack_cen_catch",'-v7.3')
        end
    else
        if cluster==1
            savefilename=sprintf('Alpha_OccipitalAll_Subj%i',subj(s));
            Alpha_SingleTrials_occ=Alpha_SingleTrials;
            ITPCalpha_nTrials_occ=ITPCalpha_nTrials;
            ITPCalpha_timevec_occ=ITPCalpha_timevec;
            ITPCalpha_jack_occ=ITPCalpha_jack;
            save(savefilename, "Alpha_SingleTrials_occ","ITPCalpha_nTrials_occ","ITPCalpha_timevec_occ","ITPCalpha_jack_occ",'-v7.3')
        else
            savefilename=sprintf('Alpha_CentralAll_Subj%i',subj(s));
            Alpha_SingleTrials_cen=Alpha_SingleTrials;
            ITPCalpha_nTrials_cen=ITPCalpha_nTrials;
            ITPCalpha_timevec_cen=ITPCalpha_timevec;
            ITPCalpha_jack_cen=ITPCalpha_jack;
            save(savefilename, "Alpha_SingleTrials_cen","ITPCalpha_nTrials_cen","ITPCalpha_timevec_cen","ITPCalpha_jack_cen",'-v7.3')
        end
    end
end
end