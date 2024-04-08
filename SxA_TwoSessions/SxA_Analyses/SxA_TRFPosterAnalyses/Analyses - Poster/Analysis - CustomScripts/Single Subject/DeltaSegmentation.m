%% Segment dela angles
subjects=17:22;

% Load Delta Phase for All Trials
occipital=input('Occipital (1) or Central (0)?: ');

for subj=subjects
    clearvars -except subjects subj occipital
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG'
    loadfilename1=sprintf('EEG_SxA_Subj%i_Results_deltaall.mat',subj);
    loadfilename2=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj);
    savefilename=sprintf('EEG_SxA_Subj%i_Results.mat',subj);
    load(loadfilename1)
    load(loadfilename1)
    load(loadfilename2,'SDATA')

    artifacts=SDATA.metadata.artifacts;
    triggers=SDATA.events.triggerChannel;
    channelnames=SDATA.info.channel_labels;
    srate=SDATA.info.sampling_rate;

    if occipital
        delta_phase_whole=delta_ITPC_whole_occ;
    else
        delta_phase_whole=delta_ITPC_whole_cen;
    end

%     %% Catch Only
%     triggercodes={1;2;3}; % Catch Target Signals per condition
%     timerange=[-1400 600];
% 
%     for c=1:size(triggercodes,1) % For conditions
%         % Segmentation
%         [Delta_phaseSingleTrials, isNotArtifact, timeVecITPC]=segmentContEEGdata(triggercodes{c} , timerange,...
%             delta_phase_whole , triggers, artifacts, srate);
% 
%         % Check how many trials are artifact-free
%         sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))
% 
%         % Remove trials with artifacts
%         artrej=input("Remove artifact trials (1=yes)? ");
%         if artrej==1
%             Delta_phaseSingleTrials=Delta_phaseSingleTrials(:,:,isNotArtifact==1);
%         end
% 
%         if occipital
%             ITPCdelta_nTrials_catch_occ(c)=size(Delta_phaseSingleTrials,3);
%             ITPCdelta_res_catch_occ{c}=circ_r(Delta_phaseSingleTrials,[], [], 3);
%             ITPCdelta_timevec_catch_occ{c}=timeVecITPC;
%         else
%             ITPCdelta_nTrials_catch_cen(c)=size(Delta_phaseSingleTrials,3);
%             ITPCdelta_res_catch_cen{c}=circ_r(Delta_phaseSingleTrials,[], [], 3);
%             ITPCdelta_timevec_catch_cen{c}=timeVecITPC;
%         end
% 
%         clear timeVecITPC Delta_phaseSingleTrials isNotArtifact artrej
%     end
% 
%     % Save
%     cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
%     if occipital
%         save(savefilename, 'ITPCdelta_res_catch_occ', 'ITPCdelta_timevec_catch_occ', 'ITPCdelta_nTrials_catch_occ',"-append");
%     else
%         save(savefilename, 'ITPCdelta_res_catch_cen', 'ITPCdelta_timevec_catch_cen', 'ITPCdelta_nTrials_catch_cen',"-append");
%     end
% 
%     %% Pre-Target Time
%     % Cut segments around (predicted) target onset. Done separately for
%     % predictable and irregular conditions, because for irregular we need to
%     % take into account that we only include the central targets, and then we
%     % need to cut around them always with the same time frame, regardless of
%     % whether the target was earlier or later.
% 
%     % Not sure this is useful, should be the same as cutting from warning
%     % signal.
% 
%     if ~ismember(subj,[14 15]) % skip for 14 and 15 because they have no target triggers
%         triggercodes1={[1 10:29];[2 210:229]}; % Target trials (only central ones for irregular)
%         triggercodes2={[3 120:129 170:179];[110:119 160:169];[130:139 180:189]};
%         timerange1=[-1400 800]; % mask offset 1 second after target was shown
%         timerange2=[-1400 800; -1320 880;-1480 720]; % making sure the segment is always the same time point, regardless of target onset
% 
%         % Rhythm and Interval
%         for c=1:size(triggercodes1,1) % For conditions
%             % Segmentation
%             [Delta_phaseSingleTrials, isNotArtifact, timeVecITPC]=segmentContEEGdata(triggercodes1{c} , timerange1,...
%                 delta_phase_whole , triggers, artifacts, srate);
% 
%             % Check how many trials are artifact-free
%             sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))
% 
%             % Remove trials with artifacts
%             artrej=input("Remove artifact trials (1=yes)? ");
%             if artrej==1
%                 Delta_phaseSingleTrials=Delta_phaseSingleTrials(:,:,isNotArtifact==1);
%             end
% 
%             if occipital
%                 ITPCdelta_nTrials_target_pred_occ(c)=size(Delta_phaseSingleTrials,3);
%                 ITPCdelta_res_target_pred_occ{c}=circ_r(Delta_phaseSingleTrials,[], [], 3);
%                 ITPCdelta_timevec_target_pred_occ{c}=timeVecITPC;
%             else
%                 ITPCdelta_nTrials_target_pred_cen(c)=size(Delta_phaseSingleTrials,3);
%                 ITPCdelta_res_target_pred_cen{c}=circ_r(Delta_phaseSingleTrials,[], [], 3);
%                 ITPCdelta_timevec_target_pred_cen{c}=timeVecITPC;
%             end
% 
%             clear timeVecITPC Delta_phaseSingleTrials isNotArtifact artrej
%         end
% 
%         % Irregular
%         for c=1:size(triggercodes2,1) % For conditions
%             % Segmentation
%             [Delta_phaseSingleTrials, isNotArtifact, timeVecITPC]=segmentContEEGdata(triggercodes2{c} , timerange2(c,:),...
%                 delta_phase_whole , triggers, artifacts, srate);
% 
%             % Check how many trials are artifact-free
%             sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))
% 
%             % Remove trials with artifacts
%             artrej=input("Remove artifact trials (1=yes)? ");
%             if artrej==1
%                 Delta_phaseSingleTrials=Delta_phaseSingleTrials(:,:,isNotArtifact==1);
%             end
% 
%             if occipital
%                 ITPCdelta_nTrials_target_irr_occ(c)=size(Delta_phaseSingleTrials,3);
%                 ITPCdelta_res_target_irr_occ{c}=circ_r(Delta_phaseSingleTrials,[], [], 3);
%                 ITPCdelta_timevec_target_irr_occ{c}=timeVecITPC;
%             else
%                 ITPCdelta_nTrials_target_irr_cen(c)=size(Delta_phaseSingleTrials,3);
%                 ITPCdelta_res_target_irr_cen{c}=circ_r(Delta_phaseSingleTrials,[], [], 3);
%                 ITPCdelta_timevec_target_irr_cen{c}=timeVecITPC;
%             end
% 
%             clear timeVecITPC Delta_phaseSingleTrials isNotArtifact artrej
%         end
% 
%         % Save
%         cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
%         if occipital
%             %save(savefilename, 'ITPCdelta_res_target_pred_occ', 'ITPCdelta_timevec_target_pred_occ', 'ITPCdelta_nTrials_target_pred_occ',"-append");
%             save(savefilename, 'ITPCdelta_res_target_irr_occ', 'ITPCdelta_timevec_target_irr_occ', 'ITPCdelta_nTrials_target_irr_occ',"-append");
%         else
%             %save(savefilename, 'ITPCdelta_res_target_pred_cen', 'ITPCdelta_timevec_target_pred_cen', 'ITPCdelta_nTrials_target_pred_cen',"-append");
%             save(savefilename, 'ITPCdelta_res_target_irr_cen', 'ITPCdelta_timevec_target_irr_cen', 'ITPCdelta_nTrials_target_irr_cen',"-append");
%         end
%     end
%% All Trials whole trial, to see ITPC across the whole trial
    triggercodes=31; % onset rhythm trial )the only one for which the length stazs mostly consistent
    timerange=[0 4000];

    for c=1:size(triggercodes,1) % For conditions
        % Segmentation
        [Delta_phaseSingleTrials, isNotArtifact, timeVecITPC]=segmentContEEGdata(triggercodes{c} , timerange,...
            delta_phase_whole , triggers, artifacts, srate);

        % Check how many trials are artifact-free
        sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))

        % Remove trials with artifacts
        artrej=input("Remove artifact trials (1=yes)? ");
        if artrej==1
            Delta_phaseSingleTrials=Delta_phaseSingleTrials(:,:,isNotArtifact==1);
        end

        if occipital
            ITPCdelta_nTrials_catch_occ(c)=size(Delta_phaseSingleTrials,3);
            ITPCdelta_res_catch_occ{c}=circ_r(Delta_phaseSingleTrials,[], [], 3);
            ITPCdelta_timevec_catch_occ{c}=timeVecITPC;
        else
            ITPCdelta_nTrials_catch_cen(c)=size(Delta_phaseSingleTrials,3);
            ITPCdelta_res_catch_cen{c}=circ_r(Delta_phaseSingleTrials,[], [], 3);
            ITPCdelta_timevec_catch_cen{c}=timeVecITPC;
        end

        clear timeVecITPC Delta_phaseSingleTrials isNotArtifact artrej
    end

%     % Save
%     cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
%     if occipital
%         save(savefilename, 'ITPCdelta_res_catch_occ', 'ITPCdelta_timevec_catch_occ', 'ITPCdelta_nTrials_catch_occ',"-append");
%     else
%         save(savefilename, 'ITPCdelta_res_catch_cen', 'ITPCdelta_timevec_catch_cen', 'ITPCdelta_nTrials_catch_cen',"-append");
%     end
end