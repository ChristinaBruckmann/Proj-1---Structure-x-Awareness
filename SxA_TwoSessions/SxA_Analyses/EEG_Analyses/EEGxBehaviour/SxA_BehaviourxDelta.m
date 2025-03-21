%% Threshold Delta Correlation
% How does the difference in perceptual threshold relate to delta ITPC
% between conditions?
clear
clc

subj=[101:103 105:106 108 110 112 113 114 116 117 118 119 121:124 126 127 129 130];
delta_TW=[800 900]; % Delta time window used for analysis, WS at 0
cluster=1; % (1-occipital, 2- central)
catchonly=0; % only catch trials?
bl_correct=1; % baseline delta itpc to theoretical basline?

if cluster==1
    elec=[25:30 62:64]; % occipital electrodes
else
    elec=[11:12 46:49]; % central electrodes
end
%% Get data
for s=1:length(subj)

    % Behavioural Data
    cd 'Y:\el-Christina\SxA\SxA_Data\Behaviour Preprocessed'
    loadname=sprintf('SxA_ResultsSubject%i_Total.mat',subj(s));
    load(loadname,'midpoints')
   
    midpoints_obj(s,1)=midpoints.rhyobj; % Objective
    midpoints_obj(s,2)=midpoints.intobj;
    midpoints_obj(s,3)=midpoints.irrobj;
    
    midpoints_subj(s,1)=midpoints.rhysubj; % Subjective
    midpoints_subj(s,2)=midpoints.intsubj;
    midpoints_subj(s,3)=midpoints.irrsubj;

    % Delta ITPC
    cd 'Y:\el-Christina\SxA\SxA_Results\Delta Results'
    loadname=sprintf('EEG_SxA_Subj%s_DeltaPhaseSingleTrials_NewFreq',subj(s));

    if cluster==1 && ~catchonly
        loadfilename=sprintf('NewFreq_OccipitalAll_Subj%i',subj(s));
        load(loadfilename,"Delta_SingleTrials_occ","ITPCdelta_timevec_occ","ITPCdelta_nTrials_occ") %(time points x electrodes x trials)
        for c=1:3
            Delta_ITPC(s,c,:,:)=circ_r(Delta_SingleTrials_occ{c},[], [], 3);
            Delta_TimeVec(s,c,:)=ITPCdelta_timevec_occ{c};
            Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_occ(c);
        end
    elseif cluster==1 && catchonly
        loadfilename=sprintf('NewFreq_OccipitalCatch_Subj%i',subj(s));
        load(loadfilename,"Delta_SingleTrials_occ_catch","ITPCdelta_timevec_occ_catch","ITPCdelta_nTrials_occ_catch") %(time points x electrodes x trials)
        for c=1:3
            Delta_ITPC(s,c,:,:)=circ_r(Delta_SingleTrials_occ_catch{c},[], [], 3);
            Delta_TimeVec(s,c,:)=ITPCdelta_timevec_occ_catch{c};
            Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_occ_catch(c);
        end
    elseif cluster==2 && ~catchonly
        loadfilename=sprintf('NewFreq_CentralAll_Subj%i',subj(s));
        load(loadfilename,"Delta_SingleTrials_cen","ITPCdelta_timevec_cen","ITPCdelta_nTrials_cen") %(time points x electrodes x trials)
        for c=1:3
            Delta_ITPC(s,c,:,:)=circ_r(Delta_SingleTrials_cen{c},[], [], 3);
            Delta_TimeVec(s,c,:)=ITPCdelta_timevec_cen{c};
            Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_cen(c);
        end
    elseif cluster==2 && catchonly
        loadfilename=sprintf('NewFreq_CentralCatch_Subj%i',subj(s));
        load(loadfilename,"Delta_SingleTrials_cen_catch","ITPCdelta_timevec_cen_catch","ITPCdelta_nTrials_cen_catch") %(time points x electrodes x trials)
        for c=1:3
            Delta_ITPC(s,c,:,:)=circ_r(Delta_SingleTrials_cen_catch{c},[], [], 3);
            Delta_TimeVec(s,c,:)=ITPCdelta_timevec_cen_catch{c};
            Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_cen_catch(c);
        end
    end
end

%% Prep and clean delta data
clearvars -except Delta_Ntrials Delta_TimeVec Delta_ITPC midpoints_subj midpoints_obj elec delta_TW bl_correct

% Merge timeVec (same for all cond and subj) and get time idx
Delta_TimeVec=squeeze(Delta_TimeVec(1,1,:));
time_idx=Delta_TimeVec>delta_TW(1)&Delta_TimeVec<delta_TW(2);

% Get Delta ITPC for elecs and time window (subj, cond, timepoints, elecs)
Delta_ITPC_clean=Delta_ITPC(:,:,time_idx,elec);
Delta_ITPC_clean=mean(Delta_ITPC_clean,4); % average across electrodes
Delta_ITPC_clean=mean(Delta_ITPC_clean,3); % average across time points

% Baseline correct to theoretical baseline
if bl_correct
    for s=1:size(Delta_ITPC_clean,1)  % for each subject
        for c=1:size(Delta_ITPC_clean,2) % for each condition
            % calculate chance angles
            angles=2*pi*rand(Delta_Ntrials(s,c),10000);
            angle_means(s,c)=mean(circ_r(angles));
            angle_std(s,c)=std(circ_r(angles));
        end
    end
    % Subtract baseline from each condition
    Delta_ITPC_clean=Delta_ITPC_clean-angle_means;
end

%% Correlate
thr_diff_rhyint_obj=midpoints_obj(:,1)-midpoints_obj(:,2);
thr_diff_rhyint_subj=midpoints_subj(:,1)-midpoints_subj(:,2);

itpc_diff_rhyint=Delta_ITPC_clean(:,1)-Delta_ITPC_clean(:,2);

% Correlation
[R_obj_subj,P_obj_subj] = corrcoef(thr_diff_rhyint_obj,thr_diff_rhyint_subj);
[R_obj_itpc,P_obj_itpc] = corrcoef(thr_diff_rhyint_obj,itpc_diff_rhyint);
[R_subj_itpc,P_subj_itpc] = corrcoef(thr_diff_rhyint_subj,itpc_diff_rhyint);

% Plot
figure("Position",[311.4000 514.6000 1504 420]); tiledlayout("flow")
nexttile;
mdl = fitlm(table(thr_diff_rhyint_subj, thr_diff_rhyint_obj)); % Fit regression line
plot(mdl,"LineWidth",3)
ylabel("Objective Threshold Difference")
xlabel("Subjective Threshold Difference")
title(sprintf("Difference Rhythm-Interval Thresholds (R=%.4f)",R_obj_subj(2)))

nexttile;
mdl = fitlm(table(thr_diff_rhyint_obj,itpc_diff_rhyint)); % Fit regression line
plot(mdl,"LineWidth",3)
ylabel("ITPC Difference")
xlabel("Objective Threshold Difference")
title(sprintf("Difference Rhythm-Interval ObjxITPC (R=%.4f)",R_obj_itpc(2)))

nexttile;
mdl = fitlm(table(thr_diff_rhyint_subj, itpc_diff_rhyint)); % Fit regression line
plot(mdl,"LineWidth",3)
ylabel("ITPC Difference")
xlabel("Subjective Threshold Difference")
title(sprintf("Difference Rhythm-Interval SubjxITPC (R=%.4f)",R_subj_itpc(2)))