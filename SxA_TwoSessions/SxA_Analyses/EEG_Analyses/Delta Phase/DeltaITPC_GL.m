%% Delta ITPC Group Level

clear
clc
subj=[17:22];
%subj=[101:103 105:108 110 112:114 116:119 121 122 124 126 127];
cd 'Z:\el-Christina\SxA\SxA_Results\Delta Results'
cluster=2; %(1 occipital, 2 central)
catchonly=0; %(1-yes)

% Time Window for Statistics
stats_tw=[700 800]; %from WS

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
SubjMean=squeeze(mean(Delta_SingleTrials(:,:,:,elec),4)); % average across electrodes
Delta=squeeze(mean(SubjMean,1)); % mean across subjects and electrodes

%% Plot
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

%% Statistics in ROI Time Window
delta_tw=SubjMean(:,:,timeVec>=stats_tw(1)&timeVec<=stats_tw(2)); % select only time window
delta_tw=mean(delta_tw,3); % average across time points

% Rhythm higher than irregular?
[h_delta(1),p_delta(1),~,stats_delta(1,:)] = ttest(delta_tw(:,1),delta_tw(:,3),"Tail","right");
% Interval higher than irregular?
[h_delta(2),p_delta(2),~,stats_delta(2,:)] = ttest(delta_tw(:,2),delta_tw(:,3),"Tail","right");
% Difference between rhythm and interval?
[h_delta(3),p_delta(3),~,stats_delta(3,:)] = ttest(delta_tw(:,1),delta_tw(:,2));

% Save
cd 'Z:\el-Christina\SxA\SxA_Results\Delta Results'
if cluster==1 && catchonly
save("DeltaGroupLevel_Occ_Catch",'Delta','subj','timeVec','nTrials')
elseif cluster==1 && ~catchonly
save("DeltaGroupLevel_Occ_All",'Delta','subj','timeVec','nTrials')
elseif cluster==2 && catchonly
save("DeltaGroupLevel_Cen_Catch",'Delta','subj','timeVec','nTrials')
elseif cluster==2 && ~catchonly
save("DeltaGroupLevel_Cen_All",'Delta','subj','timeVec','nTrials')
end