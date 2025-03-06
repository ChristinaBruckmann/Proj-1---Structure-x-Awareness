%% Alpha ITPC Group Level

clear
clc
subj=[101:103 105:108 110 112:114 116:119 121 122 124 126 127 129 130:132];
cd 'Y:\el-Christina\SxA\SxA_Results\AlphaPhaseRes'
cluster=1; %(1 occipital, 2 central)
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
        loadfilename=sprintf('Alpha_OccipitalAll_Subj%i',subj(s));
        load(loadfilename,"Alpha_SingleTrials_occ","ITPCalpha_timevec_occ","ITPCalpha_nTrials_occ") %(time points x electrodes x trials)
        for c=1:3
        Alpha_SingleTrials(s,c,:,:)=circ_r(Alpha_SingleTrials_occ{c},[], [], 3);
        Alpha_TimeVec(s,c,:)=ITPCalpha_timevec_occ{c};
        Alpha_Ntrials(s,c,:)=ITPCalpha_nTrials_occ(c);
        end
    elseif cluster==1 && catchonly
        loadfilename=sprintf('Alpha_OccipitalCatch_Subj%i',subj(s));
        load(loadfilename,"Alpha_SingleTrials_occ_catch","ITPCalpha_timevec_occ_catch","ITPCalpha_nTrials_occ_catch") %(time points x electrodes x trials)
        for c=1:3
        Alpha_SingleTrials(s,c,:,:)=circ_r(Alpha_SingleTrials_occ_catch{c},[], [], 3);
        Alpha_TimeVec(s,c,:)=ITPCalpha_timevec_occ_catch{c};
        Alpha_Ntrials(s,c,:)=ITPCalpha_nTrials_occ_catch(c);
        end
    elseif cluster==2 && ~catchonly
        loadfilename=sprintf('Alpha_CentralAll_Subj%i',subj(s));
        load(loadfilename,"Alpha_SingleTrials_cen","ITPCalpha_timevec_cen","ITPCalpha_nTrials_cen") %(time points x electrodes x trials)
        for c=1:3
        Alpha_SingleTrials(s,c,:,:)=circ_r(Alpha_SingleTrials_cen{c},[], [], 3);
        Alpha_TimeVec(s,c,:)=ITPCalpha_timevec_cen{c};
        Alpha_Ntrials(s,c,:)=ITPCalpha_nTrials_cen(c);
        end
    elseif cluster==2 && catchonly
        loadfilename=sprintf('Alpha_CentralCatch_Subj%i',subj(s));
        load(loadfilename,"Alpha_SingleTrials_cen_catch","ITPCalpha_timevec_cen_catch","ITPCalpha_nTrials_cen_catch") %(time points x electrodes x trials)
        for c=1:3
        Alpha_SingleTrials(s,c,:,:)=circ_r(Alpha_SingleTrials_cen_catch{c},[], [], 3);
        Alpha_TimeVec(s,c,:)=ITPCalpha_timevec_cen_catch{c};
        Alpha_Ntrials(s,c,:)=ITPCalpha_nTrials_cen_catch(c);
        end
    end
end

%Average
timeVec=squeeze(Alpha_TimeVec(1,1,:)); % All TV are identical
nTrials=sum(Alpha_Ntrials,1); % nTrials for chance ITPC calculation
SubjMean=squeeze(mean(Alpha_SingleTrials(:,:,:,elec),4)); % average across electrodes
Alpha=squeeze(mean(SubjMean,1)); % mean across subjects and electrodes

%% Plot
figure;
for c=1:3
    plot(timeVec,Alpha(c,:),"LineWidth",2)
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
    xline(-900)
    title('Catch Trials')
else
    xline(900)
    title('All Trials')
end
xline(0)

%% Statistics in ROI Time Window
alpha_tw=SubjMean(:,:,timeVec>=stats_tw(1)&timeVec<=stats_tw(2)); % select only time window
alpha_tw=mean(alpha_tw,3); % average across time points

% Rhythm higher than irregular?
[h_alpha(1),p_alpha(1),~,stats_alpha(1,:)] = ttest(alpha_tw(:,1),alpha_tw(:,3),"Tail","right");
% Interval higher than irregular?
[h_alpha(2),p_alpha(2),~,stats_alpha(2,:)] = ttest(alpha_tw(:,2),alpha_tw(:,3),"Tail","right");
% Difference between rhythm and interval?
[h_alpha(3),p_alpha(3),~,stats_alpha(3,:)] = ttest(alpha_tw(:,1),alpha_tw(:,2));

% Save
cd 'Z:\el-Christina\SxA\SxA_Results\AlphaPhaseRes'
if cluster==1 && catchonly
save("AlphaGroupLevel_Occ_Catch",'Alpha','subj','timeVec','nTrials')
elseif cluster==1 && ~catchonly
save("AlphaGroupLevel_Occ_All",'Alpha','subj','timeVec','nTrials')
elseif cluster==2 && catchonly
save("AlphaGroupLevel_Cen_Catch",'Alpha','subj','timeVec','nTrials')
elseif cluster==2 && ~catchonly
save("AlphaGroupLevel_Cen_All",'Alpha','subj','timeVec','nTrials')
end