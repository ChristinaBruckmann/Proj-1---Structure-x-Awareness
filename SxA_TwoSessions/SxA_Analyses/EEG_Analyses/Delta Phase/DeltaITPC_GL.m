%% Delta ITPC Group Level

clear
clc
%subj=[17:22 101:103 105:106];
subj=[101:103 105:108 113 114 117 118 119];
 cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\DeltaRes'
cluster=1; %(1 occipital, 2 central)
catchonly=0; %(1-yes)

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
 cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\DeltaRes'
if cluster==1 && catchonly
save("DeltaGroupLevel_Occ_Catch",'Delta','subj','timeVec','nTrials')
elseif cluster==1 && ~catchonly
save("DeltaGroupLevel_Occ_All",'Delta','subj','timeVec','nTrials')
elseif cluster==2 && catchonly
save("DeltaGroupLevel_Cen_Catch",'Delta','subj','timeVec','nTrials')
elseif cluster==2 && ~catchonly
save("DeltaGroupLevel_Cen_All",'Delta','subj','timeVec','nTrials')
end