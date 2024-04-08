%% Segment dela angles
clear
clc
subjects=17:22;

% Load Delta Phase for All Trials
occipital=input('Occipital (1) or Central (0)?: ');

for s=1:length(subjects)
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG'
    loadfilename1=sprintf('EEG_SxA_Subj%i_Results_deltaall.mat',subjects(s));
    loadfilename2=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subjects(s));
    savefilename=sprintf('EEG_SxA_Subj%i_Results.mat',subjects(s));
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


    triggercodes={41}; % first rhythm cue (only condition with always the same trial length
    timerange=[-800 4000];

    for c=1:size(triggercodes,1) % For conditions
        % Segmentation
        [Delta_phaseSingleTrials, isNotArtifact, timeVecITPC]=segmentContEEGdata(triggercodes{c} , timerange,...
            delta_phase_whole , triggers, artifacts, srate);

        % Check how many trials are artifact-free
        sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))

        % Remove trials with artifacts
%         artrej=input("Remove artifact trials (1=yes)? ");
%         if artrej==1
            Delta_phaseSingleTrials=Delta_phaseSingleTrials(:,:,isNotArtifact==1);
%         end

        if occipital
            ITPCdelta_nTrials_occ(c)=size(Delta_phaseSingleTrials,3);
            ITPCdelta_res_occ{c}=circ_r(Delta_phaseSingleTrials,[], [], 3);
            ITPCdelta_timevec_occ{c}=timeVecITPC;
        else
            ITPCdelta_nTrials_cen(c)=size(Delta_phaseSingleTrials,3);
            ITPCdelta_res_cen{c}=circ_r(Delta_phaseSingleTrials,[], [], 3);
            ITPCdelta_timevec_cen{c}=timeVecITPC;
        end

        clear timeVecITPC Delta_phaseSingleTrials isNotArtifact artrej
    end
    for c=1:size(triggercodes,1) % For conditions
        if occipital
        % Occipital
        ITPCdelta_res_occ_gl(s,c,:,:)=ITPCdelta_res_occ{1,c};
        ITPCdelta_timevec_occ_gl(s,c,:)=ITPCdelta_timevec_occ{1,c};
        ITPCdelta_nTrialsCond_occ_gl(s,:)=ITPCdelta_nTrials_occ;
        else
        % Central
        ITPCdelta_res_cen_gl(s,c,:,:)=ITPCdelta_res_cen{1,c};
        ITPCdelta_timevec_cen_gl(s,c,:)=ITPCdelta_timevec_cen{1,c};
        ITPCdelta_nTrialsCond_cen_gl(s,:)=ITPCdelta_nTrials_cen;
        end
    end
    ITPCdelta_res_occ=[];
    ITPCdelta_timevec_occ=[];
    ITPCdelta_nTrials_occ=[];
end

% Average
for c=1:size(triggercodes,1)
    if occipital
    gl_itpc_results_means_occ(c,:,:)=mean(ITPCdelta_res_occ_gl(:,c,:,:),1);
    else
    gl_itpc_results_means_cen(c,:,:)=mean(ITPCdelta_res_cen_gl(:,c,:,:),1);
    end
end

% Investigate amount of trials and chance angles
for c=1:size(triggercodes,1)
    % Total N of trials
    totaltrials(c)=round(mean(ITPCdelta_nTrialsCond_occ_gl(:,c),1)); % same for occ and cen

    % Change angles
    angles=2*pi*rand(totaltrials(c),5000);
    angle_mean(c)=mean(circ_r(angles));

    angle_std(c)=std(circ_r(angles));

end

% Plot
topotimeROI = [-100 -50]; % 0 is predicted target onset
warningline=[-800];
targetline=[0];
timeVector=ITPCdelta_timevec_occ_gl(1,1,:); % Same for everything
occElectrodes=[25:27 62:64]; %Po3/Po4/po7/po8/o1/o2
cenElectrodes= [11 38 46:48]; %Fz, FC1, FCz, FC2, Cz

f2=figure;
f1=figure;
for clusters=1
    deltaplotdata=[];
    deltatopodata=[];
    for c=1:size(triggercodes,1)
        if clusters==1
            deltaplotdata(c,:,:)=squeeze(gl_itpc_results_means_occ(c,:,occElectrodes));
            deltatopodata(c,:,:)=squeeze(gl_itpc_results_means_occ(c,:,:));
            titletext='ITPC Rhythm - Occipital';
            set(0, 'currentfigure', f1); 
            subplot(2,1,2)
        else
            deltaplotdata(c,:,:)=squeeze(gl_itpc_results_means_cen(c,:,cenElectrodes));
            deltatopodata(c,:,:)=squeeze(gl_itpc_results_means_cen(c,:,:));
            titletext='ITPC Delta Catch Trials - Central';
            set(0, 'currentfigure', f2); 
            subplot(2,1,2)
        end
    end

    % Time Course
    for c=1:size(triggercodes,1) % For conditions
        plot(squeeze(timeVector)', mean(deltaplotdata(c,:,:),3));
        hold on
        ylim([-0.1 0.4]);
        legend('Rhythm','Interval','Irregular')
    end
    title(titletext)
    yline(angle_mean(1),'b')
    xline(0,'label','Cue 1')
    xline(800,'label','Cue 2')
    xline(1600,'label','Cue 3')
    xline(2400,'label','WS')
    xline(3200,'label','Target')
end
