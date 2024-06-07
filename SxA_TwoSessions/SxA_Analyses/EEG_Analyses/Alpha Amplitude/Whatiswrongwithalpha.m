%% What could be wrong with our alpha signal?

%% Plot individual participants (and individual trials)
clear
clc

subj=[101:103 105 106:108 110 111 113 114 117 118 119];
TF_wavFreqs=2.^[0:1/6:5];
alpharange=8:13;
alpha_freqs=find(TF_wavFreqs>=alpharange(1)&TF_wavFreqs<=alpharange(end));
alpha_tfElectrodes=[25:30 62:64]; % Occipital only

% Load data
for s=1:length(subj)
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\TimeFreq'
    loadfilename=sprintf('EEG_SxA_Subj%i_TF_SingleTrials',subj(s));
    load(loadfilename,'alpha_Results','TF_trial_timeVec','TF_NotArtifact')
    alpha_res=alpha_Results; % Size: timepoints frequencies trials electrodes
    alpha_time_group{s}=TF_trial_timeVec; % should be the same across participants and conditions, but just in case
    alpha_ArtVectors=TF_NotArtifact;

    % Clean
    for c=1:3
    alpha_res_group{s,c}=alpha_res{c}(:,alpha_freqs,:,alpha_tfElectrodes); % Remove Artifacts (time x frequencies x trials x electrodes)
    alpha_res_group{s,c}=squeeze(mean(alpha_res_group{s,c},4));% Average across electrodes (time x frequencies x trials)
    alpha_res_group{s,c}=squeeze(mean(alpha_res_group{s,c},2)); % Average across frequencies (time x trials)
    end

end

% Plot
% Now you get TF_res_group{subject, condition}(time x trials)
for c=1:size(alpha_res_group,2) % for each condition produce one plot
    figure;
    for s=1:size(alpha_res_group,1)  % with subjects as subplots
        nexttile
        timevec=alpha_time_group{1,c};
        mean_alpha=mean(alpha_res_group{s,c},2);
        % Plot with different colours for each condition
        if c==1
            % Plot single trials
            for t=1:size(alpha_res_group{s,c},2)
                trial=plot(timevec{c},alpha_res_group{s,c}(:,t),"LineWidth",2,"Color",[0 0.4470 0.7410 0.2]);
                hold on
            end
            plot(timevec{c},mean_alpha,"LineWidth",2,"Color",[0 0 0]); % Plot mean
            xline(0)
            xline(800)
            ylim([0 20])
            xlim([-200 1500])
            title(sprintf("Rhythm Subject %i",s));
            hold off
        elseif c==2
            % Plot single trials
            for t=1:size(alpha_res_group{s,c},2)
                plot(timevec{c},alpha_res_group{s,c}(:,t),"LineWidth",2,"Color",[0.8500 0.3250 0.0980 0.2]);
                hold on
            end
            plot(timevec{c},mean_alpha,"LineWidth",2,"Color",[0 0 0]); % Plot mean
            xline(0)
            xline(800)
            ylim([0 20])
            xlim([-200 1500])
            title(sprintf("Interval Subject %i",s));
            hold off
        else
            % Plot single trials
            for t=1:size(alpha_res_group{s,c},2)
                plot(timevec{c},alpha_res_group{s,c}(:,t),"LineWidth",2,"Color",[0.9290 0.6940 0.1250 0.2]);
                hold on
            end
            plot(timevec{c},mean_alpha,"LineWidth",2,"Color",[0 0 0]); % Plot mean
            xline(0)
            xline(800)
            ylim([0 20])
            xlim([-200 1500])
            title(sprintf("Irregular Subject %i",s));
            hold off
        end
    end
end

%% Plot alpha response to mask onset
clear
clc

subj=[101:103 105];
alpharange=8:13;

% Load data
for s=1:length(subj)
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\TimeFreq'
    loadfilename=sprintf('EEG_SxA_Subj%i_AlphaResults_Mask.mat',subj(s));
    load(loadfilename,'alpha_Results','alpha_timeVecTotal')
    alpha_res=alpha_Results; % Size: timepoints frequencies trials electrodes
    alpha_time_group{s}=alpha_timeVecTotal{1}; % should be the same across participants and conditions, but just in case

    % Clean
    for c=1:3
    alpha_res_group{s,c}=alpha_res{c}(:,alpharange,:); % Take only alpha freqs
    alpha_res_group{s,c}=squeeze(mean(alpha_res_group{s,c},3));% Average across electrodes (time x frequencies x trials)
    alpha_res_group{s,c}=squeeze(mean(alpha_res_group{s,c},2)); % Average across frequencies (time x trials)
    end

end

% Plot
figure;
for s=1:length(subj)
    nexttile
    for c=1:3
        plot(alpha_time_group{s},alpha_res_group{s,c},"LineWidth",3)
        hold on
    end
    xline(0,'--','Mask Onset');
    legend ("Rhythm","Interval","Irregular")
    title(sprintf("Alpha Power Mask Onset Subject %i",s));
end

% Calculate Mean
for s=1:4
    for c=1:3
        group_alpha(s,c,:)=alpha_res_group{s,c};
    end
end
group_mean_alpha=squeeze(mean(group_alpha,1)); % average  across subjects

% Plot Mean
nexttile
for c=1:3
    plot(alpha_time_group{s},group_mean_alpha(c,:),"LineWidth",3)
    hold on
end
xline(0,'--','Mask Onset');
legend ("Rhythm","Interval","Irregular")
title("MEAN Alpha Power Mask Onset");
%% Compare first and second interval
