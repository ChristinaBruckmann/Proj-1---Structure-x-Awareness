%% What could be wrong with our alpha signal?

%% Plot individual participants (and individual trials)
clear
clc

subj=[101:103 105:108 110:114 117:119];
TF_wavFreqs=2.^[0:1/6:5];
alpharange=8:12;
alpha_freqs=find(TF_wavFreqs>=alpharange(1)&TF_wavFreqs<=alpharange(end));
alpha_tfElectrodes=[25:30 62:64]; % Occipital only

% Load data
for s=1:length(subj)
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\TimeFreq'
    loadfilename=sprintf('EEG_SxA_Subj%i_TF_SingleTrials',subj(s));
    load(loadfilename,'TF_Results_Trial','TF_trial_timeVec','TF_NotArtifact')
    alpha_res=TF_Results_Trial; % Size: timepoints frequencies trials electrodes
    alpha_time_group{s}=TF_trial_timeVec; % should be the same across participants and conditions, but just in case
    alpha_ArtVectors=TF_NotArtifact;

    % Clean
    for c=1:3
    alpha_res_group{s,c}=alpha_res{c}(:,alpha_freqs,:,alpha_tfElectrodes); % Select alpha frequencies only (time x frequencies x trials x electrodes)
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

subj=[101:103 105:108 110:114 117:119];
alpharange=8:12;

% Load data
for s=1:length(subj)
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\AlphaRes'
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
    if s==1
    legend ("Rhythm","Interval","Irregular") % only add this to one plot
    end
    title(sprintf("Alpha Amp Mask Onset Subject %i",s));
    xlim([-600 700])
end

% Calculate Mean
for s=1:length(subj)
    for c=1:3
        group_alpha(s,c,:)=alpha_res_group{s,c};
    end
end
group_mean_alpha=squeeze(mean(group_alpha,1)); % average  across subjects

% Plot Mean
figure;
%nexttile
for c=1:3
    plot(alpha_time_group{s},group_mean_alpha(c,:),"LineWidth",3)
    hold on
end
xline(0,'--','Mask Onset');
xlim([-600 700])
legend ("Rhythm","Interval","Irregular")
title("MEAN Alpha Amp Mask Onset");
%% Compare first and second interval
clear
clc

subj=[101:103 105:108 110:114 117:119];
alpharange=8:12;

% Load data
for s=1:length(subj)
    % First Interval (Cue)
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\AlphaRes'
    loadfilename1=sprintf('EEG_SxA_Subj%i_AlphaResults_FirstInterval.mat',subj(s));
    load(loadfilename1,'alpha_Results','alpha_timeVecTotal')
    alpha_res_firstinterval=alpha_Results{1}; % Size: timepoints frequencies trials electrodes
    alpha_time_group_first{s}=alpha_timeVecTotal{1}; % should be the same across participants and conditions, but just in case

    % Clean
    alpha_res_group_first{s}=alpha_res_firstinterval(:,alpharange,:); % Take only alpha freqs
    alpha_res_group_first{s}=squeeze(mean(alpha_res_group_first{s},3));% Average across electrodes (time x frequencies x trials)
    alpha_res_group_first{s}=squeeze(mean(alpha_res_group_first{s},2)); % Average across frequencies (time x trials)

     % Second Interval (Target)
    loadfilename2=sprintf('EEG_SxA_Subj%i_AlphaResults_clean.mat',subj(s));
    load(loadfilename2,'alpha_Results','alpha_timeVecTotal')
    alpha_res_secondinterval=alpha_Results{1, 2} ;  %Get only interval data // Size: timepoints frequencies trials electrodes
    alpha_time_group_second{s}=alpha_timeVecTotal{1}; % should be the same across participants and conditions, but just in case

    % Clean
    alpha_res_group_second{s}=alpha_res_secondinterval(:,alpharange,:); % Take only alpha freqs
    alpha_res_group_second{s}=squeeze(mean(alpha_res_group_second{s},3));% Average across electrodes (time x frequencies x trials)
    alpha_res_group_second{s}=squeeze(mean(alpha_res_group_second{s},2)); % Average across frequencies (time x trials)
end

% Plot
figure;
for s=1:length(subj)
    nexttile
    % Plot first interval
    plot(alpha_time_group_first{s},alpha_res_group_first{s},"LineWidth",3)
    hold on

    % Plot second interval
    plot(alpha_time_group_second{s},alpha_res_group_second{s},"LineWidth",3)

    xline(0,'--','First Stimulus');
    xline(800,'--','Second Stimulus');
    if s==1
     legend ("First Interval","Second Interval")
    end
    title(sprintf("Alpha Amp Across Intervals Subject %i",s));
    xlim([-200 1500])
end

% Calculate Mean
for s=1:length(subj)
        group_alpha_first(s,:)=alpha_res_group_first{s};
        group_alpha_second(s,:)=alpha_res_group_second{s};
end
group_mean_alpha_first=squeeze(mean(group_alpha_first,1)); % average  across subjects
group_mean_alpha_second=squeeze(mean(group_alpha_second,1)); % average  across subjects

% Plot Mean
%nexttile
figure;
plot(alpha_time_group_first{s},group_mean_alpha_first,"LineWidth",3) % Plot first interval
hold on
plot(alpha_time_group_second{s},group_mean_alpha_second,"LineWidth",3) % Plot second interval

xline(0,'--','First Stimulus');
xline(800,'--','Second Stimulus');
legend ("First Interval","Second Interval")
title("MEAN Alpha Amp Across Intervals");
xlim([-200 1500])

%% Only analyze 800ms target trials
% maybe including 800,850 and 900 trials smears the signal
% taking now only irregular 800
clear
clc

subj=[101:103 105:108 110:114 117:119];
alpharange=8:12;

% Load data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\AlphaRes\Target at 800 only'
for s=1:length(subj)
    loadfilename=sprintf('EEG_SxA_Subj%i_AlphaResults_800.mat',subj(s));
    load(loadfilename,'alpha_Results','alpha_ntrials','alpha_timeVecTotal')
    alpha_trialn(s,:)=alpha_ntrials;
    for c=1:3
        alpha_res(s,c,:,:,:)=alpha_Results{1, c}; % subj x condition x time points x frequencies x electrodes
        alpha_timevec(s,c,:)=alpha_timeVecTotal{1, c};
    end
end

% Take only alpha band and average across electrodes
alpha_res=mean(alpha_res(:,:,:,alpharange,:),5); % output:   subj x condition x time points x frequencies
alpha_res=mean(alpha_res,4); %average across frequencies

% Plot Individuals
figure;
for s=1:length(subj)
    nexttile
    for c=1:3
        plot(squeeze(alpha_timevec(s,c,:)),squeeze(alpha_res(s,c,:)),"LineWidth",3)
        hold on
    end
    xline(0,'--','Warning Signal');
    xline(800,'--','Target');
    if s==1
    legend ("Rhythm","Interval","Irregular")
    end
    title(sprintf("Alpha Amp 800ms Targets Subject %i",s));
    xlim([-200 1500])
end

% Plot Mean
alpha_res_mean=squeeze(mean(alpha_res(:,:,:),1));

figure;
for c=1:3
    plot(squeeze(alpha_timevec(1,c,:)),squeeze(alpha_res_mean(c,:)),"LineWidth",3)
    hold on
end
xline(0,'--','Warning Signal');
xline(800,'--','Target');
legend ("Rhythm","Interval","Irregular")
title("MEAN Alpha Amp Only 800ms Targets");
xlim([-200 1500])

%% Jackknifing
% Maybe one participant drives the irregular down?
clear
clc

subj=[101:103 105:108 110:114 117:119];
alpharange=8:12;

for s=1:length(subj)
    % Load Data
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\AlphaRes'
    loadfilename=sprintf('EEG_SxA_Subj%i_AlphaResults_clean.mat',subj(s));
    load(loadfilename,'alpha_Results','alpha_timeVecTotal')
    alpha_res=alpha_Results; % Size: timepoints frequencies trials electrodes
    alpha_time_group{s}=alpha_timeVecTotal{1}; % should be the same across participants and conditions, but just in case

    % Clean
    for c=1:3
    alpha_res_group_temp1(s,c,:,:,:)=alpha_res{c}(:,alpharange,:); % Take only alpha freqs
    alpha_res_group_temp2(s,c,:,:)=squeeze(mean(alpha_res_group_temp1(s,c,:,:,:),5));% Average across electrodes (time x frequencies)
    alpha_res_group(s,c,:)=squeeze(mean(alpha_res_group_temp2(s,c,:,:),4)); % Average across frequencies (time)
    end
end

figure;
for s=1:length(subj)
    % Jackknife
    curr_alpha=alpha_res_group(:,:,:);
    curr_alpha(s,:,:)=[]; % Remove Subject
    curr_alpha_mean=squeeze(mean(curr_alpha,1)); % Calculate mean of remaining subjects

    % Plot
    nexttile
    for c=1:3
        plot(alpha_time_group{s},curr_alpha_mean(c,:),"LineWidth",3)
        hold on
    end

    xline(0,'--','Warning Signal');
    xline(800,'--','Target');
    if s==1
    legend ("Rhythm","Interval","Irregular")
    end
    title(sprintf("Mean Alpha Without Subj %i",s));
    xlim([-200 1500])
end

%% Topography
% Are we selecting the correct electrodes?

clear
clc

subj=[101:103 105:108 110:114 117:119];
TF_wavFreqs=2.^[0:1/6:5];
alpharange=8:12;
alpha_freqs=find(TF_wavFreqs>=alpharange(1)&TF_wavFreqs<=alpharange(end));
timewindow=[650 800];

% Load Data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\AlphaRes'
for s=1:length(subj)
    loadfilename=sprintf('EEG_SxA_Subj%i_TF_Results',subj(s));
    load(loadfilename)
    for c=1:3
        TF_Res(s,c,:,:,:)=TF_Results{c};
        TF_Time(s,c,:)=TF_timeVecTotal{c};
    end
end

% Clean (take time window and alpha freq only, average across both)
timevec=squeeze(TF_Time(1,1,:));
TF_Res_clean=TF_Res(:,:,timevec>=timewindow(1)&timevec<=timewindow(2),alpha_freqs,:);
TF_Res_clean=squeeze(mean(TF_Res_clean,3)); % average across time window
TF_Res_clean=squeeze(mean(TF_Res_clean,3));  % average across alpha frequencies
TF_Res_clean=squeeze(mean(TF_Res_clean,1)); % average across subjects

% Plot topography for each condition
figure;
title('Alpha Amp 650-750ms pre target')
for c=1:3
    nexttile
    data_to_plot=double(TF_Res_clean(c,:));
    topoplot(data_to_plot,'head71.locs','electrodes','on','maplimits',[-0 3]); % ,'maplimits',[-0.2 0.2]
    colorbar
end

%% Maybe only block beginnings reflect true alpha amplitude before mask suppression?
% Includes Artifacts
clear
clc

subj=[101:103 105:108 110:114 117:119];
alpharange=8:12;

% Load Data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\AlphaRes'
for s=1:length(subj)
    loadfilename=sprintf('EEG_SxA_Subj%i_AlphaResults_BlockBegin.mat',subj(s));
    load(loadfilename)
    for c=1:3
    group_block(s,c,:,:,:)=alpha_Results{c}; % time x freq x elec
    group_timevec(s,c,:)=alpha_timeVecTotal{c};
    end
end

% Clean Data
group_block=squeeze(mean(group_block(:,:,:,alpharange,:),5)); % select alpha only and average across elec
group_block=squeeze(mean(group_block,4)); % average across electrodes

% Plot Individuals
figure;
for s=1:length(subj)
    nexttile
    for c=1:3
        plot(squeeze(group_timevec(s,c,:)),squeeze(group_block(s,c,:)),"LineWidth",3)
        hold on
    end
    xline(0,'--','First Trial');
    if s==1
    legend ("Rhythm","Interval","Irregular")
    end
    title(sprintf("Alpha Amp Block Beginning Subject %i",s));
    xlim([-4000 200])
end

% Plot Mean
group_block_mean=squeeze(mean(group_block(:,:,:),1));

figure;
for c=1:3
    plot(squeeze(group_timevec(1,c,:)),squeeze(group_block_mean(c,:)),"LineWidth",3)
    hold on
end

xline(0,'--','First Trial');
legend ("Rhythm","Interval","Irregular")
title("MEAN Alpha Amp Block Beginning");
xlim([-4000 200])


%% Ratio instead of absolute suppression (divide by baseline)
clear
clc

subj=[101:103 105:108 110:114 117:119];
alpharange=8:12;
baseline=[0 100]; % in relation to WS=0

% Load data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\AlphaRes'
for s=1:length(subj)
    loadfilename=sprintf('EEG_SxA_Subj%i_AlphaResults_clean.mat',subj(s));
    load(loadfilename)
    for c=1:3
    group_alpha(s,c,:,:,:)=alpha_Results{c}; % time x freq x elec
    group_timevec(s,c,:)=alpha_timeVecTotal{c};
    end
end

% Extract baseline
timevec=squeeze(group_timevec(1,2,:));
group_alpha=group_alpha(:,:,:,alpharange,:); % select alpha range only
group_alpha=squeeze(mean(group_alpha,5)); %Average across electrodes
group_alpha=squeeze(mean(group_alpha,4)); % Average across alpha freq

Alpha_Baseline=squeeze(mean(group_alpha(:,:,timevec>=baseline(1)&timevec<=baseline(2)),3)); % Average across BL time window


% Divide by baseline
for s=1:length(subj)
    for c=1:3
        Alpha_BC(s,c,:)=squeeze(group_alpha(s,c,:))/Alpha_Baseline(s,c);
    end
end

% Plot Individuals
figure;
for s=1:length(subj)
    nexttile
    for c=1:3
        plot(squeeze(group_timevec(s,c,:)),squeeze(Alpha_BC(s,c,:)),"LineWidth",3)
        hold on
    end
    xline(0,'--','Warning Signal');
    xline(800,'--','Target');
    if s==1
    legend ("Rhythm","Interval","Irregular")
    end
    title(sprintf("Alpha Ratio Subject %i",s));
    xlim([-250 1500])
end

% Plot Mean
Alpha_mean=squeeze(mean(Alpha_BC(:,:,:),1));

figure;
for c=1:3
    plot(squeeze(group_timevec(1,c,:)),squeeze(Alpha_mean(c,:)),"LineWidth",3)
    hold on
end

xline(0,'--','First Trial');
 xline(800,'--','Target');
legend ("Rhythm","Interval","Irregular")
title("MEAN Alpha Ratio");
xlim([-250 1500])
