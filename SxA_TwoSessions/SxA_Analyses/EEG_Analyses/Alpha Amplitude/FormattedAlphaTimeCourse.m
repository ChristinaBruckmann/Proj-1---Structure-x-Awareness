%% Plot the time course for alpha amplitude across trial if possible
% Requires Whole Rhythm Alpha Data for all subj
clear
clc

subj=[101   102   103   105   106   107   108   110   111 112   113   114  116  117   118   119  121   122  123 124   126   127   129   130   131   132];

%% Load Data

% Rhythm Data
for s=1:length(subj)
    % Load Data (Rhythm Trial)
    cd 'Y:\el-Christina\SxA\SxA_Results\AlphaPowerRes'
    load(sprintf("EEG_SxA_Subj%i_AlphaResults_WholeRhythmTrial.mat",subj(s)))

    % Save in GL Variable
    GL_alpha_rhy(s,:,:,:)=cell2mat(alpha_Results);
    timeVec_wholetrial=alpha_timeVecTotal{1,1}; % same across all participants

    % Load Data (Mask)
    cd 'Y:\el-Christina\SxA\SxA_Results\AlphaPowerRes'
    load(sprintf("EEG_SxA_Subj%i_AlphaResults_Mask.mat",subj(s)))

    % Save in GL Variable
    GL_alpha_rhy_mask(s,:,:,:)=alpha_Results{1,1};
    timeVec_rhy_mask=alpha_timeVecTotal{1,1}; % same across all participants
end

% Interval Data
for s=1:length(subj)
    % Load Data (Mask)
    cd 'Y:\el-Christina\SxA\SxA_Results\AlphaPowerRes'
    load(sprintf("EEG_SxA_Subj%i_AlphaResults_Mask.mat",subj(s)))

    % Save in GL Variable
    GL_alpha_int_mask(s,:,:,:)=alpha_Results{1,2};
    timeVec_int_mask=alpha_timeVecTotal{1,1}; % same across all participants

    % Load Data (First Interval)
    cd 'Y:\el-Christina\SxA\SxA_Results\AlphaPowerRes'
    load(sprintf("EEG_SxA_Subj%i_AlphaResults_FirstInterval.mat",subj(s)))

    % Save in GL Variable
    GL_alpha_firstint(s,:,:,:)=alpha_Results{1,1};
    timeVec_firstint=alpha_timeVecTotal{1,1}; % same across all participants

     % Load Data (Target Interval)
    cd 'Y:\el-Christina\SxA\SxA_Results\AlphaPowerRes'
    load(sprintf("EEG_SxA_Subj%i_AlphaResults_clean.mat",subj(s)))

    % Save in GL Variable
    GL_alpha_tarint(s,:,:,:)=alpha_Results{1,2};
    timeVec_tarint=alpha_timeVecTotal{1,1}; % same across all participants
end


%% Average (in: subj x tp x freq x elec // out: 1xtp)

% Rhythm
GL_alpha_rhy=squeeze(mean(GL_alpha_rhy,4)); % across electrodes
GL_alpha_rhy=squeeze(mean(GL_alpha_rhy,3)); % across frequencies
GL_alpha_rhy=squeeze(mean(GL_alpha_rhy,1)); % across subjects 

GL_alpha_rhy_mask=squeeze(mean(GL_alpha_rhy_mask,4)); % across electrodes
GL_alpha_rhy_mask=squeeze(mean(GL_alpha_rhy_mask,3)); % across frequencies
GL_alpha_rhy_mask=squeeze(mean(GL_alpha_rhy_mask,1)); % across subjects 

% Interval
GL_alpha_int_mask=squeeze(mean(GL_alpha_int_mask,4)); % across electrodes
GL_alpha_int_mask=squeeze(mean(GL_alpha_int_mask,3)); % across frequencies
GL_alpha_int_mask=squeeze(mean(GL_alpha_int_mask,1)); % across subjects 

GL_alpha_firstint=squeeze(mean(GL_alpha_firstint,4)); % across electrodes
GL_alpha_firstint=squeeze(mean(GL_alpha_firstint,3)); % across frequencies
GL_alpha_firstint=squeeze(mean(GL_alpha_firstint,1)); % across subjects 

GL_alpha_tarint=squeeze(mean(GL_alpha_tarint,4)); % across electrodes
GL_alpha_tarint=squeeze(mean(GL_alpha_tarint,3)); % across frequencies
GL_alpha_tarint=squeeze(mean(GL_alpha_tarint,1)); % across subjects 

% Overall Y-Axis Limits
ymin=min([min(GL_alpha_rhy), min(GL_alpha_rhy_mask),min(GL_alpha_firstint), min(GL_alpha_int_mask),min(GL_alpha_tarint)]);
ymax=max([max(GL_alpha_rhy), max(GL_alpha_rhy_mask),max(GL_alpha_firstint), max(GL_alpha_int_mask),max(GL_alpha_tarint)]);
ymin=ymin-0.1;
ymax=ymax+0.1;


%% Make Combined Vectors

% Calculate Empty (Jittered) times and convert to NaN
jitter1=nan(1,300);
int_length=length([GL_alpha_firstint GL_alpha_tarint]);
rhy_length=length(GL_alpha_rhy);
jitter2=nan(1,rhy_length-int_length);

% Combine All
rhythm_total=[GL_alpha_rhy_mask jitter1 GL_alpha_rhy];
interval_total=[GL_alpha_int_mask jitter1 GL_alpha_firstint jitter2 GL_alpha_tarint];
rhy_total_timeVec=[timeVec_rhy_mask jitter1 timeVec_wholetrial];
int_total_timeVec=[timeVec_int_mask jitter1 timeVec_firstint jitter2 timeVec_tarint];


% Calculate Rhythm Events
[~,mask_onset_idx_rhy] = min(timeVec_rhy_mask-0, [], ComparisonMethod = "abs"); % find mask onset time for rhythm
mask_onset_time_rhy=zeros(1,length(timeVec_rhy_mask));
mask_onset_time_rhy(mask_onset_idx_rhy)=1; % mark the moment of mask onset time

[~,firststim_idx_rhy] = min(timeVec_wholetrial-0, [], ComparisonMethod = "abs"); % find first stimulus onset
[~,secondstim_idx_rhy] = min(timeVec_wholetrial-900, [], ComparisonMethod = "abs"); % find second stimulus onset
[~,thirdstim_idx_rhy] = min(timeVec_wholetrial-1800, [], ComparisonMethod = "abs"); % find third stimulus onset
[~,ws_idx_rhy] = min(timeVec_wholetrial-2700, [], ComparisonMethod = "abs"); % find WS onset
[~,target_idx_rhy] = min(timeVec_wholetrial-3600, [], ComparisonMethod = "abs"); % find target onset
rhy_trial_events=zeros(1,length(timeVec_wholetrial));
rhy_trial_events([firststim_idx_rhy secondstim_idx_rhy thirdstim_idx_rhy ws_idx_rhy target_idx_rhy])=1; % mark the moment of mask onset time

jitter1_events=zeros(1,length(jitter1));

rhy_events_total=[mask_onset_time_rhy  jitter1_events rhy_trial_events]; % complete rhythm event vector!
rhy_events_idx=find((rhy_events_total==1));


% Calculate Interval Events
[~,mask_onset_idx_int] = min(timeVec_int_mask-0, [], ComparisonMethod = "abs"); % find mask onset time for rhythm
mask_onset_time_int=zeros(1,length(timeVec_int_mask));
mask_onset_time_int(mask_onset_idx_int)=1; % mark the moment of mask onset time

[~,firststim_idx] = min(timeVec_firstint-0, [], ComparisonMethod = "abs"); % find first stimulus onset
[~,secondstim_idx] = min(timeVec_firstint-900, [], ComparisonMethod = "abs");% find second stimulus onset
firstint_stim_time=zeros(1,length(timeVec_firstint));
firstint_stim_time([firststim_idx secondstim_idx])=1; % mark the moment of mask onset time

[~,ws_idx] = min(timeVec_tarint-0, [], ComparisonMethod = "abs"); % find first stimulus onset
[~,target_idx] = min(timeVec_tarint-900, [], ComparisonMethod = "abs");% find second stimulus onset
tarint_stim_time=zeros(1,length(timeVec_tarint));
tarint_stim_time([ws_idx target_idx])=1; % mark the moment of mask onset time

jitter2_events=zeros(1,length(jitter2));

int_events_total=[mask_onset_time_int  jitter1_events firstint_stim_time jitter2_events tarint_stim_time]; % complete interval event vector!
int_events_idx=find((int_events_total==1));


%% Plot
figure; 

% Separate Plots
subplot(3,1,1);
plot(1:length(rhythm_total),rhythm_total,"LineWidth",3)
xlim([0 length(rhythm_total)])
box("off")
xline(rhy_events_idx(1),"Label","Mask Onset","LabelHorizontalAlignment","left","LabelVerticalAlignment","bottom")
xline(rhy_events_idx(2),"Label","First Cue","LabelHorizontalAlignment","left","LabelVerticalAlignment","bottom")
xline(rhy_events_idx(3),"Label","Second Cue","LabelHorizontalAlignment","left","LabelVerticalAlignment","bottom")
xline(rhy_events_idx(4),"Label","Third Cue","LabelHorizontalAlignment","left","LabelVerticalAlignment","bottom")
xline(rhy_events_idx(5),"Label","Warning Signal","LabelHorizontalAlignment","left","LabelVerticalAlignment","bottom")
xline(rhy_events_idx(6),"Label","Target","LineStyle","--","LabelHorizontalAlignment","left")
title("Rhythm")

subplot(3,1,2);
plot(1:length(interval_total),interval_total,"LineWidth",3,'Color',[0.85,0.33,0.10])
xlim([0 length(interval_total)])
box("off")
xline(int_events_idx(1),"Label","Mask Onset","LabelHorizontalAlignment","left","LabelVerticalAlignment","bottom")
xline(int_events_idx(2),"Label","First Cue","LabelHorizontalAlignment","left","LabelVerticalAlignment","bottom")
xline(int_events_idx(3),"Label","Second Cue","LabelHorizontalAlignment","left","LabelVerticalAlignment","bottom")
xline(int_events_idx(4),"Label","Warning Signal","LabelHorizontalAlignment","left","LabelVerticalAlignment","bottom")
xline(int_events_idx(5),"Label","Target","LineStyle","--","LabelHorizontalAlignment","left")
title("Interval")

% One Plot
subplot(3,1,3);
plot(1:length(rhythm_total),rhythm_total,"LineWidth",3)
hold on
plot(1:length(interval_total),interval_total,"LineWidth",3)
xlim([0 length(interval_total)])
box("off")
xline(rhy_events_idx(1),"Label","Mask Onset","LabelHorizontalAlignment","left","LabelVerticalAlignment","bottom")
xline(rhy_events_idx(2),"Label","First Cue","LabelHorizontalAlignment","left","LabelVerticalAlignment","bottom")
xline(rhy_events_idx(3),"Label","Second Cue","LabelHorizontalAlignment","left","LabelVerticalAlignment","bottom")
xline(rhy_events_idx(4),"Label","Third Cue","LabelHorizontalAlignment","left","LabelVerticalAlignment","bottom")
xline(rhy_events_idx(5),"Label","Warning Signal","LabelHorizontalAlignment","left","LabelVerticalAlignment","bottom")
xline(rhy_events_idx(6),"Label","Target","LineStyle","--","LabelHorizontalAlignment","left")
title("Rhythm and Interval")
%% Plot Rhythm Time Course

% Plot Mask Onset
figure; tiledlayout('flow')
nexttile;
plot(timeVec_rhy_mask,GL_alpha_rhy_mask,"LineWidth",2)
xlim([min(timeVec_rhy_mask) max(timeVec_rhy_mask)])
ylim([ymin ymax])
xline(0,"LineWidth",1,"Label","Mask Onset","LabelVerticalAlignment","bottom")
box("off")

 
%Plot Whole Rhythm Trial
nexttile; plot(timeVec_wholetrial,GL_alpha_rhy,"LineWidth",2)
xlim([min(timeVec_wholetrial) max(timeVec_wholetrial)])
xline(0,"LineWidth",1,"Label","First Cue","LabelVerticalAlignment","bottom")
xline(900,"LineWidth",1,"Label","Second Cue","LabelVerticalAlignment","bottom")
xline(1800,"LineWidth",1,"Label","Third Cue","LabelVerticalAlignment","bottom")
xline(2700,"LineWidth",1,"Label","Warning Signal","LabelVerticalAlignment","bottom")
xline(3600,'--',"LineWidth",1,"Label","Target")
ylim([ymin ymax])
box("off")

%% Plot Interval Time Course

% Plot Mask Onset
figure; tiledlayout('flow')
nexttile;
plot(timeVec_int_mask,GL_alpha_int_mask,"LineWidth",2)
xlim([min(timeVec_int_mask) max(timeVec_int_mask)])
ylim([ymin ymax])
xline(0,"LineWidth",1,"Label","Mask Onset","LabelVerticalAlignment","bottom")
box("off")
 
% Plot first interval
nexttile; plot(timeVec_firstint,GL_alpha_firstint,"LineWidth",2)
xlim([min(timeVec_firstint) max(timeVec_firstint)])
xline(0,"LineWidth",1,"Label","First Cue","LabelVerticalAlignment","bottom")
xline(900,"LineWidth",1,"Label","Second Cue","LabelVerticalAlignment","bottom")
ylim([ymin ymax])
box("off")

% Plot target interval
nexttile; plot(timeVec_tarint,GL_alpha_tarint,"LineWidth",2)
xlim([min(timeVec_tarint) max(timeVec_tarint)])
xline(0,"LineWidth",1,"Label","Warning Signal","LabelVerticalAlignment","bottom")
xline(900,"LineWidth",1,"Label","Target","LabelVerticalAlignment","bottom")
ylim([ymin ymax])
box("off")