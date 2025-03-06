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
    timeVec_mask=alpha_timeVecTotal{1,1}; % same across all participants
end

% Interval Data
for s=1:length(subj)
    % Load Data (Mask)
    cd 'Y:\el-Christina\SxA\SxA_Results\AlphaPowerRes'
    load(sprintf("EEG_SxA_Subj%i_AlphaResults_Mask.mat",subj(s)))

    % Save in GL Variable
    GL_alpha_int_mask(s,:,:,:)=alpha_Results{1,2};
    timeVec_mask=alpha_timeVecTotal{1,1}; % same across all participants

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
%% Plot Rhythm Time Course

% Plot Mask Onset
figure; tiledlayout('flow')
nexttile;
plot(timeVec_mask,GL_alpha_rhy_mask,"LineWidth",2)
xlim([min(timeVec_mask) max(timeVec_mask)])
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
plot(timeVec_mask,GL_alpha_int_mask,"LineWidth",2)
xlim([min(timeVec_mask) max(timeVec_mask)])
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