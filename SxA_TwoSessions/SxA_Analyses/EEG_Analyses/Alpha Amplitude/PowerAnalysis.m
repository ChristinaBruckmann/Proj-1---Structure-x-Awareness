% Power Calculation Alpha
clear
clc
subj=[101   102   103   105   106   107   108   110   111 112   113   114  116  117   118   119  121   122  123 124   126   127   129   130   131   132];
roi=[850 900]; % WS at 0, target at 900

% Load Data
for s=1:length(subj)
    % Load Data (Target Interval)
    cd 'Y:\el-Christina\SxA\SxA_Results\AlphaPowerRes'
    load(sprintf("EEG_SxA_Subj%i_AlphaResults_clean.mat",subj(s)))

    % Save in GL Variable
    for c=1:3 % for each condition
    GL_alpha(s,c,:,:,:)=alpha_Results{1,c};
    timevec=alpha_timeVecTotal{1,1}; % same across all participants
    end
end

% Get roi data & average across electrodes and frequencies
GL_alpha=squeeze(mean(GL_alpha,5));
GL_alpha=squeeze(mean(GL_alpha,4));
GL_alpha= GL_alpha(:,:,timevec>roi(1)&timevec<roi(2)); % select time window
GL_alpha=squeeze(mean(GL_alpha,3)); % average across time window

%% Cohens d (av) for target -500ms
GL_mean=mean(GL_alpha,1); % mean across subjects (for each condition in dim 2)
GL_std=std(GL_alpha,1); % (std across subjects for each condition in dim 2)
fprintf("\n<strong>Summary Statistics: </strong>")
fprintf("\nRhythm Mean %.3f (std: %.3f)",GL_mean(1),GL_std(1));
fprintf("\nInterval Mean %.3f (std: %.3f)",GL_mean(2),GL_std(2));
fprintf("\nIrregular Mean %.3f (std: %.3f)\n",GL_mean(3),GL_std(3));

% Calculate paired-samples cohens d with averaged STDs
pooled_std_rhyint=sqrt((GL_std(1)^2+GL_std(2)^2))/2;
pooled_std_rhyirr=sqrt((GL_std(1)^2+GL_std(3)^2))/2;
pooled_std_intirr=sqrt((GL_std(2)^2+GL_std(3)^2))/2;

cd_rhyint=abs(((GL_mean(1)-GL_mean(2))/pooled_std_rhyint));
cd_rhyirr=abs(((GL_mean(1)-GL_mean(3))/pooled_std_rhyirr));
cd_intirr=abs(((GL_mean(2)-GL_mean(3))/pooled_std_intirr));

fprintf("\n<strong>Cohen's d (av): </strong>")
fprintf("\nRhythm - Interval %.3f ",cd_rhyint);
fprintf("\nRhythm - Irregular %.3f",cd_rhyirr);
fprintf("\nInterval - Irregular %.3f\n",cd_intirr);
%% Calculate sample size needed for 90% power
required_SS_rhyint=round(((1.96 + 1.282)^2/cd_rhyint)+(1.96^2)/2);
required_SS_rhyirr=round(((1.96 + 1.282)^2/cd_rhyirr)+(1.96^2)/2);
required_SS_intirr=round(((1.96 + 1.282)^2/cd_intirr)+(1.96^2)/2);

fprintf("\n<strong>Required Sample Size (90 Percent Power): </strong>")
fprintf("\nRhythm - Interval %i ",required_SS_rhyint);
fprintf("\nRhythm - Irregular %i",required_SS_rhyirr);
fprintf("\nInterval - Irregular %i\n",required_SS_intirr);