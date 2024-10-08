%% Cross_Trial_Logregression_GL
% Aggregates the regression coefficients which are the output from the SS Log-Reg and associated R scripts
% First predictor is the contrast, the second one is the power at the current time point, outcome is subj and obj perception
% Done for each frequency and each time point
% (misspelled 'Rhythm' as 'rhyhm' when saving, keeping the change here)
clear
clc
subj=[101:103, 105:108, 111, 113, 114];
plotsubj=1;
timep=[250 1200]; % (data aligned to warning signal, target at 800)

cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\SubjObj\Routput'
%% Load Data
% Data Structure (Subject, Obj(1)-Sub(2),Rhy(1)-Int(2), ContrastCoeff(1)-PowerCoeff(2))
for s=1:length(subj)

    % Objective Rhythm
    loadfilename1=sprintf('Subj%i_ObjectiveRhyhmAlpha_contrast.txt',subj(s));
    loadfilename2=sprintf('Subj%i_ObjectiveRhyhmAlpha_power.txt',subj(s));
    logregcoeff(s,1,1,1,:,:)=readmatrix(loadfilename1); % read and save contrast coefficients
    logregcoeff(s,1,1,2,:,:)=readmatrix(loadfilename2); % read and save power coefficients

    % Objective Interval
    loadfilename1=sprintf('Subj%i_ObjectiveIntervalAlpha_contrast.txt',subj(s));
    loadfilename2=sprintf('Subj%i_ObjectiveIntervalAlpha_power.txt',subj(s));
    logregcoeff(s,1,2,1,:,:)=readmatrix(loadfilename1); % read and save contrast coefficients
    logregcoeff(s,1,2,2,:,:)=readmatrix(loadfilename2); % read and save power coefficients

    % Subjective Rhythm
    loadfilename1=sprintf('Subj%i_SubjectiveRhyhmAlpha_contrast.txt',subj(s));
    loadfilename2=sprintf('Subj%i_SubjectiveRhyhmAlpha_power.txt',subj(s));
    logregcoeff(s,2,1,1,:,:)=readmatrix(loadfilename1); % read and save contrast coefficients
    logregcoeff(s,2,1,2,:,:)=readmatrix(loadfilename2); % read and save power coefficients

    % Subjective Interval
    loadfilename1=sprintf('Subj%i_SubjectiveIntervalAlpha_contrast.txt',subj(s));
    loadfilename2=sprintf('Subj%i_SubjectiveIntervalAlpha_power.txt',subj(s));
    logregcoeff(s,2,2,1,:,:)=readmatrix(loadfilename1); % read and save contrast coefficients
    logregcoeff(s,2,2,2,:,:)=readmatrix(loadfilename2); % read and save power coefficients

    % Plot Individual Subjects
    if plotsubj
        % Objective
        figure; title(sprintf("Objective Response Subject %i",subj(s)))
        subplot(2,2,1)
        imagesc(timep(1):timep(2),1:40, squeeze(logregcoeff(s,1,1,1,:,:))') % Plot Contrast Coefficients - Rhythm
        axis xy; colorbar; hold on
        line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 2)
        title(sprintf('Objective Contrast Coefficients - Rhythm - Subject %i',subj(s)))

        subplot(2,2,2)
        imagesc(timep(1):timep(2),1:40, squeeze(logregcoeff(s,1,1,2,:,:))') % Plot Power Coefficients - Rhythm
        axis xy; colorbar; hold on
        line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 2)
        title('Objective Power Coefficients - Rhythm')

        subplot(2,2,3)
        imagesc(timep(1):timep(2),1:40, squeeze(logregcoeff(s,1,2,1,:,:))') % Plot Contrast Coefficients - Interval
        axis xy; colorbar; hold on
        line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 2)
        title('Objective Contrast Coefficients - Interval')

        subplot(2,2,4)
        imagesc(timep(1):timep(2),1:40, squeeze(logregcoeff(s,1,2,2,:,:))') % Plot Power Coefficients - Interval
        axis xy; colorbar; hold on
        line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 2)
        title('Objective Power Coefficients - Interval')

        % Subjective
        figure; title(sprintf("Subjective Response Subject %i",subj(s)))
        subplot(2,2,1)
        imagesc(timep(1):timep(2),1:40, squeeze(logregcoeff(s,2,1,1,:,:))') % Plot Contrast Coefficients - Rhythm
        axis xy; colorbar; hold on
        line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 2)
        title(sprintf('Subjective Contrast Coefficients - Rhythm - Subject %i',subj(s)))

        subplot(2,2,2)
        imagesc(timep(1):timep(2),1:40, squeeze(logregcoeff(s,2,1,2,:,:))') % Plot Power Coefficients - Rhythm
        axis xy; colorbar; hold on
        line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 2)
        title('Subjective Power Coefficients - Rhythm')

        subplot(2,2,3)
        imagesc(timep(1):timep(2),1:40, squeeze(logregcoeff(s,2,2,1,:,:))') % Plot Contrast Coefficients - Interval
        axis xy; colorbar; hold on
        line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 2)
        title('Subjective Contrast Coefficients - Interval')

        subplot(2,2,4)
        imagesc(timep(1):timep(2),1:40, squeeze(logregcoeff(s,2,2,2,:,:))') % Plot Power Coefficients - Interval
        axis xy; colorbar; hold on
        line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 2)
        title('Subjective Power Coefficients - Interval')
    end
end

clearvars loadfilename1 loadfilename2 s

%% GL Analysis
% Average across subjects
logregcoeff_mean=mean(logregcoeff,1);

% Average across conditions (for now, as both are predictive)
logregcoeff_mean=mean(logregcoeff_mean,3);

% Plot raw
for objsubj=[1:2]
    figure('Position', [100 100 1200 500]);
    subplot(1,2,1)
%     heatmap(squeeze(logregcoeff_mean(1,objsubj,1,1,:,:))','GridVisible','off','ColorLimits',[-0.5 1.2]) % Plot Contrast Coefficients
 imagesc(timep(1):timep(2),1:40, squeeze(logregcoeff_mean(1,objsubj,1,1,:,:))') % Plot Contrast Coefficients
 axis xy; colorbar; hold on
 line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 4)
    if objsubj==1
        title('Objective Contrast Coefficients - Predictive Conditions')
    else
        title('Subjective Contrast Coefficients - Predictive Conditions')
    end
    subplot(1,2,2)
%     heatmap(squeeze(logregcoeff_mean(1,objsubj,1,2,:,:))','GridVisible','off','ColorLimits',[-0.5 1.2]) % Plot Power Coefficients
    imagesc(timep(1):timep(2),1:40, squeeze(logregcoeff_mean(1,objsubj,1,2,:,:))') % Plot Power Coefficients
    axis xy; colorbar; hold on
 line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 4)
    if objsubj==1
        title('Objective Power Coefficients - Predictive Conditions')
    else
        title('Subjective Power Coefficients - Predictive Conditions')
    end
end


% Plot t values
for objsubj=[1:2]
    figure('Position', [100 100 1200 500]);
    subplot(1,2,1)
    %     heatmap(squeeze(logregcoeff_mean(1,objsubj,1,1,:,:))','GridVisible','off','ColorLimits',[-0.5 1.2]) % Plot Contrast Coefficients
    mean_for_plot=squeeze(logregcoeff_mean(1,objsubj,1,1,:,:))';
    se_for_plot=squeeze(std(logregcoeff(:,objsubj,1,1,:,:),[],1))'/sqrt(size(logregcoeff,1));
    imagesc(timep(1):timep(2),1:40, mean_for_plot./se_for_plot) % Plot Power Coefficients
    axis xy; colorbar; hold on
    line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 4)
    if objsubj==1
        title('Objective Contrast Coefficients - Predictive Conditions')
    else
        title('Subjective Contrast Coefficients - Predictive Conditions')
    end
    subplot(1,2,2)
    %     heatmap(squeeze(logregcoeff_mean(1,objsubj,1,2,:,:))','GridVisible','off','ColorLimits',[-0.5 1.2]) % Plot Power Coefficients
    mean_for_plot=squeeze(logregcoeff_mean(1,objsubj,1,2,:,:))';
    se_for_plot=squeeze(std(logregcoeff(:,objsubj,1,2,:,:),[],1))'/sqrt(size(logregcoeff,1));
    imagesc(timep(1):timep(2),1:40, mean_for_plot./se_for_plot) % Plot Power Coefficients
    axis xy; colorbar; hold on
    line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 4)
    if objsubj==1
        title('Objective Power Coefficients - Predictive Conditions')
    else
        title('Subjective Power Coefficients - Predictive Conditions')
    end
end

clearvars objsubj
%% Cluster Based Permutation Test
logregcoeff_partmean=mean(logregcoeff,1); % Average across participants
% Objective (averaged across predictive conditions)
objectivedata_contrast=squeeze(mean(logregcoeff_partmean(:,1,:,1,:,:),3)); % Average across predictive conditions, create separate matrices for contrast and power coeff
objectivedata_power=squeeze(mean(logregcoeff_partmean(:,1,:,2,:,:),3)); 

% Subjective (averaged across predictive conditions)
subjectivedata_contrast=squeeze(mean(logregcoeff_partmean(:,2,:,1,:,:),3)); 
subjectivedata_power=squeeze(mean(logregcoeff_partmean(:,2,:,2,:,:),3)); 

% This currently tests objective vs. subjective (because I have that data here)
% Other contrasts that might make more sense are predictive vs. non-predictive or correct vs. incorrect, including subjects as random factors
[clustersTrue_con, trueT_P_con, maxSumPermDistribution_con]=clusterBasedPermTest(objectivedata_contrast,subjectivedata_contrast, 1);
[clustersTrue_pow, trueT_P_pow, maxSumPermDistribution_pow]=clusterBasedPermTest(objectivedata_power,subjectivedata_contrast, 1);

% Outputs -
% clusters: a clusterNumber X 4 matrix with true detected clusters (in non-shuffled
% data). columns: 1: first sample of cluster; 2: number of samples in
% cluster; 3: tSum statistic of cluster; 4: p-value of the cluster relative
% to the permutation distribution.
% trueT_P: a 2 X tN matrix with point by point t and p values (rows 1 and 2).
% maxSumPermDistribution: the constructed null distribution of maxTsum values.
