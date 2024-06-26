%% Baseline Test Cross_Trial_Logregression_GL
% Similar to the Cross_Trial_Logregression_GL script, but optimized for examining different baselines
% Only Interval Objective to reduce amount of data to process
clear
clc
subj=[101, 102, 103, 105, 111, 114];
plotsubj=1;
timep=[250 1200]; % (data aligned to warning signal, target at 800)

cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\SubjObj\BaselineTesting'
%% Load Data
% Data Structure (Subject, Baseline(1) - Uncorrected(2) ,Divided (1) - Subtracted (2), zscore(1) - raw(2), contrast(1)-power(2))
% Divided/Subtracted not appliccable to uncorrected, just assigned 1

cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\SubjObj\BaselineTesting'
for s=1:length(subj)

    % Baseline - Divided - Zscored
    loadfilename1=sprintf('bc_dv_z_Subj%i_SubjectiveIntervalAlpha_power.txt',subj(s));
    loadfilename2=sprintf('bc_dv_z_Subj%i_SubjectiveIntervalAlpha_contrast.txt',subj(s));
    logregcoeff(s,1,1,1,1,:,:)=readmatrix(loadfilename1); % read and save contrast coefficients
    logregcoeff(s,1,1,1,2,:,:)=readmatrix(loadfilename2); % read and save power coefficients

    % Baseline - Subtracted - Zscored
    loadfilename1=sprintf('bc_sb_z_Subj%i_SubjectiveIntervalAlpha_power.txt',subj(s));
    loadfilename2=sprintf('bc_sb_z_Subj%i_SubjectiveIntervalAlpha_contrast.txt',subj(s));
    logregcoeff(s,1,2,1,1,:,:)=readmatrix(loadfilename1); % read and save contrast coefficients
    logregcoeff(s,1,2,1,2,:,:)=readmatrix(loadfilename2); % read and save power coefficients

    % Baseline - Divided
    loadfilename1=sprintf('bc_dv_Subj%i_SubjectiveIntervalAlpha_power.txt',subj(s));
    loadfilename2=sprintf('bc_dv_Subj%i_SubjectiveIntervalAlpha_contrast.txt',subj(s));
    logregcoeff(s,1,1,2,1,:,:)=readmatrix(loadfilename1); % read and save contrast coefficients
    logregcoeff(s,1,1,2,2,:,:)=readmatrix(loadfilename2); % read and save power coefficients

    %Baseline - Subtracted
    loadfilename1=sprintf('bc_sb_Subj%i_SubjectiveIntervalAlpha_power.txt',subj(s));
    loadfilename2=sprintf('bc_sb_Subj%i_SubjectiveIntervalAlpha_contrast.txt',subj(s));
    logregcoeff(s,1,2,2,1,:,:)=readmatrix(loadfilename1); % read and save contrast coefficients
    logregcoeff(s,1,2,2,2,:,:)=readmatrix(loadfilename2); % read and save power coefficients

    % Uncorrected - Zscored
    loadfilename1=sprintf('uc_dv_z_Subj%i_SubjectiveIntervalAlpha_power.txt',subj(s));
    loadfilename2=sprintf('uc_dv_z_Subj%i_SubjectiveIntervalAlpha_contrast.txt',subj(s));
    logregcoeff(s,2,1,1,1,:,:)=readmatrix(loadfilename1); % read and save contrast coefficients
    logregcoeff(s,2,1,1,2,:,:)=readmatrix(loadfilename2); % read and save power coefficients

    % Uncorrected
    loadfilename1=sprintf('uc_dv_Subj%i_SubjectiveIntervalAlpha_power.txt',subj(s));
    loadfilename2=sprintf('uc_dv_Subj%i_SubjectiveIntervalAlpha_contrast.txt',subj(s));
    logregcoeff(s,2,1,2,1,:,:)=readmatrix(loadfilename1); % read and save contrast coefficients
    logregcoeff(s,2,1,2,2,:,:)=readmatrix(loadfilename2); % read and save power coefficients +

end
%% GL Analysis
% Average across subjects
logregcoeff_mean=mean(logregcoeff,1);


% Plot raw
for weights=[1:2] % once for contrast coefficients, once for power coefficients (two separate figures)
    figure('Position', [100 100 1200 500]);
    counter=1;

    if weights==1 % Title
        title1='Contrast Coefficients ';
    else
        title1='Power Coefficients ';
    end

    for baseline=[1:2]
        if baseline==1 % only look at subtracted/divided for baseline corrected data
            title2='- Baseline Corrected '; %title
            for divided=[1:2] % once for divided once for subtracted baseline

                if divided==1 %title
                    title3='- Division ';
                else
                    title3='- Subtraction ';
                end

                for zscored=[1:2]
                    if zscored==1 %title
                        title4='- zscored';
                    else
                        title4='- raw';
                    end

                    subplot(2,3,counter)
                    imagesc(timep(1):timep(2),2.^[0:1/6:5], squeeze(logregcoeff_mean(1,baseline,divided,zscored,weights,:,:))') % Plot Contrast Coefficients
                    axis xy; colorbar; hold on
                    line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 4)
                    title(strcat(title1,title2,title3,title4))
                    % Choose Title
                    counter=counter+1;
                end
            end
        else
            title2='- Uncorrected';
            title3=' ';
            for zscored=[1:2] % if uncorrected, only loop over z scored vs. raw
                if zscored==1 %title
                    title4='- zscored';
                else
                    title4='- raw';
                end
                subplot(2,3,counter)
                imagesc(timep(1):timep(2),2.^[0:1/6:5], squeeze(logregcoeff_mean(1,baseline,1,zscored,weights,:,:))') % Plot Contrast Coefficients
                axis xy; colorbar; hold on
                line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 4)
                title(strcat(title1,title2,title3,title4))
                counter=counter+1;
            end
        end
    end
end


% % Plot t values
% for baseline=[1:2]
%     figure('Position', [100 100 1200 500]);
%     subplot(1,2,1)
%     %     heatmap(squeeze(logregcoeff_mean(1,objsubj,1,1,:,:))','GridVisible','off','ColorLimits',[-0.5 1.2]) % Plot Contrast Coefficients
%     mean_for_plot=squeeze(logregcoeff_mean(1,baseline,1,1,:,:))';
%     se_for_plot=squeeze(std(logregcoeff(:,baseline,1,1,:,:),[],1))'/sqrt(size(logregcoeff,1));
%     imagesc(timep(1):timep(2),1:40, mean_for_plot./se_for_plot) % Plot Power Coefficients
%     axis xy; colorbar; hold on
%     line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 4)
%     if baseline==1
%         title('Objective Contrast Coefficients - Predictive Conditions')
%     else
%         title('Subjective Contrast Coefficients - Predictive Conditions')
%     end
%     subplot(1,2,2)
%     %     heatmap(squeeze(logregcoeff_mean(1,objsubj,1,2,:,:))','GridVisible','off','ColorLimits',[-0.5 1.2]) % Plot Power Coefficients
%     mean_for_plot=squeeze(logregcoeff_mean(1,baseline,1,2,:,:))';
%     se_for_plot=squeeze(std(logregcoeff(:,baseline,1,2,:,:),[],1))'/sqrt(size(logregcoeff,1));
%     imagesc(timep(1):timep(2),1:40, mean_for_plot./se_for_plot) % Plot Power Coefficients
%     axis xy; colorbar; hold on
%     line([800 800], [0 40], 'Color', [1 1 1], 'LineWidth', 4)
%     if baseline==1
%         title('Objective Power Coefficients - Predictive Conditions')
%     else
%         title('Subjective Power Coefficients - Predictive Conditions')
%     end
% end