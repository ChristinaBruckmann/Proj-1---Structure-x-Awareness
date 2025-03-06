% Plot Individual Participants Alpha Log Reg
clear 
clc

% Load data 
cd 'Y:\el-Christina\SxA\SxA_Results\LogRegResults\Full Analysis'
load('EEG_SxA_LogRes_GL_FullAnalysis.mat')

% Average across electrodes (output: subj x cond x freq x time point x params)
obj_mean=squeeze(mean(GL_res_obj,3)); % Input: subj x cond x elec x freq x time point x params, output averaged across elecs
subj_mean=squeeze(mean(GL_res_subj,3));

% Remove first subj (for some reason just 0 - check this!)
obj_mean=obj_mean(2:end,:,:,:,:);
subj_mean=subj_mean(2:end,:,:,:,:);

clear GL_res_obj GL_res_subj 

%% Plot Slope Coefficients for Contrast and Power
% Parameters: intercept, contrast, power, interaction (?)
for s=1:size(obj_mean,1)
scale_colouraxes=1; % 0 - independent axes, 1 - scaled within row, 2 - scaled across all subplots of figure, 3 - scaled across all figures and subplots
zscoring=0; % z-score coefficients (0 - no, 1 - within each subplot, 2- across each row, 3- across all subplots of one figure, 4 - across all plots) - since this changes the main variable, when changing this, run script from top!

% Select data for subject
obj_mean_ind=squeeze(obj_mean(s,:,:,:,:));
subj_mean_ind=squeeze(subj_mean(s,:,:,:,:));
% Zscore - 3 and 4 dont seem to work, no idea why
switch zscoring
    case 1
        obj_mean_ind=zscore(obj_mean_ind,[],[2 3]); % zscore across frequencies and time points, not conditions or parameters
        subj_mean_ind=zscore(subj_mean_ind,[],[2 3]);
    case 2
        obj_mean_ind=zscore(obj_mean_ind,[],[1 2 3]); % zscore across frequencies, time points and conditions, but not parameters
        subj_mean_ind=zscore(subj_mean_ind,[],[1 2 3]);
    case 3
        obj_mean_ind(:,:,:,2:3)=zscore(obj_mean_ind(:,:,:,2:3),[],'all'); % zscore across all data, but still separate for subj and obj (only 2 and 3 because otherwise it will scale also across the intercept values etc and distort the results)
        subj_mean_ind(:,:,:,2:3)=zscore(subj_mean_ind(:,:,:,2:3),[],'all');
    case 4
        combined(1,:,:,:,:)=obj_mean_ind; % combine to get zscores across both data
        combined(2,:,:,:,:)=subj_mean_ind;
        combined(:,:,:,:,2:3)=zscore(combined(:,:,:,:,2:3),[],'all'); % zscore across all data, but still separate for subj and obj (only 2 and 3 for last dim because otherwise it will scale also across the intercept values etc and distort the results)
        obj_mean_ind=squeeze(combined(1,:,:,:,:)); % separate again for plotting
        subj_mean_ind=squeeze(combined(2,:,:,:,:));
end

% Objective
figure('Position', [244.2000 157 1.5208e+03 759.2000]);
figtitle=sprintf('Objective Subject %i',s);
sgtitle(figtitle); % Main title
for c=1:3 % For each condition
    subplot(2,3,c)
    data_to_plot=squeeze(obj_mean_ind(c,:,:,2)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(obj_mean_ind,2), data_to_plot) % Plot
    axis xy; hold on
    line([900 900], [0 40], 'Color', [0 0 0], 'LineWidth', 2) % Mark Target Onset
    line([0 0], [0 40], 'Color', [0 0 0], 'LineWidth', 2) % Mark WS Onset
    if c == 1
        title('Rhythm'); % Subtitle
    elseif c == 2
        title('Interval');
    else
        title('Irregular');
    end
    colorbar;
    switch scale_colouraxes
        case 1
            caxis([min(obj_mean_ind(:,:,:,2),[],'all') max(obj_mean_ind(:,:,:,2),[],'all')]) % scale for current parameter
        case 2
            caxis([min(obj_mean_ind(:,:,:,2:3),[],'all') max(obj_mean_ind(:,:,:,2:3),[],'all')]) % scale for all relevant coefficients
        case 3
            globalmin=min([min(obj_mean_ind(:,:,:,2:3),[],'all') min(subj_mean_ind(:,:,:,2:3),[],'all')]);
            globalmax=max([max(obj_mean_ind(:,:,:,2:3),[],'all') max(subj_mean_ind(:,:,:,2:3),[],'all')]);
            caxis([globalmin globalmax]) % Scale across obj and subj
    end
end

for c=1:3
    subplot(2,3,c+3)
    data_to_plot=squeeze(obj_mean_ind(c,:,:,3)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(obj_mean_ind,2), data_to_plot) % Plot Power Coefficients
    axis xy; colorbar; hold on
    line([900 900], [0 40], 'Color', [0 0 0], 'LineWidth', 2)
     line([0 0], [0 40], 'Color', [0 0 0], 'LineWidth', 2) % Mark WS Onset
     if c == 1
        title('Rhythm'); % Subtitle
    elseif c == 2
        title('Interval');
    else
        title('Irregular');
     end
     switch scale_colouraxes % Set colour axes
         case 1
             caxis([min(obj_mean_ind(:,:,:,3),[],'all') max(obj_mean_ind(:,:,:,3),[],'all')]) % scale for current parameter
         case 2
             caxis([min(obj_mean_ind(:,:,:,2:3),[],'all') max(obj_mean_ind(:,:,:,2:3),[],'all')]) % scale for all parameters
         case 3
             globalmin=min([min(obj_mean_ind(:,:,:,2:3),[],'all') min(subj_mean_ind(:,:,:,2:3),[],'all')]);
             globalmax=max([max(obj_mean_ind(:,:,:,2:3),[],'all') max(subj_mean_ind(:,:,:,2:3),[],'all')]);
             caxis([globalmin globalmax]) % Scale across obj and subj
     end
end

annotation('textbox', [0.01, 0.7, 0.1, 0.1], 'String', 'Contrast Coefficients', ... % Add Row Labels
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
annotation('textbox', [0.01, 0.25, 0.1, 0.1], 'String', 'Power Coefficients', ...
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);


% Subjective 
figure('Position', [244.2000 157 1.5208e+03 759.2000]);
figtitle=sprintf('Subjective Subject %i',s);
sgtitle(figtitle); % Main title
for c=1:3 % For each condition
    subplot(2,3,c)
    data_to_plot=squeeze(subj_mean_ind(c,:,:,2)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(subj_mean_ind,2), data_to_plot) % Plot 
    axis xy; colorbar; hold on
    line([900 900], [0 40], 'Color', [0 0 0], 'LineWidth', 2) % Mark Target Onset
     line([0 0], [0 40], 'Color', [0 0 0], 'LineWidth', 2) % Mark WS Onset
    if c == 1
        title('Rhythm'); % Subtitle
    elseif c == 2
        title('Interval');
    else
        title('Irregular');
    end

     switch scale_colouraxes
        case 1
            caxis([min(subj_mean_ind(:,:,:,2),[],'all') max(subj_mean_ind(:,:,:,2),[],'all')]) % scale for current parameter
        case 2
            caxis([min(subj_mean_ind(:,:,:,2:3),[],'all') max(subj_mean_ind(:,:,:,2:3),[],'all')]) % scale for all parameters
        case 3
            globalmin=min([min(obj_mean_ind(:,:,:,2:3),[],'all') min(subj_mean_ind(:,:,:,2:3),[],'all')]);
            globalmax=max([max(obj_mean_ind(:,:,:,2:3),[],'all') max(subj_mean_ind(:,:,:,2:3),[],'all')]);
            caxis([globalmin globalmax]) % Scale across obj and subj
     end

end

for c=1:3
    subplot(2,3,c+3)
    data_to_plot=squeeze(subj_mean_ind(c,:,:,3)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(subj_mean_ind,2), data_to_plot) % Plot Power Coefficients
    axis xy; colorbar; hold on
    line([900 900], [0 40], 'Color', [0 0 0], 'LineWidth', 2)
     line([0 0], [0 40], 'Color', [0 0 0], 'LineWidth', 2) % Mark WS Onset
    if c == 1
        title('Rhythm'); % Subtitle
    elseif c == 2
        title('Interval');
    else
        title('Irregular');
    end

    switch scale_colouraxes
        case 1
            caxis([min(subj_mean_ind(:,:,:,3),[],'all') max(subj_mean_ind(:,:,:,3),[],'all')]) % scale for current parameter
        case 2
            caxis([min(subj_mean_ind(:,:,:,2:3),[],'all') max(subj_mean_ind(:,:,:,2:3),[],'all')]) % scale for all parameters
        case 3
            globalmin=min([min(obj_mean_ind(:,:,:,2:3),[],'all') min(subj_mean_ind(:,:,:,2:3),[],'all')]);
            globalmax=max([max(obj_mean_ind(:,:,:,2:3),[],'all') max(subj_mean_ind(:,:,:,2:3),[],'all')]);
            caxis([globalmin globalmax]) % Scale across obj and subj
    end
end

annotation('textbox', [0.01, 0.7, 0.1, 0.1], 'String', 'Contrast Coefficients', ... % Add Row Labels
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
annotation('textbox', [0.01, 0.25, 0.1, 0.1], 'String', 'Power Coefficients', ...
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
end