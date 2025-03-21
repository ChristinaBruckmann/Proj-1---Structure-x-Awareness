% Group Level Analysis for Alpha Log Reg
clear 
clc

% Load data 
cd 'Y:\el-Christina\SxA\SxA_Results\LogRegResults\Raw'
load('EEG_SxA_LogRes_GL_Raw.mat')

% Average across participants
obj_mean=squeeze(mean(GL_res_obj,1)); % Input: subj x cond x elec x freq x time point x params, output averaged across participants
subj_mean=squeeze(mean(GL_res_subj,1));

% Average across electrodes
obj_mean=squeeze(mean(obj_mean,2)); % Input: cond x elec x freq x time point x params, output averaged across elecs (cond x freqs x tp x params)
subj_mean=squeeze(mean(subj_mean,2));


%% Plot Slope Coefficients for Contrast and Power
% Parameters: intercept, contrast, power, interaction (?)
scale_colouraxes=1; % 0 - independent axes, 1 - scaled within row, 2 - scaled across all subplots of figure, 3 - scaled across all figures and subplots
zscoring=0; % z-score coefficients (0 - no, 1 - within each subplot, 2- across each row, 3- across all subplots of one figure, 4 - across all plots) - since this changes the main variable, when changing this, run script from top!

% Zscore - 3 and 4 dont seem to work, no idea why
switch zscoring
    case 1
        obj_mean=zscore(obj_mean,0,[2 3]); % zscore across frequencies and time points, not conditions or parameters
        subj_mean=zscore(subj_mean,0,[2 3]);
    case 2
        obj_mean=zscore(obj_mean,0,[1 2 3]); % zscore across frequencies, time points and conditions, but not parameters
        subj_mean=zscore(subj_mean,0,[1 2 3]);
    case 3
        obj_mean(:,:,:,2:3)=zscore(obj_mean(:,:,:,2:3),0,'all'); % zscore across all data, but still separate for subj and obj (only 2 and 3 because otherwise it will scale also across the intercept values etc and distort the results)
        subj_mean(:,:,:,2:3)=zscore(subj_mean(:,:,:,2:3),0,'all');
    case 4
        combined(1,:,:,:,:)=obj_mean; % combine to get zscores across both data
        combined(2,:,:,:,:)=subj_mean;
        combined(:,:,:,:,2:3)=zscore(combined(:,:,:,:,2:3),0,'all'); % zscore across all data, but still separate for subj and obj (only 2 and 3 for last dim because otherwise it will scale also across the intercept values etc and distort the results)
        obj_mean=squeeze(combined(1,:,:,:,:)); % separate again for plotting
        subj_mean=squeeze(combined(2,:,:,:,:));
end

% Objective
figure('Position', [244.2000 157 1.5208e+03 759.2000]);
sgtitle('Objective Performance'); % Main title
for c=1:3 % For each condition
    subplot(2,3,c)
    data_to_plot=squeeze(obj_mean(c,:,:,2)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(obj_mean,2), data_to_plot) % Plot
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
            caxis([min(obj_mean(:,:,:,2),[],'all') max(obj_mean(:,:,:,2),[],'all')]) % scale for current parameter
        case 2
            caxis([min(obj_mean(:,:,:,2:3),[],'all') max(obj_mean(:,:,:,2:3),[],'all')]) % scale for all relevant coefficients
        case 3
            globalmin=min([min(obj_mean(:,:,:,2:3),[],'all') min(subj_mean(:,:,:,2:3),[],'all')]);
            globalmax=max([max(obj_mean(:,:,:,2:3),[],'all') max(subj_mean(:,:,:,2:3),[],'all')]);
            caxis([globalmin globalmax]) % Scale across obj and subj
    end
end

for c=1:3
    subplot(2,3,c+3)
    data_to_plot=squeeze(obj_mean(c,:,:,3)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(obj_mean,2), data_to_plot) % Plot Power Coefficients
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
             caxis([min(obj_mean(:,:,:,3),[],'all') max(obj_mean(:,:,:,3),[],'all')]) % scale for current parameter
         case 2
             caxis([min(obj_mean(:,:,:,2:3),[],'all') max(obj_mean(:,:,:,2:3),[],'all')]) % scale for all parameters
         case 3
             globalmin=min([min(obj_mean(:,:,:,2:3),[],'all') min(subj_mean(:,:,:,2:3),[],'all')]);
             globalmax=max([max(obj_mean(:,:,:,2:3),[],'all') max(subj_mean(:,:,:,2:3),[],'all')]);
             caxis([globalmin globalmax]) % Scale across obj and subj
     end
end

annotation('textbox', [0.01, 0.7, 0.1, 0.1], 'String', 'Contrast Coefficients', ... % Add Row Labels
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
annotation('textbox', [0.01, 0.25, 0.1, 0.1], 'String', 'Power Coefficients', ...
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);


% Subjective 
figure('Position', [244.2000 157 1.5208e+03 759.2000]);
sgtitle('Subjective Perception'); % Main title
for c=1:3 % For each condition
    subplot(2,3,c)
    data_to_plot=squeeze(subj_mean(c,:,:,2)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(subj_mean,2), data_to_plot) % Plot 
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
            caxis([min(subj_mean(:,:,:,2),[],'all') max(subj_mean(:,:,:,2),[],'all')]) % scale for current parameter
        case 2
            caxis([min(subj_mean(:,:,:,2:3),[],'all') max(subj_mean(:,:,:,2:3),[],'all')]) % scale for all parameters
        case 3
            globalmin=min([min(obj_mean(:,:,:,2:3),[],'all') min(subj_mean(:,:,:,2:3),[],'all')]);
            globalmax=max([max(obj_mean(:,:,:,2:3),[],'all') max(subj_mean(:,:,:,2:3),[],'all')]);
            caxis([globalmin globalmax]) % Scale across obj and subj
     end

end

for c=1:3
    subplot(2,3,c+3)
    data_to_plot=squeeze(subj_mean(c,:,:,3)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(subj_mean,2), data_to_plot) % Plot Power Coefficients
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
            caxis([min(subj_mean(:,:,:,3),[],'all') max(subj_mean(:,:,:,3),[],'all')]) % scale for current parameter
        case 2
            caxis([min(subj_mean(:,:,:,2:3),[],'all') max(subj_mean(:,:,:,2:3),[],'all')]) % scale for all parameters
        case 3
            globalmin=min([min(obj_mean(:,:,:,2:3),[],'all') min(subj_mean(:,:,:,2:3),[],'all')]);
            globalmax=max([max(obj_mean(:,:,:,2:3),[],'all') max(subj_mean(:,:,:,2:3),[],'all')]);
            caxis([globalmin globalmax]) % Scale across obj and subj
    end
end

annotation('textbox', [0.01, 0.7, 0.1, 0.1], 'String', 'Contrast Coefficients', ... % Add Row Labels
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
annotation('textbox', [0.01, 0.25, 0.1, 0.1], 'String', 'Power Coefficients', ...
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);

%% T-Test Maps (Each Condition Against 0)
% For each time point (not statistically valid of course, but just as a proxy of whether cluster based permutation is even worth doing
scale_colouraxes=1; % 0 - independent axes, 1 - scaled within row, 2 - scaled across all subplots of figure, 3 - scaled across all figures and subplots
% Average Across Electrodes
T_data_obj=squeeze(mean(GL_res_obj,3));
T_data_subj=squeeze(mean(GL_res_subj,3));

% T-Test Calculations

mean_obj=squeeze(mean(T_data_obj,1)); % Mean across subjects for each data point
mean_subj=squeeze(mean(T_data_subj,1));
std_obj=squeeze(std(T_data_obj,[],1)); % STD across subjects for each data point
std_subj=squeeze(std(T_data_subj,[],1));

t_obj=mean_obj./(std_obj-size(GL_res_obj,1)); % get t values 
t_subj=mean_subj./(std_subj-size(GL_res_subj,1));

p_obj=tcdf(t_obj,(size(GL_res_obj,1)-1));
p_subj=tcdf(t_subj,(size(GL_res_subj,1)-1));

% Plot
% Objective
figure('Position', [244.2000 157 1.5208e+03 759.2000]);
sgtitle('Objective Performance p-values'); % Main title
for c=1:3 % For each condition
    subplot(2,3,c)
    data_to_plot=squeeze(p_obj(c,:,:,2)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(p_obj,2), data_to_plot) % Plot
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
            caxis([min(p_obj(:,:,:,2),[],'all') max(p_obj(:,:,:,2),[],'all')]) % scale for current parameter
        case 2
            caxis([min(p_obj(:,:,:,2:3),[],'all') max(p_obj(:,:,:,2:3),[],'all')]) % scale for all relevant coefficients
        case 3
            globalmin=min([min(p_obj(:,:,:,2:3),[],'all') min(p_obj(:,:,:,2:3),[],'all')]);
            globalmax=max([max(p_obj(:,:,:,2:3),[],'all') max(p_obj(:,:,:,2:3),[],'all')]);
            caxis([globalmin globalmax]) % Scale across obj and subj
    end
end

for c=1:3
    subplot(2,3,c+3)
    data_to_plot=squeeze(p_obj(c,:,:,3)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(p_obj,2), data_to_plot) % Plot Power Coefficients
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
             caxis([min(p_obj(:,:,:,3),[],'all') max(p_obj(:,:,:,3),[],'all')]) % scale for current parameter
         case 2
             caxis([min(p_obj(:,:,:,2:3),[],'all') max(p_obj(:,:,:,2:3),[],'all')]) % scale for all parameters
         case 3
             globalmin=min([min(p_obj(:,:,:,2:3),[],'all') min(p_subj(:,:,:,2:3),[],'all')]);
             globalmax=max([max(p_obj(:,:,:,2:3),[],'all') max(p_subj(:,:,:,2:3),[],'all')]);
             caxis([globalmin globalmax]) % Scale across obj and subj
     end
end

annotation('textbox', [0.01, 0.7, 0.1, 0.1], 'String', 'Contrast Coefficients', ... % Add Row Labels
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
annotation('textbox', [0.01, 0.25, 0.1, 0.1], 'String', 'Power Coefficients', ...
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);


% Subjective 
figure('Position', [244.2000 157 1.5208e+03 759.2000]);
sgtitle('Subjective Perception p-values'); % Main title
for c=1:3 % For each condition
    subplot(2,3,c)
    data_to_plot=squeeze(p_subj(c,:,:,2)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(p_subj,2), data_to_plot) % Plot 
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
            caxis([min(p_subj(:,:,:,2),[],'all') max(p_subj(:,:,:,2),[],'all')]) % scale for current parameter
        case 2
            caxis([min(p_subj(:,:,:,2:3),[],'all') max(p_subj(:,:,:,2:3),[],'all')]) % scale for all parameters
        case 3
            globalmin=min([min(obj_mean(:,:,:,2:3),[],'all') min(p_subj(:,:,:,2:3),[],'all')]);
            globalmax=max([max(obj_mean(:,:,:,2:3),[],'all') max(p_subj(:,:,:,2:3),[],'all')]);
            caxis([globalmin globalmax]) % Scale across obj and subj
     end

end

for c=1:3
    subplot(2,3,c+3)
    data_to_plot=squeeze(p_subj(c,:,:,3)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(p_subj,2), data_to_plot) % Plot Power Coefficients
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
            caxis([min(p_subj(:,:,:,3),[],'all') max(p_subj(:,:,:,3),[],'all')]) % scale for current parameter
        case 2
            caxis([min(p_subj(:,:,:,2:3),[],'all') max(p_subj(:,:,:,2:3),[],'all')]) % scale for all parameters
        case 3
            globalmin=min([min(p_obj(:,:,:,2:3),[],'all') min(p_subj(:,:,:,2:3),[],'all')]);
            globalmax=max([max(p_obj(:,:,:,2:3),[],'all') max(p_subj(:,:,:,2:3),[],'all')]);
            caxis([globalmin globalmax]) % Scale across obj and subj
    end
end

annotation('textbox', [0.01, 0.7, 0.1, 0.1], 'String', 'Contrast Coefficients', ... % Add Row Labels
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
annotation('textbox', [0.01, 0.25, 0.1, 0.1], 'String', 'Power Coefficients', ...
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);

