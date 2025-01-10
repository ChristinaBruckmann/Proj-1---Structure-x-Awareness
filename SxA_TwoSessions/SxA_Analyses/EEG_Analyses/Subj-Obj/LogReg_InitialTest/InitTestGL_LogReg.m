% Group Level Analysis for Alpha Log Reg
clear 
clc
%subj=[102 105:108 110 112:114 116:119 121 122 124 126 127 129 130]; 

% Load data 
cd 'Y:\el-Christina\SxA\SxA_Results\LogRegResults'
load('GL_Res_LogReg.mat') 

obj_mean=squeeze(mean(res_obj,1)); % Input: subj x cond x freq x time point x params, output averaged across participants
subj_mean=squeeze(mean(res_subj,1));

clear obj_GL subj_GL obj_GL_mat subj_GL_mat f s tp

%% Plot Slope Coefficients for Contrast and Power
scale_colouraxes=1; % 0 - independent axes, 1 - scaled within row, 2 - scaled across all subplots of figure, 3 - scaled across all figures and subplots

% Objective
figure('Position', [2007.4000, 723.4000, 1767.2000, 754.4000]);
sgtitle('Objective Performance'); % Main title
for c=1:3 % For each condition
    subplot(2,3,c)
    data_to_plot=squeeze(obj_mean(c,:,:,2)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(obj_mean,2), data_to_plot) % Plot
    axis xy; hold on
    line([900 900], [0 40], 'Color', [0 0 0], 'LineWidth', 4) % Mark Target Onset
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
            caxis([min(obj_mean,'all') max(obj_mean,'all')]) % scale for all parameters
        case 3
            globalmin=min([min(obj_mean,[],'all') min(subj_mean,[],'all')]);
            globalmax=max([max(obj_mean,[],'all') max(subj_mean,[],'all')]);
            caxis([globalmin globalmax]) % Scale across obj and subj
    end
end

for c=1:3
    subplot(2,3,c+3)
    data_to_plot=squeeze(obj_mean(c,:,:,3)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(obj_mean,2), data_to_plot) % Plot Power Coefficients
    axis xy; colorbar; hold on
    line([900 900], [0 40], 'Color', [0 0 0], 'LineWidth', 4)
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
             caxis([min(obj_mean,[],'all') max(obj_mean,[],'all')]) % scale for all parameters
         case 3
             globalmin=min([min(obj_mean,[],'all') min(subj_mean,[],'all')]);
             globalmax=max([max(obj_mean,[],'all') max(subj_mean,[],'all')]);
             caxis([globalmin globalmax]) % Scale across obj and subj
     end
end

annotation('textbox', [0.01, 0.7, 0.1, 0.1], 'String', 'Contrast Coefficients', ... % Add Row Labels
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
annotation('textbox', [0.01, 0.25, 0.1, 0.1], 'String', 'Power Coefficients', ...
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);


% Subjective 
figure('Position', [2007.4000, 723.4000, 1767.2000, 754.4000]);
sgtitle('Subjective Perception'); % Main title
for c=1:3 % For each condition
    subplot(2,3,c)
    data_to_plot=squeeze(subj_mean(c,:,:,2)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(subj_mean,2), data_to_plot) % Plot 
    axis xy; colorbar; hold on
    line([900 900], [0 40], 'Color', [0 0 0], 'LineWidth', 4) % Mark Target Onset
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
            caxis([min(subj_mean,[],'all') max(subj_mean,[],'all')]) % scale for all parameters
        case 3
            globalmin=min([min(obj_mean,[],'all') min(subj_mean,[],'all')]);
            globalmax=max([max(obj_mean,[],'all') max(subj_mean,[],'all')]);
            caxis([globalmin globalmax]) % Scale across obj and subj
     end

end

for c=1:3
    subplot(2,3,c+3)
    data_to_plot=squeeze(subj_mean(c,:,:,3)); % Choose condition and parameter to plot
    imagesc(timeVec,1:size(subj_mean,2), data_to_plot) % Plot Power Coefficients
    axis xy; colorbar; hold on
    line([900 900], [0 40], 'Color', [0 0 0], 'LineWidth', 4)
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
            caxis([min(subj_mean,[],'all') max(subj_mean,[],'all')]) % scale for all parameters
        case 3
            globalmin=min([min(obj_mean,[],'all') min(subj_mean,[],'all')]);
            globalmax=max([max(obj_mean,[],'all') max(subj_mean,[],'all')]);
            caxis([globalmin globalmax]) % Scale across obj and subj
    end
end

annotation('textbox', [0.01, 0.7, 0.1, 0.1], 'String', 'Contrast Coefficients', ... % Add Row Labels
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
annotation('textbox', [0.01, 0.25, 0.1, 0.1], 'String', 'Power Coefficients', ...
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 12);
