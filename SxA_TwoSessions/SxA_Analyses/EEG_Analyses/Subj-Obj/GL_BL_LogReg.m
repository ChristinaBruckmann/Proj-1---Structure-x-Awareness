% Group Level Analysis for Alpha Log Reg
% For baseline data (excludes -100 to WS because of uninterpretable results)
clear 
clc

% Load data 
cd 'Y:\el-Christina\SxA\SxA_Results\LogRegResults\Baseline'
load('EEG_SxA_LogRes_GL_Baseline.mat')

% Average across participants
obj_mean=squeeze(mean(GL_res_obj,1)); % Input: subj x cond x elec x freq x time point x params, output averaged across participants
subj_mean=squeeze(mean(GL_res_subj,1));

% Average across electrodes
obj_mean=squeeze(mean(obj_mean,2)); % Input: cond x elec x freq x time point x params, output averaged across elecs (cond x freqs x tp x params)
subj_mean=squeeze(mean(subj_mean,2));


%% Plot Slope Coefficients for Contrast and Power
% Parameters: intercept, contrast, power, interaction (?)
scale_colouraxes=1; % 0 - independent axes, 1 - scaled within row, 2 - scaled across all subplots of figure, 3 - scaled across all figures and subplots
zscoring=0; % z-score coefficients (0 - no, 1 - within each subplot )across time points but never freqs, 2- across each row, 3- across all subplots of one figure, 4 - across all plots) - since this changes the main variable, when changing this, run script from top!
tw=[50 1]; % min  and max of time window in ms (from WS)
plotfreqs=5:30; % which frequencies to plot (lower ones also seem to be affected by the baseline)

% Calculate TW and reduce data size to selected time window
[~,minidx]=min(abs(timeVec-tw(1)));
[~,maxidx]=max(abs(timeVec-tw(2)));
obj_mean=obj_mean(:,:,minidx:maxidx,:);
subj_mean=subj_mean(:,:,minidx:maxidx,:);
redtimeVec=timeVec(minidx:maxidx);

% Select frequencies of interest
freqidx=1:size(obj_mean,2);
freqidx=find(freqidx<min(plotfreqs) | freqidx>max(plotfreqs)); %which frequencies to exclude
obj_mean(:,freqidx,:,:)=NaN;
subj_mean(:,freqidx,:,:)=NaN;

% Zscore
switch zscoring
    case 1
        % All this complication because zscore cannot handle NaNs
        obj_mean_z= obj_mean(:,plotfreqs,:,:); % select only frequencies of interest
        subj_mean_z=subj_mean(:,plotfreqs,:,:);
        obj_mean_z=zscore(obj_mean_z,0,[3]); % zscore time points, not conditions or parameters
        subj_mean_z=zscore(subj_mean_z,0,[3]);
        obj_mean(:,plotfreqs,:,:)=obj_mean_z; % add z scored frequencies of interest to the excluded NaNs
        subj_mean(:,plotfreqs,:,:)=subj_mean_z;
    case 2
        obj_mean_z= obj_mean(:,plotfreqs,:,:); % select only frequencies of interest
        subj_mean_z=subj_mean(:,plotfreqs,:,:);
        obj_mean_z=zscore(obj_mean_z,0,[1 3]); % zscore time points and conditions, not parameters
        subj_mean_z=zscore(subj_mean_z,0,[1 3]);
        obj_mean(:,plotfreqs,:,:)=obj_mean_z; % add z scored frequencies of interest to the excluded NaNs
        subj_mean(:,plotfreqs,:,:)=subj_mean_z;
    case 3
        obj_mean_z= obj_mean(:,plotfreqs,:,:); % select only frequencies of interest
        subj_mean_z=subj_mean(:,plotfreqs,:,:);
        obj_mean_z=zscore(obj_mean_z,0,[1 3]); % zscore across time points, conditions and parameters
        subj_mean_z=zscore(subj_mean_z,0,[1 3 4]);
        obj_mean(:,plotfreqs,:,:)=obj_mean_z; % add z scored frequencies of interest to the excluded NaNs
        subj_mean(:,plotfreqs,:,:)=subj_mean_z;
    case 4
        combined(1,:,:,:,:)=obj_mean; % combine to get zscores across both data
        combined(2,:,:,:,:)=subj_mean;
        combined_z=combined(:,:,plotfreqs,:,:); % select only frequencies of interest
        combined_z(:,:,:,:,2:3)=zscore(combined_z(:,:,:,:,2:3),[],'all'); % zscore across all data (only 2 and 3 for last dim because otherwise it will scale also across the intercept values etc and distort the results)
        combined(:,:,plotfreqs,:,:)=combined_z;
        obj_mean=squeeze(combined(1,:,:,:,:)); % separate again for plotting
        subj_mean=squeeze(combined(2,:,:,:,:));
end

% Objective
figure('Position', [244.2000 157 1.5208e+03 759.2000]);
sgtitle('Objective Performance'); % Main title
for c=1:3 % For each condition
    subplot(2,3,c)
    data_to_plot=squeeze(obj_mean(c,:,:,2)); % Choose condition and parameter to plot
    imagesc(redtimeVec,1:size(obj_mean,2), data_to_plot) % Plot
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
    imagesc(redtimeVec,1:size(obj_mean,2), data_to_plot) % Plot Power Coefficients
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
    imagesc(redtimeVec,1:size(subj_mean,2), data_to_plot) % Plot 
    axis xy; colorbar;hold on
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
    imagesc(redtimeVec,1:size(subj_mean,2), data_to_plot) % Plot Power Coefficients
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

tw=[50 1]; % min  and max of time window in ms (from WS)
plotfreqs=5:30; % which frequencies to plot (lower ones also seem to be affected by the baseline)

% Average Across Electrodes
T_data_obj=squeeze(mean(GL_res_obj,3));
T_data_subj=squeeze(mean(GL_res_subj,3));

% Calculate TW and reduce data size to selected time window
[~,minidx]=min(abs(timeVec-tw(1)));
[~,maxidx]=max(abs(timeVec-tw(2)));
T_data_obj=T_data_obj(:,:,:,minidx:maxidx,:);
T_data_subj=T_data_subj(:,:,:,minidx:maxidx,:);
redtimeVec=timeVec(minidx:maxidx);

% Select frequencies of interest
freqidx=1:size(T_data_obj,2);
freqidx=find(freqidx<min(plotfreqs) | freqidx>max(plotfreqs)); %which frequencies to exclude
T_data_obj(:,:,freqidx,:,:)=NaN;
T_data_subj(:,:,freqidx,:,:)=NaN;

% T-Test Calculations

mean_obj=squeeze(mean(T_data_obj,1,"omitnan")); % Mean across subjects for each data point
mean_subj=squeeze(mean(T_data_subj,1,"omitnan"));
std_obj=squeeze(std(T_data_obj,[],1,"omitnan")); % STD across subjects for each data point
std_subj=squeeze(std(T_data_subj,[],1,"omitnan"));

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
    imagesc(redtimeVec,1:size(p_subj,2), data_to_plot) % Plot 
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
    imagesc(redtimeVec,1:size(p_subj,2), data_to_plot) % Plot Power Coefficients
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

