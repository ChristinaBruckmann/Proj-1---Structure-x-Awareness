%% SxA TF Group Level Average and Plot

clear
clc
basec=1; % baseline correction?

%% Load Individual Subjects and Average
subj=[101:103 105:108 110 112:114 116:119 121 122 124 126 127 129 130];

% Load Data
for s=1:length(subj)
    cd 'Y:\el-Christina\SxA\SxA_Results\New TF Results'
    loadfilename=sprintf('EEG_SxA_Subj%i_TF_SingleTrials.mat',subj(s));
    load(loadfilename,'TF_trial_timeVec');
    load(loadfilename,'TF_Results_Trial'); % three conditions: time points x frequencies x trials x electrodes

    % Store
    for c=1:3
        TF_temp=squeeze(mean(TF_Results_Trial{1, c},3)); % Average across trials (saves memory)
        tf_results(s,c,:,:,:)=TF_temp; % output: subject x condition x time points x frequencies x electrodes
        gl_tf_timeVec(s,c,:)=TF_trial_timeVec{1, c}; % output: subject x condition x time points
        clear TF_temp
    end
end

% Average
 gl_tf_res_means=squeeze(mean(tf_results,1)); % average across subjects (output: con x tp x freq x elec)

% Save group level results for faster loading
cd 'Y:\el-Christina\SxA\SxA_Results\New TF Results'
save("GL_TF_Res","gl_tf_res_means","gl_tf_timeVec","subj")

%% Plot TF 
clear
clc

% Params
freqs=1:40;
elecs=[25:30 62:64]; % occipital

% Load Data
cd 'Y:\el-Christina\SxA\SxA_Results\New TF Results'
load("GL_TF_Res","gl_tf_res_means","gl_tf_timeVec")
timeVec=squeeze(gl_tf_timeVec(1,1,:))'; % same time vec for all subj and conditions anyway

% Plot
TF_plot=figure;
for c=1:size(gl_tf_res_means,1) % For conditions
    %minmaxscale= [min(mean(datatoplot,3),[],"all") max(mean(datatoplot,3),[],"all")];
    %totalminmax=[-3 3];
    subplot(3,1,c);imagesc(timeVec, freqs, squeeze(mean(gl_tf_res_means(c,:,freqs,elecs),4))', [1 6]); axis xy; hold all
    xline(0,'w--','Warning Signal');
    xline(800,'w--','Predicted Target');
    colorbar();
    xlabel('Time (ms)', 'FontWeight','bold', 'FontSize', 10);
    ylabel('Frequency (Hz)', 'FontWeight','bold', 'FontSize', 10);
    if c==1
        title('Rhythm - Pre-Target Frequency Amplitudes');
    elseif c==2
        title('Interval - Pre-Target Frequency Amplitudes');
    elseif c==3
        title('Irregular - Pre-Target Frequency Amplitudes');
    end
end

%% Plot TF Better
clear
clc

% Params
freqs=1:40;
elecs=[25:30 62:64]; % occipital

% Load Data
cd 'Y:\el-Christina\SxA\SxA_Results\New TF Results'
load("GL_TF_Res","gl_tf_res_means","gl_tf_timeVec")
timeVec=squeeze(gl_tf_timeVec(1,1,:))'; % same time vec for all subj and conditions anyway

% Plot
TF_plot=figure;
for c=1:size(gl_tf_res_means,1) % For conditions
    %minmaxscale= [min(mean(datatoplot,3),[],"all") max(mean(datatoplot,3),[],"all")];
    %totalminmax=[-3 3];
    subplot(1,3,c);imagesc(timeVec, freqs, squeeze(mean(gl_tf_res_means(c,:,freqs,elecs),4))', [1 6]); axis xy; hold all
    xline(0,'w--','Warning Signal');
    xline(800,'w--','Predicted Target');
    colorbar();
    xlabel('Time (ms)', 'FontWeight','bold', 'FontSize', 10);
    ylabel('Frequency (Hz)', 'FontWeight','bold', 'FontSize', 10);
    if c==1
        title('Rhythm - Pre-Target Frequency Amplitudes');
    elseif c==2
        title('Interval - Pre-Target Frequency Amplitudes');
    elseif c==3
        title('Irregular - Pre-Target Frequency Amplitudes');
    end
end
%% Plot Topography

clear
clc

% Params
freqs=8:12;
time_ROI=[750 850]; % from WS
bc=1; % Baseline correct?
bl_win=[0 100]; % Baseline Correct to which time window?

% Load Data
cd 'Y:\el-Christina\SxA\SxA_Results\New TF Results'
load("GL_TF_Res","gl_tf_res_means","gl_tf_timeVec")
timeVec=squeeze(gl_tf_timeVec(1,1,:))'; % same time vec for all subj and conditions anyway
timewindow=timeVec>=time_ROI(1) & timeVec<=time_ROI(2);

% Baseline Correct for each electrode (useful to see suppression)
if bc
baseline=mean(gl_tf_res_means(:,timeVec>=bl_win(1)&timeVec<=bl_win(2),:,:),2);
gl_tf_res_means=gl_tf_res_means-baseline;
end

% Plot
plot_topo_data=squeeze(mean(gl_tf_res_means(:,timewindow,freqs,:),3)); % average across frequencies
plot_topo_data=double(squeeze(mean(plot_topo_data,2)))'; % average across time points in ROI window
minVal=min(plot_topo_data,[],"all");
maxVal=max(plot_topo_data,[],"all");
figure('Position', [113.80, 186.60, 1320.00, 492.00]); 
for c=1:3
    % topography
    subplot(1,4,c)
    topoplot(plot_topo_data(:,c),'head64.locs','electrodes','on','style','map','shading','interp','maplimits',[minVal,maxVal],'whitebk','on');
    colorbar('off'); % Turn off the color bar
    if c==1
        title('Rhythm','FontSize',12)
    elseif c==2
        title('Interval','FontSize',12)
    elseif c==3
        title('Irregular','FontSize',12)
    end
end

% Make Nice
sgtitle("Alpha Band (8-12Hz) Suppression Before Target Onset")
subplot(1,4,4);
c = colorbar; % Display the color bar
caxis([minVal,maxVal]); % Set the same color axis for all topoplots
axis off; % Turn off the axis in the color bar subplot
c.Position = [0.8, 0.35, 0.02, 0.4]; % position nicely
% Get the default ticks and show only every second tick
defaultTicks = c.Ticks; % Get the default tick values
c.Ticks = defaultTicks(1:2:end); % Keep every second tick
c.FontSize = 12; % Change to your desired font size

%exportgraphics(gcf,'AlphaTopo.jpg','ContentType','vector');
%% Plot Mean TR of Predictable
% average_Tf_pred=squeeze(mean(squeeze(mean(gl_tf_res_means(1:2,:,:,:),4)),1));
% figure;
% imagesc(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, gl_alpha_wavFreqs(1).alpha_wavFreqs, average_Tf_pred'); axis xy; hold all
% xline(0,'r--','Warning Signal');
% xline(800,'b--','Predicted Target');
% colorbar();
% xlabel('Time (ms)', 'FontWeight','bold', 'FontSize', 10);
% ylabel('Frequency (Hz)', 'FontWeight','bold', 'FontSize', 10);
% %% Plot each conditions without alpha
% for c=1:3
% TF_noalpha(c,:,:)=squeeze(mean(gl_tf_res_means(c,:,:,:),4));
% end
% TF_noalpha(:,:,7:12)=NaN;
% for c=1:3
% imAlpha=ones(size(squeeze(TF_noalpha(c,:,:))));
% imAlpha(isnan(squeeze(TF_noalpha(c,:,:))))=0;
% imagesc(squeeze(TF_noalpha(c,:,:))','AlphaData',imAlpha');
% set(gca,'color',0*[1 1 1]);
% end
% figure;
% for c=1:height(TF_noalpha) % For conditions
%     %minmaxscale= [min(mean(datatoplot,3),[],"all") max(mean(datatoplot,3),[],"all")];
%     %totalminmax=[-3 3];
%     subplot(3,1,c);imagesc(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, gl_alpha_wavFreqs(1).alpha_wavFreqs, squeeze(TF_noalpha(c,:,:))','AlphaData',imAlpha'); axis xy; hold all
%     xline(0,'r--','Warning Signal');
%     xline(800,'b--','Predicted Target');
%     ax = gca;
%     ax.CLim = [1 3];
%     colorbar();
%     xlabel('Time (ms)', 'FontWeight','bold', 'FontSize', 10);
%     ylabel('Frequency (Hz)', 'FontWeight','bold', 'FontSize', 10);
%     if c==1
%         title('Rhythm - Pre-Target Frequency Amplitudes');
%     elseif c==2
%         title('Interval - Pre-Target Frequency Amplitudes');
%     elseif c==3
%         title('Irregular - Pre-Target Frequency Amplitudes');
%     end
% end
% %% Plot Differences between Conditions
% for c=1:3
% Allconditions_TF(c,:,:)=squeeze(mean(gl_tf_res_means(c,:,:,:),4));
% end
% % Baseline Correct
% if basec
%     if ~catchonly
%         baselinetp=[0 100]; % baseline time points in relation to WS (WS at time point 0)
%     else
%         baselinetp=[-800 -700]; % baseline time points in relation to Target (Target at 0)
%     end
%     for c=1:height(Allconditions_TF)
%         Allconditions_TF(c,:,:)=baselineCorrectSegmentedData(squeeze(Allconditions_TF(c,:,:)), gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, baselinetp);
%     end
% else
% end
% % Rhythm - Irregular
% rhy_min_irr=Allconditions_TF(1,:,:)-Allconditions_TF(3,:,:);
% % Interval - Irregular
% int_min_irr=Allconditions_TF(2,:,:)-Allconditions_TF(3,:,:);
% % Rhythm - Irregular
% rhy_min_int=Allconditions_TF(1,:,:)-Allconditions_TF(2,:,:);
% figure;
% subplot(3,1,1);imagesc(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, gl_alpha_wavFreqs(1).alpha_wavFreqs, squeeze(rhy_min_irr)'); axis xy; hold all
% title('Difference Rhythm - Irregular');
% xline(0,'r--','Warning Signal');