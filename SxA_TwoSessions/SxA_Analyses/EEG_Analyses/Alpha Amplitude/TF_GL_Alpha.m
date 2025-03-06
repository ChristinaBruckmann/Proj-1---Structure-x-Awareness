%% GL Alpha Analysis with (correct) Error Bands
% SxA Group Level Alpha Analysis
% Plots also variance between subjects
% "Raw" alpha plots as well as differences between conditions
% Requires data to already be free of artifacts.

clear
clc

subj=[101:103 105:108 110:114 116:119 121:124 126 127 129:132];
alpharange=8:12;
basec=1;
catchonly=0;
var=1; % variance plot?
perm_cluster=1; % Cluster Based Permutation?


% Time Window for Statistics
stats_tw=[800 900]; %from WS
shade=1; % Shade stats window?

%% Load Data and Average
for s=1:length(subj)
    cd 'Y:\el-Christina\SxA\SxA_Results\AlphaPowerRes'
    if catchonly
        loadfilename=sprintf('EEG_SxA_Subj%i_AlphaResults_Catch.mat',subj(s));
    else
        loadfilename=sprintf('EEG_SxA_Subj%i_AlphaResults_clean.mat',subj(s));
    end
    gl_tf_timeVec(s)=load(loadfilename,'alpha_timeVecTotal');
    gl_alpha_wavFreqs(s)=load(loadfilename,'alpha_wavFreqs');
    subj_res=load(loadfilename,'alpha_Results'); 
    ntrials(s)=load(loadfilename,'alpha_ntrials');

    % array instead of struct
    for c=1:3
    tf_results(s,c,:,:,:)=cell2mat(subj_res.alpha_Results(c)); % subject x condition x time point x freq x electrodes
    end
end

% Average
alpha_results(:,:,:,:)=squeeze(mean(tf_results(:,:,:,alpharange,:),4)); % Select alpha-band only and average across alpha-frequencies (output: subj x cond x time point x electrodes)
alpha_results_avg(:,:,:)=squeeze(mean(alpha_results(:,:,:,:),4)); % average across electrodes

% Baseline Correct for each subject (?)
if basec
    if ~catchonly
        baselinetp=[0 100]; % baseline time points in relation to WS (WS at time point 0)
    else
        baselinetp=[-800 -700]; % baseline time points in relation to Target (Target at 0)
    end
    for s=1:size(alpha_results_avg,1)
        for c=1:3
            alpha_results_avg(s,c,:)=baselineCorrectSegmentedData(squeeze(alpha_results_avg(s,c,:)), gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, baselinetp);
        end
    end
end

if perm_cluster
    % Average across subjects
    Data_Rhy=squeeze(mean(alpha_results_avg(:,1,:),1));
    Data_Int=squeeze(mean(alpha_results_avg(:,2,:),1));
    Data_Irr=squeeze(mean(alpha_results_avg(:,3,:),1));

    % Rhy - Irr
    [Rhy_Clusters]=clusterBasedPermTest(Data_Rhy, Data_Irr, 1);

    % Int - Irr
    [Int_Clusters]=clusterBasedPermTest(Data_Int, Data_Irr, 1) ;

    % Rhy - Int
    [RhyInt_Clusters]=clusterBasedPermTest(Data_Rhy, Data_Int, 1);
end

%% Plot
% Error Bars between Subjects for each condition separately
% Plot with error bars for each condition
figure("Position",[362.6000 137.8000 781.6000 550.4000]);
colourvec=[0.00,0.45,0.74; 0.85,0.33,0.10;0.93,0.69,0.13];

if var
    for c=1:3
        data_to_plot=squeeze(alpha_results_avg(:,c,:)); % select condition data
        varplot(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1},data_to_plot','linewidth',2,"color",colourvec(c,:));
        hold on;
    end
else
    data_to_plot=squeeze(mean(alpha_results_avg(:,:,:),1));
    plot(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, data_to_plot)
end

% Make Nice
xline(0,'--','Warning Signal','LabelVerticalAlignment','bottom',"FontSize",12);
xline(900,'--','Predicted Target',"FontSize",12);
xlabel("time (ms)","FontSize",14)
ylabel("Band-Limited Alpha Amplitude (8-12 Hz)","FontSize",14)
box off
ax=gca;
ax.FontSize = 12;
xlim([min(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}), max(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1})])
xticks=[min(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}):200:max(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1})];

if basec
    title("Alpha Suppression - Baseline Corrected","FontSize",16)
else
    title("Alpha Suppression - Uncorrected","FontSize",16)
end

% Add Clusters
if perm_cluster
    for cl=1:size(Rhy_Clusters,1) % Rhythm-Irregular
        startcl=timeVec(Rhy_Clusters(cl,1)); % Get start point of cluster
        endcl= startcl+Rhy_Clusters(cl,2);% Get end point of cluster
        line([startcl endcl], [ylimits(1)+0.07 ylimits(1)+0.07], 'Color',[0 0.4470 0.7410], 'LineWidth', 2)
    end

    for cl=1:size(Int_Clusters,1) % Interval - Irregular
        startcl=timeVec(Int_Clusters(cl,1)); % Get start point of cluster
        endcl= startcl+Int_Clusters(cl,2);% Get end point of cluster
        line([startcl endcl], [ylimits(1)+0.06 ylimits(1)+0.06], 'Color',[0.8500 0.3250 0.0980], 'LineWidth', 2)
    end

    for cl=1:size(RhyInt_Clusters,1) %  Rhythm - Interval
        startcl=timeVec(RhyInt_Clusters(cl,1)); % Get start point of cluster
        endcl= startcl+RhyInt_Clusters(cl,2);% Get end point of cluster
        line([startcl endcl], [ylimits(1)+0.05 ylimits(1)+0.05], 'Color',[0 0 0], 'LineWidth', 2)
    end
end

% ROI Window Shading
if shade
     y_limits = ylim;
     patch([stats_tw(1) stats_tw(2) stats_tw(2) stats_tw(1)], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], 'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
legend(['Rhythm','',"Interval",'','Irregular',''])
 %% Plot without error bars 
% figure;
% data_to_plot=squeeze(mean(alpha_results_avg(:,:,:),1));
% createfigure_alpha_novar(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, data_to_plot)
% %% Difference Bars
% % Calculate differences between conditions for each participant
% 
% for s=1:size(alpha_results_avg,1)
%     % Difference rhythm to irregular
%     diff_rhy_irr(s,:)=alpha_results_avg(s,1,:)-alpha_results_avg(s,3,:); % rhy - irr
%     % Difference interval to irregular
%     diff_int_irr(s,:)=alpha_results_avg(s,2,:)-alpha_results_avg(s,3,:); % int - irr
% end
% 
% % plot with error bars
% figure;
% colourvec=[0.00,0.45,0.74; 0.85,0.33,0.10];
% 
% varplot(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1},diff_rhy_irr','linewidth',2,"color",colourvec(1,:)); % diff rhythm-irr
% hold on;
% varplot(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1},diff_int_irr','linewidth',2,"color",colourvec(2,:)); % diff interval-irr
% 
% xline(0,'--','Warning Signal');
% xline(800,'--','Predicted Target');
% yline(0)
% xlabel("time (ms)")
% ylabel("Alpha Amplitude (8-12 Hz) Difference")
% legend(['Rhythm','',"Interval",''])
% if basec
%     title("Difference in Alpha Suppression - Baseline Corrected")
% else
%     title("Difference in Alpha Suppression - Uncorrected")
% end

%% Statistics in ROI Time Window
timevec=gl_tf_timeVec(1).alpha_timeVecTotal{1, 1};
alpha_tw=alpha_results_avg(:,:,timevec>=stats_tw(1)&timevec<=stats_tw(2)); % select only time window
alpha_tw=mean(alpha_tw,3); % average across time points


% Rhythm better than irregular?
[h_alpha(1),p_alpha(1),~,stats_alpha(1,:)] = ttest(alpha_tw(:,1),alpha_tw(:,3));
% Interval better than irregular?
[h_alpha(2),p_alpha(2),~,stats_alpha(2,:)] = ttest(alpha_tw(:,2),alpha_tw(:,3));
% Difference between rhythm and interval?
[h_alpha(3),p_alpha(3),~,stats_alpha(3,:)] = ttest(alpha_tw(:,1),alpha_tw(:,2));


% %% Figure Functions
% function createfigure_alpha_novar(X1, YMatrix1)
% %CREATEFIGURE(X1, YMatrix1)
% %  X1:  vector of x data
% %  YMATRIX1:  matrix of y data
% 
% %  Auto-generated by MATLAB on 23-May-2024 10:27:59
% 
% % Create figure
% figure1 = figure;
% 
% % Create axes
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');
% 
% % Create multiple lines using matrix input to plot
% plot(X1,YMatrix1,'LineWidth',3,'Parent',axes1);
% 
% % set(plot1(1),'DisplayName','Rhythm');
% % set(plot1(2),'DisplayName','Interval');
% % set(plot1(3),'DisplayName','Irregular');
% 
% % Create xline
% xline(0,'Parent',axes1,'Alpha',1,'LineStyle','--','LabelOrientation','horizontal','LineWidth',1,'FontSize',15,'Label','Warning Signal');
% 
% % Create matlab.graphics.interactor.ListOfPointsHighlight
% %matlab.graphics.interactor.ListOfPointsHighlight('Visible','off','VertexData',zeros(3,0));
% 
% % Create xline
% xline(900,'Parent',axes1,'Alpha',1,'LineStyle','--','LabelOrientation','horizontal','LineWidth',1,'FontSize',15,'Label','Predicted Target');
% 
% % % Create matlab.graphics.interactor.ListOfPointsHighlight
% % matlab.graphics.interactor.ListOfPointsHighlight('Visible','off','VertexData',zeros(3,0));
% 
% % Create title
% title('Alpha Power between WS and Target');
% 
% xlim(axes1,[-400 1500]);
% hold(axes1,'off');
% % Set the remaining axes properties
% set(axes1,'FontSize',15,'XTick',[0 400 800 1200 1600]);
% xlabel("time(ms)")
% ylabel("Band-Limited Amplitude (µs)")
% % Create legend
% legend1 = legend(['Rhythm', 'Interval',"Irregular","",""]);
% set(legend1,'Position',[0.660288949743352 0.582990712020913 0.158190580302694 0.132720591180465]);
% end