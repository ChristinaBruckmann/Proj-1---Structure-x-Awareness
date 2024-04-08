%% Delta ITPC Group Level
clear
clc

catchonly=0;
ITIonly=1;

if catchonly
    subj=17:22;
else
    subj=[14 15 17:22];
end

% Load Data
for s=1:length(subj)
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
    loadfilename=sprintf('EEG_SxA_Subj%i_Results.mat',subj(s));
    for c=1:3
    if catchonly
        ITPCdelta_elec_gl(s,:)=cell2mat(struct2cell(load(loadfilename,'ITPCdelta_elec_catch')));
        load(loadfilename,'ITPCdelta_res_catch');
        ITPCdelta_res_gl(s,c,:,:)=ITPCdelta_res_catch{1,c};
        load(loadfilename,'ITPCdelta_timevec_catch');
        ITPCdelta_timevec_gl(s,c,:)=ITPCdelta_timevec_catch{1,c};
        ITPCdelta_wavFreqs_gl(s)=load(loadfilename,'ITPCdelta_wavFreqs_catch');
        ITPCdelta_nTrialsCond_gl(s,:)=cell2mat(struct2cell(load(loadfilename,'ITPCdelta_nTrialsCond_catch')));
    elseif ITIonly
        ITPCdelta_elec_gl(s,:)=cell2mat(struct2cell(load(loadfilename,'ITPCdelta_elec_ITI')));
        load(loadfilename,'ITPCdelta_res_ITI');
        ITPCdelta_res_gl(s,c,:,:)=ITPCdelta_res_ITI{1,c};
        load(loadfilename,'ITPCdelta_timevec_ITI');
        ITPCdelta_timevec_gl(s,c,:)=ITPCdelta_timevec_ITI{1,c};
        ITPCdelta_wavFreqs_gl(s)=load(loadfilename,'ITPCdelta_wavFreqs_ITI');
        ITPCdelta_nTrialsCond_gl(s,:)=cell2mat(struct2cell(load(loadfilename,'ITPCdelta_nTrialsCond_ITI')));
    else
        ITPCdelta_elec_gl(s,:)=cell2mat(struct2cell(load(loadfilename,'ITPCdelta_elec')));
        load(loadfilename,'ITPCdelta_res');
        ITPCdelta_res_gl(s,c,:,:)=ITPCdelta_res{1,c};
        load(loadfilename,'ITPCdelta_timevec');
        ITPCdelta_timevec_gl(s,c,:)=ITPCdelta_timevec{1,c};
        ITPCdelta_wavFreqs_gl(s)=load(loadfilename,'ITPCdelta_wavFreqs');
        ITPCdelta_nTrialsCond_gl(s,:)=cell2mat(struct2cell(load(loadfilename,'ITPCdelta_nTrialsCond')));
    end
    end
end

% Average and Save
% for c=1:3 % for each condition
%     for s=1:length(subj) % for each subject
%         gl_itpc_results(s,:,:)=ITPCdelta_res_gl(s,c);
%     end
%     gl_itpc_results_means(c,:,:)=mean(gl_itpc_results(:,:,:),1);
% end

for c=1:3
    gl_itpc_results_means(c,:,:)=mean(ITPCdelta_res_gl(:,c,:,:),1);
end

% Investigate amount of trials and chance angles
for c=1:3 
    % Total N of trials
    totaltrials(c)=round(mean(ITPCdelta_nTrialsCond_gl(:,c),1));

    % Change angles
    angles=2*pi*rand(totaltrials(c),1000);
    angle_mean(c)=mean(circ_r(angles));
    angle_std(c)=std(circ_r(angles));
end

%% check ITPC time course and topography
if catchonly
    timeROI = [-100 -50];
    baseRange = [-800 -700];
    warningline=[-800];
    targetline=[0];
elseif ITIonly
    timeROI = [-900 -100];
else
    timeROI = [700 750];
    baseRange = [0 100];
    warningline=[0];
    targetline=[800];
end

occelonly=input('Display occipital electrodes only? (Yes-1)');
if occelonly
    dispElectrodes = [25:30 62:64]; % Occipital only
else
    dispElectrodes = ITPCdelta_elec_gl(1);
end

timeVector=ITPCdelta_timevec_gl(1,1,:);

if ~ITIonly
baselinecorr=input('Baseline correct? (Yes-1)');
% for c=1:3
%     if baselinecorr
%         deltaplotdata(c,:,:)=baselineCorrectSegmentedData(squeeze(gl_itpc_results_means(c,:,dispElectrodes)), timeVector, baseRange);
%         deltatopodata(c,:,:)=baselineCorrectSegmentedData(squeeze(gl_itpc_results_means(c,:,:)), timeVector, baseRange);
%     else
%         deltaplotdata(c,:,:)=squeeze(gl_itpc_results_means(c,:,dispElectrodes));
%         deltatopodata(c,:,:)=squeeze(gl_itpc_results_means(c,:,:));
%     end
% end
else
    baselinecorr=0;
end

% Correct to ITI Mean ITPC
if baselinecorr
    % Load ITI Data
    for s=1:length(subj)
        cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
        loadfilename=sprintf('EEG_SxA_Subj%i_Results.mat',subj(s));
        for c=1:3
            ITPCdelta_elec_gl_ITI(s,:)=cell2mat(struct2cell(load(loadfilename,'ITPCdelta_elec_ITI')));
            load(loadfilename,'ITPCdelta_res_ITI');
            ITPCdelta_res_gl_ITI(s,c,:,:)=ITPCdelta_res_ITI{1,c};
            load(loadfilename,'ITPCdelta_timevec_ITI');
            ITPCdelta_timevec_gl_ITI(s,c,:)=ITPCdelta_timevec_ITI{1,c};
            ITPCdelta_wavFreqs_gl_ITI(s)=load(loadfilename,'ITPCdelta_wavFreqs_ITI');
            ITPCdelta_nTrialsCond_gl_ITI(s,:)=cell2mat(struct2cell(load(loadfilename,'ITPCdelta_nTrialsCond_ITI')));
        end
    end

    % Calculate Mean ITI ITPC Across Time Points
    for c=1:3
    mean_ITI_ITPC(c,:)=mean(squeeze(mean(ITPCdelta_res_gl_ITI(s,c,:,:),1)),1);
    end

    % Subtract from each electrode at each time point of ITPC data
    for c=1:3
    deltaplotdata(c,:,:)=minus(squeeze(gl_itpc_results_means(c,:,dispElectrodes)),mean_ITI_ITPC(c,dispElectrodes));
    deltatopodata(c,:,:)=minus(squeeze(gl_itpc_results_means(c,:,:)),mean_ITI_ITPC(c,:));
    end

else
    for c=1:3
    deltaplotdata(c,:,:)=squeeze(gl_itpc_results_means(c,:,dispElectrodes));
    deltatopodata(c,:,:)=squeeze(gl_itpc_results_means(c,:,:));
    end
end

% mean time course across occipital electrodes
figure;
for c=1:3 % For conditions
    subplot(1,3,c)
    plot(squeeze(timeVector)', mean(deltaplotdata(c,:,:),3));
    ylim([-0.1 0.4])
    if ~ITIonly
        xline(warningline,'LineStyle', '-', 'Color', [0 0 0])
        xline(targetline,'LineStyle', '- -', 'Color', [0 0 0])
    end

    if ~baselinecorr
    yline(angle_mean(c))
    end

    if c==1
        title('ITPC Rhythm')
    elseif c==2
        title('ITPC Interval')
    elseif c==3
        title('ITPC Irregular')
    end
end

%In one graph
figure;
for c=1:3 % For conditions
    plot(squeeze(timeVector)', mean(deltaplotdata(c,:,:),3));
    hold on
    ylim([-0.1 0.4]);
    legend('Rhythm','Interval','Irregular')
end
title('ITPC Delta')

if ~ITIonly
xline(warningline,'LineStyle', '-', 'Color', [0 0 0])
xline(targetline,'LineStyle', '- -', 'Color', [0 0 0])
end

if ~baselinecorr
    yline(angle_mean(1),'b')
    yline(angle_mean(2),'r')
    yline(angle_mean(3),'y')
end



% Topography
figure;
for c=1:3 % For conditions
    % topography
    subplot(1,3,c)
    curr_topo_data=squeeze(mean(deltatopodata(c,timeVector>=timeROI(1) & timeVector<=timeROI(2), :),2));
    topoplot(curr_topo_data, 'head71.locs'); % ,'maplimits',[-0.2 0.2]
    colorbar
    if c==1
        title('Rhythm')
    elseif c==2
        title('Interval')
    elseif c==3
        title('Irregular')
    end
end

if ITIonly
    createfigure(squeeze(timeVector)', mean(deltaplotdata(c,:,:),3))
end

function createfigure(X1, YMatrix1)
%CREATEFIGURE(X1, YMatrix1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 28-Sep-2023 18:42:32

% Create figure
figure;

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',2);
set(plot1(1),'DisplayName','Rhythm',...
    'Color',[0 0.450980392156863 0.741176470588235]);
set(plot1(2),'DisplayName','Interval',...
    'Color',[0.850980392156863 0.329411764705882 0.101960784313725]);
set(plot1(3),'DisplayName','Irregular');

% Create yline
yline(0.0826204512909749,'Color',[0 0.450980392156863 0.741176470588235],...
    'LineStyle','-.',...
    'LineWidth',2,...
    'FontSize',14,...
    'Label',{'Expected ITPC for Random Phase Angles'});

% Create matlab.graphics.interactor.ListOfPointsHighlight
matlab.graphics.interactor.ListOfPointsHighlight('Visible','off',...
    'VertexData',zeros(3,0));

% Create matlab.graphics.interactor.ListOfPointsHighlight
matlab.graphics.interactor.ListOfPointsHighlight('Visible','off',...
    'VertexData',zeros(3,0));

% Create yline
yline(0.0844975036744117,...
    'Color',[0.850980392156863 0.329411764705882 0.101960784313725],...
    'LineStyle','-.',...
    'LineWidth',2);

% Create matlab.graphics.interactor.ListOfPointsHighlight
matlab.graphics.interactor.ListOfPointsHighlight('Visible','off',...
    'VertexData',zeros(3,0));

% Create yline
yline(0.0662905166420247,...
    'Color',[0.929411764705882 0.690196078431373 0.129411764705882],...
    'LineStyle','-.',...
    'LineWidth',2);

% Create matlab.graphics.interactor.ListOfPointsHighlight
matlab.graphics.interactor.ListOfPointsHighlight('Visible','off',...
    'VertexData',zeros(3,0));

% Create ylabel
ylabel('Intertrial Phase Coherence (0.88-1.77 Hz)');

% Create xlabel
xlabel('Time (ms)');

% Create title
title({'Delta Band (0.88-1.77 Hz) Intertrial Phase Coherence at Trial Beginning ','(Occipital)'});

% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[-0.1 0.4]);
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontSize',14,'YTickLabel',...
    {'-0.1','','0','','0.1','','0.2','','0.3','','0.4'});
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'FontSize',13);

end