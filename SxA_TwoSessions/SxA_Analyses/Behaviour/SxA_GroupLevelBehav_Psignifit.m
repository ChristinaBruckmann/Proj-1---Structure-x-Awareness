%% Group Average Behavioural
clear
clc

subjects=[101:103 105:106 108 110 111 112 113 114 116 117 118 119 121:124 126 127 129 130 131 132 133];
cd 'Y:\el-Christina\SxA\SxA_Data\Behaviour Preprocessed'

% Load Data And Merge
dataforpsignifit_total=zeros(2,3,10,3); %obj/subj, condition, levels, variables

for s=1:length(subjects)
loadname=sprintf('SxA_ResultsSubject%i_Total.mat',subjects(s));
load(loadname,'dataforpsignifit','midpoints')
% Objective
dataforpsignifit_total(1,1,:,:)=squeeze(dataforpsignifit_total(1,1,:,:)) + dataforpsignifit{1, 1}.obj;
dataforpsignifit_total(1,2,:,:)=squeeze(dataforpsignifit_total(1,2,:,:)) + dataforpsignifit{1, 2}.obj;
dataforpsignifit_total(1,3,:,:)=squeeze(dataforpsignifit_total(1,3,:,:)) + dataforpsignifit{1, 3}.obj;
% Subjective
dataforpsignifit_total(2,1,:,:)=squeeze(dataforpsignifit_total(2,1,:,:)) + dataforpsignifit{1, 1}.subj;
dataforpsignifit_total(2,2,:,:)=squeeze(dataforpsignifit_total(2,2,:,:)) + dataforpsignifit{1, 2}.subj;
dataforpsignifit_total(2,3,:,:)=squeeze(dataforpsignifit_total(2,3,:,:)) + dataforpsignifit{1, 3}.subj;

% Also load individual psigifit results to obtain thresholds
% Objective
midpoints_obj(s,1)=midpoints.rhyobj;
midpoints_obj(s,2)=midpoints.intobj;
midpoints_obj(s,3)=midpoints.irrobj;
% Subjective 
midpoints_subj(s,1)=midpoints.rhysubj;
midpoints_subj(s,2)=midpoints.intsubj;
midpoints_subj(s,3)=midpoints.irrsubj;
end

% Fix levels (imported weirdly)
dataforpsignifit_total(1,1,:,1)=[1:10];
dataforpsignifit_total(1,2,:,1)=[1:10];
dataforpsignifit_total(1,3,:,1)=[1:10];

dataforpsignifit_total(2,1,:,1)=[1:10];
dataforpsignifit_total(2,2,:,1)=[1:10];
dataforpsignifit_total(2,3,:,1)=[1:10];

% Run Psignifit 
[psignifitsresults,midpoints]=fitpsychcurves(dataforpsignifit_total);

%% Plot Results
plotpsychometric(psignifitsresults,midpoints)

%Save
% csvwrite('midpoints_obj.csv',midpoints_obj)
% csvwrite('midpoints_subj.csv',midpoints_subj)

 %% Statistical Analysis
% Perform repeated measures ANOVA 

% Specify the factor levels (conditions)
factorNames = {'Condition'};
factorLevels = {'Condition1', 'Condition2', 'Condition3'};% Assuming 1 - rhythm, 2 - interval , 3- irregular
%factorLevels = [1 2 3]; % Assuming 1 - rhythm, 2 - interval , 3- irregular

% Create a table from the data matrix (each column is a variable)
dataTable_obj = array2table(midpoints_obj, 'VariableNames', {'Condition1', 'Condition2', 'Condition3'});
dataTable_subj = array2table(midpoints_subj, 'VariableNames', {'Condition1', 'Condition2', 'Condition3'});

% Fit the repeated measures ANOVA model
rmModel_obj = fitrm(dataTable_obj, 'Condition1-Condition3 ~ 1', 'WithinDesign', factorLevels);
rmModel_subj = fitrm(dataTable_subj, 'Condition1-Condition3 ~ 1', 'WithinDesign', factorLevels);

% Conduct repeated measures ANOVA
ranovaResults_obj = ranova(rmModel_obj);
ranovaResults_subj = ranova(rmModel_subj);

% Display the ANOVA table
disp('OBJECTIVE Repeated Measures ANOVA Results:');
disp(ranovaResults_obj);
disp('SUBJECTIVE Repeated Measures ANOVA Results:');
disp(ranovaResults_subj);


%% T-Tests for Planned Contrasts (one tailed)
% Rhythm better than irregular?
[h_obj(1),p_obj(1),~,stats_obj(1,:)] = ttest(midpoints_obj(:,1),midpoints_obj(:,3),"Tail","left");
[h_subj(1),p_subj(1),~,stats_subj(1,:)] = ttest(midpoints_subj(:,1),midpoints_subj(:,3),"Tail","left");
% Interval better than irregular?
[h_obj(2),p_obj(2),~,stats_obj(2,:)] = ttest(midpoints_obj(:,2),midpoints_obj(:,3),"Tail","left");
[h_subj(2),p_subj(2),~,stats_subj(2,:)] = ttest(midpoints_subj(:,2),midpoints_subj(:,3),"Tail","left");
% Difference between rhythm and interval?
[h_obj(3),p_obj(3),~,stats_obj(3,:)] = ttest(midpoints_obj(:,1),midpoints_obj(:,2),"Tail","left");
[h_subj(3),p_subj(3),~,stats_subj(3,:)] = ttest(midpoints_subj(:,1),midpoints_subj(:,2),"Tail","left");

%% Run Psignifit
function [psignifitsresults,midpoints]=fitpsychcurves(dataforpsignifit)

% Prepare struct with fitting options
fitting_options_obj=struct;
fitting_options_subj=struct;

fitting_options_obj.expType = 'nAFC';
fitting_options_obj.expN = 2;
fitting_options_subj.expType = 'YesNo';

fitting_options_obj.sigmoidName = 'logistic';
fitting_options_subj.sigmoidName = 'logistic';

% % Simple Fit
% fitting_options_obj.fixedPars = [NaN;NaN;0;NaN;0]; % threshold, width, lapse rate, guess rate,eta
% fitting_options_subj.fixedPars = [NaN;NaN;0;NaN;0]; % threshold, width, lapse rate, guess rate,eta

% Fit Curves for each condition - data including outliers
for conditions=1:3

    % Objective
    psignifitsresults{conditions}.obj = psignifit(squeeze(dataforpsignifit(1, conditions,:,:)),fitting_options_obj);

    % Subjective
    psignifitsresults{conditions}.subj = psignifit(squeeze(dataforpsignifit(2, conditions,:,:)),fitting_options_subj);

end

% Extract and compare Midpoints
midpoints.rhyobj=psignifitsresults{1}.obj.Fit(1);
midpoints.intobj=psignifitsresults{2}.obj.Fit(1);
midpoints.irrobj=psignifitsresults{3}.obj.Fit(1);

midpoints.rhysubj=psignifitsresults{1}.subj.Fit(1);
midpoints.intsubj=psignifitsresults{2}.subj.Fit(1);
midpoints.irrsubj=psignifitsresults{3}.subj.Fit(1);

% How does int and rhythm shift compared to irr
midpoints.objdiff.int=midpoints.intobj-midpoints.irrobj;
midpoints.objdiff.rhy=midpoints.rhyobj-midpoints.irrobj;

midpoints.subjdiff.int=midpoints.intsubj-midpoints.irrsubj;
midpoints.subjdiff.rhy=midpoints.rhysubj-midpoints.irrsubj;
end
%% Plot Psychometric
function []=plotpsychometric(psignifitsresults,midpoints)
%% Plot
% Plot Curves
plotOptions1.lineColor = [0.00,0.45,0.74];
plotOptions1.dataColor = [0.00,0.45,0.74];
plotOptions1.CIthresh = true;  
plotOptions1.dataSize=1;
plotOptions1.lineWidth = 3;

plotOptions1.labelSize      = 15;                   % font size labels
plotOptions1.fontSize       = 10;                   % font size numbers
plotOptions1.plotAsymptote  = false;                 % plot Asympotes 
plotOptions2.lineColor = [0.85,0.33,0.10];
plotOptions2.dataColor = [0.85,0.33,0.10];
plotOptions2.lineWidth = 3;
plotOptions2.dataSize=1;
plotOptions2.CIthresh = true;  
plotOptions2.plotAsymptote  = false;                 % plot Asympotes 

plotOptions3.lineColor = [0.93,0.69,0.13];
plotOptions3.dataColor = [0.93,0.69,0.13];
plotOptions3.dataSize=1;
plotOptions3.lineWidth = 3;
plotOptions3.CIthresh = true;  
plotOptions3.plotAsymptote  = false;                 % plot Asympotes 
plotOptions3.xlabel='Stimulus Intensity';
plotOptions3.ylabel='Proporgdfgftion Correct';

fig=figure;
fig.Position=[328.2000,  152.2000,  762.4000,  611.2000];
subplot(2,1,1); [hline]=plotPsych(psignifitsresults{1}.obj,plotOptions1);
hold on
subplot(2,1,1); [hline2]=plotPsych(psignifitsresults{2}.obj,plotOptions2);
subplot(2,1,1); [hline3]=plotPsych(psignifitsresults{3}.obj,plotOptions3);

% Customize
title('Objective Performance',"FontSize",15)
leg1=legend([hline,hline2,hline3],'Rhythm','Interval','Irregular');
leg1.Position=[0.75, 0.6267, 0.1, 0.1];
ylabel('Proportion Correct',"FontSize",14)
xlabel("Stimulus Level","FontSize",14)
% axesHandles = findall(fig, 'Type', 'axes');
% set(axesHandles, 'FontSize', 12);
hold off

plotOptions1.ylabel='Proportion Seen';

subplot(2,1,2); [hline]=plotPsych(psignifitsresults{1}.subj,plotOptions1);
hold on
subplot(2,1,2); [hline2]=plotPsych(psignifitsresults{2}.subj,plotOptions2);
subplot(2,1,2); [hline3]=plotPsych(psignifitsresults{3}.subj,plotOptions3);

% Customize
title('Subjective Perception',"FontSize",15)
%legend([hline,hline2,hline3],'Rhythm','Interval','Irregular',"Font Size",11)
leg2=legend([hline, hline2, hline3], 'Rhythm', 'Interval', 'Irregular', 'FontSize', 11);
leg2.Position=[0.75, 0.16    0.1    0.1];
xlabel("Stimulus Level","FontSize",14)
ylabel("Proportion 'Seen'","FontSize",14)

% Whole Figure Properties
axesHandles = findall(fig, 'Type', 'axes');
set(axesHandles, 'FontSize', 12);
lines = findall(fig, 'Type', 'Line'); % Find all line objects in the figure
% Initialize an array to hold handles of straight lines
straightLines = [];

% Loop through each line to check its data
for i = 1:length(lines)
    xData = get(lines(i), 'XData');
    yData = get(lines(i), 'YData');
    
    % Check for straight lines: 
    % (can either be vertical or horizontal)
    if all(diff(xData) == 0) || all(diff(yData) == 0) % Vertical or Horizontal
        straightLines = [straightLines; lines(i)]; % Store the handle
    end
end

% Display or modify straight lines as needed
set(straightLines, 'LineWidth', 2); % Example: Change width of straight linesset(lines, 'LineWidth', 2); % Set the line width of all lines to 2

% Plot Midpoints
figure;
plot(1,midpoints.objdiff.rhy,'.r','MarkerSize',20)
hold on
plot(1, midpoints.objdiff.int,'.g','MarkerSize',20)
plot(2,midpoints.subjdiff.rhy,'.g','MarkerSize',20)
plot(2,midpoints.subjdiff.int,'.r','MarkerSize',20)
xlim([0 3])
ylim([-3 3])
xticks(1:2)
xticklabels({'Objective','Subjective'})
title('Midpoint Difference')
yline(0,'-','Irregular Midpoint');
legend('Rhythm','Interval')
end

