%% Group Average Behavioural
clear
clc

subjects=[1:13];
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_Pilot\Official Pilot\Results\Advisory Board Results'

% Load Data And Merge
dataforpsignifit_total=zeros(2,3,10,3); %obj/subj, condition, levels, variables

for s=1:length(subjects)
loadname=sprintf('SxAPilot_ResultsSubject%i.mat',subjects(s));
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

% Plot Results
plotpsychometric(psignifitsresults,midpoints)

 %% Statistical Analysis
% % Perform repeated measures ANOVA 
% 
% % Specify the factor levels (conditions)
% factorNames = {'Condition'};
% factorLevels = {'Condition1', 'Condition2', 'Condition3'};% Assuming 1 - rhythm, 2 - interval , 3- irregular
% %factorLevels = [1 2 3]; % Assuming 1 - rhythm, 2 - interval , 3- irregular
% 
% % Create a table from the data matrix (each column is a variable)
% dataTable_obj = array2table(midpoints_obj, 'VariableNames', {'Condition1', 'Condition2', 'Condition3'});
% dataTable_subj = array2table(midpoints_subj, 'VariableNames', {'Condition1', 'Condition2', 'Condition3'});
% 
% % Fit the repeated measures ANOVA model
% rmModel_obj = fitrm(dataTable_obj, 'Condition1-Condition3 ~ 1', 'WithinDesign', factorLevels);
% rmModel_subj = fitrm(dataTable_subj, 'Condition1-Condition3 ~ 1', 'WithinDesign', factorLevels);
% 
% % Conduct repeated measures ANOVA
% ranovaResults_obj = ranova(rmModel_obj);
% ranovaResults_subj = ranova(rmModel_subj);
% 
% % Display the ANOVA table
% disp('OBJECTIVE Repeated Measures ANOVA Results:');
% disp(ranovaResults_obj);
% disp('SUBJECTIVE Repeated Measures ANOVA Results:');
% disp(ranovaResults_subj);
% 
% % Post-hoc comparisons (if needed) using multcompare function
% % Specify the comparison method (e.g., 'tukey-kramer', 'bonferroni', etc.)
% % Replace 'ComparisonMethod' with your preferred method
% posthocResults_obj = multcompare(rmModel_obj);
% posthocResults_subj = multcompare(rmModel_subj, 'Condition', 'By', 'Condition', 'ComparisonType', 'tukey-kramer');
% 
% 
% % Display the post-hoc comparisons table
% disp('OBJECTIVE Post-hoc Comparisons:');
% disp(posthocResults_obj);
% 
% disp('SUBJECTIVE Post-hoc Comparisons:');
% disp(posthocResults_subj);

%% T-Tests for Planned Contrasts
% Rhythm better than irregular?
[h_obj(1),p_obj(1)] = ttest(midpoints_obj(:,1),midpoints_obj(:,3));
[h_subj(1),p_subj(1)] = ttest(midpoints_subj(:,1),midpoints_subj(:,3));
% Interval better than irregular?
[h_obj(2),p_obj(2)] = ttest(midpoints_obj(:,2),midpoints_obj(:,3));
[h_subj(2),p_subj(2)] = ttest(midpoints_subj(:,2),midpoints_subj(:,3));
% Difference between rhythm and interval?
[h_obj(3),p_obj(3)] = ttest(midpoints_obj(:,1),midpoints_obj(:,2));
[h_subj(3),p_subj(3)] = ttest(midpoints_subj(:,1),midpoints_subj(:,2));
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
plotOptions1.dataSize=2;
plotOptions1.lineWidth = 2;
plotOptions2.lineColor = [0.85,0.33,0.10];
plotOptions2.dataColor = [0.85,0.33,0.10];
plotOptions2.lineWidth = 2;
plotOptions2.dataSize=2;
plotOptions2.CIthresh = true;  
plotOptions3.lineColor = [0.93,0.69,0.13];
plotOptions3.dataColor = [0.93,0.69,0.13];
plotOptions3.dataSize=2;
plotOptions3.lineWidth = 2;
plotOptions3.CIthresh = true;  

figure;
subplot(2,1,1); [hline]=plotPsych(psignifitsresults{1}.obj,plotOptions1);
hold on
subplot(2,1,1); [hline2]=plotPsych(psignifitsresults{2}.obj,plotOptions2);
subplot(2,1,1); [hline3]=plotPsych(psignifitsresults{3}.obj,plotOptions3);
title('Objective All Trials')
legend([hline,hline2,hline3],'Rhythm','Interval','Irregular')
hold off

subplot(2,1,2); [hline]=plotPsych(psignifitsresults{1}.subj,plotOptions1);
hold on
subplot(2,1,2); [hline2]=plotPsych(psignifitsresults{2}.subj,plotOptions2);
subplot(2,1,2); [hline3]=plotPsych(psignifitsresults{3}.subj,plotOptions3);
title('Subjective All Trials')
legend([hline,hline2,hline3],'Rhythm','Interval','Irregular')

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

