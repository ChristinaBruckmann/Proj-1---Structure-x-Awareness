% Fitting
% Data needs to be arranged in a nx3 matrix (stimulus level | nCorrect |
% ntotal)
% Need to be sure psignifit is 

clear
clc
totalsubn=13;
for subj_n=1:totalsubn
clearvars -except subj_n totalsubjn
%% Load Data
filename=sprintf('ResultsSubjectOfficial%i',subj_n);
load(filename)
%% Prepare struct with fitting options
fitting_options_obj=struct;
fitting_options_subj=struct;

fitting_options_obj.expType = 'nAFC';
fitting_options_obj.expN = 2;
fitting_options_subj.expType = 'YesNo';

fitting_options_obj.sigmoidName = 'logistic';
fitting_options_subj.sigmoidName = 'logistic';

%% Fit Curves for each condition - data including outliers

% Objective

psignifitsresults.irrobj = psignifit(dataforpsignifit.irrobj,fitting_options_obj);
psignifitsresults.intobj = psignifit(dataforpsignifit.intobj,fitting_options_obj);
psignifitsresults.rhyobj = psignifit(dataforpsignifit.rhyobj,fitting_options_obj);

% Subjective
psignifitsresults.irrsubj = psignifit(dataforpsignifit.irrsubj,fitting_options_subj);
psignifitsresults.intsubj = psignifit(dataforpsignifit.intsubj,fitting_options_subj);
psignifitsresults.rhysubj = psignifit(dataforpsignifit.rhysubj,fitting_options_subj);

%% Fit Curves for each condition reduced - Data without outliers

% Objective
psignifitsresults.irrobj_red = psignifit(red_dataforpsignifit.irrobj,fitting_options_obj);
psignifitsresults.intobj_red = psignifit(red_dataforpsignifit.intobj,fitting_options_obj);
psignifitsresults.rhyobj_red = psignifit(red_dataforpsignifit.rhyobj,fitting_options_obj);

% Subjective
psignifitsresults.irrsubj_red = psignifit(red_dataforpsignifit.irrsubj,fitting_options_subj);
psignifitsresults.intsubj_red = psignifit(red_dataforpsignifit.intsubj,fitting_options_subj);
psignifitsresults.rhysubj_red = psignifit(red_dataforpsignifit.rhysubj,fitting_options_subj);

%% Plot
plotOptions1.lineColor = [0,0,0];
plotOptions1.dataColor = [0,0,0];
plotOptions1.CIthresh = true;  
plotOptions1.dataSize=70;
plotOptions1.lineWidth = 1.5;
plotOptions2.lineColor = [1,0,0];
plotOptions2.dataColor = [1,0,0];
plotOptions2.lineWidth = 1.5;
plotOptions2.dataSize=50;
plotOptions2.CIthresh = true;  
plotOptions3.lineColor = [0,0.7,0.5];
plotOptions3.dataColor = [0,0.7,0.5];
plotOptions3.dataSize=30;
plotOptions3.lineWidth = 1.5;
plotOptions3.CIthresh = true;  

figure;
subplot(2,1,1); [hline]=plotPsych(psignifitsresults.irrobj,plotOptions1);
hold on
subplot(2,1,1); [hline2]=plotPsych(psignifitsresults.intobj,plotOptions2);
subplot(2,1,1); [hline3]=plotPsych(psignifitsresults.rhyobj,plotOptions3);
title('Objective All Trials')
legend([hline,hline2,hline3],'Irregular','Interval','Rhythm')
hold off

subplot(2,1,2); [hline]=plotPsych(psignifitsresults.irrsubj,plotOptions1);
hold on
subplot(2,1,2); [hline2]=plotPsych(psignifitsresults.intsubj,plotOptions2);
subplot(2,1,2); [hline3]=plotPsych(psignifitsresults.rhysubj,plotOptions3);
title('Subjective All Trials')
legend([hline,hline2,hline3],'Irregular','Interval','Rhythm')

figure;
subplot(2,1,1); [hline]=plotPsych(psignifitsresults.irrobj_red,plotOptions1);
hold on
subplot(2,1,1); [hline2]=plotPsych(psignifitsresults.intobj_red,plotOptions2);
subplot(2,1,1); [hline3]=plotPsych(psignifitsresults.rhyobj_red,plotOptions3);
title('Objective Without RT Outliers')
legend([hline,hline2,hline3],'Irregular','Interval','Rhythm')

subplot(2,1,2); [hline,hdata]=plotPsych(psignifitsresults.irrsubj_red,plotOptions1);
hold on
subplot(2,1,2); [hline2,hdata2]=plotPsych(psignifitsresults.intsubj_red,plotOptions2);
subplot(2,1,2); [hline3,hdata3]=plotPsych(psignifitsresults.rhysubj_red,plotOptions3);
title('Subjective Without RT Outliers')
legend([hline,hline2,hline3],'Irregular','Interval','Rhythm')


% % Get goodness-of-fit parameter psignifitsresults.irrobj.devianceResiduals
% for x=1:10
% irrobj_pred_all(x)=
% intobj_pred_all(x)=
% rhyobj_pred_all(x)=
% irrsubj_pred_all(x)=
% intsubj_pred_all(x)=
% rhysubj_pred_all(x)=
% 
% irrobj_pred_red(x)=
% intobj_pred_red(x)=
% rhyobj_pred_red(x)=
% irrsubj_pred_red(x)=
% intsubj_pred_red(x)=
% rhysubj_pred_red(x)=
% end
% 
% r_square_irrobj_all=1-
% r_square_intobj_all=1-
% r_square_rhyobj_all=1-
% r_square_irrsubj_all=1-
% r_square_intsubj_all=1-
% r_square_rhysubj_all=1-
% 
% % Get goodness-of-fit parameter
% r_square_irrobj_red=1-((psignifitsresults.irrobj.devianceResiduals^2)/
% r_square_intobj_red=1-
% r_square_rhyobj_red=1-
% r_square_irrsubj_red=1-
% r_square_intsubj_red=1-
% r_square_rhysubj_red=1-
% 
% % Improvement of GOF
% improvement(1)=r_square_irrobj_all

% Extract and compare Midpoints
midpoint.irrobj=psignifitsresults.irrobj.Fit(1);
midpoint.intobj=psignifitsresults.intobj.Fit(1);
midpoint.rhyobj=psignifitsresults.rhyobj.Fit(1);

midpoint.irrsubj=psignifitsresults.irrsubj.Fit(1);
midpoint.intsubj=psignifitsresults.intsubj.Fit(1);
midpoint.rhysubj=psignifitsresults.rhysubj.Fit(1);

% How does int and rhythm shift compared to irr
midpoint.objdiff.int=midpoint.intobj-midpoint.irrobj;
midpoint.objdiff.rhy=midpoint.rhyobj-midpoint.irrobj;

%% Save
save(filename, 'psignifitsresults' ,'midpoint', "-append")
end