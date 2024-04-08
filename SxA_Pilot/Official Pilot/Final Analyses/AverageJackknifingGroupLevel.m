% Calculate Group Level Result based on jackknifing values
% Averages across all jackknifing results from all participants and plots
% curves
%%
clear
clc

subj=1:13;

load('gl_dataforpsignifit_jk.mat')

% Merge Data for Psignifit Objective per condition
rhyobj_all=zeros(size(gl_dataforpsignifit_jk{1, 1}.rhyobj));
intobj_all=rhyobj_all;
irrobj_all=rhyobj_all;

for s=subj
jackgroup.rhyobj_all=rhyobj_all+ gl_dataforpsignifit_jk{1, s}.rhyobj;
jackgroup.intobj_all=rhyobj_all+ gl_dataforpsignifit_jk{1, s}.intobj;  
jackgroup.irrobj_all=rhyobj_all+ gl_dataforpsignifit_jk{1, s}.irrobj;  
end

save('average_jackknife_grouplevel_data','jackgroup')

clear
clc

%% Fit and Plot Curve

load('average_jackknife_grouplevel_data')

% Prepare struct with fitting options
fitting_options_obj=struct;

fitting_options_obj.expType = 'nAFC';
fitting_options_obj.expN = 2;

fitting_options_obj.sigmoidName = 'logistic';

% Fit Objective Curves

grouplevel_res.avgcurves.irrobj = psignifit(jackgroup.irrobj_all,fitting_options_obj);
grouplevel_res.avgcurves.intobj = psignifit(jackgroup.intobj_all,fitting_options_obj);
grouplevel_res.avgcurves.rhyobj = psignifit(jackgroup.rhyobj_all,fitting_options_obj);

% Plot
plotOptions1.lineColor = [0,0,0];
plotOptions1.dataColor = [0,0,0];
plotOptions1.CIthresh = true;
plotOptions1.dataSize=2;
plotOptions1.lineWidth = 1.5;
plotOptions2.lineColor = [1,0,0];
plotOptions2.dataColor = [1,0,0];
plotOptions2.lineWidth = 1.5;
plotOptions2.dataSize=2;
plotOptions2.CIthresh = true;
plotOptions3.lineColor = [0,0.7,0.5];
plotOptions3.dataColor = [0,0.7,0.5];
plotOptions3.dataSize=2;
plotOptions3.lineWidth = 1.5;
plotOptions3.CIthresh = true;

figure(1);
[hline]=plotPsych(grouplevel_res.avgcurves.irrobj,plotOptions1);
hold on
[hline2]=plotPsych(grouplevel_res.avgcurves.intobj,plotOptions2);
[hline3]=plotPsych(grouplevel_res.avgcurves.rhyobj,plotOptions3);
title('Objective Average Curves')
legend([hline,hline2,hline3],'Irregular','Interval','Rhythm')
hold off

save('grouplevel_res','grouplevel_res')