% Load
load('grouplevel_res.mat')

%% Prepare struct with fitting options
fitting_options_obj=struct;

fitting_options_obj.expType = 'nAFC';
fitting_options_obj.expN = 2;

fitting_options_obj.sigmoidName = 'logistic';

% Fit Objective Curves

psignifitsresults.irrobj = psignifit(grouplevel_res.dataforpsignifit.irrobj,fitting_options_obj);
psignifitsresults.intobj = psignifit(grouplevel_res.dataforpsignifit.intobj,fitting_options_obj);
psignifitsresults.rhyobj = psignifit(grouplevel_res.dataforpsignifit.rhyobj,fitting_options_obj);


% Plot
plotOptions1.lineColor = [1.00,0.75,0.00];
plotOptions1.dataColor = [1.00,0.75,0.00];
plotOptions1.CIthresh = false;
plotOptions1.dataSize=2;
plotOptions1.lineWidth = 2;
plotOptions2.lineColor = [0.85,0.33,0.10];
plotOptions2.dataColor = [0.85,0.33,0.10];
plotOptions2.lineWidth = 2;
plotOptions2.dataSize=2;
plotOptions2.CIthresh = false;
plotOptions3.lineColor = [0.00,0.45,0.74];
plotOptions3.dataColor = [0.00,0.45,0.74];
plotOptions3.dataSize=2;
plotOptions3.lineWidth = 2;
plotOptions3.CIthresh = false;

figure;
[hline]=plotPsych(psignifitsresults.irrobj,plotOptions1);
hold on
[hline2]=plotPsych(psignifitsresults.intobj,plotOptions2);
[hline3]=plotPsych(psignifitsresults.rhyobj,plotOptions3);
if subj==1
title('Objective')
end
figure(3)
legend([hline3, hline2,hline],'Rhythm','Interval','Irregular')

hold off