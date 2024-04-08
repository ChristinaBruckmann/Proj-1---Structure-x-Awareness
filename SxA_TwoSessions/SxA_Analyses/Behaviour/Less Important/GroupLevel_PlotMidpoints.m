%SxA_GroupLevelBehav_Thresholds

clear
clc

subjects=[14 15 17:22];

% Load Data And Merge

for subjn=1:length(subjects)

loadname=sprintf('SxA_ResultsSubject%i_Total.mat',subjects(subjn));
load(loadname)

midpoint_objrhy(subjn)=psignifitsresults{1, 1}.obj.Fit(1);
midpoint_objint(subjn)=psignifitsresults{1, 2}.obj.Fit(1);
midpoint_objirr(subjn)=psignifitsresults{1, 3}.obj.Fit(1);

midpoint_subjrhy(subjn)=psignifitsresults{1, 1}.subj.Fit(1);
midpoint_subjint(subjn)=psignifitsresults{1, 2}.subj.Fit(1);
midpoint_subjirr(subjn)=psignifitsresults{1, 3}.subj.Fit(1);

% midpointshift_ss_rhyobj(subjn)=midpoint_objrhy(subjn)-midpoint_objirr(subjn); % rhythm compared to irregular
% midpointshift_ss_intobj(subjn)=midpoint_objint(subjn)-midpoint_objirr(subjn); % interval compared to irregular
% midpointshift_ss_rhysubj(subjn)=midpoint_subjrhy(subjn)-midpoint_subjirr(subjn); % rhythm compared to irregular
% midpointshift_ss_intsubj(subjn)=midpoint_subjint(subjn)-midpoint_subjirr(subjn); % interval compared to irregular

end

%% Plot all midpoints per Condition
figure;
plot(midpoint_objrhy,1,'o','MarkerFaceColor',[0.00,0.45,0.74])
hold on
plot(mean(midpoint_objrhy),1,'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerSize',5)
plot(midpoint_objint,1,'o','MarkerFaceColor',[0.85,0.33,0.10])
plot(mean(midpoint_objint),1,'o','MarkerFaceColor',[0.85,0.33,0.10],'MarkerSize',5)
plot(midpoint_objirr,1,'o','MarkerFaceColor',[0.93,0.69,0.13])
plot(mean(midpoint_objirr),1,'o','MarkerFaceColor',[0.93,0.69,0.13],'MarkerSize',5)




%% Model Fit Group Level

% Calculate difference between conditions:
midpointshift_gl_rhyobj=mean(midpointshift_ss_rhyobj); % rhythm compared to irregular
midpointshift_gl_intobj=mean(midpointshift_ss_intobj); % interval compared to irregular
midpointshift_gl_rhysubj=mean(midpointshift_ss_rhysubj); % rhythm compared to irregular
midpointshift_gl_intsubj=mean(midpointshift_ss_intsubj); % interval compared to irregular

% % Plot individual midpoint shift
% midpointshift_ss_rhyobj(subjn)
% midpointshift_ss_intobj(subjn)
% midpointshift_ss_rhysubj(subjn)
% midpointshift_ss_intsubj(subjn)


% Plot difference in midpoints
f5=figure;

plot(1,midpointshift_gl_rhyobj,'o','LineWidth',2,'Color','b')
hold on
plot(1,midpointshift_gl_intobj,'o','LineWidth',2,'Color','r')
plot(2,midpointshift_gl_rhysubj,'o','LineWidth',2,'Color','b')
plot(2,midpointshift_gl_intsubj,'o','LineWidth',2,'Color','r')
title('Difference in Midpoint of Psychometric Curve - Accuracy')
ylabel('Midpoint Difference')
xlim([0 3])
xticks([1:2])
xlabel(['Objective' 'Subjective'])
legend('Rhythm-Irregular','Interval-Irregular')
hold off