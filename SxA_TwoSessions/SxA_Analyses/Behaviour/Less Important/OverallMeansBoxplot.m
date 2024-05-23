% Group Average Behavioural
clear
clc

cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\Behavioural Preprocessed'
subjects=[101:103 105: 106 108 113 114 117 118 119];

% Load Data And Merge

for s=1:size(subjects,2)
loadname=sprintf('SxA_ResultsSubject%i_Total.mat',subjects(s));
load(loadname,'means')
% Objective
objmeans(s,:)=means.objgrandmeanpercondition;
% Subjective
subjmeans(s,:)=means.subjgrandmeanpercondition;
end

% Plot
figure;
subplot(2,1,1)
b=bar(mean(objmeans,1),'Facecolor', 'flat');
hold on
% plot(1,rhythm_survey,'ob')
% plot(2,interval_survey,'or')
xticklabels({'Rhythm', 'Interval','Irregular'})
ylabel('% Correct')
ylim([0.5 1.05])
title('Objective Performance')
b.CData(1,:) = [0.00,0.45,0.74];
b.CData(2,:) = [0.85,0.33,0.10];
b.CData(3,:) = [0.93,0.69,0.13];
b.EdgeColor = [1 1 1];
b.FaceAlpha = [0.5];

for i=1:size(objmeans,1)
    plot([1,2,3],[objmeans(i,1),objmeans(i,2),objmeans(i,3)],'ro', "LineWidth", 0.7, "Color", 'k')
hold on
end

% Plot Subjective
subplot(2,1,2)
b=bar(mean(subjmeans,1),'Facecolor', 'flat');
hold on
% plot(1,rhythm_survey,'ob')
% plot(2,interval_survey,'or')
xticklabels({'Rhythm', 'Interval','Irregular'})
ylabel('% Reported as Seen')
ylim([0.0 1.05])
title('Subjective Perception')
b.CData(1,:) = [0.00,0.45,0.74];
b.CData(2,:) = [0.85,0.33,0.10];
b.CData(3,:) = [0.93,0.69,0.13];
b.EdgeColor = [1 1 1];
b.FaceAlpha = [0.5];

for i=1:size(objmeans,1)
    plot([1,2,3],[subjmeans(i,1),subjmeans(i,2),subjmeans(i,3)],'ro', "LineWidth", 0.7, "Color", 'k')
hold on
end