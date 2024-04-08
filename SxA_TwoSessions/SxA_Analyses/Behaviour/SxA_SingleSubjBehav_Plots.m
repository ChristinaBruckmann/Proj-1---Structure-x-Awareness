%% SxA Single Subject Plots

clc
clear
subj_n=input('Subject Number? ');
session_n=input('Session(1/2) or Merged (3)? ');

if ismember(session_n,[1,2])
    loadname=sprintf('SxA_ResultsSubject%i_Session%i',subj_n,session_n);
elseif session_n==3
    loadname=sprintf('SxA_ResultsSubject%i_Total',subj_n);
end 

load(loadname)


% Plot Means
figure;
idx_rhy=find(means.blockmeansobj(:,1)==1);
idx_int=find(means.blockmeansobj(:,1)==2);
idx_irr=find(means.blockmeansobj(:,1)==3);
plot(idx_rhy,means.blockmeansobj(idx_rhy,2),'-o', 'Color', 'b');
hold on
plot(idx_int,means.blockmeansobj(idx_int,2),'-o', 'Color', 'r');
plot(idx_irr,means.blockmeansobj(idx_irr,2),'-o', 'Color', 'g');
legend('Rhythm','Interval','Irregular')
title('Accuracy across blocks')
xlabel('Block')
ylabel('Correct in %')
hold off

% Accuracy per PAS Rating
figure;
plot(0:3, means.accuracyperPAS,'o','LineWidth',2)
title('Subjective Perception vs. Accuracy')
xlabel('PAS Level')
ylabel('Correct in %')
xlim([-1 4])
ylim([0.4 1.1])
xticks(0:3)

% Accuracy per binary PAS Rating
figure;
plot([0 1],means.accuracybinaryPAS,'o','LineWidth',2)
title('Subjective Perception Binary vs. Accuracy')
xlabel('Visibility')
ylabel('Correct in %')
xlim([-1 2])
ylim([0.4 1.1])
xticks([0 1])

% Plot Psychometric Curves
%% Plot Psychometric
% Plot Curves
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
legend('Rhythm','Interval');
