% Group Level Analysis
clear
clc

totalsubj=3;

% Load all data
for subjn=1:totalsubj
    filename1=sprintf('ResultsSubject%i.mat',subjn);
    ssresults(subjn)=load(filename1);
end

%% Plot results for all subjects
% These are across all trials (including later excluded irregular trials)

%% Sanity Check: Plot graded PAS vs. Accuracy
f1=figure;
subplot(1,2,1);
for subjn=1:totalsubj
plot(1:4,ssresults(subjn).objmeanperPASrating,'-o','LineWidth',2)
title('Subjective Perception vs. Accuracy')
xlabel('PAS Level')
ylabel('Correct in %')
xlim([0 5])
ylim([0.4 1.1])
xticks(1:4)
hold on
end

% Group Average
subplot(1,2,2);
for subjn=1:totalsubj
objPAStemp(subjn,:)=ssresults(subjn).objmeanperPASrating;
end
objmeanPASgroup=mean(objPAStemp,1);
plot(1:4,objmeanPASgroup,'-o','LineWidth',2)
title('Subjective Perception vs. Accuracy Group Level')
xlabel('PAS Level')
ylabel('Correct in %')
xlim([0 5])
ylim([0.4 1.1])
xticks(1:4)
%% d' and criterion



%% These are after eliminating extra irregular trials
%% Means for each condition

f3=figure;
subplot(1,2,1)
for subjn=1:totalsubj
plot(1:3,ssresults(subjn).objgrandmeanpercondition,'-o','LineWidth',2)
title('Accuracy per Condition')
xlabel('Condition')
ylabel('Correct in %')
xlim([0 4])
ylim([0 1.1])
xticks([1 2 3])
hold on
end

subplot(1,2,2)
for subjn=1:totalsubj
plot(1:3,ssresults(subjn).subjgrandmeanpercondition,'-o','LineWidth',2)
title('Subjective Perception per Condition')
xlabel('Condition')
ylabel('Subjective Perception')
xlim([0 4])
ylim([0 1.1])
xticks([1 2 3])
hold on
end