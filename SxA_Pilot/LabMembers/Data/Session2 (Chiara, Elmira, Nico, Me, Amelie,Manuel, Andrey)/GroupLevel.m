% Group Level Analysis
clear
clc

totalsubj=6;
remove3=1; % Exclude subj3 for all group averages? (no proper curve fit and generally weird results
remove6=1;

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
xticklabels({'0' '1' '2' '3'})
hold on
end

% Group Average
subplot(1,2,2);
for subjn=1:totalsubj
objPAStemp(subjn,:)=ssresults(subjn).objmeanperPASrating;
end

if remove3
% remove subj 6
objPAStemp(3,:)=[];
end

if remove6
% remove subj 6
objPAStemp(6,:)=[];
end

objmeanPASgroup=mean(objPAStemp,1);
plot(1:4,objmeanPASgroup,'-o','LineWidth',2)
title('Subjective Perception vs. Accuracy Group Level')
xlabel('PAS Level')
ylabel('Correct in %')
xlim([0 5])
ylim([0.4 1.1])
xticks(1:4)
xticklabels({'0' '1' '2' '3'})

%% %%%%%%%%%%%% These are after eliminating extra irregular trials %%%%%%%%%%%%%%
%% Means for each condition (averaged across contrasts)

f3=figure;
subplot(2,2,1)
for subjn=1:totalsubj
plot(1:3,ssresults(subjn).objgrandmeanpercondition,'-o','LineWidth',2)
title('Accuracy per Condition')
ylabel('Correct in %')
xlim([0 4])
ylim([0 1.1])
xticks([1 2 3])
xticklabels({'Rhythm', 'Interval', 'Irregular'})
hold on
end

subplot(2,2,3)
for subjn=1:totalsubj
plot(1:3,ssresults(subjn).subjgrandmeanpercondition,'-o','LineWidth',2)
title('Subjective Perception per Condition')
ylabel('Subjective Perception')
xlim([0 4])
ylim([0 1.1])
xticks([1 2 3])
xticklabels({'Rhythm', 'Interval', 'Irregular'})
hold on
end

% Group Average averaged across all contrasts

for subjn=1:totalsubj %objective
objpercondtemp(subjn,:)=ssresults(subjn).objgrandmeanpercondition;
end

% % remove subj 3
if remove3
 objpercondtemp(3,:)=[];
end
if remove6
 objpercondtemp(6,:)=[];
end

subplot(2,2,2);
objmeanpercondgroup=mean(objpercondtemp,1);
plot(1:3,objmeanpercondgroup,'-o','LineWidth',2)
title('Accuracy per Condition - Group Level')
ylabel('Accuracy')
xlim([0 4])
ylim([0 1.1])
xticks([1 2 3])
xticklabels({'Rhythm', 'Interval', 'Irregular'})

for subjn=1:totalsubj % subjective
subjpercondtemp(subjn,:)=ssresults(subjn).subjgrandmeanpercondition;
end

% % remove subj 3
if remove3
 subjpercondtemp(3,:)=[];
end
if remove6
 subjpercondtemp(6,:)=[];
end

subjmeanpercondgroup=mean(subjpercondtemp,1);


subplot(2,2,4);
plot(1:3,subjmeanpercondgroup,'-o','LineWidth',2)
title('Subjective Perception per Condition - Group Level')
ylabel('Subjective Perception')
xlim([0 4])
ylim([0 1.1])
xticks([1 2 3])
xticklabels({'Rhythm', 'Interval', 'Irregular'})

%% Mean for each contrast for each condition

for subjn=1:totalsubj
objmeanpercont_rhy(subjn,:)=ssresults(subjn).objmeancondcontrast(1,:);
objmeanpercont_int(subjn,:)=ssresults(subjn).objmeancondcontrast(2,:);
objmeanpercont_irr(subjn,:)=ssresults(subjn).objmeancondcontrast(3,:);

subjmeanpercont_rhy(subjn,:)=ssresults(subjn).subjmeancondcontrast(1,:);
subjmeanpercont_int(subjn,:)=ssresults(subjn).subjmeancondcontrast(2,:);
subjmeanpercont_irr(subjn,:)=ssresults(subjn).subjmeancondcontrast(3,:);
end

% Averaged across subjects
group_objmeanpercont_rhy=mean(objmeanpercont_rhy,1);
group_objmeanpercont_int=mean(objmeanpercont_int,1);
group_objmeanpercont_irr=mean(objmeanpercont_irr,1);

group_subjmeanpercont_rhy=mean(subjmeanpercont_rhy,1);
group_subjmeanpercont_int=mean(subjmeanpercont_int,1);
group_subjmeanpercont_irr=mean(subjmeanpercont_irr,1);

f4=figure;
subplot(1,2,1)
plot(1:10,group_objmeanpercont_rhy,'-o','LineWidth',2)
title('Accuracy per Contrast Level - Group')
ylabel('% Correct')
xlim([0 11])
ylim([0 1.1])
xticks([1:10])
xlabel({'Contrasts'})
hold on
plot(1:10,group_objmeanpercont_int,'-o','LineWidth',2)
plot(1:10,group_objmeanpercont_irr,'-o','LineWidth',2)
legend('Rhythm','Interval','Irregular')
hold off

subplot(1,2,2)
plot(1:10,group_subjmeanpercont_rhy,'-o','LineWidth',2)
title('Subjective Perception per Contrast Level - Group')
ylabel('Subjective Perception')
xlim([0 11])
ylim([0 1.1])
xticks([1:10])
xlabel({'Contrasts'})
hold on
plot(1:10,group_subjmeanpercont_int,'-o','LineWidth',2)
plot(1:10,group_subjmeanpercont_irr,'-o','LineWidth',2)
legend('Rhythm','Interval','Irregular')
hold off

%% Model Fit Group Level

% Calculate difference between conditions:
for subjn=1:totalsubj % irregular minus interval/rhythm
    midpoint_objirr(subjn)=ssresults(subjn).modelfit.objirr.b;
    midpoint_objint(subjn)=ssresults(subjn).modelfit.subjint.b;
    midpoint_objrhy(subjn)=ssresults(subjn).modelfit.objrhy.b;

    midpoint_subjirr(subjn)=ssresults(subjn).modelfit.subjirr.b;
    midpoint_subjint(subjn)=ssresults(subjn).modelfit.subjint.b;
    midpoint_subjrhy(subjn)=ssresults(subjn).modelfit.subjrhy.b;

    midpointshift_ss_rhyobj(subjn)=midpoint_objrhy(subjn)-midpoint_objirr(subjn); % rhythm compared to irregular
    midpointshift_ss_intobj(subjn)=midpoint_objint(subjn)-midpoint_objirr(subjn); % interval compared to irregular
    midpointshift_ss_rhysubj(subjn)=midpoint_subjrhy(subjn)-midpoint_subjirr(subjn); % rhythm compared to irregular
    midpointshift_ss_intsubj(subjn)=midpoint_subjint(subjn)-midpoint_subjirr(subjn); % interval compared to irregular
end

if remove3
% remove subj 3
midpointshift_ss_rhyobj(3)=[];
midpointshift_ss_intobj(3)=[];
midpointshift_ss_rhysubj(3)=[];
midpointshift_ss_intsubj(3)=[];
end

if remove63
% remove subj 3
midpointshift_ss_rhyobj(6)=[];
midpointshift_ss_intobj(6)=[];
midpointshift_ss_rhysubj(6)=[];
midpointshift_ss_intsubj(6)=[];
end

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

