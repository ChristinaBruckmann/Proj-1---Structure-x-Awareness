totalsubj=13;
nblocks=10;
%% Plot Performance Across Blocks
for subj=1:totalsubj
% Load
clearvars -except totalsubj subj nblocks
filename=sprintf('ResultsSubjectOfficial%i',subj);
load(filename,'means')

% Plot Means Per Block
idx_rhy=find(means.blockmeansobj(:,1)==1);
idx_int=find(means.blockmeansobj(:,1)==2);
idx_irr=find(means.blockmeansobj(:,1)==3);
p1(subj)=subplot(5,3,subj); plot(1:nblocks,means.blockmeansobj(:,2))
plot(idx_rhy,means.blockmeansobj(idx_rhy,2),'-o', 'Color', 'b');
hold on
plot(idx_int,means.blockmeansobj(idx_int,2),'-o', 'Color', 'r');
plot(idx_irr,means.blockmeansobj(idx_irr,2),'-o', 'Color', 'g');
currtitle=sprintf('Subject %i',subj);
title(currtitle)
xlabel('Block')
ylabel('Performance')
ylim([0.4 1])
if subj==totalsubj
legend('Rhythm','Interval','Irregular')
end
hold off

end
%% Plot Curves
for subj=1:totalsubj
% Load
clearvars -except totalsubj subj
filename=sprintf('ResultsSubjectOfficial%i',subj);
load(filename)
% Plot Curves
plotOptions1.lineColor = [0,0,0];
plotOptions1.dataColor = [0,0,0];
plotOptions1.CIthresh = true;  
plotOptions1.dataSize=10;
plotOptions1.lineWidth = 1.5;
plotOptions1.labelSize = 10;
plotOptions2.lineColor = [1,0,0];
plotOptions2.dataColor = [1,0,0];
plotOptions2.lineWidth = 1.5;
plotOptions2.dataSize=10;
plotOptions2.CIthresh = true;  
plotOptions2.labelSize = 10;
plotOptions3.lineColor = [0,0.7,0.5];
plotOptions3.dataColor = [0,0.7,0.5];
plotOptions3.dataSize=10;
plotOptions3.lineWidth = 1.5;
plotOptions3.CIthresh = true;
plotOptions3.labelSize = 10;

figure(1);
subplot(5,3,subj); [hline]=plotPsych(psignifitsresults.irrobj,plotOptions1);
hold on
subplot(5,3,subj); [hline2]=plotPsych(psignifitsresults.intobj,plotOptions2);
subplot(5,3,subj); [hline3]=plotPsych(psignifitsresults.rhyobj,plotOptions3);
if subj==1
title('Objective')
end
if subj==totalsubj
legend([hline,hline2,hline3],'Irregular','Interval','Rhythm')
end
hold off

plotOptions1.yLabel='Visibility';
plotOptions2.yLabel='Visibility';
plotOptions3.yLabel='Visibility';

figure(2);
subplot(5,3,subj); [hline]=plotPsych(psignifitsresults.irrsubj,plotOptions1);
hold on
subplot(5,3,subj); [hline2]=plotPsych(psignifitsresults.intsubj,plotOptions2);
subplot(5,3,subj); [hline3]=plotPsych(psignifitsresults.rhysubj,plotOptions3);
if subj==1
title('Subjective')
end
if subj==totalsubj
legend([hline,hline2,hline3],'Irregular','Interval','Rhythm')
end
hold off
end

%% Calculate and Plot Midpoints Across Participants
clearvars -except totalsubj subj
for subj=1:totalsubj
% Load
filename=sprintf('ResultsSubjectOfficial%i',subj);
load(filename)
midpoint.irrobj=psignifitsresults.irrobj.Fit(1);
midpoint.intobj=psignifitsresults.intobj.Fit(1);
midpoint.rhyobj=psignifitsresults.rhyobj.Fit(1);

midpoint.irrsubj=psignifitsresults.irrsubj.Fit(1);
midpoint.intsubj=psignifitsresults.intsubj.Fit(1);
midpoint.rhysubj=psignifitsresults.rhysubj.Fit(1);

% How does int and rhythm shift compared to irr
midpoint.objdiff.int=midpoint.intobj-midpoint.irrobj;
midpoint.objdiff.rhy=midpoint.rhyobj-midpoint.irrobj;

midpoint.subjdiff.int=midpoint.intsubj-midpoint.irrsubj;
midpoint.subjdiff.rhy=midpoint.rhysubj-midpoint.irrsubj;

midpointsall.objdiff.int(subj)=midpoint.objdiff.int;
midpointsall.objdiff.rhy(subj)=midpoint.objdiff.rhy;

midpointsall.subjdiff.int(subj)=midpoint.subjdiff.int;
midpointsall.subjdiff.rhy(subj)=midpoint.subjdiff.rhy;

save(filename, 'midpoint', "-append")
end

% Calculate mean
meanmidpointdiff.objint=mean(midpointsall.objdiff.int);
meanmidpointdiff.objrhy=mean(midpointsall.objdiff.rhy);

meanmidpointdiff.subjint=mean(midpointsall.subjdiff.int);
meanmidpointdiff.subjrhy=mean(midpointsall.subjdiff.rhy);
figure;

plot(1, meanmidpointdiff.objint,'o','color','g')
hold on
plot(1,meanmidpointdiff.objrhy,'o','color','r')
plot(2,meanmidpointdiff.subjint,'o','Color','g')
plot(2,meanmidpointdiff.subjrhy,'o','Color','r')
xlim([0 4])
xticks(1:2)