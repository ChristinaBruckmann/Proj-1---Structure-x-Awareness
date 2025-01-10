%% Single Subject Analysis (Structure x Awareness)

function SxA_SingleSubjBehav_Analysis(subj_n,session_n,plots)
cd 'Y:\el-Christina\SxA\SxA_Data\Raw\Behaviour Raw\Second Batch (Post-TRF)'

% subj_n=input('Subject Number? ');
% session_n=input('Session(1/2) or Merged (3)? ');
% plots=input('Generate Plots (0/1)? '); 

if ismember(session_n,[1,2])
    loadname=sprintf('SxA1.1_s%i-session%i.mat',subj_n,session_n);
    savename=sprintf('SxA_ResultsSubject%i_Session%i',subj_n,session_n);
elseif session_n==3
    loadname=sprintf('SxA1.1_s%i-merged.mat',subj_n);
    savename=sprintf('SxA_ResultsSubject%i_Total',subj_n);
end 

load(loadname)
%% Clean Up Data File & Basic Checks
cd 'Y:\el-Christina\SxA\SxA_Data\Behaviour Preprocessed'
% Check frame rate:
meanfrrate=mean(subresults.data{:,"Frame Rate"});

if meanfrrate==60
    sprintf('Frame rate check passed. Mean frame rate across all trials is 60hz.')
else
    sprintf('ATTENTION. Frame rate check FAILED. Mean frame rate across all trials is %i.',meanfrrate)
end

% Add variable: Session
session_number(1:height(subresults.data),1)=session_n;
subresults.data=addvars(subresults.data,session_number,'After',"Subject",'NewVariableNames',"Session");

% Change missing variables to NaN
subresults.data{:,"Visibility Response"}(subresults.data{:,"Visibility Response"}==99)=NaN;
subresults.data{:,"Reaction Time"}(subresults.data{:,"Reaction Time"}==9999)=NaN;
subresults.data{:,"Orientation Reseponse"}(subresults.data{:,"Orientation Reseponse"}==99)=NaN;
subresults.data{:,"Correct/Incorrect"}(subresults.data{:,"Correct/Incorrect"}==99)=NaN;

% Add variable: Categorize PAS responses into binary data

set(0, 'DefaultFigurePosition', [430 230 560 420])

BinaryVisibility=zeros(1,height(subresults.data));

for ivis=1:height(subresults.data)
    if subresults.data{ivis,"Visibility Response"}== 0 % invisible
        BinaryVisibility(ivis)=0;
    elseif subresults.data{ivis,"Visibility Response"}==99 || isnan(subresults.data{ivis,"Visibility Response"}) % missing
        BinaryVisibility(ivis)=NaN;
    else % visible
        BinaryVisibility(ivis)=1;
    end
end

subresults.data=addvars(subresults.data,BinaryVisibility','After',"Visibility Response",'NewVariableNames','Binary Visibility');

% Create Total List of Trials In Correct Order and add to data file
countrhy=1;
countirr=1;
countint=1;

for icon=1:height(subresults.data)
    if subresults.data{icon,"Condition"} == 3 % irregular
        Alltrials(icon,:)=subresults.trialmatrix3(countirr,:);
        countirr=countirr+1;
    elseif subresults.data{icon,"Condition"}==1 % rhythm
        Alltrialstemp=subresults.trialmatrix1(countrhy,:);
        Alltrialstemp(:,3)=NaN;
        Alltrials(icon,:)= Alltrialstemp;
        countrhy=countrhy+1;
    elseif subresults.data{icon,"Condition"}==2 %interval 
        Alltrialstemp=subresults.trialmatrix2(countint,:);
        Alltrialstemp(:,3)=NaN;
        Alltrials(icon,:)= Alltrialstemp;
        countint=countint+1;
    end
end

subresults.data=addvars(subresults.data,Alltrials(:,1),'After',"Trial",'NewVariableNames','Contrast Level');
subresults.data=addvars(subresults.data,Alltrials(:,2),'After',"Trial",'NewVariableNames','Orientation');
subresults.data=addvars(subresults.data,Alltrials(:,3),'After',"Trial",'NewVariableNames','Irregular Target Time');
save(savename, 'subresults')


% Remove Additional Irregular Trials before Analysis
% Change (30.10.23): Only for Session 1 and Merged, because for Session 2 we
% need the info for all trials for the EEG.
if ~session_n==2
idxaddirr1=find((subresults.data{:,"Irregular Target Time"} == 1) & (subresults.data{:,"Contrast Level"} ~= 0)); % Do not remove catch trials which were assigned random target times
idxaddirr2=find(subresults.data{:,"Irregular Target Time"} == 5 & subresults.data{:,"Contrast Level"} ~= 0); 

alldataclean=subresults.data;
idxaddirrall=sort([idxaddirr1; idxaddirr2]);
alldataclean(idxaddirrall,:)=[];
else
alldataclean=subresults.data;
end
save(savename, 'alldataclean',"-append")

%% Calculate and save all the means (per condition, per contrast, per rating etc.)
nblocks=max(subresults.data{:,'Block'});
means=calculate_means(alldataclean,nblocks);
save(savename,'means',"-append")

% Plot Means Per Block
if plots
plotmeans(means)
end

%% Fit Psychometric Curves
% Prepare Data for Psignifit
[dataforpsignifit]=prepforpsignifit(alldataclean);
save(savename, 'dataforpsignifit',"-append")

% Run Psignifit
[psignifitsresults,midpoints]=fitpsychcurves(dataforpsignifit,plots);
save(savename, 'psignifitsresults',"midpoints","-append")

% Plot Psignifit Results
if plots
plotpsychometric(psignifitsresults,midpoints)
end
%% Function: Calculate Means
function [means]=calculate_means(alldataclean,nblocks)
%%%%% Means for all conditions

conditions=1:3;

% Logical Indices
indxirr=find(alldataclean{:,"Condition"} == 3);
indxint=find(alldataclean{:,"Condition"} == 2);
indxrhy=find(alldataclean{:,"Condition"} == 1);

nanidx=find(alldataclean{:,"Correct/Incorrect"}==99);
alldataclean{nanidx,"Correct/Incorrect"}=NaN;

% Objective Mean
means.objgrandmeanpercondition(1)=mean(alldataclean{indxrhy,"Correct/Incorrect"},"omitnan");
means.objgrandmeanpercondition(2)=mean(alldataclean{indxint,"Correct/Incorrect"},"omitnan");
means.objgrandmeanpercondition(3)=mean(alldataclean{indxirr,"Correct/Incorrect"},"omitnan");

% Subjective Mean
means.subjgrandmeanpercondition(1)=mean(alldataclean{indxrhy,"Binary Visibility"},"omitnan");
means.subjgrandmeanpercondition(2)=mean(alldataclean{indxint,"Binary Visibility"},"omitnan");
means.subjgrandmeanpercondition(3)=mean(alldataclean{indxirr,"Binary Visibility"},"omitnan");


% Calculate Objective Performance for each Contrast and Condition

fullresultsobj=[];

% Make table samples x contrasts (objective)
for idxcond=1:3 % for each condition
    conditiondata=alldataclean(find(alldataclean{:,"Condition"}==idxcond),:); % get trials
    results=[];
    for idxcont=1:10 % find each contrast
        contrastindex=find(conditiondata{:,"Contrast Level"}==idxcont);
        results(1:length(conditiondata{contrastindex,"Correct/Incorrect"}),idxcont)=conditiondata{contrastindex,"Correct/Incorrect"};
        results(length(conditiondata{contrastindex,"Correct/Incorrect"}):end,idxcont)=NaN;
    end
    fullresultsobj{idxcond}=results; % calculate mean across trials for each contrast
end

% Mean Objective

for idxcond=1:3
    means.objmeancondcontrast(idxcond,:)=mean(fullresultsobj{idxcond},1,'omitnan');
end

%%%% Calculate Subjective Performance for each Contrast and Condition

fullresultssubj=[];

% Make table samples x contrasts (subjective)
for idxcond=1:3 % for each condition
    conditiondata=alldataclean(find(alldataclean{:,"Condition"}==idxcond),:); % get trials
    results=[];
    for idxcont=1:10 % find each contrast
        contrastindex=find(conditiondata{:,"Contrast Level"}==idxcont);
        results(1:length(conditiondata{contrastindex,"Binary Visibility"}),idxcont)=conditiondata{contrastindex,"Binary Visibility"};
        results(length(conditiondata{contrastindex,"Binary Visibility"}):end,idxcont)=NaN;
    end
    fullresultssubj{idxcond}=results; % calculate mean across trials for each contrast
end

% Mean Subjective
means.subjmeancondcontrast=[];
for idxcond=1:3
    means.subjmeancondcontrast(idxcond,:)=mean(fullresultssubj{idxcond},1,'omitnan');
end

%%% Calculate Means Per Block
for block=1:nblocks

% Logical Indices
blockidx=find(alldataclean{:,"Block"} == block);

% Objective Mean
means.blockmeansobj(block,2)=mean(alldataclean{blockidx,"Correct/Incorrect"},'omitnan');
means.blockmeansobj(block,1)=mean(alldataclean{blockidx,"Condition"});

% Subjective Mean
means.blockmeanssubj(block,2)=mean(alldataclean{blockidx,"Binary Visibility"},'omitnan');
means.blockmeansobj(block,1)=mean(alldataclean{blockidx,"Condition"});

end

% Mean for different visibility ratings
idx0=find(alldataclean{:,"Binary Visibility"}==0);
idx1=find(alldataclean{:,"Binary Visibility"}==1);

means.accuracybinaryPAS(1)=mean(alldataclean{idx0,"Correct/Incorrect"});
means.accuracybinaryPAS(2)=mean(alldataclean{idx1,"Correct/Incorrect"});

idx0=find(alldataclean{:,"Visibility Response"}==0);
idx1=find(alldataclean{:,"Visibility Response"}==1);
idx2=find(alldataclean{:,"Visibility Response"}==2);
idx3=find(alldataclean{:,"Visibility Response"}==3);

means.accuracyperPAS(1)=mean(alldataclean{idx0,"Correct/Incorrect"});
means.accuracyperPAS(2)=mean(alldataclean{idx1,"Correct/Incorrect"});
means.accuracyperPAS(3)=mean(alldataclean{idx2,"Correct/Incorrect"});
means.accuracyperPAS(4)=mean(alldataclean{idx3,"Correct/Incorrect"});

end

%% Plot Means
function []=plotmeans(means)
% Accuracy across blocks
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
end
%% Prepare Data for Psignifit
function [dataforpsignifit]=prepforpsignifit(alldataclean)
% Calculate Mean for Each Contrast and Prepare Psignifit
for condition=1:3
    idxcondition=find(alldataclean{:,"Condition"} == condition);
    currentdata=alldataclean(idxcondition,:);
    for ncont=1:10
        idxcont=find(currentdata{:,"Contrast Level"} == ncont);

        % Objective
        tempdataforpsignifit.obj(ncont,3)=length(idxcont);
        tempdataforpsignifit.obj(ncont,2)=length(find(currentdata{idxcont,"Correct/Incorrect"} == 1));
        tempdataforpsignifit.obj(ncont,1)=ncont;
        
        % Subjective
        tempdataforpsignifit.subj(ncont,3)=length(idxcont);
        tempdataforpsignifit.subj(ncont,2)=length(find(currentdata{idxcont,"Binary Visibility"} == 1));
        tempdataforpsignifit.subj(ncont,1)=ncont;
    end

   dataforpsignifit{condition}=tempdataforpsignifit;
end
end
%% Run Psignifit
function [psignifitsresults,midpoints]=fitpsychcurves(dataforpsignifit,plots)
% Prepare struct with fitting options
fitting_options_obj=struct;
fitting_options_subj=struct;

fitting_options_obj.expType = 'nAFC';
fitting_options_obj.expN = 2;
fitting_options_subj.expType = 'YesNo';

fitting_options_obj.sigmoidName = 'logistic';
fitting_options_subj.sigmoidName = 'logistic';

% Fit Curves for each condition - data including outliers
for conditions=1:3

    % Objective
    psignifitsresults{conditions}.obj = psignifit(dataforpsignifit{1, conditions}.obj,fitting_options_obj);

    % Subjective
    psignifitsresults{conditions}.subj = psignifit(dataforpsignifit{1, conditions}.subj,fitting_options_subj);

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
legend('Rhythm','Interval')
end
end