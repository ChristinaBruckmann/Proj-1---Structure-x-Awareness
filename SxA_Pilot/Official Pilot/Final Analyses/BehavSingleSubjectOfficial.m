%% Behavioural Analysis Structure X Awareness Data

clc
clear

load PilotOfficialSxA1.1_s11.mat
filename1=sprintf('ResultsSubjectOfficial%i',subresults.data{1,"Subject"});

plots=1; % create plots?
newfitmethod=1;
%% Clean Up Data File & Basic Checks

% Check frame rate:
meanfrrate=mean(subresults.data{:,"Frame Rate"});

if meanfrrate==60
    sprintf('Frame rate check passed. Mean frame rate across all trials is 60hz.')
else
    sprintf('ATTENTION. Frame rate check FAILED. Mean frame rate across all trials is %i.',meanfrrate)
end

% Add variable: Categorize PAS responses into binary data

set(0, 'DefaultFigurePosition', [430 230 560 420])

BinaryVisibility=zeros(1,height(subresults.data));

for ivis=1:height(subresults.data)
    if subresults.data{ivis,"Visibility Response"}== 0 % invisible
        BinaryVisibility(ivis)=0;
    elseif subresults.data{ivis,"Visibility Response"}==99 % missing
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
%% Sanity Check: Plot graded PAS vs. Accuracy

PASratings=0:3;

idx0=find(subresults.data{:,"Visibility Response"}==0);
idx1=find(subresults.data{:,"Visibility Response"}==1);
idx2=find(subresults.data{:,"Visibility Response"}==2);
idx3=find(subresults.data{:,"Visibility Response"}==3);

objmeanperPASrating(1)=mean(subresults.data{idx0,"Correct/Incorrect"});
objmeanperPASrating(2)=mean(subresults.data{idx1,"Correct/Incorrect"});
objmeanperPASrating(3)=mean(subresults.data{idx2,"Correct/Incorrect"});
objmeanperPASrating(4)=mean(subresults.data{idx3,"Correct/Incorrect"});

if plots
f1=figure;
plot(PASratings,objmeanperPASrating,'o','LineWidth',2)
title('Subjective Perception vs. Accuracy')
xlabel('PAS Level')
ylabel('Correct in %')
xlim([-1 4])
ylim([0.4 1.1])
xticks(PASratings)
end
%% Plot Binary PAS vs. Accuracy

idx0=find(subresults.data{:,"Binary Visibility"}==0);
idx1=find(subresults.data{:,"Binary Visibility"}==1);

binobjsubjmean(1)=mean(subresults.data{idx0,"Correct/Incorrect"});
binobjsubjmean(2)=mean(subresults.data{idx1,"Correct/Incorrect"});

if plots
f2=figure;
plot([0 1],binobjsubjmean,'o','LineWidth',2)
title('Subjective Perception Binary vs. Accuracy')
xlabel('Visibility')
ylabel('Correct in %')
xlim([-1 2])
ylim([0.4 1.1])
xticks([0 1])
end
%% Check for orientation bias (and remove?) - across conditions?

% Hit rate & false alarm rate (vertical=0, horizontal=1)
idxallver=find(subresults.data{:,"Orientation"}==0); % all trial indices for vertical
idxallhor=find(subresults.data{:,"Orientation"}==1); % all trial indices for horizontal

totalver = length(idxallver); % N total vertical targets
totalhor = length(idxallhor); % N total horizontal targets

ansver=subresults.data{idxallver,"Correct/Incorrect"}; % response evaluation to all vertical targets (0 incorr and 1 corr)
anshor=subresults.data{idxallhor,"Correct/Incorrect"}; % response evaluation to all horizontal targets (0 incorr and 1 corr)

HIThor = length(find(anshor==1)); % hit rate horizontal (given the target was horizontal, you said horizontal)
FAver = length(find(anshor==0)); % false alarm vertical (given the target was horizontal, you said vertical)

HITver= length(find(ansver==1)); % hit rate vertical
FAhor = length(find(ansver==0));  % false alarm horizontal (given the target was vertical, answer was horizontal)

SDTresults.hitrateprobver=HITver/totalver;
SDTresults.hitrateprobhor=HIThor/totalhor;
SDTresults.falsealarmver=FAver/totalhor;
SDTresults.falsealarmhor=FAhor/totalver;

% Criterion
SDTresults.criterionhor=-(norminv(SDTresults.hitrateprobhor) + norminv(SDTresults.falsealarmhor)) / 2;
SDTresults.criterionver=-(norminv(SDTresults.hitrateprobver) + norminv(SDTresults.falsealarmver)) / 2;

% d' for all conditions combined
SDTresults.dprimevert=norminv(SDTresults.falsealarmver)-norminv(SDTresults.hitrateprobver);
SDTresults.dprimehor=norminv(SDTresults.falsealarmhor)-norminv(SDTresults.hitrateprobhor);

%% Remove Additional Irregular Trials before Analysis

idxaddirr1=find(subresults.data{:,"Irregular Target Time"} == 1);
idxaddirr2=find(subresults.data{:,"Irregular Target Time"} == 5);

alldataclean=subresults.data;
idxaddirrall=sort([idxaddirr1; idxaddirr2]);
alldataclean(idxaddirrall,:)=[];

%% Means for all conditions

conditions=1:3;

% Logical Indices
indxirr=find(alldataclean{:,"Condition"} == 3);
indxint=find(alldataclean{:,"Condition"} == 2);
indxrhy=find(alldataclean{:,"Condition"} == 1);


% Objective Mean
objgrandmeanpercondition(1)=mean(alldataclean{indxrhy,"Correct/Incorrect"});
objgrandmeanpercondition(2)=mean(alldataclean{indxint,"Correct/Incorrect"});
objgrandmeanpercondition(3)=mean(alldataclean{indxirr,"Correct/Incorrect"});

% Subjective Mean
subjgrandmeanpercondition(1)=mean(alldataclean{indxrhy,"Binary Visibility"});
subjgrandmeanpercondition(2)=mean(alldataclean{indxint,"Binary Visibility"});
subjgrandmeanpercondition(3)=mean(alldataclean{indxirr,"Binary Visibility"});

if plots
f3=figure;
plot(conditions,objgrandmeanpercondition,'o','LineWidth',2)
title('Accuracy per Condition')
xlabel('Condition')
ylabel('Correct in %')
xlim([0 4])
ylim([0 1.1])
xticks([1 2 3])


f4=figure;
plot(conditions,subjgrandmeanpercondition,'o','LineWidth',2)
title('Subjective Perception per Condition')
xlabel('Condition')
ylabel('Subjective Perception')
xlim([0 4])
ylim([0 1.1])
xticks([1 2 3])
end
%% Calculate Objective Performance for each Contrast and Condition

% N Contrasts
nContrasts=10;

fullresultsobj=[];

% Make table samples x contrasts (objective)
for idxcond=1:3 % for each condition
    conditiondata=alldataclean(find(alldataclean{:,"Condition"}==idxcond),:); % get trials
    results=[];
    for idxcont=1:10 % find each contrast
        contrastindex=find(conditiondata{:,"Contrast Level"}==idxcont);
        results(:,idxcont)=conditiondata{contrastindex,"Correct/Incorrect"};
    end
    fullresultsobj{idxcond}=results; % calculate mean across trials for each contrast
end

% Mean Objective

for idxcond=1:3
    objmeancondcontrast(idxcond,:)=mean(fullresultsobj{idxcond},1,'omitnan');
end

%% Calculate Subjective Performance for each Contrast and Condition

% Model Function
fullresultssubj=[];

% Make table samples x contrasts (subjective)
for idxcond=1:3 % for each condition
    conditiondata=alldataclean(find(alldataclean{:,"Condition"}==idxcond),:); % get trials
    results=[];
    for idxcont=1:10 % find each contrast
        contrastindex=find(conditiondata{:,"Contrast Level"}==idxcont);
        results(:,idxcont)=conditiondata{contrastindex,"Binary Visibility"};
    end
    fullresultssubj{idxcond}=results; % calculate mean across trials for each contrast
end

% Mean Subjective
subjmeancondcontrast=[];
for idxcond=1:3
    subjmeancondcontrast(idxcond,:)=mean(fullresultssubj{idxcond},1,'omitnan');
end

%% Preparing Data for Psignifit
% Data needs to be arranged in a nx3 matrix (stimulus level | nCorrect | ntotal)

% Prepare data structs
dataforpsignifit.irrobj=(1:10)'; % Stimulus Levels
dataforpsignifit.irrobj(:,3)=12; % nTotal
dataforpsignifit.intobj=dataforpsignifit.irrobj; % copy for other conditions
dataforpsignifit.rhyobj=dataforpsignifit.irrobj;

dataforpsignifit.irrsubj=(1:10)'; % Stimulus Levels
dataforpsignifit.irrsubj(:,3)=12; % nTotal
dataforpsignifit.intsubj=dataforpsignifit.irrsubj; % copy for other conditions
dataforpsignifit.rhysubj=dataforpsignifit.irrsubj;

% Get Correct/incorrect and binary subj response for each condition 
indxirr=find(alldataclean{:,"Condition"} == 3);
indxint=find(alldataclean{:,"Condition"} == 2);
indxrhy=find(alldataclean{:,"Condition"} == 1);

irreval=alldataclean{indxirr,["Contrast Level" "Correct/Incorrect" "Binary Visibility"]};
inteval=alldataclean{indxint,["Contrast Level" "Correct/Incorrect" "Binary Visibility"]};
rhyeval=alldataclean{indxrhy,["Contrast Level" "Correct/Incorrect" "Binary Visibility"]};

% Calculate nCorrect per condition and Stimulus Level

% Get Correct/incorrect and binary subj response for each stimulus Level
for stimlev=1:10
indxtempirr=find(irreval(:,1) == stimlev);
indxtempint=find(inteval(:,1) == stimlev);
indxtemprhy=find(rhyeval(:,1) == stimlev);

irreval_currlevobj=irreval(indxtempirr,2);
inteval_currlevobj=inteval(indxtempint,2);
rhyeval_currlevobj=rhyeval(indxtemprhy,2);

irreval_currlevsubj=irreval(indxtempirr,3);
inteval_currlevsubj=inteval(indxtempint,3);
rhyeval_currlevsubj=rhyeval(indxtemprhy,3);

dataforpsignifit.irrobj(stimlev,2)=sum(irreval_currlevobj);
dataforpsignifit.intobj(stimlev,2)=sum(inteval_currlevobj);
dataforpsignifit.rhyobj(stimlev,2)=sum(rhyeval_currlevobj);

dataforpsignifit.irrsubj(stimlev,2)=sum(irreval_currlevsubj);
dataforpsignifit.intsubj(stimlev,2)=sum(inteval_currlevsubj);
dataforpsignifit.rhysubj(stimlev,2)=sum(rhyeval_currlevsubj);
end


% Put all into one struct
save(filename1, 'alldataclean', 'dataforpsignifit', 'objmeanperPASrating','objgrandmeanpercondition','subjgrandmeanpercondition','subjmeancondcontrast','objmeancondcontrast','subresults','SDTresults')
