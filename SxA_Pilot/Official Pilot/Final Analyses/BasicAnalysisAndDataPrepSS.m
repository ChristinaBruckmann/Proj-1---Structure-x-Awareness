%% Full Analysis Script Single Subject - From Raw Data to Curves and Stats
clc
clear

totalsubn=13;
nblocks=10;
for subj_n=1:totalsubn
% Load Data
filename_toload=sprintf('PilotOfficialSxA1.1_s%i.mat',subj_n);
load(filename_toload)
filename=sprintf('ResultsSubjectOfficial%i',subresults.data{1,"Subject"});

%% DATA PREP
% Check frame rate:
meanfrrate=mean(subresults.data{:,"Frame Rate"});

if meanfrrate==60
    sprintf('Frame rate check passed. Mean frame rate across all trials is 60hz.')
else
    sprintf('ATTENTION. Frame rate check FAILED. Mean frame rate across all trials is %i.',meanfrrate)
end

% Clean Data 1 (Binary PAS, Integrate Trial Matrix into Results File)
subresults=clean_1(subresults);
save(filename,'subresults')

% Remove Additional Irregular Trials before Analysis

idxaddirr1=find(subresults.data{:,"Irregular Target Time"} == 1);
idxaddirr2=find(subresults.data{:,"Irregular Target Time"} == 5);

alldataclean=subresults.data;
idxaddirrall=sort([idxaddirr1; idxaddirr2]);
alldataclean(idxaddirrall,:)=[];

save(filename, 'alldataclean',"-append")

%% Calculate and save all the means (per condition, per contrast, per rating etc.)
means=calculate_means(alldataclean,nblocks);
% Plot Means Per Block
idx_rhy=find(means.blockmeansobj(:,1)==1);
idx_int=find(means.blockmeansobj(:,1)==2);
idx_irr=find(means.blockmeansobj(:,1)==3);
p1(subj_n)=subplot(5,3,subj_n); plot(1:nblocks,means.blockmeansobj(:,2))
plot(idx_rhy,means.blockmeansobj(idx_rhy,2),'-o', 'Color', 'b');
hold on
plot(idx_int,means.blockmeansobj(idx_int,2),'-o', 'Color', 'r');
plot(idx_irr,means.blockmeansobj(idx_irr,2),'-o', 'Color', 'g');
legend('Rhythm','Interval','Irregular')
hold off
save(filename,'means',"-append")
%% Remove Outliers Based on RT (automatically saves new data file and give the amount of rejected trials as output)
n_outliers=Remove_RT_Outliers(alldataclean,filename);
%% Prepare Data for Psignifit
[dataforpsignifit]=psignifit_prep(alldataclean);
save(filename, 'dataforpsignifit',"-append")
end

linkaxes(p1,'xy')
%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function: Prepare Data for Psignifit
% Data needs to be arranged in a nx3 matrix (stimulus level | nCorrect | ntotal)

function [dataforpsignifit]=psignifit_prep(alldataclean)
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
end
%% Function: Calculate all the Means
function [means]=calculate_means(alldataclean,nblocks)
%%%%% Means for all conditions

conditions=1:3;

% Logical Indices
indxirr=find(alldataclean{:,"Condition"} == 3);
indxint=find(alldataclean{:,"Condition"} == 2);
indxrhy=find(alldataclean{:,"Condition"} == 1);


% Objective Mean
means.objgrandmeanpercondition(1)=mean(alldataclean{indxrhy,"Correct/Incorrect"});
means.objgrandmeanpercondition(2)=mean(alldataclean{indxint,"Correct/Incorrect"});
means.objgrandmeanpercondition(3)=mean(alldataclean{indxirr,"Correct/Incorrect"});

% Subjective Mean
means.subjgrandmeanpercondition(1)=mean(alldataclean{indxrhy,"Binary Visibility"});
means.subjgrandmeanpercondition(2)=mean(alldataclean{indxint,"Binary Visibility"});
means.subjgrandmeanpercondition(3)=mean(alldataclean{indxirr,"Binary Visibility"});


%%%%% Calculate Objective Performance for each Contrast and Condition

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
        results(:,idxcont)=conditiondata{contrastindex,"Binary Visibility"};
    end
    fullresultssubj{idxcond}=results; % calculate mean across trials for each contrast
end

% Mean Subjective
means.subjmeancondcontrast=[];
for idxcond=1:3
    means.subjmeancondcontrast(idxcond,:)=mean(fullresultssubj{idxcond},1,'omitnan');
end

%%%%%% Calculate Means Per Block
for block=1:nblocks

% Logical Indices
blockidx=find(alldataclean{:,"Block"} == block);

% Objective Mean
means.blockmeansobj(block,2)=mean(alldataclean{blockidx,"Correct/Incorrect"});
means.blockmeansobj(block,1)=mean(alldataclean{blockidx,"Condition"});

% Subjective Mean
means.blockmeanssubj(block,2)=mean(alldataclean{blockidx,"Binary Visibility"});
means.blockmeansobj(block,1)=mean(alldataclean{blockidx,"Condition"});

end

end
%% Function: Clean and Sort Data 1
function [subresults]=clean_1(subresults)
% Add variable: Categorize PAS responses into binary data
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
end
%% Function: Remove Outliers Based on RT
function [nrejected]=Remove_RT_Outliers(alldataclean,filename)

    % Get all reaction times and standardize
    allRTs=alldataclean{:,"Reaction Time"};
    RT_zscore=zscore(allRTs); % zscore

    % Outliers based on STD
    outlieridx_std=find(RT_zscore>3 | RT_zscore<(-3));% Index above and below 3
    nrejected=length(outlieridx_std); % Calculate total number of rejected trials
    %figtitle1=sprintf("Participant: %i, Rejected trials: %i",part,nrejected);
    %figure(1); sp1(part)=subplot(3,4,part); histogram(RT_zscore); xline(3); title(figtitle1); 

%     % Determine if they are all from one condition
%     outlierconditions=alldataclean{outlieridx_std,"Condition"};
%     figtitle2=sprintf("Participant: %i, Condition of Outliers",part);
%     figure(2); sp2(part)=subplot(3,4,part); histogram(outlierconditions); title(figtitle2);
% 
%     % Determine if they are all from one contrast
%     outliercontrasts=alldataclean{outlieridx_std,"Contrast Level"};
%     figtitle3=sprintf("Participant: %i, Contrast of Outliers",part);
%     figure(3); sp3(part)=subplot(3,4,part); histogram(outliercontrasts); title(figtitle3);
% 
%     % Determine if all of them are from the same visibility level
%     outliercontrasts=alldataclean{outlieridx_std,"Visibility Response"};
%     figtitle4=sprintf("Participant: %i, Visibility of Outliers",part);
%     figure(4); sp4(part)=subplot(3,4,part); histogram(outliercontrasts); title(figtitle4);
% 
%     % Determine if all of them are correct/incorrect
%     outliercontrasts=alldataclean{outlieridx_std,"Correct/Incorrect"};
%     figtitle5=sprintf("Participant: %i, Accuracy of Outliers",part);
%     figure(5); sp5(part)=subplot(3,4,part); histogram(outliercontrasts); title(figtitle5);

    % Remove Outliers
    tempalldata=alldataclean;
    tempalldata(outlieridx_std,:)=[];
    alldataclean_outliersrejec=tempalldata;

    %%% Prepare cleaned data for psignifit
    % Data needs to be arranged in a nx3 matrix (stimulus level | nCorrect | ntotal)
    
    % Prepare data structs
    red_dataforpsignifit.irrobj=(1:10)'; % Stimulus Levels
    red_dataforpsignifit.intobj=red_dataforpsignifit.irrobj; % copy for other conditions
    red_dataforpsignifit.rhyobj=red_dataforpsignifit.irrobj;
    
    red_dataforpsignifit.irrsubj=(1:10)'; % Stimulus Levels
    red_dataforpsignifit.intsubj=red_dataforpsignifit.irrsubj; % copy for other conditions
    red_dataforpsignifit.rhysubj=red_dataforpsignifit.irrsubj;
    
    % Get Correct/incorrect and binary subj response for each condition 
    indxirr=find(alldataclean_outliersrejec{:,"Condition"} == 3);
    indxint=find(alldataclean_outliersrejec{:,"Condition"} == 2);
    indxrhy=find(alldataclean_outliersrejec{:,"Condition"} == 1);
    
    irreval=alldataclean_outliersrejec{indxirr,["Contrast Level" "Correct/Incorrect" "Binary Visibility"]};
    inteval=alldataclean_outliersrejec{indxint,["Contrast Level" "Correct/Incorrect" "Binary Visibility"]};
    rhyeval=alldataclean_outliersrejec{indxrhy,["Contrast Level" "Correct/Incorrect" "Binary Visibility"]};
    
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
    
    red_dataforpsignifit.irrobj(stimlev,2)=sum(irreval_currlevobj);
    red_dataforpsignifit.intobj(stimlev,2)=sum(inteval_currlevobj);
    red_dataforpsignifit.rhyobj(stimlev,2)=sum(rhyeval_currlevobj);
    
    red_dataforpsignifit.irrsubj(stimlev,2)=sum(irreval_currlevsubj);
    red_dataforpsignifit.intsubj(stimlev,2)=sum(inteval_currlevsubj);
    red_dataforpsignifit.rhysubj(stimlev,2)=sum(rhyeval_currlevsubj);
    
    % Get the total n of trials at this level for each condition
    red_dataforpsignifit.irrobj(stimlev,3)=length(find(alldataclean_outliersrejec{indxirr,"Contrast Level"} == stimlev));
    red_dataforpsignifit.intobj(stimlev,3)=length(find(alldataclean_outliersrejec{indxint,"Contrast Level"} == stimlev));
    red_dataforpsignifit.rhyobj(stimlev,3)=length(find(alldataclean_outliersrejec{indxrhy,"Contrast Level"} == stimlev));

    red_dataforpsignifit.irrsubj(stimlev,3)=red_dataforpsignifit.irrobj(stimlev,3);
    red_dataforpsignifit.intsubj(stimlev,3)=red_dataforpsignifit.intobj(stimlev,3);
    red_dataforpsignifit.rhysubj(stimlev,3)=red_dataforpsignifit.rhyobj(stimlev,3);
    end


    % Save 
    save(filename, 'alldataclean_outliersrejec', 'red_dataforpsignifit' ,"-append")
end