% Get data in blocks for psignifit
%% Trying the Psignifit Toolbox for Fitting Psychometric Curves
% Data needs to be arranged in a nx3 matrix (stimulus level | nCorrect | ntotal)
clear
clc

part=1;

resfilename=sprintf('ResultsSubjectOfficial%i.mat',part);
    load(resfilename)

% Prepare data structs
expdataforpsignifit.irrobj=(1:10)'; % Stimulus Levels
expdataforpsignifit.irrobj(:,3)=4; % nTotal
expdataforpsignifit.intobj=expdataforpsignifit.irrobj; % copy for other conditions
expdataforpsignifit.rhyobj=expdataforpsignifit.irrobj;
expdataforpsignifit.irrobj(:,3)=3;  % 4 blocks times 3 reps

expdataforpsignifit.irrsubj=(1:10)'; % Stimulus Levels
expdataforpsignifit.irrsubj(:,3)=4; % nTotal
expdataforpsignifit.intsubj=expdataforpsignifit.irrsubj; % copy for other conditions
expdataforpsignifit.rhysubj=expdataforpsignifit.irrsubj;
expdataforpsignifit.irrsubj(:,3)=3; % 4 blocks times 3 reps

% % Get Correct/Incorrect and binary subj response for each condition 
% indxirr=find(alldataclean{:,"Condition"} == 3);
% indxint=find(alldataclean{:,"Condition"} == 2);
% indxrhy=find(alldataclean{:,"Condition"} == 1);
% 
% irreval=alldataclean{indxirr,["Block" "Contrast Level" "Correct/Incorrect" "Binary Visibility"]};
% inteval=alldataclean{indxint,["Block" "Contrast Level" "Correct/Incorrect" "Binary Visibility"]};
% rhyeval=alldataclean{indxrhy,["Block" "Contrast Level" "Correct/Incorrect" "Binary Visibility"]};

% Calculate nCorrect per condition and Stimulus Level for each block

% Get Correct/incorrect and binary subj response for each stimulus Level
for nblock=1:10
    % Extract all trials from that block
    blockidx=find(alldataclean{:,"Block"} == nblock);
    allblockdata=alldataclean{blockidx,["Condition" "Contrast Level" "Correct/Incorrect" "Binary Visibility"]};
        for stimlev=1:10
        indxtemp=find(allblockdata(:,2) == stimlev);
        
        eval_currlevobj=allblockdata(indxtemp,3);
        
        eval_currlevsubj=allblockdata(indxtemp,4);
        
        blockres_obj(stimlev,1)=sum(eval_currlevobj);
        blockres_subj(stimlev,1)=sum(eval_currlevsubj);
        end

        if allblockdata(1,1)==1
            expdataforpsignifit.rhyobj(:,2)=[expdataforpsignifit.rhyobj(:,2);blockres_obj];
            expdataforpsignifit.rhysubj(2)=[expdataforpsignifit.rhysubj(:,2);blockres_subj];
        elseif allblockdata(1,1)==2
            expdataforpsignifit.intobj(2)=[expdataforpsignifit.intobj(2);blockres_obj];
            expdataforpsignifit.intsubj(2)=[expdataforpsignifit.intsubj(2);blockres_subj];
        elseif allblockdata(1,1)==3
            expdataforpsignifit.irrobj(2)=[expdataforpsignifit.irrobj(2);blockres_obj];
            expdataforpsignifit.irrsubj(2)=[expdataforpsignifit.irrsubj(2);blockres_subj];
        end
        
        expdataforpsignifit.irrobj(stimlevel,2)=sum(irreval_currlevobj);
        expdataforpsignifit.intobj(stimlevel,2)=sum(inteval_currlevobj);
        expdataforpsignifit.rhyobj(stimlevel,2)=sum(rhyeval_currlevobj);
        
        expdataforpsignifit.irrsubj(stimlevel,2)=sum(irreval_currlevsubj);
        expdataforpsignifit.intsubj(stimlevel,2)=sum(inteval_currlevsubj);
        expdataforpsignifit.rhysubj(stimlevel,2)=sum(rhyeval_currlevsubj);
end