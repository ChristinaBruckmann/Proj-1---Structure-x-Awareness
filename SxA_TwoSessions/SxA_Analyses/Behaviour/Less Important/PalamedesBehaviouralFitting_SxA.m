%% SxA - Hierarchical Logistic Regression for Behavioural Data with Palamedes
% Christina Bruckmann last updated 25.03.24
clear
clc


% Fitting curves to individual subjcts
subj=[101:103 105 106 108 110:113];
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data\Behavioural Preprocessed'

for s=1:length(subj)
    % Load Data (Pal takes input similarly to psignifit, so we will use that matrix)
    loadfilename=sprintf('SxA_ResultsSubject%i_Total.mat',subj(s));
    load(loadfilename,'dataforpsignifit')

    % Model Input
    parmguess_obj=[1,0.5,0.5,0]; % threshold, slope, guess rate (fixed at 0.5), lapsereate
    parmguess_subj=[1,0.5,0,0]; % threshold, slope, guess rate (fixed at 0.5), lapsereate
    
    % Need to define reasonable areas where the function should search
    sub_searchGrid.alpha = [0:.01:9];    %threshold between 0 an 10
    sub_searchGrid.beta = 10.^[-1:.01:2]; %slope
    sub_searchGrid.gamma = [0:.01:.1]; % guess rate
    sub_searchGrid.lambda = [0:.01:.06]; % lapse rate

    obj_searchGrid=sub_searchGrid;
    obj_searchGrid.gamma = 0.45:0.01:0.55;
    PF=@PAL_Logistic;

    for c=1:3 % separate for all conditions, because there is no reason to assume that any parameters are the same across conditions

        % Prepare for Fitting Subj
        levels_sub=dataforpsignifit{1, c}.subj(:,1)'; % contrast levels
        ncorr_sub=dataforpsignifit{1, c}.subj(:,2)'; % n trials correct
        ntotal_sub=dataforpsignifit{1, c}.subj(:,3)'; % n trials total

        % Prepare for Fitting  Obj
        levels_obj=dataforpsignifit{1, c}.obj(:,1)'; % contrast levels
        ncorr_obj=dataforpsignifit{1, c}.obj(:,2)'; % n trials correct
        ntotal_obj=dataforpsignifit{1, c}.obj(:,3)'; % n trials total

        % Fit model simultaneously to all conditions (within subj or obj, not across obviously)
        params_sub(s,c,:)=PAL_PFML_Fit(levels_sub, ncorr_sub, ntotal_sub, sub_searchGrid, [1,1,0,1], PF);
        params_obj(s,c,:)=PAL_PFML_Fit(levels_obj, ncorr_obj, ntotal_obj, obj_searchGrid, [1,1,0,1], PF);

%         % Plot
%         StimLevelsFineGrain=[min(levels_sub):max(levels_sub)./1000:max(levels_sub)];
%         ProportionCorrectModel = PF(params_sub(8,3,:),StimLevelsFineGrain);
%         ProportionCorrectObserved=ncorr_sub./ntotal_sub; 
% 
%         figure('name','Maximum Likelihood Psychometric Function Fitting');
%         axes
%         hold on
%         plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',4);
%         plot(levels_sub,ProportionCorrectObserved,'k.','markersize',40);
%         set(gca, 'fontsize',16);
%         set(gca, 'Xtick',levels_sub);
%         axis([min(levels_sub) max(levels_sub) 0 1]);
%         xlabel('Stimulus Intensity');
%         ylabel('proportion correct');

        % Estimate goodness of fit
        [Dev_sub(s,c,:), pDev_sub(s,c,:), DevSim_sub(s,c,:), converged_sub(s,c,:)] = PAL_PFML_GoodnessOfFit(levels_sub, ncorr_sub, ntotal_sub, params_sub(s,c,:), [1,1,0,1], 1000, PF);
        [Dev_obj(s,c,:), pDev_obj(s,c,:), DevSim_obj(s,c,:), converged_obj(s,c,:)] = PAL_PFML_GoodnessOfFit(levels_obj, ncorr_obj, ntotal_obj, params_obj(s,c,:), [1,1,0,1], 1000, PF);
    end
    % Save parameters
end

%% Within subjects RM anova on slope and threhsold

% Only select subjects wiht a good enough fit of all curves
s_GoodFitIdx=all(pDev_sub>=0.05,2); % 0.05 as a cut off
o_GoodFitIdx=all(pDev_obj>=0.05,2);

anova_params_sub=params_sub(s_GoodFitIdx,:,:);
anova_params_obj=params_obj(o_GoodFitIdx,:,:);

