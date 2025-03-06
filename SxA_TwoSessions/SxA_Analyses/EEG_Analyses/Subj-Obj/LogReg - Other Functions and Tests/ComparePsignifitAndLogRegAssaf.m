% Compare Results for MultipleLogReg and Psignifit

clear 
clc

%% Load and extract multlogreg thresh
% Load
cd 'Y:\el-Christina\SxA\SxA_Results\LogRegResults\ContrastOnly'
load('Y:\el-Christina\SxA\SxA_Results\LogRegResults\ContrastOnly\EEG_SxA_LogRes_GL_Contrast_Simple.mat')

% Remove 107 (because of corrupted behavioural data)
rem_idx=subj==107;
res_obj(rem_idx,:,:)=[];
res_subj(rem_idx,:,:)=[];
subj(rem_idx)=[];

% Extract threshold value
res_obj=squeeze(res_obj); % remove empty dimensions (out: subj x cond x params)
res_subj=squeeze(res_subj);

res_obj_slope=res_obj(:,:,2); % get slope for contrast
res_subj_slope=res_subj(:,:,2);

res_obj_inter=res_obj(:,:,1); % get intercept
res_subj_inter=res_subj(:,:,1);

res_obj_thr=-(res_obj_inter./res_obj_slope); % get threshold 
res_subj_thr=-(res_subj_inter./res_subj_slope);

clear res_obj_inter res_subj_inter res_obj_slope res_subj_slope res_subj res_obj
%% Load and extract psignifit

res_obj_all(:,:,1)=res_obj_thr; % dimensions: subjects, conditions, estimation method (1-assaf, 2-psignifit)
res_subj_all(:,:,1)=res_subj_thr; % dimensions: subjects, conditions, estimation method (1-assaf, 2-psignifit)
for s=1:length(subj)
    cd 'Y:\el-Christina\SxA\SxA_Data\Behaviour Preprocessed'
    % Load Data
    loadfilename=sprintf('SxA_ResultsSubject%i_Total.mat',subj(s));
    load(loadfilename,'midpoints')
    
    % Add to Array
    res_obj_all(s,1,2)=midpoints.rhyobj;
    res_obj_all(s,2,2)=midpoints.intobj;
    res_obj_all(s,3,2)=midpoints.irrobj;
    res_subj_all(s,1,2)=midpoints.rhysubj;
    res_subj_all(s,2,2)=midpoints.intsubj;
    res_subj_all(s,3,2)=midpoints.irrsubj;
end

%% Calculate Difference Between Method 1 and 2
obj_diff=res_obj_all(:,:,1)-res_obj_all(:,:,2);
subj_diff=res_subj_all(:,:,1)-res_subj_all(:,:,2);

overall_meandiff=mean([mean(obj_diff,'all') mean(subj_diff,'all')]);
meandiff_subj=mean([mean(obj_diff,2) mean(subj_diff,2)],2);
mean_abs_diff=mean(abs(obj_diff),'all')+mean(abs(obj_diff),'all');

%% Calculate Correlation between thresholds
corr1=corr(res_obj_all(:,1,1),res_obj_all(:,1,2));
corr2=corr(res_obj_all(:,2,1),res_obj_all(:,2,2));
corr3=corr(res_obj_all(:,3,1),res_obj_all(:,3,2));
mean_corr2=mean([corr1 corr2 corr3]);
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Analyses\EEG_Analyses\Subj-Obj\LogReg - Other Functions and Tests'