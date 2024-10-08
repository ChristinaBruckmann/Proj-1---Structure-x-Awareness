%% SxA Group Level Psignifit-Fit
% Determine which participants to exclude based on bad fit values.
% Bad fit defined as p-value of deviance (results.Fit(3)) smaller than 0.05
% in cases of overdispersion, deviance might incorrectly signal a bad fit.
% This code corrects for overdispersion.

% Parameters
subjects=[101:103 105:106 108 110 111 112 113 114 116 117 118 119 121:124 126 127 129 130 131 132];
bad_fit_cutoff=0.05;
corr_overdisp=1; % correct for overdispersion?


cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\Behavioural Preprocessed'

% Load all the data and save fit value
for s=1:length(subjects)
    loadfilename=sprintf('SxA_ResultsSubject%i_Total.mat',subjects(s));
    load(loadfilename, 'psignifitsresults')
    fit_value(s,:)=[subjects(s),psignifitsresults{1, 1}.obj.Fit(3),psignifitsresults{1, 1}.subj.Fit(3)]; % subject number, objective curve fit, subjective curve fit.
end

% Extract the dispersion-corrected deviance (deviance / dispersion)
deviance = result.Fit(2);          % Extract the deviance
dispersion = result.Fit(4);        % Extract the dispersion parameter
adjusted_deviance = deviance / dispersion;

% Extract the adjusted p-value
adjusted_p_value = result.Fit(5);  % Extract the p-value adjusted for dispersion

% Display results
fprintf('Deviance: %.2f\n', deviance);
fprintf('Dispersion: %.2f\n', dispersion);
fprintf('Adjusted Deviance: %.2f\n', adjusted_deviance);
fprintf('Adjusted p-value: %.4f\n', adjusted_p_value);

% Assess the fit based on adjusted p-value
if adjusted_p_value > 0.05
    disp('Good fit (after accounting for overdispersion).');
else
    disp('Poor fit (even after accounting for overdispersion).');
end

% Determine bad subjects
bad_fit_obj=subjects(fit_value(:,2)<bad_fit_cutoff);
bad_fit_subj=subjects(fit_value(:,3)<bad_fit_cutoff);

% Print results to command window
disp('<strong>Bad Objective Fit:</strong>')
cprintf(bad_fit_obj)
disp('<strong>Bad Subjective Fit:</strong>')