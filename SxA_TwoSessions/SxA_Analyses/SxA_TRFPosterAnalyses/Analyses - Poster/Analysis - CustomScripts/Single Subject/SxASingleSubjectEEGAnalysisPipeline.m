%% SxA Single Subject EEG Analysis Pipeline 
clear
clc

subj=input("Subject Number? ");

% Sanity Check - Evoked Response to Mask Onset
sxa_erpsanity_ss(subj, 1)
sanitypassed=input('Continue? (1-yes) ');
if ~sanitypassed==1
    return
end
disp('Sanity Check Completed')

% CNV
sxa_cnv_ss(subj,1)

% Alpha Amplitude (WS to predicted target)
sxa_TF_alpaamp(subj,1) % Run TF Analysis, save results, and create plots

% Delta Phase
sxa_deltaphase_ss(subj,1) % Run delta phase analysis, save results, and create plots

disp('Analysis Complete - Results Saved')