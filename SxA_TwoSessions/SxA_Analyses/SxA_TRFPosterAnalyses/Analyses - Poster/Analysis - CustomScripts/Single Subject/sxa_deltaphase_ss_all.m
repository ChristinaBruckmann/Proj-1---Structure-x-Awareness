%% Extract ITPC from whole data
function []=sxa_deltaphase_ss_all(subj)

disp('Starting Delta Phase Analysis')

% Load Data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG'

loadfilename=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj);
savefilename=sprintf('EEG_SxA_Subj%i_Results_deltaall.mat',subj);
load(loadfilename)

% Parameters
ITPCdelta_wavFreqs=[1.25*2^-0.75 1.25*2^0.75];
srate=SDATA.info.sampling_rate;
ITPCdelta_elec= 1:64;

% Re-Reference for different clusters
data_occ=referenceContEEGdata(SDATA.data,71); %nose

if subj==19 % has a noisy mastoid, always reference to nose
data_cen=data_occ;
else
data_cen=referenceContEEGdata(SDATA.data,[69 71]); % mastoids
end

% band-pass filter 0.5-2 Hz butterworth, 24dB/octave
bpFilteredData_occ = bandPassFilter(min(ITPCdelta_wavFreqs),max(ITPCdelta_wavFreqs),data_occ(:,ITPCdelta_elec),srate);
bpFilteredData_cen = bandPassFilter(min(ITPCdelta_wavFreqs),max(ITPCdelta_wavFreqs),data_cen(:,ITPCdelta_elec),srate);

% Hilbert Transform
hil_occ = hilbert(bpFilteredData_occ);
hil_cen = hilbert(bpFilteredData_cen);

% Extract Phase
delta_ITPC_whole_occ = angle(hil_occ); % phase
delta_ITPC_whole_cen = angle(hil_cen); % phase

cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
save(savefilename, 'delta_ITPC_whole_occ','delta_ITPC_whole_cen','-v7.3')
end