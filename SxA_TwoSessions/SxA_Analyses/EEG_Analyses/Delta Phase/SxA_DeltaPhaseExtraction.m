%% Re-Analysis of Delta Phase from scratch to make sure there are no bugs in the first script. Only occipital for now.
% Addition: Individual Trial data gets saved to that they can all be
% plotted and compared.

function SxA_DeltaPhaseExtraction(subj)

%% Extract ITPC from whole data (with different reference depending on cluster)

for s=1:length(subj)
disp('Starting Delta Phase Analysis')

% Load Data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Preprocessed'

loadfilename=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj(s));
savefilename=sprintf('EEG_SxA_Subj%i_DeltaPhaseSingleTrials.mat',subj(s));
load(loadfilename)

% Parameters
ITPCdelta_wavFreqs=[1.25*2^-0.75 1.25*2^0.75];
srate=SDATA.info.sampling_rate;
ITPCdelta_elec= 1:64;

% Re-Reference for different clusters
if subj(s)==101 % has a noisy nose, always refernce to mastoids
    data_occ=referenceContEEGdata(SDATA.data,[69 71]);
else
    data_occ=referenceContEEGdata(SDATA.data,71); %nose
end

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
delta_phase_whole_occ = angle(hil_occ); % phase
delta_phase_whole_cen = angle(hil_cen); % phase

cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\DeltaRes'
save(savefilename, 'delta_phase_whole_occ','delta_phase_whole_cen','-v7.3')
end
end