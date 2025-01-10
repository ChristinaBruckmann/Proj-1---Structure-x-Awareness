%% Analysis of Alpha Phase
% Addition: Individual Trial data gets saved to that they can all be plotted and compared.

function SxA_AlphaPhaseExtraction(subj)

%% Extract ITPC from whole data (with different reference depending on cluster)

for s=1:length(subj)
disp('Starting Alpha Phase Analysis')

% Load Data
cd 'D:\'

loadfilename=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj(s));
savefilename=sprintf('EEG_SxA_Subj%i_AlphaPhaseSingleTrials.mat',subj(s));
load(loadfilename)

% Parameters
ITPCalpha_wavFreqs=[10*2^-0.3 10*2^0.3]; %confirm that this makes sense
srate=SDATA.info.sampling_rate;
ITPCalpha_elec= 1:64;

% Re-Reference for different clusters
if ismember(subj(s),[101,108,121]) % has a noisy nose, always refernce to mastoids
    data_occ=referenceContEEGdata(SDATA.data,[69 71]);
else
    data_occ=referenceContEEGdata(SDATA.data,71); %nose
end

% if subj==19 % has a noisy mastoid, always reference to nose
%     data_cen=data_occ;
% else
%     data_cen=referenceContEEGdata(SDATA.data,[69 71]); % mastoids
% end

% band-pass filter 0.5-2 Hz butterworth, 24dB/octave
bpFilteredData_occ = bandPassFilter(min(ITPCalpha_wavFreqs),max(ITPCalpha_wavFreqs),data_occ(:,ITPCalpha_elec),srate);
%bpFilteredData_cen = bandPassFilter(min(ITPCalpha_wavFreqs),max(ITPCalpha_wavFreqs),data_cen(:,ITPCalpha_elec),srate);

% Hilbert Transform
hil_occ = hilbert(bpFilteredData_occ);
%hil_cen = hilbert(bpFilteredData_cen);


% Extract Phase
alpha_phase_whole_occ = angle(hil_occ); % phase
%alpha_phase_whole_cen = angle(hil_cen); % phase

cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\AlphaPhaseRes'
save(savefilename, 'alpha_phase_whole_occ','-v7.3')
end
end