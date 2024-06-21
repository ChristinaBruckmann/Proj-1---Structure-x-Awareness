function []=sxa_TF_alpaamp(subj,plots,catchonly)
%% TF Analysis - Alpha Amplitude
% Input: Subject Number, generate plots (1/0)
irrtartimes=[3 4 5]; % when irregular targets can appear
basec=0; %baselinecorrection?
maskonly=0; % alpha power around mask onsest only?
firstinterval=0; % compare alpha power of first interval to second interval
block_begin=1; % only the beginning 5 seconds of each block?
power=0; % power==1, amplitude==0
artrej=1; % reject artifact trials before saving? (automaticallt disabled for block_begin to identify the first trial of a block

if maskonly || firstinterval || block_begin% catch only not relevant for mask onset
    catchonly=0;
end

if maskonly && firstinterval || maskonly && block_begin || firstinterval && block_begin
    error('Incompatible time-windows to analyze. Set only one of them to 1.')
end

% clc
% clear 
disp('Starting Alpha Amplitude Analysis')

% % Parameters
if catchonly
    triggercodes={1;2;3};
    timerange=[-1000 700];
elseif maskonly
    triggercodes={31;32;33};
    timerange=[-600 700];
elseif firstinterval
    triggercodes={42};
    timerange=[-200 1500];
elseif block_begin
    triggercodes={31;32;33}; % no code for start of block, but this takes all trials beginnings, later selects the first trial of each block
    timerange=[-4000 200]; %back five seconds from the first trial of each block
else
    triggercodes={71;72;73}; % Warning Signals per condition
    timerange=[-200 1500];
end

% Parameters (mask onset alpha)


paddingLength=500; % in ms
alpha_wavFreqs=1:40; % Range taken from Elmira, Assaf uses 1:30 in plos bio
alpharange=8:12;

% Load Data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Preprocessed'
% subj=input("Subject Number? ");
%occelonly=input("Display occipital electrodes only (1-yes, 0-all)? ");
occelonly=1;
loadfilename=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj);
if catchonly
    savefilename=sprintf('EEG_SxA_Subj%i_AlphaResults_Catch.mat',subj);
elseif maskonly
     savefilename=sprintf('EEG_SxA_Subj%i_AlphaResults_Mask.mat',subj);
elseif firstinterval
    savefilename=sprintf('EEG_SxA_Subj%i_AlphaResults_FirstInterval.mat',subj);
elseif block_begin
    artrej=0; %do not automatically reject artifacts to be able to count the block beginning
    savefilename=sprintf('EEG_SxA_Subj%i_AlphaResults_BlockBegin.mat',subj);
elseif artrej
    savefilename=sprintf('EEG_SxA_Subj%i_AlphaResults_clean.mat',subj);
else
     savefilename=sprintf('EEG_SxA_Subj%i_AlphaResults_raw.mat',subj);
end
load(loadfilename)

% Extract Info and Data from File
data=SDATA.data;
artifacts=SDATA.metadata.artifacts;
triggers=SDATA.events.triggerChannel;
channelnames=SDATA.info.channel_labels;
srate=SDATA.info.sampling_rate;
if occelonly
alpha_tfElectrodes=[25:30 62:64]; % Occipital only
else
alpha_tfElectrodes=1:64; % All electrodes
end

% Segmentation
for c=1:size(triggercodes,1) % For conditions
    [segmentedData, isNotArtifact, timeVec]=segmentContEEGdata(triggercodes{c}, timerange+[-paddingLength paddingLength], data, triggers, artifacts, srate);


        % Get index for irregular trials without target in time window if not catch trials or mask only
        if c==3 && ~catchonly && ~maskonly && ~block_begin
            cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\Behavioural Preprocessed'
            load(sprintf('SxA_ResultsSubject%i_Session2.mat',subj),'alldataclean')
            idx_notar=ismember(alldataclean{(alldataclean{:,'Condition'}==c),'Irregular Target Time'},irrtartimes);
        else
            idx_notar=ones(size(segmentedData,3));
        end

        % Check how many trials are artifact-free
        sprintf('Proportion of artifact-free trials: %.2f', mean(isNotArtifact))

        % Remove trials with artifacts
        %artrej=input("Remove artifact trials (1=yes)? ");
        if artrej==1 && c==3 && ~catchonly && ~maskonly && ~block_begin
            segmentedData=segmentedData(:,:,isNotArtifact & idx_notar);
        elseif artrej==1
            segmentedData=segmentedData(:,:,isNotArtifact==1);
        end

        % If Block-Begin, search for only the first trial beginning of each block
        if block_begin
            if c==1 || c==2 % block_length of 52, three blocks total
                block_begin_idx=[1 53 105];
                block_vector=zeros(1,156)'; % 156 Trials total
                block_vector(block_begin_idx)=1;
            else % block_length of 60, four blocks total
                block_begin_idx=[1 61 121 181];
                block_vector=zeros(1,240)'; % 240 Trials total
                block_vector(block_begin_idx)=1;
            end


           % Select only Block_Begin Trials without Artifacts
           clean_firsttrials=isNotArtifact==1&block_vector==1;
           if sum(clean_firsttrials)==0
               disp('No artifact-free block beginnings')
           else
               fprintf('%i clean block beginnings\n',sum(clean_firsttrials))
           end
           block_begin_notartifact{c}=isNotArtifact(block_begin_idx);
           segmentedData=segmentedData(:,:,block_vector==1);
        end

    % TF Analysis
    tic
    parfor el=1:length(alpha_tfElectrodes)
        [wvlt_amp, ~] = morletwave(alpha_wavFreqs, 12, squeeze(segmentedData(:,alpha_tfElectrodes(el),:))', srate, 0, 'waitbar', 'on'); % time points x frequnencies x trials
        inducedMat=squeeze(mean(wvlt_amp,3)); % average over trials
        condResults(:,:,el)=inducedMat'; % time points x frequencies x electrodes
    end
    toc

    % Convert to power
    if power
        condResults=condResults^2;
    end

    % Remove Padding
    condResults=condResults(timeVec>=timerange(1) & timeVec<=timerange(2) ,:,:);
    timeVec=timeVec(timeVec>=timerange(1) & timeVec<=timerange(2));
    alpha_Results{c}=condResults; % Time Points, Frequencies, Electrodes
    alpha_timeVecTotal{c}=timeVec;
    alpha_ntrials{c}=size(segmentedData,3);
    clear condResults segmentedData timeVec
end
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\AlphaRes'
if ~block_begin
save(savefilename, 'alpha_Results', 'alpha_timeVecTotal', 'alpha_wavFreqs', 'alpha_tfElectrodes','alpha_ntrials',"power",'irrtartimes');
else
    save(savefilename, 'alpha_Results', 'alpha_timeVecTotal', 'alpha_wavFreqs', 'alpha_tfElectrodes','alpha_ntrials',"power",'irrtartimes','block_begin_notartifact');
end
disp('Alpha Amplitude Results Saved')
%% Plot TF
if plots
    % Baseline Correct?
    % basec=input("Baseline Correction (1-yes)? ");
    if ~catchonly
        baselinetp=[0 100]; % baseline time points in relation to WS (WS at time point 0)
    else
        baselinetp=[-800 -700]; % baseline time points in relation to Target (Target at 0)
    end


    for c=1:size(triggercodes,1) % For conditions

        if basec
            %     [~,startbaseidx]=(min(abs(timeVec-baselinetp(1))));
            %     [~,endbaseidx]=(min(abs(timeVec-baselinetp(2))));
            %     baseline=mean(tfResults(startbaseidx:endbaseidx,:,:), 1);
            %     datatoplot=tfResults-baseline;
            datatoplot=baselineCorrectSegmentedData(alpha_Results{c}, alpha_timeVecTotal{c}, baselinetp);
        else
            datatoplot=alpha_Results{c};
        end

        %minmaxscale= [min(mean(datatoplot,3),[],"all") max(mean(datatoplot,3),[],"all")];
        minmaxscale=[-3 3];

        figure;
        imagesc(alpha_timeVecTotal{c}, alpha_wavFreqs, squeeze(mean(datatoplot,3))',minmaxscale); axis xy; hold all
        if catchonly
            xline(-800,'r--','Warning Signal');
            xline(0,'b--','Predicted Target');
        elseif maskonly
            xline(0,'r--','Mask Onset');
        else
            xline(0,'r--','Warning Signal');
            xline(800,'b--','Predicted Target');
        end

        colorbar();
        xlabel('Time (ms)', 'FontWeight','bold', 'FontSize', 10);
        ylabel('Frequency (Hz)', 'FontWeight','bold', 'FontSize', 10);
        if c==1
            heading=sprintf('Rhythm - Pre-Target Alpha Amplitude (%i trials)',alpha_ntrials{c});
        elseif c==2
            heading=sprintf('Interval - Pre-Target Alpha Amplitude (%i trials)',alpha_ntrials{c});
        elseif c==3
            heading=sprintf('Irregular - Pre-Target Alpha Amplitude (%i trials)',alpha_ntrials{c});
        end
        title(heading)
    end

    figure;
    for c=1:size(triggercodes,1) % For conditions
        if basec
            %     [~,startbaseidx]=(min(abs(timeVec-baselinetp(1))));
            %     [~,endbaseidx]=(min(abs(timeVec-baselinetp(2))));
            %     baseline=mean(tfResults(startbaseidx:endbaseidx,:,:), 1);
            %     datatoplot=tfResults-baseline;
            datatoplot=baselineCorrectSegmentedData(alpha_Results{c}, alpha_timeVecTotal{c}, baselinetp);
        else
            datatoplot=alpha_Results{c};
        end
        % Plot Alpha Only (8-12 hz)
        plot(alpha_timeVecTotal{c},mean(mean(datatoplot(:,alpharange,:),3),2), 'LineWidth', 2)
        % alphaamp=tfResults(:,alpharange,:);
        % meanelec=mean(alphaamp,3);
        % varplot(mean(alphaamp,3))
        legend('Rhythm','Interval', 'Irregular')
        title("Alpha Power between WS and Target")
        hold on;
    end
    if catchonly
        xline(-800,'r--','Warning Signal');
        xline(0,'b--','Predicted Target');
    elseif maskonly
        xline(0,'r--','Mask Onset');
    else
        xline(0,'r--','Warning Signal');
        xline(800,'b--','Predicted Target');
    end
end
end