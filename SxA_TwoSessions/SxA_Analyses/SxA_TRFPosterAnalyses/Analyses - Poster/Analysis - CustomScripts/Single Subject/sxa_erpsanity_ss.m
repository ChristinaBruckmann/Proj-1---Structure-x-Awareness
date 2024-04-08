function []=sxa_erpsanity_ss(subj, plots)
% Sanity Check: Calculate ERP to onset of noise background Single Subject
% Wihout need for outside functions. Attention, might give issues regarding
% artifact rejection and trial splitting

% clear
% clc
disp('Starting Sanity Check - ERP at Mask Onset')
% Load data
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG'
%subj=input("Subject Number? ");
loadfilename=sprintf('EEG_SxA_Subj%i_Session2_pp.mat',subj);
load(loadfilename)

% Parameters
onsettriggers=[31 32 33];
timerange=[-.2 .7]; % in seconds

% Extract Info and Data from File
data=SDATA.data(SDATA.metadata.artifacts==0,:);
triggers=SDATA.events.triggerChannel(SDATA.metadata.artifacts==0,:);
channelnames=SDATA.info.channel_labels;
srate=SDATA.info.sampling_rate;
timepoints=(0:length(data)-1)/srate; % time points in seconds

% Find trigger moments in EEG data
onsetidx=find(ismember(triggers,onsettriggers)); % find trigger idx

for n_onset=1:length(onsetidx)
currtime=timepoints(onsetidx(n_onset)); % find time point of current onset
minustime=currtime+timerange(1); % calculate time points of range
plustime=currtime+timerange(2);
[~,minusidx]=min(abs(timepoints-minustime)); % find closest idx to this time point
[~,plusidx]=min(abs(timepoints-plustime));
timeindices(n_onset,:)=[minusidx; onsetidx(n_onset); plusidx];
end

% Calculate ERP
ERPs=[];
for elecidx=1:width(data)
    for trial=1:length(onsetidx)
    evokedresp(:,trial)=data(timeindices(trial,1):timeindices(trial,3),elecidx);
    end
    ERPs(:,elecidx)=mean(evokedresp,2);
end

% Plot
if plots
    timeinms=timerange(1)*1000:(timerange(2)*1000+22);
    for elecidx=1:width(data)
        subplot(10,8,elecidx); plot(timeinms,ERPs(:,elecidx))
        title(channelnames{elecidx,1})
        xline(0,'--')
    end

    % Select region of interest (some occipital electrodes)
    ROIelecname=["o1","iz","oz"];

    roiidx=contains(channelnames,ROIelecname,'IgnoreCase',true);

    % Average across trials and across ROI electrodes
    allroiidx=find(roiidx==1);
    for elecidx=1:length(ROIelecname)
        currelec=allroiidx(elecidx);
        for trial=1:length(onsetidx)
            evokedresp(:,trial)=data(timeindices(trial,1):timeindices(trial,3),currelec);
        end
        ROIERPs(:,elecidx)=mean(evokedresp,2);
    end

    ROImeanERP=mean(ROIERPs,2);

    % Plot
    timeinms=timerange(1)*1000:(timerange(2)*1000+22);
    figure;
    for elecidx=1:length(ROIelecname)
        subplot(2,round(length(ROIelecname)/2),elecidx); plot(timeinms,ROIERPs(:,elecidx))
        title(channelnames{allroiidx(elecidx),1})
        xline(0,'--')
    end

    figure; plot(timeinms,ROImeanERP)
    title(sprintf('Average Activity from ROI channels %s', strcat(channelnames{allroiidx,1})))
    xline(0,'--')
end
end