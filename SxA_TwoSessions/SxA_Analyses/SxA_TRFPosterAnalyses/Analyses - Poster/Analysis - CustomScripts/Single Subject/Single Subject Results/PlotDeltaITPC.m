% Plot Single Subject Results
clear
clc

subj=[14 15 17:22];

%% Delta Phase
for s=subj

    clear -except subj
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
    loadfilename=sprintf('EEG_SxA_Subj%i_Results.mat',s);

    load(loadfilename,'ITPCdelta_res','ITPCdelta_timevec','ITPCdelta_elec')

    % check ITPC time course and topography
    timeROI = [600 750];
    baseRange = [0 100];
    occelonly = 1;
    if occelonly
        dispElectrodes = [25:30 62:64]; % Occipital only
    else
        dispElectrodes = ITPCdelta_elec;
    end

    % mean time course across occipital electrodes
    figure;
    for c=1:size(ITPCdelta_res,2) % For conditions
        subplot(1,3,c)
        plot(ITPCdelta_timevec{c}, squeeze(mean(ITPCdelta_res{c}(:,dispElectrodes), 2)));
        ylim([0 0.5])
        line([0 0],ylim,'LineStyle', '-', 'Color', [0 0 0])
        line([700 700],ylim,'LineStyle', '- -', 'Color', [0 0 0])
        if c==1
            title('ITPC Rhythm')
        elseif c==2
            title('ITPC Interval')
        elseif c==3
            title('ITPC Irregular')
        end
    end

    figure;
    for c=1:size(ITPCdelta_res,2) % For conditions
        % topography (in relation to baseline)
        subplot(1,3,c)
        bp_phase= baselineCorrectSegmentedData(ITPCdelta_res{c}, ITPCdelta_timevec{c}, baseRange);
        topoplot(squeeze(mean(bp_phase(ITPCdelta_timevec{c}>=timeROI(1) & ITPCdelta_timevec{c}<=timeROI(2), :))), 'head71.locs');
        colorbar
        if c==1
            title('Rhythm')
        elseif c==2
            title('Interval')
        elseif c==3
            title('Irregular')
        end
    end
end
