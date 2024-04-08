% Plot Single Subject Results
clear
clc

subj=[14 15];
%% CNV
for s=subj

    clear -except subj
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
    loadfilename=sprintf('EEG_SxA_Subj%i_Results.mat',s);

    load(loadfilename,'CNV_electodes','CNV_timeVec','CNV_perElec','CNV_Mean','CNV_allroiidx','CNV_channelnames')

    % Baseline correct
    basec=1;
    baselinetp=[0 100]; %baseline time points in relation to WS (WS at time point 0)

        if basec
            %     [~,startbaseidx]=(min(abs(timeVec-baselinetp(1))));
            %     [~,endbaseidx]=(min(abs(timeVec-baselinetp(2))));
            %     baseline=mean(tfResults(startbaseidx:endbaseidx,:,:), 1);
            %     datatoplot=tfResults-baseline;
            datatoplot=baselineCorrectSegmentedData(CNV_Mean, CNV_timeVec, baselinetp);
        else
            datatoplot=CNV_Mean;
        end


    % Plot
    figure;
    for elecidx=1:length(CNV_electodes)
        subplot(2,round(length(CNV_electodes)/2),elecidx); plot(CNV_timeVec,CNV_perElec(:,elecidx))
        title(CNV_channelnames{CNV_allroiidx(elecidx),1})
        x1=xline(0,'r--','Warning Signal');
        x2=xline(800,'b--','Predicted Target');
        x1.FontSize = 8;
        x2.FontSize = 8;
    end

    figure; plot(CNV_timeVec,datatoplot)
    title(sprintf('Average Activity from ROI channels %s', strcat(CNV_channelnames{CNV_allroiidx,1})))
    xline(0,'r--','Warning Signal');
    xline(800,'b--','Predicted Target');


end