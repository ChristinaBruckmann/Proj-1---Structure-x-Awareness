% Plot Single Subject Results
clear
clc

subj=[22];
%% TF and Alpha Amplitude
% load
for s=subj
clear -except subj
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
loadfilename=sprintf('EEG_SxA_Subj%i_Results.mat',s);
load(loadfilename,'alpha_timeVecTotal','alpha_wavFreqs','alpha_Results')

% Baseline correct
basec=0;
baselinetp=[0 100]; %baseline time points in relation to WS (WS at time point 0)
for c=1:width(alpha_Results)
    if basec
        %     [~,startbaseidx]=(min(abs(timeVec-baselinetp(1))));
        %     [~,endbaseidx]=(min(abs(timeVec-baselinetp(2))));
        %     baseline=mean(tfResults(startbaseidx:endbaseidx,:,:), 1);
        %     datatoplot=tfResults-baseline;
        datatoplot{c}=baselineCorrectSegmentedData(alpha_Results{c}, alpha_timeVecTotal{c}, baselinetp);
    else
        datatoplot{c}=alpha_Results{c};
    end
alpharange=8:12;
minval(c)=min(mean(datatoplot{c},3),[],"all");
maxval(c)=max(mean(datatoplot{c},3),[],"all");
end

% Plot TF 
figure;
totalminmax=[min(minval) max(maxval)];

for c=1:width(alpha_Results) % For conditions


    %minmaxscale= [min(mean(datatoplot,3),[],"all") max(mean(datatoplot,3),[],"all")];
    %totalminmax=[-3 3];

    subplot(3,1,c);imagesc(alpha_timeVecTotal{c}, alpha_wavFreqs, squeeze(mean(datatoplot{c},3))',totalminmax); axis xy; hold all
    xline(0,'r--','Warning Signal');
    xline(800,'b--','Predicted Target');
    colorbar();
    xlabel('Time (ms)', 'FontWeight','bold', 'FontSize', 10);
    ylabel('Frequency (Hz)', 'FontWeight','bold', 'FontSize', 10);
    if c==1
        title('Rhythm - Pre-Target Frequency Amplitudes');
    elseif c==2
        title('Interval - Pre-Target Frequency Amplitudes');
    elseif c==3
        title('Irregular - Pre-Target Frequency Amplitudes');
    end
end

figure;
for c=1:width(alpha_Results) % For conditions
    % Plot Alpha Only (8-12 hz)
    plot(alpha_timeVecTotal{c},mean(mean(datatoplot{c}(:,alpharange,:),3),2), 'LineWidth', 2)
    % alphaamp=tfResults(:,alpharange,:);
    % meanelec=mean(alphaamp,3);
    % varplot(mean(alphaamp,3))
    legend('Rhythm','Interval', 'Irregular')
    title("Alpha Power between WS and Target")
    hold on;
end
xline(0,'r--','Warning Signal');
xline(800,'b--','Predicted Target');

% Alpha Variance Plot
figure;
for c=1:3
varplot(alpha_timeVecTotal{c},squeeze(mean(datatoplot{c}(:,alpharange,:),3)),'linewidth',2)
hold on
end
xline(0,'r--','Warning Signal');
xline(800,'b--','Predicted Target');
end