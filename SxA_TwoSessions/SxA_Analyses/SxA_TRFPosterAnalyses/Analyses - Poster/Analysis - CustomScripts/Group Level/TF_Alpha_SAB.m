%% Time Frequency and Alpha Amplitude Group level for SAB Report
clear
clc

subj=[14 15];

% Load Data
for s=1:length(subj)
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
loadfilename=sprintf('EEG_SxA_Subj%i_Results_forfigure.mat',subj(s));
gl_tf_timeVec(s)=load(loadfilename,'alpha_timeVecTotal');
gl_alpha_wavFreqs(s)=load(loadfilename,'alpha_wavFreqs');
tf_results(s)=load(loadfilename,'alpha_Results');
end

% Average and Save
    for s=1:length(subj) % for each subject
        gl_tf_results(s,:,:,:)=tf_results(s).alpha_Results;
    end
    gl_tf_results_means(:,:,:)=mean(gl_tf_results(:,:,:,:),1);

% Plot TF 
figure;



    %minmaxscale= [min(mean(datatoplot,3),[],"all") max(mean(datatoplot,3),[],"all")];
    %totalminmax=[-3 3];

    imagesc(gl_tf_timeVec(1).alpha_timeVecTotal, gl_alpha_wavFreqs(1).alpha_wavFreqs, squeeze(mean(gl_tf_results_means(:,:,:),3))'); axis xy; hold all
    xline(0,'r--','Warning Signal');
    xline(800,'g--','Predicted Target');
    colorbar();
    xlabel('Time (ms)', 'FontWeight','bold', 'FontSize', 10);
    ylabel('Frequency (Hz)', 'FontWeight','bold', 'FontSize', 10);
    title('Irregular - Pre-Target Frequency Amplitudes');



% Plot Alpha Only
alpharange=8:12;

% Baseline Correct
baselinetp=[0 100]; %baseline time points in relation to WS (WS at time point 0)
for c=1:height(gl_tf_results_means)
    datatoplot(c,:,:,:)=baselineCorrectSegmentedData(squeeze(gl_tf_results_means(c,:,:,:)), gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, baselinetp);
end

figure;
for c=1:height(gl_tf_results_means) % For conditions
    % Plot Alpha Only (8-12 hz)
    plot(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, mean(mean(datatoplot(c,:,alpharange,:),4),3), 'LineWidth', 2);
    %plot(alpha_timeVecTotal{c},mean(mean(datatoplot{c}(:,alpharange,:),3),2), 'LineWidth', 2)
    % alphaamp=tfResults(:,alpharange,:);
    % meanelec=mean(alphaamp,3);
    % varplot(mean(alphaamp,3))
    legend('Rhythm','Interval', 'Irregular')
    title("Alpha Power between WS and Target")
    hold on;
end
xline(0,'r--','Warning Signal');
xline(800,'b--','Predicted Target');


% Variance Plot
figure;
for c=1:3
varplot(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1},squeeze(mean(datatoplot(c,:,alpharange,:),4)),'linewidth',2)
hold on
end
xline(0,'r--','Warning Signal');
xline(800,'b--','Predicted Target');