%% Time Frequency Group level
clear
clc

catchonly=0;
basec=1; % baseline correction?
%subj=[14 15 17:22];
subj=[101:103 105 106:108 111 113 114 117 118 119];
%subj=[17:22 101:103 105 106:108 113 114];

% Load Data
for s=1:length(subj)
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\AlphaRes'
    if catchonly
        loadfilename=sprintf('EEG_SxA_Subj%i_AlphaResults_Catch.mat',subj(s));
    else
        loadfilename=sprintf('EEG_SxA_Subj%i_AlphaResults.mat',subj(s));
    end
    gl_tf_timeVec(s)=load(loadfilename,'alpha_timeVecTotal');
    gl_alpha_wavFreqs(s)=load(loadfilename,'alpha_wavFreqs');
    tf_results(s)=load(loadfilename,'alpha_Results');
    ntrials(s)=load(loadfilename,'alpha_ntrials');
end

% Average and Save
for c=1:3 % for each condition
    for s=1:length(subj) % for each subject
        gl_tf_results(s,:,:,:)=tf_results(s).alpha_Results{c};
        alphatrials(s,c)=ntrials(s).alpha_ntrials{c};
    end
    gl_tf_results_means(c,:,:,:)=mean(gl_tf_results(:,:,:,:),1);
end

% Plot TF 
figure;

for c=1:height(gl_tf_results_means) % For conditions

    %minmaxscale= [min(mean(datatoplot,3),[],"all") max(mean(datatoplot,3),[],"all")];
    %totalminmax=[-3 3];

    subplot(3,1,c);imagesc(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, gl_alpha_wavFreqs(1).alpha_wavFreqs, squeeze(mean(gl_tf_results_means(c,:,:,:),4))', [1 6]); axis xy; hold all
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

%% Plot Mean TF of Predictable
average_Tf_pred=squeeze(mean(squeeze(mean(gl_tf_results_means(1:2,:,:,:),4)),1));
if basec
    if ~catchonly
        baselinetp=[0 100]; % baseline time points in relation to WS (WS at time point 0)
    else
        baselinetp=[-800 -700]; % baseline time points in relation to Target (Target at 0)
    end
    average_Tf_pred=baselineCorrectSegmentedData(average_Tf_pred, gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, baselinetp);
end
figure;
imagesc(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, gl_alpha_wavFreqs(1).alpha_wavFreqs, average_Tf_pred'); axis xy; hold all
xline(800,'Alpha',1,'Color',[0 0 0],'FontWeight','bold','LineStyle','--','LineWidth',1,'FontSize',15,'Label','Predicted Target');
xline(0,'Alpha',1,'Color',[0 0 0],'FontWeight','bold','LineStyle','--','LineWidth',1,'FontSize',15,'Label','Warning Signal');
colorbar();
ylabel('Frequency (Hz)','FontWeight','bold');
xlabel('Time (ms)','FontWeight','bold');
fig=gca;
set(fig,'FontSize',15,'Layer','top','XTick',[0 400 800 1200],'YTick',[10 20 30 40]);
%% Plot each conditions without alpha
for c=1:3
TF_noalpha(c,:,:)=squeeze(mean(gl_tf_results_means(c,:,:,:),4));
end

TF_noalpha(:,:,7:12)=NaN;

for c=1:3
imAlpha=ones(size(squeeze(TF_noalpha(c,:,:))));
imAlpha(isnan(squeeze(TF_noalpha(c,:,:))))=0;
imagesc(squeeze(TF_noalpha(c,:,:))','AlphaData',imAlpha');
set(gca,'color',0*[1 1 1]);
end

figure;
for c=1:height(TF_noalpha) % For conditions

    %minmaxscale= [min(mean(datatoplot,3),[],"all") max(mean(datatoplot,3),[],"all")];
    %totalminmax=[-3 3];

    subplot(3,1,c);imagesc(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, gl_alpha_wavFreqs(1).alpha_wavFreqs, squeeze(TF_noalpha(c,:,:))','AlphaData',imAlpha'); axis xy; hold all
    xline(0,'r--','Warning Signal');
    xline(800,'b--','Predicted Target');
    ax = gca;
    ax.CLim = [1 3];
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

%% Plot Differences between Conditions
for c=1:3
Allconditions_TF(c,:,:)=squeeze(mean(gl_tf_results_means(c,:,:,:),4));
end

% Baseline Correct
if basec
    if ~catchonly
        baselinetp=[0 100]; % baseline time points in relation to WS (WS at time point 0)
    else
        baselinetp=[-800 -700]; % baseline time points in relation to Target (Target at 0)
    end

    for c=1:height(Allconditions_TF)
        Allconditions_TF(c,:,:)=baselineCorrectSegmentedData(squeeze(Allconditions_TF(c,:,:)), gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, baselinetp);
    end
else
end

% Rhythm - Irregular
rhy_min_irr=Allconditions_TF(1,:,:)-Allconditions_TF(3,:,:);
% Interval - Irregular
int_min_irr=Allconditions_TF(2,:,:)-Allconditions_TF(3,:,:);
% Rhythm - Irregular
rhy_min_int=Allconditions_TF(1,:,:)-Allconditions_TF(2,:,:);

figure;
subplot(3,1,1);imagesc(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, gl_alpha_wavFreqs(1).alpha_wavFreqs, squeeze(rhy_min_irr)'); axis xy; hold all
title('Difference Rhythm - Irregular');
xline(0,'r--','Warning Signal');
xline(800,'b--','Predicted Target');
ax = gca;
ax.CLim = [-0.5 0.2];
colorbar();
xlabel('Time (ms)', 'FontWeight','bold', 'FontSize', 10);
ylabel('Frequency (Hz)', 'FontWeight','bold', 'FontSize', 10);
subplot(3,1,2);imagesc(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, gl_alpha_wavFreqs(1).alpha_wavFreqs, squeeze(int_min_irr)'); axis xy; hold all
title('Difference Interval - Irregular');
xline(0,'r--','Warning Signal');
xline(800,'b--','Predicted Target');
ax = gca;
ax.CLim = [-0.5 0.2];
colorbar();
xlabel('Time (ms)', 'FontWeight','bold', 'FontSize', 10);
ylabel('Frequency (Hz)', 'FontWeight','bold', 'FontSize', 10);
subplot(3,1,3);imagesc(gl_tf_timeVec(1).alpha_timeVecTotal{1, 1}, gl_alpha_wavFreqs(1).alpha_wavFreqs, squeeze(rhy_min_int)'); axis xy; hold all
title('Difference Rhythm - Interval');
xline(0,'r--','Warning Signal');
xline(800,'b--','Predicted Target');
ax = gca;
ax.CLim = [-0.5 0.2];
colorbar();
xlabel('Time (ms)', 'FontWeight','bold', 'FontSize', 10);
ylabel('Frequency (Hz)', 'FontWeight','bold', 'FontSize', 10);


