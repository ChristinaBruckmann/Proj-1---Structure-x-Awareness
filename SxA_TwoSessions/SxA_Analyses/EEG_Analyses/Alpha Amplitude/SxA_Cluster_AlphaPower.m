%% SxA Alpha Power Cluster Based Permutation (Group Level)
% Currently two tailed, for one tailed ask Assaf if you can interfere with his code

subj=[101:103 105:108 110:114 116:119 121 122 124 126 127 129:132];
catchonly=0; % only catch trials?
freq=8:12;

%% Load Data
for s=1:length(subj)
    cd 'Y:\el-Christina\SxA\SxA_Results\AlphaPowerRes'
    if catchonly
        loadfilename=sprintf('EEG_SxA_Subj%i_AlphaResults_Catch.mat',subj(s));
    else
        loadfilename=sprintf('EEG_SxA_Subj%i_AlphaResults_clean.mat',subj(s));
    end
    gl_tf_timeVec(s)=load(loadfilename,'alpha_timeVecTotal');
    gl_alpha_wavFreqs(s)=load(loadfilename,'alpha_wavFreqs');
    subj_res=load(loadfilename,'alpha_Results'); 
    ntrials(s)=load(loadfilename,'alpha_ntrials');

    % Save as array and separate by condition

    Data_Rhy(s,:,:,:)=cell2mat(subj_res.alpha_Results(1)); % subject x condition x time point x freq x electrodes
    Data_Int(s,:,:,:)=cell2mat(subj_res.alpha_Results(2)); % subject x condition x time point x freq x electrodes
    Data_Irr(s,:,:,:)=cell2mat(subj_res.alpha_Results(3)); % subject x condition x time point x freq x electrodes
end

% Average across electrodes and frequencies
Data_Rhy=squeeze(mean(mean(Data_Rhy(:,:,freq,:),4),3));
Data_Int=squeeze(mean(mean(Data_Int(:,:,freq,:),4),3));
Data_Irr=squeeze(mean(mean(Data_Irr(:,:,freq,:),4),3));

% Rhy - Irr
[Rhy_Clusters]=clusterBasedPermTest(Data_Rhy, Data_Irr, 1);

% Int - Irr
[Int_Clusters]=clusterBasedPermTest(Data_Int, Data_Irr, 1) ;

% Rhy - Int
[RhyInt_Clusters]=clusterBasedPermTest(Data_Rhy, Data_Int, 1);


%% Plot
timeVec=gl_tf_timeVec(1).alpha_timeVecTotal{1, 1};

% Plot Alpha
figure("Position",[548.2000  257.8000  758.4000  580.8000]);
plot(timeVec,mean(Data_Rhy,1),'LineWidth',2)
hold on
plot(timeVec,mean(Data_Int,1),'LineWidth',2)
plot(timeVec,mean(Data_Irr,1),'LineWidth',2)

ylimits = ylim;

% Add Clusters
for cl=1:size(Rhy_Clusters,1) % Rhythm-Irregular
    startcl=timeVec(Rhy_Clusters(cl,1)); % Get start point of cluster
    endcl= startcl+Rhy_Clusters(cl,2);% Get end point of cluster
    line([startcl endcl], [ylimits(1)+0.07 ylimits(1)+0.07], 'Color',[0 0.4470 0.7410], 'LineWidth', 2)
end

for cl=1:size(Int_Clusters,1) % Interval - Irregular
    startcl=timeVec(Int_Clusters(cl,1)); % Get start point of cluster
    endcl= startcl+Int_Clusters(cl,2);% Get end point of cluster
    line([startcl endcl], [ylimits(1)+0.06 ylimits(1)+0.06], 'Color',[0.8500 0.3250 0.0980], 'LineWidth', 2)
end

for cl=1:size(RhyInt_Clusters,1) %  Rhythm - Interval
    startcl=timeVec(RhyInt_Clusters(cl,1)); % Get start point of cluster
    endcl= startcl+RhyInt_Clusters(cl,2);% Get end point of cluster
    line([startcl endcl], [ylimits(1)+0.05 ylimits(1)+0.05], 'Color',[0 0 0], 'LineWidth', 2)
end


% Make nice
xline(0,'LineStyle','--','LabelOrientation','horizontal','LineWidth',1,'FontSize',15,'Label','Warning Signal');
xline(900,'LineStyle','--','LabelOrientation','horizontal','LineWidth',1,'FontSize',15,'Label','Predicted Target');
title('Alpha Power between WS and Target');
xlabel("time(ms)")
ylabel("band Limited Amplitude (µs)")
legend1 = legend(['Rhythm', 'Interval',"Irregular","",""]);
set(legend1,'Position',[0.660288949743352 0.582990712020913 0.158190580302694 0.132720591180465]);