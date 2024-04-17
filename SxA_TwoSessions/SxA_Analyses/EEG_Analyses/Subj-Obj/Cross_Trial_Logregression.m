%% Subj-Obj Power Cross Trial Regression Analysis
% For one participant, fit a logistic curve for each time point and each frequency within the alpha range
% single Trial TF data has to have alrady been calculated
clear
clc

% Parameters
subj=[101 102 103 105 106 107 108 111 113 114];
freqrange=[1:40];
electrodes=[25:30 62:64]; % occipital
timep=[250 1200]; % (data aligned to warning signal, target at 800)
cond=[1 2]; % which conditions (rhythm, interval, irregular)

% R or matlab?
matlab=0; 
%% Load Data

% Load TF Data
for s=subj
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\TimeFreq'
    loadfilename=sprintf('EEG_SxA_Subj%i_TF_SingleTrials',s);
    load(loadfilename,'TF_Results_Trial','TF_trial_timeVec','TF_NotArtifact')
    savefilename=sprintf('EEG_SxA_Subj%i_Results_SubjObj.mat',s);
    TF_res=TF_Results_Trial; % Size: frequencies timepoints trials electrodes
    TF_time=TF_trial_timeVec;
    TF_ArtVectors=TF_NotArtifact;

    % Convert to power and average across electrodes
    for c=cond % for each condition (only predictable for now)
        % Select Data
        data=TF_res{1, c};
        timevec=TF_time{1, c};
        idx=timevec>=timep(1)&timevec<=timep(2);
        data=mean(data(idx,freqrange,:,electrodes),4); % average across electrodes (time, freq,trials,elecs)
        % Convert to Power
        data_power{c}=squeeze(abs(data).^2);
    end

    % Load Behavioural Data
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\Behavioural Preprocessed'
    loadfilename=sprintf('SxA_ResultsSubject%i_Session2',s);
    load(loadfilename,'alldataclean')
    behaviour=alldataclean;
    clear alldataclean

    % For each subject, create a table with contrast level, subj, obj for each
    % trial (artifacts removed)
    for c=cond
        logidx=logical(behaviour{:,'Condition'}==c);
        obj_resp=behaviour{logidx,'Correct/Incorrect'};
        subj_resp=behaviour{logidx,'Binary Visibility'};
        contrast=behaviour{logidx,'Contrast Level'};
        all_behav{c}=[contrast(logical(TF_ArtVectors{c})) obj_resp(logical(TF_ArtVectors{c}))  subj_resp(logical(TF_ArtVectors{c}))]; % save trial list with artifacts removed

        % Remove Catch Trials
        catchidx=all_behav{c}(:,1)~=0; % select non-catch trials only
        all_behav{c}=all_behav{c}(catchidx,:); % remove catch trials
        data_power{c}=data_power{c}(:,:,catchidx);

        % For irreular, remove trials with target before 800
        if c==3
            targettime=behaviour{logidx,'Irregular Target Time'};
            targettime=targettime(logical(TF_ArtVectors{c})); % remove artifacts
            targettime=targettime(catchidx,:); % remove catch trials
            targetidx=ismember(targettime,[3,4,5]); % find trials with target at 800ms or later
            all_behav{c}=all_behav{c}(targetidx,:);
            data_power{c}=data_power{c}(:,:,targetidx); % remove early target trials from data
        end
    end

    cd  'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\SubjObj'
    savefilename2=sprintf("RInput_LogReg_S%i",s);
    save(savefilename2,'all_behav','data_power')
end
    %% Cross-Trial Regression
    if matlab
        % For each condition
        for c=cond
            currobj=all_behav{c}(:,2);
            currsubj=all_behav{c}(:,3);
            contrasts=all_behav{c}(:,1);
            % At each frequency
            for freq= 1:size(data_power{c},1)
                % At each time point
                for tp=1:size(data_power{c},2)
                    % Select relevant data
                    currpow=squeeze(data_power{c}(freq,tp,:));

                    % run logistic regression
                    [bob,~,statsob]=mnrfit([contrasts, currpow, currpow.*contrasts], currobj+1); % +1 needed as this function doesnt take 0 as category value
                    [bsu,~,statssu]=mnrfit([contrasts, currpow, currpow.*contrasts], currsubj+1);

                    % Save weights and stats
                    tp_bob(tp,:)=bob;
                    tp_bsu(tp,:)=bsu;
                    tp_statsob(tp,:)=statsob;
                    tp_statssu(tp,:)=statssu;
                end
                % Save weights and stats
                freq_bob(freq,:,:)=tp_bob;
                freq_bsu(freq,:,:)=tp_bsu;
                freq_statsob(freq,:,:)=tp_statsob;
                freq_statssu(freq,:,:)=tp_statssu;
            end
            % Save weights and stats
            total_bob(c,:,:,:)=freq_bob; % (condition, frequency, timepoint,variables)
            total_bsu(c,:,:,:)=freq_bsu;
            total_pob(c,:,:,:)=freq_statsob.p;
            total_psu(c,:,:,:)=freq_statssu.p;
        end
        % Save
        cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data Analysis\SxA_EEG_Analyses_Current\Results\SubjObj'
        save(savefilename,'total_bob','total_bsu','total_statsob','total_statssu','cond','electrodes','freqrange','timep')
    else

        % Export to R

        % Run R

        % Save
        cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\SxA_EEG_Analyses_Current\Results\SubjObj'
        save(savefilename,'total_bob','total_bsu','total_statsob','total_statssu','cond','electrodes','freqrange','timep')
    end
end

%% Plot
% For each condition
subj=[18:22];
plot_var=3; %intercept,contrast,power,interaction
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data Analysis\SxA_EEG_Analyses_Current\Results\SubjObj'
for s=subj
    % Load
    loadfilename=sprintf('EEG_SxA_Subj%i_Results_SubjObj.mat',s);
    load(loadfilename)
    figure;
    for c=1:2

        % Title
        if c==1
            obtitle='Objective Rhythm';
            sutitle='Subjective Rhythm';
        elseif c==2
            obtitle='Objective Interval';
            sutitle='Subjective Interval';
        else
            obtitle='Objective Irregular';
            sutitle='Subjective Irregular';
        end

        % Plot Power Weights

        subplot(2,2,c);title(obtitle);imagesc(flip(squeeze(total_bob(c,:,:,plot_var)))); colorbar;
        caxis([-1 1]);
        xline(800-timep(1));
        subplot(2,2,c+2);title(sutitle);imagesc(flip(squeeze(total_bsu(c,:,:,plot_var)))); colorbar;
        xline(800-timep(1));
        caxis([-1 1]);
    end
    sgtitle(sprintf('Alpha Power Weights Subj %i',s))
end

%         % Extract P values and average
%         for freq=1:size(test,1)
%             clear test10
%             clear test11
%             for tp=1:size(test,2)
%                 test10(tp,:)=squeeze(tempstatsob(freq, tp).p);
%                 test11(tp,:)=squeeze(tempstatssu(freq, tp).p);
%             end
%             p_values_ob(freq,:,:)=test10;
%             p_values_su(freq,:,:)=test11;
%         end
%
%         p_val_ob=squeeze(mean(p_values_ob,1));
%         p_val_su=squeeze(mean(p_values_su,1));
%
%         % Plot p-values obj
%         figure;
%         for weight=1:size(p_val_ob,2) % for each weight
%             subplot(1,4,weight)
%             plot(1:size(p_val_ob,1),p_val_ob(:,weight))
%             if weight==1
%                 title('P-Value Intercept')
%             elseif weight==2
%                 title('P-Value Contrast')
%             elseif weight==3
%                 title('P-Value Alpha Power')
%             elseif weight==4
%                 title('P-Value Interaction Contrast-Alpha')
%             end
%         end