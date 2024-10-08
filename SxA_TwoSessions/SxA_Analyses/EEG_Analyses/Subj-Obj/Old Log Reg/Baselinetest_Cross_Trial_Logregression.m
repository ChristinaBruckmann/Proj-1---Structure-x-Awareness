%% Subj-Obj Power Cross Trial Regression Analysis
% This code is similar to Cross-Trial_Logregression, but calculates and saves several baseline and zscore options per participant for later comparison
clear
clc

% Parameters
subj=[101 102 103 105 114];
%freqrange=[1:40];
freqrange=[1:31];
electrodes=[25:30 62:64]; % occipital
timep=[250 1200]; % (data aligned to warning signal, target at 800), make sure the values you choose are within the time limits you chose for the TF analysis
cond=[1 2 3]; % which conditions (rhythm, interval, irregular)

% R or matlab?
matlab=0; 

% baselinecorr=1; % Correct for pre-trial power (needs to be extracted with the TF function first)
% subtract=1; % subtract (1) or divide(0) by baseline? (only aplied if baselinecorr=1)
% calczscore=1; % z-score contrasts  and power values?
%% Load Data
% Load TF Data
for s=subj
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\TimeFreq'
    loadfilename=sprintf('EEG_SxA_Subj%i_TF_SingleTrials',s);
    load(loadfilename,'TF_Results_Trial','TF_trial_timeVec','TF_NotArtifact')
    savefilename=sprintf('EEG_SxA_Subj%i_Results_SubjObj.mat',s);
    TF_res_WS=TF_Results_Trial; % Size: timepoints frequencies trials electrodes
    TF_time_WS=TF_trial_timeVec;
    TF_ArtVectors_WS=TF_NotArtifact;

    % Cycle through the different options and save all of them
    for baselinecorr=[0,1] %once corrected, once uncorrected
         for calczscore=[0,1] % one time z-scored, one time without
             for subtract= [0,1] % one time baseline subtracted once, divided
                if baselinecorr
                    clearvars TF_Results_Trial TF_trial_timeVec TF_NotArtifact

                    % Load Baseline Data
                    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\TimeFreq'
                    loadfilename=sprintf('EEG_SxA_Subj%i_TF_SingleTrials_BL',s);
                    load(loadfilename,'TF_Results_Trial','TF_trial_timeVec','TF_NotArtifact')
                    TF_res_BL=TF_Results_Trial; % Dimensions: frequencies x timepoints x trials x electrodes
                    TF_time_BL=TF_trial_timeVec;
                    TF_ArtVectors_BL=TF_NotArtifact;
                    clearvars TF_Results_Trial TF_trial_timeVec TF_NotArtifact

                    % Merge Artifact Vectors
                    for c=cond
                        TF_ArtVector{c}=TF_ArtVectors_WS{c} & TF_ArtVectors_BL{c}; % A trial needs to be artifact-free (1) in both BL and WS period. If one of these has an artifact (0) it will be rejected.
                    end

                else
                    TF_ArtVector=TF_ArtVectors_WS; % If we do not baseline correct, the overall TF artifact vector is equivalent to the artifact vector of the WS data.
                end

                % Convert to power and average across electrodes
                for c=cond % for each condition
                    if baselinecorr
                        % Select data
                        data=TF_res_WS{1, c};
                        base_data=mean(TF_res_BL{1, c},1); % calculate the mean power during baseline for each trial, frequency and electrode, average across time points
                        % Baseline correct
                        if subtract
                            data=data-base_data; % for each trial, subtract baseline mean from each time point
                        else % divide
                            data=data./base_data; % or, for each trial, divide amp at each time point by the baseline mean at that point
                        end
                    else
                        % Select Data
                        data=TF_res_WS{1, c};
                    end

                    timevec=TF_time_WS{1, c};
                    idx=timevec>=timep(1)&timevec<=timep(2); % Select time-points
                    data=mean(data(idx,freqrange,logical(TF_ArtVector{c}),electrodes),4); % average across electrodes (time, freq,trials,elecs), remove artifact trials at this point
                    fprintf('Percent Artifact Free Trials: %.4f \n',mean(TF_ArtVector{c}));
                    data_power{c}=squeeze(abs(data).^2); % Convert to Power
                end

                % Load Behavioural Data
                cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\Behavioural Preprocessed'
                loadfilename=sprintf('SxA_ResultsSubject%i_Session2',s);
                load(loadfilename,'alldataclean')
                behaviour=alldataclean;
                clear alldataclean

                % For each subject, create a table with contrast level, subj, obj for each trial (artifacts removed)
                for c=cond
                    logidx=logical(behaviour{:,'Condition'}==c);
                    obj_resp=behaviour{logidx,'Correct/Incorrect'};
                    subj_resp=behaviour{logidx,'Binary Visibility'};
                    contrast=behaviour{logidx,'Contrast Level'};
                    all_behav{c}=[contrast(logical(TF_ArtVector{c})) obj_resp(logical(TF_ArtVector{c}))  subj_resp(logical(TF_ArtVector{c}))]; % save trial list with artifacts removed

                    % Remove Catch Trials
                    catchidx=all_behav{c}(:,1)~=0; % select non-catch trials only
                    all_behav{c}=all_behav{c}(catchidx,:); % remove catch trials
                    data_power{c}=data_power{c}(:,:,catchidx);

                    % For irreular, remove trials with target before 800
                    if c==3
                        targettime=behaviour{logidx,'Irregular Target Time'};
                        targettime=targettime(logical(TF_ArtVector{c})); % remove artifacts
                        targettime=targettime(catchidx,:); % remove catch trials
                        targetidx=ismember(targettime,[3,4,5]); % find trials with target at 800ms or later
                        all_behav{c}=all_behav{c}(targetidx,:);
                        data_power{c}=data_power{c}(:,:,targetidx); % remove early target trials from data
                    end
                end

                % Z-Score
                if calczscore
                    for c=cond
                        data_power{c}=zscore(data_power{c},0,'all');
                        all_behav{1, c}(:,1)=zscore(all_behav{1, c}(:,1),0,'all');
                    end
                end

                cd  'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\SubjObj'

                % Choose file name
                if baselinecorr 
                    if subtract
                        if calczscore
                            savefilename2=sprintf("RInput_LogReg_S%i_bc_sb_z",s); % bc(baselinecorrect), sb(subtract), z(scored)
                        else
                            savefilename2=sprintf("RInput_LogReg_S%i_bc_sb",s); % bc(baselinecorrect), sb(subtract)
                        end
                    else
                        if calczscore
                            savefilename2=sprintf("RInput_LogReg_S%i_bc_dv_z",s); % bc(baselinecorrect), dv(divide), z(scored)
                        else
                            savefilename2=sprintf("RInput_LogReg_S%i_bc_dv",s); % bc(baselinecorrect), dv(divide)
                        end
                    end
                else
                    if subtract
                        if calczscore
                            savefilename2=sprintf("RInput_LogReg_S%i_uc_sb_z",s); % uc(uncorrected), sb(subtract), z(scored)
                        else
                            savefilename2=sprintf("RInput_LogReg_S%i_uc_sb",s); %  uc(uncorrected), sb(subtract)
                        end
                    else
                        if calczscore
                            savefilename2=sprintf("RInput_LogReg_S%i_uc_dv_z",s); %  uc(uncorrected), dv(divide), z(scored)
                        else
                            savefilename2=sprintf("RInput_LogReg_S%i_uc_dv",s); % uc(uncorrected), dv(divide)
                        end
                    end
                end
                save(savefilename2,'all_behav','data_power')
            end

        end
    end
end
%% %% Cross-Trial Regression
% if matlab
%     % For each condition
%     for c=cond
%         currobj=all_behav{c}(:,2);
%         currsubj=all_behav{c}(:,3);
%         contrasts=all_behav{c}(:,1);
%         % At each frequency
%         for freq= 1:size(data_power{c},1)
%             % At each time point
%             for tp=1:size(data_power{c},2)
%                 % Select relevant data
%                 currpow=squeeze(data_power{c}(freq,tp,:));
% 
%                 % run logistic regression
%                 [bob,~,statsob]=mnrfit([contrasts, currpow, currpow.*contrasts], currobj+1); % +1 needed as this function doesnt take 0 as category value
%                 [bsu,~,statssu]=mnrfit([contrasts, currpow, currpow.*contrasts], currsubj+1);
% 
%                 % Save weights and stats
%                 tp_bob(tp,:)=bob;
%                 tp_bsu(tp,:)=bsu;
%                 tp_statsob(tp,:)=statsob;
%                 tp_statssu(tp,:)=statssu;
%             end
%             % Save weights and stats
%             freq_bob(freq,:,:)=tp_bob;
%             freq_bsu(freq,:,:)=tp_bsu;
%             freq_statsob(freq,:,:)=tp_statsob;
%             freq_statssu(freq,:,:)=tp_statssu;
%         end
%         % Save weights and stats
%         total_bob(c,:,:,:)=freq_bob; % (condition, frequency, timepoint,variables)
%         total_bsu(c,:,:,:)=freq_bsu;
%         total_pob(c,:,:,:)=freq_statsob.p;
%         total_psu(c,:,:,:)=freq_statssu.p;
%     end
%     % Save
%     cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data Analysis\SxA_EEG_Analyses_Current\Results\SubjObj'
%     save(savefilename,'total_bob','total_bsu','total_statsob','total_statssu','cond','electrodes','freqrange','timep')
% else
%     % Export to R
%     disp("Data for R Prepared. Continue Analysis in R now.")
% end
% 
% %% Plot
% % For each condition
% subj=[18:22];
% plot_var=3; %intercept,contrast,power,interaction
% cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\Data Analysis\SxA_EEG_Analyses_Current\Results\SubjObj'
% for s=subj
%     % Load
%     loadfilename=sprintf('EEG_SxA_Subj%i_Results_SubjObj.mat',s);
%     load(loadfilename)
%     figure;
%     for c=1:2
% 
%         % Title
%         if c==1
%             obtitle='Objective Rhythm';
%             sutitle='Subjective Rhythm';
%         elseif c==2
%             obtitle='Objective Interval';
%             sutitle='Subjective Interval';
%         else
%             obtitle='Objective Irregular';
%             sutitle='Subjective Irregular';
%         end
% 
%         % Plot Power Weights
% 
%         subplot(2,2,c);title(obtitle);imagesc(flip(squeeze(total_bob(c,:,:,plot_var)))); colorbar;
%         caxis([-1 1]);
%         xline(800-timep(1));
%         subplot(2,2,c+2);title(sutitle);imagesc(flip(squeeze(total_bsu(c,:,:,plot_var)))); colorbar;
%         xline(800-timep(1));
%         caxis([-1 1]);
%     end
%     sgtitle(sprintf('Alpha Power Weights Subj %i',s))
% end
% 
% %         % Extract P values and average
% %         for freq=1:size(test,1)
% %             clear test10
% %             clear test11
% %             for tp=1:size(test,2)
% %                 test10(tp,:)=squeeze(tempstatsob(freq, tp).p);
% %                 test11(tp,:)=squeeze(tempstatssu(freq, tp).p);
% %             end
% %             p_values_ob(freq,:,:)=test10;
% %             p_values_su(freq,:,:)=test11;
% %         end
% %
% %         p_val_ob=squeeze(mean(p_values_ob,1));
% %         p_val_su=squeeze(mean(p_values_su,1));
% %
% %         % Plot p-values obj
% %         figure;
% %         for weight=1:size(p_val_ob,2) % for each weight
% %             subplot(1,4,weight)
% %             plot(1:size(p_val_ob,1),p_val_ob(:,weight))
% %             if weight==1
% %                 title('P-Value Intercept')
% %             elseif weight==2
% %                 title('P-Value Contrast')
% %             elseif weight==3
% %                 title('P-Value Alpha Power')
% %             elseif weight==4
% %                 title('P-Value Interaction Contrast-Alpha')
% %             end
% %         end