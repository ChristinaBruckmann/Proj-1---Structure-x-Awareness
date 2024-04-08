%% Plot Delta ITPC Time Course
% Plots Delta time course for both catch trials only, all trials, occipital
% and central clusters with the possibility to split subjects into batches
clear
clc

batches=[2]; % 3 is combined
subj1=[17:22];
subj2=[101:103 105:108 113 114];
subj3=[17:22 101:103 105:108 113 114];
clusters=[1 2]; % (1-occipital, 2- central)
trialtypes= [1 2]; % (1-all trials, 2 - catch trials)
plotindividuals=[0]; % plot also individual participants
%% Dir
cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\SxA_TwoSessions\SxA_Data\EEG Results\DeltaRes'

%% Delta Average
for cbatch=1:length(batches)

    figure; t = tiledlayout('flow');
    % Define (change this, this is buggy for the wrong input)
    if batches(cbatch)==1
        subj=subj1;
    elseif batches(cbatch)==2
        subj=subj2;
    elseif batches(cbatch)==3
        subj=subj3;
    end

    for ccluster=clusters
        % Define
        if ccluster==1
            elec=[25:30 62:64]; % occipital electrodes
        else
            elec=[11:12 46:49]; % central electrodes
        end

        for trialt=trialtypes
            % Define
            if trialt==1
                catchonly=0;
            else
                catchonly=1;
            end
            % Load
            for s=1:length(subj)
                if ccluster==1 && ~catchonly
                    loadfilename=sprintf('OccipitalAll_Subj%i',subj(s));
                    load(loadfilename,"Delta_SingleTrials_occ","ITPCdelta_timevec_occ","ITPCdelta_nTrials_occ") %(time points x electrodes x trials)
                    for c=1:3
                        Delta_SingleTrials(s,c,:,:)=circ_r(Delta_SingleTrials_occ{c},[], [], 3);
                        Delta_TimeVec(s,c,:)=ITPCdelta_timevec_occ{c};
                        Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_occ(c);
                    end
                    subtitle='All Trials Occipital';
                elseif ccluster==1 && catchonly
                    loadfilename=sprintf('OccipitalCatch_Subj%i',subj(s));
                    load(loadfilename,"Delta_SingleTrials_occ_catch","ITPCdelta_timevec_occ_catch","ITPCdelta_nTrials_occ_catch") %(time points x electrodes x trials)
                    for c=1:3
                        Delta_SingleTrials(s,c,:,:)=circ_r(Delta_SingleTrials_occ_catch{c},[], [], 3);
                        Delta_TimeVec(s,c,:)=ITPCdelta_timevec_occ_catch{c};
                        Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_occ_catch(c);
                    end
                    subtitle='Catch Occipital';
                elseif ccluster==2 && ~catchonly
                    loadfilename=sprintf('CentralAll_Subj%i',subj(s));
                    load(loadfilename,"Delta_SingleTrials_cen","ITPCdelta_timevec_cen","ITPCdelta_nTrials_cen") %(time points x electrodes x trials)
                    for c=1:3
                        Delta_SingleTrials(s,c,:,:)=circ_r(Delta_SingleTrials_cen{c},[], [], 3);
                        Delta_TimeVec(s,c,:)=ITPCdelta_timevec_cen{c};
                        Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_cen(c);
                    end
                    subtitle='All Trials Central';
                elseif ccluster==2 && catchonly
                    loadfilename=sprintf('CentralCatch_Subj%i',subj(s));
                    load(loadfilename,"Delta_SingleTrials_cen_catch","ITPCdelta_timevec_cen_catch","ITPCdelta_nTrials_cen_catch") %(time points x electrodes x trials)
                    for c=1:3
                        Delta_SingleTrials(s,c,:,:)=circ_r(Delta_SingleTrials_cen_catch{c},[], [], 3);
                        Delta_TimeVec(s,c,:)=ITPCdelta_timevec_cen_catch{c};
                        Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_cen_catch(c);
                    end
                    subtitle='Catch Central';
                end
            end

            % Average
            timeVec=squeeze(Delta_TimeVec(1,1,:)); % All TV are identical
            nTrials=mean(Delta_Ntrials,1); % nTrials for chance ITPC calculation

            for c=1:3
                SubjMean=squeeze(mean(Delta_SingleTrials(:,c,:,elec),1));
                Delta(c,:)=squeeze(mean(SubjMean,2)); % mean across subjects and electrodes

                % Total N of trials
                totaltrials(c)=round(mean(nTrials(c),1));

                % Chance angles
                angles=2*pi*rand(totaltrials(c),10000);
                angle_mean(c)=mean(circ_r(angles));
                angle_std(c)=std(circ_r(angles));

                % Average for each participant across trials and participants
                for s=1:length(subj)
                PartMean(s,c,:)=squeeze(mean(Delta_SingleTrials(s,c,:,elec),4));
                end
            end


            % Plot
            nexttile
            for c=1:3

                if c==1
                    colourvec1=[0.00,0.45,0.74];
                    colourvec2=[0.75,0.87,0.96];
                elseif c==2
                    colourvec1=[0.85,0.33,0.10];
                    colourvec2=[0.95,0.65,0.71];
                elseif c==3
                    colourvec1=[0.93,0.69,0.13];
                    colourvec2=[1.00,0.93,0.75];
                end

                plot(timeVec,Delta(c,:),"LineWidth",2,'Color',colourvec1)
                ylim([0 0.6])
                hold on

                if plotindividuals
                    for s=1:length(subj)
                        plot(timeVec,squeeze(PartMean(s,c,:)),"LineWidth",1,'Color',colourvec2)
                    end
                end
                yline(angle_mean(c))
            end

            if catchonly
                xlim([-1500 400])
                xticks([-1200:400:400])
                xline(-800)
                title('Catch Trials')
            else
                xlim([-700 1200])
                xticks([-400:400:1200])
                xticklabels({'-1200','-800','-400','0','400','800'})
                xline(800)
                title('All Trials')
            end
            xline(0)
            yticks([0:0.2:0.6])
            box off
            title(subtitle)
            if cbatch==3
                headline=sprintf('Delta ITPC All Subjects');
            else
                headline=sprintf('Delta ITPC Batch %i',cbatch);
            end
            title(t,headline)
            clearvars -except t subj elec subj1 subj2 subj3 clusters catchonly trialtypes cbatch trialt ccluster batches plotindividuals
        end
    end
end
