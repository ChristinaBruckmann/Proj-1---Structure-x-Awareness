%% New Plot Delta

%% Plot Delta ITPC Time Course
% Plots Delta time course for both catch trials only, all trials, occipital
% and central clusters with the possibility to split subjects into batches
clear
clc

avg_done=0; % Average already computed (saves time) or compute here?
%subj=[101:103 105:108 110 112:114 116:119 121 122 123 124 126 127 129 130 131 132];
subj=[102   103   105   106   107   108   110   112   113   114   117   118   119  121   122   124   126   127   129   130];
clusters=[1]; % (1-occipital, 2- central)
trialtypes= [1]; % (1-all trials, 2 - catch trials)
plotindividuals=[0]; % plot also individual participants
varplotting=1; % add error bars to plot?

% Time Window for Statistics
stats_tw=[800 900]; %from WS
shade=0; % shade the statistics window?
clust_perm=1; % plot cluster based perm results?

%% Dir
cd ' Y:\el-Christina\SxA\SxA_Results\New Delta Results'

%% Delta Average
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

        % Load Data
        if ~avg_done % If average has already been computed and saved (saves time!)
            for s=1:length(subj)
                if ccluster==1 && ~catchonly
                    loadfilename=sprintf('NewFreq_OccipitalAll_Subj%i',subj(s));
                    load(loadfilename,"Delta_SingleTrials_occ","ITPCdelta_timevec_occ","ITPCdelta_nTrials_occ") %(time points x electrodes x trials)
                    for c=1:3
                        Delta_ITPC(s,c,:,:)=circ_r(Delta_SingleTrials_occ{c},[], [], 3);
                        Delta_TimeVec(s,c,:)=ITPCdelta_timevec_occ{c};
                        Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_occ(c);
                    end
                    subtitle='All Trials Occipital';
                elseif ccluster==1 && catchonly
                    loadfilename=sprintf('NewFreq_OccipitalCatch_Subj%i',subj(s));
                    load(loadfilename,"Delta_SingleTrials_occ_catch","ITPCdelta_timevec_occ_catch","ITPCdelta_nTrials_occ_catch") %(time points x electrodes x trials)
                    for c=1:3
                        Delta_ITPC(s,c,:,:)=circ_r(Delta_SingleTrials_occ_catch{c},[], [], 3);
                        Delta_TimeVec(s,c,:)=ITPCdelta_timevec_occ_catch{c};
                        Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_occ_catch(c);
                    end
                    subtitle='Catch Occipital';
                elseif ccluster==2 && ~catchonly
                    loadfilename=sprintf('NewFreq_CentralAll_Subj%i',subj(s));
                    load(loadfilename,"Delta_SingleTrials_cen","ITPCdelta_timevec_cen","ITPCdelta_nTrials_cen") %(time points x electrodes x trials)
                    for c=1:3
                        Delta_ITPC(s,c,:,:)=circ_r(Delta_SingleTrials_cen{c},[], [], 3);
                        Delta_TimeVec(s,c,:)=ITPCdelta_timevec_cen{c};
                        Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_cen(c);
                    end
                    subtitle='All Trials Central';
                elseif ccluster==2 && catchonly
                    loadfilename=sprintf('NewFreq_CentralCatch_Subj%i',subj(s));
                    load(loadfilename,"Delta_SingleTrials_cen_catch","ITPCdelta_timevec_cen_catch","ITPCdelta_nTrials_cen_catch") %(time points x electrodes x trials)
                    for c=1:3
                        Delta_ITPC(s,c,:,:)=circ_r(Delta_SingleTrials_cen_catch{c},[], [], 3);
                        Delta_TimeVec(s,c,:)=ITPCdelta_timevec_cen_catch{c};
                        Delta_Ntrials(s,c,:)=ITPCdelta_nTrials_cen_catch(c);
                    end
                    subtitle='Catch Central';
                end
            end

            % Time Vec
            timeVec=squeeze(Delta_TimeVec(1,1,:)); % All TV are identical

            % For all subjects and conditions calculate the chance ITPC (for later BL correction)

            for part=1:length(Delta_Ntrials)
                for con=1:width(Delta_Ntrials)
                    % Chance angles
                    angles=2*pi*rand(Delta_Ntrials(part,con),10000);
                    angle_means(part,con)=mean(circ_r(angles));
                    angle_std(part,con)=std(circ_r(angles));
                end


                % Subtract baseline from each condition
                Delta_ITPC=Delta_ITPC-angle_means;
            end

            % Average for each participant across electrodes
            %PartMean=squeeze(mean(Delta_ITPC(:,:,:,elec),4));
            PartMean=Delta_ITPC;

                % Save
                if ccluster==1 && ~catchonly
                    savefilename="GL_Delta_Occ_All";
                    save(savefilename,"PartMean","Delta_Ntrials","angle_means","angle_std","timeVec")
                elseif ccluster==1 && catchonly
                    savefilename="GL_Delta_Occ_Catch";
                    save(savefilename,"PartMean","Delta_Ntrials","angle_means","angle_std","timeVec")
                elseif ccluster==2 && ~catchonly
                    savefilename="GL_Delta_Cen_All";
                    save(savefilename,"PartMean","Delta_Ntrials","angle_means","angle_std","timeVec")
                elseif ccluster==2 && catchonly
                    savefilename="GL_Delta_Cen_Catch";
                    save(savefilename,"PartMean","Delta_Ntrials","angle_means","angle_std","timeVec")
                end
        else
            % If average has already been computed, just load this
            if ccluster==1 && ~catchonly
                load("GL_Delta_Occ_All")
            elseif ccluster==1 && catchonly
                load("GL_Delta_Occ_Catch")
            elseif ccluster==2 && ~catchonly
                load("GL_Delta_Cen_All")
            elseif ccluster==2 && catchonly
                load("GL_Delta_Cen_Catch")
            end
        end
    end
end
%% Plot
cd ' Y:\el-Christina\SxA\SxA_Results\New Delta Results'
figure("Position",[269.8000  145.0000  778.2000  617.0000]);
t = tiledlayout('flow');
for ccluster=clusters

      % Define
    if ccluster==1
        elec=[25:30 62:64]; % occipital electrodes
    else
        elec=[11:12 46:49]; % central electrodes
    end

    for trialt=trialtypes
        if trialt==1
            catchonly=0;
        else
            catchonly=1;
        end

        % If average has already been computed, just load this
        if ccluster==1 && ~catchonly
            load("GL_Delta_Occ_All")
        elseif ccluster==1 && catchonly
            load("GL_Delta_Occ_Catch")
        elseif ccluster==2 && ~catchonly
            load("GL_Delta_Cen_All")
        elseif ccluster==2 && catchonly
            load("GL_Delta_Cen_Catch")
        end

        % Cluster Based Permutation
        if clust_perm
            if trialt==1
                cluster_data=squeeze(mean(PartMean(:,:,:,elec),4)); %average across electodes
            else
                cluster_data=Part_Delta;
            end
            % Rhy - Irr
            [Rhy_Clusters]=clusterBasedPermTest(squeeze(cluster_data(:,1,:)), squeeze(cluster_data(:,3,:)), 1);

            % Int - Irr
            [Int_Clusters]=clusterBasedPermTest(squeeze(cluster_data(:,2,:)), squeeze(cluster_data(:,3,:)), 1) ;

            % Rhy - Int
            [RhyInt_Clusters]=clusterBasedPermTest(squeeze(cluster_data(:,1,:)), squeeze(cluster_data(:,2,:)), 1);

        end

        % Average
        if trialt==1  % Correct way
            % Average Across Participants
            Delta=squeeze(mean(PartMean(:,:,:,elec),1));

            % Average Across Electrodes
            Delta=squeeze(mean(Delta,3));
            Part_Delta=squeeze(mean(PartMean(:,:,:,elec),4));
        else % apparently did not save electrodes separately for catch trials, thus, do this:
            Delta=squeeze(mean(PartMean(:,:,:),1));  % Average Across Participants
            Part_Delta=PartMean;
        end

        nexttile
        for c=1:3

            % Choose Colours
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

            % Plot (with variances?)
            if varplotting
                varplot(timeVec,squeeze(Part_Delta(:,c,:))',"LineWidth",2,'Color',colourvec1)
            else
                plot(timeVec,Delta(c,:),"LineWidth",2,'Color',colourvec1)
            end

            ylim([0 0.6])
            hold on

            % Plot Individual Means
            if plotindividuals
                for s=1:length(subj)
                    plot(timeVec,squeeze(PartMean(s,c,:)),"LineWidth",1,'Color',colourvec2)
                end
            end
        end

        % Plot Significant Clusters
        if clust_perm
            ylimits = ylim;
            for cl=1:size(Rhy_Clusters,1) % Rhythm-Irregular
                startcl=timeVec(Rhy_Clusters(cl,1)); % Get start point of cluster
                endcl= startcl+Rhy_Clusters(cl,2);% Get end point of cluster
                line([startcl endcl],[ylimits(2)-0.08 ylimits(2)-0.08], 'Color',[0 0.4470 0.7410], 'LineWidth', 2)
            end

            for cl=1:size(Int_Clusters,1) % Interval - Irregular
                startcl=timeVec(Int_Clusters(cl,1)); % Get start point of cluster
                endcl= startcl+Int_Clusters(cl,2);% Get end point of cluster
                line([startcl endcl],[ylimits(2)-0.1 ylimits(2)-0.1], 'Color',[0.8500 0.3250 0.0980], 'LineWidth', 2)
            end

            for cl=1:size(RhyInt_Clusters,1) %  Rhythm - Interval
                startcl=timeVec(RhyInt_Clusters(cl,1)); % Get start point of cluster
                endcl= startcl+RhyInt_Clusters(cl,2);% Get end point of cluster
                line([startcl endcl], [ylimits(2)-0.12 ylimits(2)-0.12], 'Color',[0 0 0], 'LineWidth', 2)
            end
        end

        % Choose
        if catchonly
            xlim([-1500 400])
            xticks([-1200:400:400])
            xline(-900)
           % line(0)
             if ccluster==1
                title('Occipital Cluster')
            else
                title('Central Cluster')
            end
        else
            xlim([-700 1200])
            xticks([-400:400:1200])
            xline(900,"-","Target","FontSize",12,"LabelOrientation","horizontal")
            xline(0,"-", "Warning Signal","FontSize",12,"LabelOrientation","horizontal")
            if ccluster==1
                title('Occipital Cluster')
            else
                title('Central Cluster')
            end
        end
       % xline(0)
        yticks([0:0.2:0.6])
        ylabel('Delta (0.66-1.87 Hz) ITPC')
        xlabel('Time (ms)')
        ax=gca;
        ax.FontSize = 12;
        box off
        % title(subtitle)

        if shade % Shades the statistics window
            y_limits = ylim;
            patch([stats_tw(1) stats_tw(2) stats_tw(2) stats_tw(1)], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], 'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end

        headline=sprintf('Delta Band (0.66-1.87 Hz) Intertrial Phase Coherence');

        title(t,headline)
        %clearvars -except t subj elec subj1 subj2 subj3 clusters catchonly trialtypes cbatch trialt ccluster batches plotindividuals varplotting stats_tw shade timeVec
    end
end
%% Topography
cd ' Y:\el-Christina\SxA\SxA_Results\New Delta Results'
    figure; t = tiledlayout('flow');

    for trialt=trialtypes
        % Define
        if trialt==1
            catchonly=0;
            timeROI=[750 850]; % target at 900
        else
            catchonly=1;
            timeROI=[-200 -50]; % target at 0
        end

       % Load Group Data
        if ccluster==1 && ~catchonly
            load("GL_Delta_Occ_All")
        elseif ccluster==1 && catchonly
            load("GL_Delta_Occ_Catch")
        elseif ccluster==2 && ~catchonly
            load("GL_Delta_Cen_All")
        elseif ccluster==2 && catchonly
            load("GL_Delta_Cen_Catch")
        end

        % Average Across Participants
        Delta=squeeze(mean(PartMean,1));

        for c=1:3
            GL_Mean(c,:,:)=squeeze(mean(Delta(c,:,:),1)); % average across participants (input: subj x cond x time points x electrodes | output: condition x time points x electrodes)
            vectimeROI=timeVec>=timeROI(1) & timeVec<=timeROI(2);
            GL_mean_tw(c,:)=squeeze(mean(GL_Mean(c,vectimeROI,:),2)); % average across time window (output: condition x electrodes)
        end


        % Plot
        maxVal=max(GL_mean_tw,[],"all");
        figure('Position', [113.80, 186.60, 1320.00, 492.00]); % Example of setting the size;
        for c=1:3
            % topography
            subplot(1,4,c)
            curr_topo_data=GL_mean_tw(c,:);
            topoplot(curr_topo_data,'head64.locs','electrodes','on','style','map','shading','interp','maplimits',[0  maxVal]); 
            colorbar('off'); % Turn off the color bar
            if c==1
                title('Rhythm','FontSize',12)
            elseif c==2
                title('Interval','FontSize',12)
            elseif c==3
                title('Irregular','FontSize',12)
            end
        end

        % Make Nice
        sgtitle("Delta ITPC Topography before Target Onset")
        subplot(1,4,4);
        c = colorbar; % Display the color bar
        caxis([0,maxVal]); % Set the same color axis for all topoplots
        axis off; % Turn off the axis in the color bar subplot
        c.Position = [0.8, 0.35, 0.02, 0.4]; % position nicely
        % Get the default ticks and show only every second tick
        defaultTicks = c.Ticks; % Get the default tick values
        c.Ticks = defaultTicks(1:2:end); % Keep every second tick
        c.FontSize = 12; % Change to your desired font size

        clearvars -except t subj elec subj1 subj2 subj3 clusters catchonly trialtypes cbatch trialt ccluster batches plotindividuals stats_tw
    end

%% Statistics in ROI Time Window
cd ' Y:\el-Christina\SxA\SxA_Results\New Delta Results'
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
            stats_tw=[800 900]; % target at 900
        else
            catchonly=1;
            stats_tw=[-100 -0]; % target at 0
        end

        % Load Group Data
        if ccluster==1 && ~catchonly
            load("GL_Delta_Occ_All")
        elseif ccluster==1 && catchonly
            load("GL_Delta_Occ_Catch")
        elseif ccluster==2 && ~catchonly
            load("GL_Delta_Cen_All")
        elseif ccluster==2 && catchonly
            load("GL_Delta_Cen_Catch")
        end

        delta_tw=PartMean(:,:,timeVec>=stats_tw(1)&timeVec<=stats_tw(2),:); % select only time window
        delta_tw=squeeze(mean(delta_tw,4)); % average across electrodes
        delta_tw=squeeze(mean(delta_tw,3)); % average across time points

        % Rhythm higher than irregular?
        [h_delta(1),p_delta(1),~,stats_delta(1,:)] = ttest(delta_tw(:,1),delta_tw(:,3),"Tail","right");
        % Interval higher than irregular?
        [h_delta(2),p_delta(2),~,stats_delta(2,:)] = ttest(delta_tw(:,2),delta_tw(:,3),"Tail","right");
        % Difference between rhythm and interval?
        [h_delta(3),p_delta(3),~,stats_delta(3,:)] = ttest(delta_tw(:,1),delta_tw(:,2));
    end
end