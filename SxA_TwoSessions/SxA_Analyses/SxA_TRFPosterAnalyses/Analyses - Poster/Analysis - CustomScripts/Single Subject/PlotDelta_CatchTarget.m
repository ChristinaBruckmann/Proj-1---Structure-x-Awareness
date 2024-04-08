% Delta Plots for Targets (once only catch, one with target)

clear
clc

subj=17:22;

ITPCdelta_elec=64;

% Load Data Catch
for s=1:length(subj)
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
    loadfilename=sprintf('EEG_SxA_Subj%i_Results.mat',subj(s));
    for c=1:3
        % Occipital
        load(loadfilename,'ITPCdelta_res_catch_occ');
        ITPCdelta_res_occ_gl(s,c,:,:)=ITPCdelta_res_catch_occ{1,c};
        load(loadfilename,'ITPCdelta_timevec_catch_occ');
        ITPCdelta_timevec_occ_gl(s,c,:)=ITPCdelta_timevec_catch_occ{1,c};
        ITPCdelta_nTrialsCond_occ_gl(s,:)=cell2mat(struct2cell(load(loadfilename,'ITPCdelta_nTrials_catch_occ')));

        % Central
        load(loadfilename,'ITPCdelta_res_catch_cen');
        ITPCdelta_res_cen_gl(s,c,:,:)=ITPCdelta_res_catch_cen{1,c};
        load(loadfilename,'ITPCdelta_timevec_catch_cen');
        ITPCdelta_timevec_cen_gl(s,c,:)=ITPCdelta_timevec_catch_cen{1,c};
        ITPCdelta_nTrialsCond_cen_gl(s,:)=cell2mat(struct2cell(load(loadfilename,'ITPCdelta_nTrials_catch_cen')));
    end
end

% Average
for c=1:3
    gl_itpc_results_means_occ(c,:,:)=mean(ITPCdelta_res_occ_gl(:,c,:,:),1);
    gl_itpc_results_means_cen(c,:,:)=mean(ITPCdelta_res_cen_gl(:,c,:,:),1);
end

% Investigate amount of trials and chance angles
for c=1:3 
    % Total N of trials
    totaltrials(c)=round(mean(ITPCdelta_nTrialsCond_occ_gl(:,c),1)); % same for occ and cen

    % Change angles
    angles=2*pi*rand(totaltrials(c),5000);
    angle_mean(c)=mean(circ_r(angles));

    angle_std(c)=std(circ_r(angles));

end

% Plot
topotimeROI = [-100 -50]; % 0 is predicted target onset
warningline=[-800];
targetline=[0];
timeVector=ITPCdelta_timevec_cen_gl(1,1,:); % Same for everything
occElectrodes=[25:27 62:64]; %Po3/Po4/po7/po8/o1/o2
cenElectrodes= [11 38 46:48]; %Fz, FC1, FCz, FC2, Cz

f2=figure;
f1=figure;
for clusters=1:2
    deltaplotdata=[];
    deltatopodata=[];
    for c=1:3
        if clusters==1
            deltaplotdata(c,:,:)=squeeze(gl_itpc_results_means_occ(c,:,occElectrodes));
            deltatopodata(c,:,:)=squeeze(gl_itpc_results_means_occ(c,:,:));
            titletext='ITPC Delta Catch Trials - Occipital';
            set(0, 'currentfigure', f1); 
            subplot(2,1,2)
        else
            deltaplotdata(c,:,:)=squeeze(gl_itpc_results_means_cen(c,:,cenElectrodes));
            deltatopodata(c,:,:)=squeeze(gl_itpc_results_means_cen(c,:,:));
            titletext='ITPC Delta Catch Trials - Central';
            set(0, 'currentfigure', f2); 
            subplot(2,1,2)
        end
    end

    % Time Course
    for c=1:3 % For conditions
        plot(squeeze(timeVector)', mean(deltaplotdata(c,:,:),3));
        hold on
        ylim([-0.1 0.4]);
        legend('Rhythm','Interval','Irregular')
    end
    title(titletext)
    xline(warningline,'LineStyle', '-', 'Color', [0 0 0])
    xline(targetline,'LineStyle', '- -', 'Color', [0 0 0])
    yline(angle_mean(1),'b')
   yline(angle_mean(2),'r')
   yline(angle_mean(3),'y')
   xlim([-1400 600])
end



%% ALL TARGETS

clearvars -except 'f1' 'f2'
clc

subj=17:22;


ITPCdelta_elec=64;

% Load Data All Predictable Targets
for s=1:length(subj)
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
    loadfilename=sprintf('EEG_SxA_Subj%i_Results.mat',subj(s));
    for c=1:2
        % Occipital
        load(loadfilename,'ITPCdelta_res_target_pred_occ');
        ITPCdelta_res_occ_pred_gl(s,c,:,:)=ITPCdelta_res_target_pred_occ{1,c};
        load(loadfilename,'ITPCdelta_timevec_target_pred_occ');
        ITPCdelta_timevec_occ_pred_gl(s,c,:)=ITPCdelta_timevec_target_pred_occ{1,c};
        ITPCdelta_nTrialsCond_occ_pred_gl(s,:)=cell2mat(struct2cell(load(loadfilename,'ITPCdelta_nTrials_target_pred_occ')));

        % Central
        load(loadfilename,'ITPCdelta_res_target_pred_cen');
        ITPCdelta_res_cen_pred_gl(s,c,:,:)=ITPCdelta_res_target_pred_cen{1,c};
        load(loadfilename,'ITPCdelta_timevec_target_pred_cen');
        ITPCdelta_timevec_cen_pred_gl(s,c,:)=ITPCdelta_timevec_target_pred_cen{1,c};
        ITPCdelta_nTrialsCond_cen_pred_gl(s,:)=cell2mat(struct2cell(load(loadfilename,'ITPCdelta_nTrials_target_pred_cen')));
    end
end

% Load Data All Irregular Targets
for s=1:length(subj)
    cd 'C:\Users\cbruckmann\Documents\PhD Projects\Proj1 - StructurexAwareness\Two Session Version - Code and Data\Data Analysis\EEG\Analysis - CustomScripts\Single Subject\Single Subject Results'
    loadfilename=sprintf('EEG_SxA_Subj%i_Results.mat',subj(s));
    for c=1:3
        % Occipital
        load(loadfilename,'ITPCdelta_res_target_irr_occ');
        ITPCdelta_res_occ_irr_gl(s,c,:,:)=ITPCdelta_res_target_irr_occ{1,c};
        load(loadfilename,'ITPCdelta_timevec_target_irr_occ');
        ITPCdelta_timevec_occ_irr_gl(s,c,:)=ITPCdelta_timevec_target_irr_occ{1,c};
        ITPCdelta_nTrialsCond_occ_irr_gl(s,:)=cell2mat(struct2cell(load(loadfilename,'ITPCdelta_nTrials_target_irr_occ')));

        % Central
        load(loadfilename,'ITPCdelta_res_target_irr_cen');
        ITPCdelta_res_cen_irr_gl(s,c,:,:)=ITPCdelta_res_target_irr_cen{1,c};
        load(loadfilename,'ITPCdelta_timevec_target_irr_cen');
        ITPCdelta_timevec_cen_irr_gl(s,c,:)=ITPCdelta_timevec_target_irr_cen{1,c};
        ITPCdelta_nTrialsCond_cen_irr_gl(s,:)=cell2mat(struct2cell(load(loadfilename,'ITPCdelta_nTrials_target_irr_cen')));
    end
end

% Average Predictable
for c=1:2
    gl_itpc_results_means_occ(c,:,:)=mean(ITPCdelta_res_occ_pred_gl(:,c,:,:),1);
    gl_itpc_results_means_cen(c,:,:)=mean(ITPCdelta_res_cen_pred_gl(:,c,:,:),1);
end

% Average Irregular and Add
gl_itpc_results_means_occ(3,:,:)=squeeze(mean(mean(ITPCdelta_res_occ_irr_gl,2),1));
gl_itpc_results_means_cen(3,:,:)=squeeze(mean(mean(ITPCdelta_res_cen_irr_gl,2),1));

% Investigate amount of trials and chance angles predictable
for c=1:2
    % Total N of trials
    totaltrials(c)=round(mean(ITPCdelta_nTrialsCond_occ_pred_gl(:,c),1)); % same for occ and cen
    % Change angles
    angles=2*pi*rand(totaltrials(c),1000);
    angle_mean(c)=mean(circ_r(angles));
    angle_std(c)=std(circ_r(angles));
end
% Add irregular
totaltrials(3)=round(mean(sum(ITPCdelta_nTrialsCond_cen_irr_gl,2)));
angles=2*pi*rand(totaltrials(3),5000);
angle_mean(3)=mean(circ_r(angles));
angle_std(3)=std(circ_r(angles));

% Plot Target Trials
topotimeROI = [-100 -50]; % 0 is predicted target onset
warningline=[-800];
targetline=[0];
timeVector=ITPCdelta_timevec_cen_pred_gl(1,1,:); % Same for everything
occElectrodes=[25:27 62:64]; %Po3/Po4/po7/po8/o1/o2
cenElectrodes= [11 38 46:48]; %Fz, FC1, FCz, FC2, Cz

for clusters=1:2
    deltaplotdata=[];
    deltatopodata=[];

    subplot(2,1,1)
    for c=1:3
        if clusters==1
            deltaplotdata(c,:,:)=squeeze(gl_itpc_results_means_occ(c,:,occElectrodes));
            deltatopodata(c,:,:)=squeeze(gl_itpc_results_means_occ(c,:,:));
            titletext='ITPC All Targets - Occipital';
            set(0, 'currentfigure', f1);
            subplot(2,1,1)
        else
            deltaplotdata(c,:,:)=squeeze(gl_itpc_results_means_cen(c,:,cenElectrodes));
            deltatopodata(c,:,:)=squeeze(gl_itpc_results_means_cen(c,:,:));
            titletext='ITPC All Targets - Central';
            set(0, 'currentfigure', f2);
            subplot(2,1,1)
        end
    end

    % Time Course
    for c=1:3 % For conditions
        plot(squeeze(timeVector)', mean(deltaplotdata(c,:,:),3));
        hold on
        ylim([-0.1 0.4]);
        legend('Rhythm','Interval','Irregular')
    end
    title(titletext)
    xline(warningline,'LineStyle', '-', 'Color', [0 0 0])
    xline(targetline,'LineStyle', '- -', 'Color', [0 0 0])
    yline(angle_mean(1),'b')
    yline(angle_mean(2),'r')
    yline(angle_mean(3),'y')
    xlim([-1400 600])
end



%% Plot Delta Topography (No baseline)

timeROI = [-600 -50];
timeVector=ITPCdelta_timevec_cen_pred_gl(1,1,:); % Same for everything


figure;subplot(1,3,c)
figuretitle=sprintf('Delta ITPC Time Course (%ims to %ims)',timeROI);
for c=1:3 % For conditions
    deltatopodata=[];
    % topography (in relation to baseline)
    subplot(1,3,c)
    deltatopodata(:,:)=squeeze(gl_itpc_results_means_occ(c,:,:));
    deltatopodata(:,65:71)=0;
    topoplot(squeeze(mean(deltatopodata(squeeze(timeVector)'>=timeROI(1) & squeeze(timeVector)'<=timeROI(2), :))),'head64.locs','electrodes','on','style','map','shading','interp','maplimits',[0 0.3]);
    colorbar
    if c==1
        title('Rhythm')
    elseif c==2
        title('Interval')
    elseif c==3
        title('Irregular')
    end
end
sgtitle(figuretitle)
%save('topoplotissue','timeVector','gl_itpc_results_means_occ')