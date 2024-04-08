% Obj-Subj Cluster-Based Permutation Analysis
% Input: Subject Results from Cross_Trial_Logregression
clear
clc
subj=[18:22];
calc_var=3; % intercept,contrast,power,interaction

% Load b-matrices
for s=1:length(subj)
    loadfilename=sprintf('EEG_SxA_Subj%i_Results_SubjObj.mat',subj(s));
    load(loadfilename)

    % Average conditions
    tempmeanob=squeeze(mean(total_bob,1)); % original dimensions (conditions x freq x timepoints x calc_var)
    tempmeansu=squeeze(mean(total_bsu,1)); % original dimensions (conditions x freq x timepoints x calc_var)

    % Store
    b_mat_ob(s,:,:)=tempmeanob(:,:,calc_var); % only chose relevant b (calc_var) and add all subj
    b_mat_su(s,:,:)=tempmeansu(:,:,calc_var); % output dim (subj x freq x tp)
end

% Average across subjects and plot (Raw average, no permutation)
b_avg_ob=squeeze(mean(b_mat_ob,1));
b_avg_su=squeeze(mean(b_mat_su,1));

figure; imagesc(flipud(b_avg_ob)); colorbar; 
caxis([-1 1]);
xline(800-timep(1));
set(gca,'YDir','normal')
title('Objective Obtained');

figure; imagesc(flipud(b_avg_su)); colorbar;
xline(800-timep(1));
caxis([-1 1]);
set(gca,'YDir','normal')
title('Subjective Obtained');

% Permutation
nperm=100;
mult_opt=[-1 1];
for currp=1:nperm
    for s=1:(length(subj))
    % Obj
    mult=mult_opt(randi(2)); % Choose randomly whether pos or neg sign
    curr_ob(s,:,:)=b_mat_ob(s,:,:)*mult;

    % Subj
    mult=mult_opt(randi(2));
    curr_su(s,:,:)=b_mat_su(s,:,:)*mult;
    end
    
    % Average Across Subjects and Store
    perm_ob(currp,:,:)=mean(curr_su,1);
    perm_su(currp,:,:)=mean(curr_su,1);
end

% Plot average across perms
figure; imagesc(flipud(squeeze(mean(perm_ob,1)))); colorbar; 
caxis([-1 1]);
xline(800-timep(1));
set(gca,'YDir','normal')
title('Objective Permutation Mean');

figure; imagesc(flipud(squeeze(mean(perm_su,1)))); colorbar;
xline(800-timep(1));
caxis([-1 1]);
set(gca,'YDir','normal')
title('Subjective Permutation Mean');