% Merge LogReg results from Cluster
clear
clc

subj=[101:103 105:108 110 112:114 117:119 121 122 124 126 127 129 130:132];

for s=1:length(subj)
    % Load
    load(sprintf('EEG_SxA_LogRes_GL_BL_Subj%i.mat',subj(s)))

    % Save
    GL_res_obj(s,:,:,:,:,:)=res_obj(end,:,:,:,:,:);
    GL_res_subj(s,:,:,:,:,:)=res_subj(end,:,:,:,:,:);
end

% Save Group Level
save("EEG_SxA_LogRes_GL_Baseline","GL_res_subj","GL_res_obj","subj","timeVec","elecs","freqs")