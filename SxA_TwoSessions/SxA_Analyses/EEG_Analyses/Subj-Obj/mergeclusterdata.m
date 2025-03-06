%% Merge Cluster Data
cd 'Y:\el-Christina\SxA\SxA_Results\LogRegResults\Full Analysis\BaselineCorrected'

fulldatasubj=[108, 110, 114, 119, 122, 123, 126, 132];
GL_res_obj=zeros(1,3,9,30,59,5); % subj x cond x elecs x freqs x time points x params
GL_res_subj=zeros(1,3,9,30,59,5); 
% Load Data
for s=1:length(fulldatasubj)
    loadname=sprintf('EEG_SxA_LogRes_GL_BL_Subj%i.mat',fulldatasubj(s));
    load(loadname)
    GL_BL_res_obj=cat(1,GL_res_obj,res_obj);
    GL_BL_res_subj=cat(1,GL_res_subj,res_subj);
end

save('EEG_SxA_LogRes_GL_BL','GL_BL_res_obj','GL_BL_res_subj','timeVec')
