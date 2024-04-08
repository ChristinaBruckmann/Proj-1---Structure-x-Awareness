% What do I need?
% From 3: First block of data, all the background info. First 50 trials of trial matrix 3.
% From 33: Second block of data, first 40 trials of trial matrix 1.
% from 31: The rest of the trials and all other parts of the trial matrices

% Merge results:
load PilotProjectSxA1.1_s31.mat
datarest=subresults.data;
load PilotProjectSxA1.1_s33.mat
datablock2=subresults.data;
load PilotProjectSxA1.1_s3.mat
subresults.data=[subresults.data; datablock2; datarest];

% Fix subject number
subresults.data.Subject(:)=3;

save PilotProjectSxA1.1_s3.mat

% Fix trial matrix:

load PilotProjectSxA1.1_s31.mat
trialmatrix3part2=subresults.trialmatrix3(51:end,:);
load PilotProjectSxA1.1_s3.mat
subresults.trialmatrix3=[subresults.trialmatrix3(1:50,:); trialmatrix3part2];

save PilotProjectSxA1.1_s3.mat

load PilotProjectSxA1.1_s33.mat
trialmatrix1firstblock=subresults.trialmatrix1(1:40,:);
load PilotProjectSxA1.1_s31.mat
trialmatrixrest=subresults.trialmatrix1(41:end,:);
load PilotProjectSxA1.1_s3.mat
subresults.trialmatrix1=[trialmatrix1firstblock; trialmatrixrest];

save PilotProjectSxA1.1_s3.mat

load PilotProjectSxA1.1_s31.mat
trialmatrix2all=subresults.trialmatrix2;
load PilotProjectSxA1.1_s3.mat
subresults.trialmatrix2=trialmatrix2all;

save PilotProjectSxA1.1_s3.mat
