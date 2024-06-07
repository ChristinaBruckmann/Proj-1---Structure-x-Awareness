% Multiple Logistic Regression Simulation

%% fixed values of aData
clear
clc
close all

contrastVals=1:10;
contrastThreshold=5;
aDataVals=-5:4;
b=[0,1,1,5]; % chosen weights (intercept, contrast, alpha)

contrastData=sort(repmat(contrastVals',10,1)); % controlled factor (contrast levels)
aData=repmat(aDataVals', 10,1);
pC=1./(1+exp(-1*(b(1)+b(2)*(contrastData-contrastThreshold)+b(3)*aData))); % Percent correct

figure; surf(contrastVals,aDataVals,reshape(pC,10,10)); % either surf or mesh, whichever you like more


%% random values of aData - fit regression model to simluated data
clear
clc
close all

contrastVals=1:10;
trialsPerContrastVal=100;
contrastThreshold=5;
noiseStd=1;

b=[0,1,1]; % chosen weights (intercept, contrast, alpha)

contrastData=repmat(contrastVals',trialsPerContrastVal,1); % controlled factor (contrast levels)
aData=2*randn(size(contrastData)); % random factor (e.g. spontaneous alpha)
pC=1./(1+exp(-1*(b(1)+b(2)*(contrastData-contrastThreshold)+b(3)*aData))); % Percent correct

% create simulated binary data based on pC
isCorrect=rand(size(pC))<pC; % the higher pC, the more likely to have a correct response
% isCorrect=rand(size(pC))<pC+noiseStd*randn(size(pC)); % this version allows to add more noise when choosing correct/wrong for a given pC 

% run logistic regression
nominalResp=isCorrect+1; % needed as this function doesnt take 0 as category value
[b,~,stats]=mnrfit([contrastData, aData], nominalResp);
b %coefficients
stats.t %t vals
%% 

