%% develop matlab solution for logistic regression on binary data with non-zero asymptote and compare to R

%% test matlab's fitglm against R 
clear
clc

nObs=1000;
bSim=[-5, 10, -5];
noiseStd=0.5;

xData=zscore(rand(nObs,(length(bSim)-1)));

pC=1./(1+exp(-1*(bSim(1)+bSim(2)*xData(:,1)+bSim(3)*xData(:,2)+noiseStd*randn(nObs,1))));
resp=pC>0.5;

mdl = fitglm(xData, resp, "Distribution", "binomial")

dataForR=[resp xData];
save('SimDataforR.mat', 'dataForR', '-v6')

%% custom code logistic regression vs. fitglm - binned 1 predictor
clear 
clc

bSim=[-2, 0.5];
xData=kron([1:10],ones(1,1000))';

noiseStd=1.65;
% noiseStd=3;

pC=1./(1+exp(-1*([ones(size(xData,1),1), xData]*bSim'+noiseStd*randn(length(xData),1))));
resp=pC>0.5;

tic; mdl = fitglm(xData, resp, "Distribution", "binomial"); toc

% plot
meansBins=mean(reshape(resp,1000,10));
xForPlot=(min(xData)-10):0.01:(max(xData)+10);
predictedP=1./(1+exp(-1*(mdl.Coefficients{1,1}+mdl.Coefficients{2,1}*xForPlot)));

figure; hold on
plot(xForPlot, predictedP)
scatter(unique(xData),meansBins)

midpoint=-mdl.Coefficients{1,1}/mdl.Coefficients{2,1};
slope=mdl.Coefficients{2,1};

title(['midpoint=' num2str(midpoint) ', slope=' num2str(slope)])

% using custom function
tic; [coeffs, resSorted, errSorted]=logRegression(xData,resp,0); toc

fitGlmCoeffs=mdl.Coefficients{:,1}'
customFuncCoeffs=coeffs

%% custom code logistic regression vs. fitglm - n continuous predictors
clear
clc

nPredictors=3;
bSim=[-2, 0.5 1 -0.5]; % length needs to b nPredictors+1
xData=randn(10000,nPredictors);

noiseStd=1.65;
noiseStd=3;

pC=1./(1+exp(-1*([ones(size(xData,1),1), xData]*bSim'+noiseStd*randn(length(xData),1))));
resp=pC>0.5;

tic; mdl = fitglm(xData, resp, "Distribution", "binomial"); toc

% using custom function
tic; [coeffs, resSorted, errSorted]=logRegression(xData,resp,0); toc

fitGlmCoeffs=mdl.Coefficients{:,1}'
customFuncCoeffs=coeffs

clear 
clc

bSim=[-2, 0.5];
xData=kron([1:10],ones(1,1000))';

noiseStd=1.65;
% noiseStd=3;

pC=1./(1+exp(-1*([ones(size(xData,1),1), xData]*bSim'+noiseStd*randn(length(xData),1))));
resp=pC>0.5;

tic; mdl = fitglm(xData, resp, "Distribution", "binomial"); toc

% plot
meansBins=mean(reshape(resp,1000,10));
xForPlot=(min(xData)-10):0.01:(max(xData)+10);
predictedP=1./(1+exp(-1*(mdl.Coefficients{1,1}+mdl.Coefficients{2,1}*xForPlot)));

figure; hold on
plot(xForPlot, predictedP)
scatter(unique(xData),meansBins)

midpoint=-mdl.Coefficients{1,1}/mdl.Coefficients{2,1};
slope=mdl.Coefficients{2,1};

title(['midpoint=' num2str(midpoint) ', slope=' num2str(slope)])

% using custom function
tic; [coeffs, resSorted, errSorted]=logRegression(xData,resp,0); toc

fitGlmCoeffs=mdl.Coefficients{:,1}'
customFuncCoeffs=coeffs


%% custom code logistic regression vs. fitting logistic function to bin means - binned 1 predictor, chance level = 50%
% cant compare with fitglm as it doesnt accept non-zero asymptote
clear
clc

chanceLev=0.5;
bSim=[-8, 2];
%bSim=[-2, 0.5];

xData=kron([1:0.5:10],ones(1,1000))';

noiseStd=1.65;
% noiseStd=4;

pC=chanceLev+(1-chanceLev)./(1+exp(-1*([ones(size(xData,1),1), xData]*bSim'+noiseStd*randn(length(xData),1))));
resp=rand(size(pC))<pC;

meansBins=mean(reshape(resp,1000,length(unique(xData))));

[bestFitParams, resFitSorted]=fitPsychometric(unique(xData), meansBins');

% using custom function
tic; [coeffs, resSorted, errSorted]=logRegression(xData,resp,0.5); toc

fitToMeansCoeffs=bestFitParams(1:2)
customFuncCoeffs=[-coeffs(1)/coeffs(2) coeffs(2)]

%%
