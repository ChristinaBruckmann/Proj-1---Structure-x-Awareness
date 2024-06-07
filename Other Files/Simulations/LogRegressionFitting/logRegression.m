function [bestFitParams, resFitSorted, resErrSorted]=logRegression(xData, binaryResp, chanceRate)

% fit logistic regression model to binary response data

% params for fitting
nStartingPoints=3; % the higher the number the longer it runs, but less chances for local minima
sl0_range=[-0.5 0 0.5]';


% define range for starting point for each model param (to avoid local minimum)
xRange=max(xData)-min(xData);
sh0_range=linspace(min(min(xData))-max(xRange), max(max(xData))+max(xRange), nStartingPoints)';

startingPointsSlopes=createStartPoints(sl0_range, size(xData,2));
startingPointsAll=[kron(sh0_range, ones(size(startingPointsSlopes,1),1)) repmat(startingPointsSlopes, length(sh0_range),1)];

options.MaxFunEvals= 2000;
options.MaxIter= 2000;

errFunToMinimize = @(params)err_logisticReg_binary(params, xData, binaryResp, chanceRate);

% run fitting algorithm with different starting points
resFit=zeros(size(startingPointsAll,1), size(startingPointsAll,2)*2);
resErr=zeros(size(startingPointsAll,1), 1);

for startPoint=1:size(startingPointsAll,1)
    [estParamsCurr, errCurr]=fminsearch(errFunToMinimize, startingPointsAll(startPoint,:), options);
    resFit(startPoint,:)=[startingPointsAll(startPoint,:) estParamsCurr];
    resErr(startPoint)=errCurr;
end

% sort according to smallest error
[resErrSorted,i]=sort(resErr);
resFitSorted=resFit(i,:);

bestFitParams=resFitSorted(1, end-size(xData,2):end); % average across best fits for stability


function LL=err_logisticReg_binary(params, xData, binaryResp, chanceRate)

% error of logistic regression model, calculated using ML estimate given binomial distribution for each observation
pC=chanceRate+(1-chanceRate)./(1+exp(-1*([ones(size(xData,1),1), xData]*params')));

pC(binaryResp==0)=1-pC(binaryResp==0);
pC(pC==0)=eps;

LL=-sum(log(pC));



function startPointMat=createStartPoints(startPointRange, nParams)

% creates a matrix of all start point combinations across n predictors
if size(startPointRange,1)==1
    startPointRange=startPointRange';
end

nSP=length(startPointRange);

startPointMat=[];
for n=1:nParams
    startPointMat=[startPointMat kron(repmat(startPointRange,nSP^(n-1),1), ones(nSP^(nParams-n),1))];    
end

%%

