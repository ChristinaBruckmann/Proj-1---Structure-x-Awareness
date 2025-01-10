function [bestFitParams, resFitSorted, resErrSorted, startingPointsSorted]=fitPsychometric_binaryData(xData, binaryAccData, asymptoteBounds)

% fitPsychometric - fit psychometric function to probability data calculated at fixed predictor levels.
% uses fmincon optimization to minimize least square error, with constraints on parameters

% Inputs -
% xData:          n-by-p matrix of x values of psychometric function, n=observations, p=predictors
% binaryAccData:  n-by-1 vector of binary accuracy values for each level of x, must be 0 or 1
% asymptoteBounds:lower and upper Bounds for lower and upper Asymptotes: [LBLA UBLA LBUA UBUA]. In most
%                   cases LBLA=0 and UBUA=1. If lBlA=uBlA or lBuA=uBuA, these asymptotes will be fixed to
%                   this value.
%
% Outputs -
% bestFitParams:        best estimate of intercept, slopes, and asymptotes
% resFitSorted:         estimated parameters from all starting points used for fitting, sorted by increasing fit error
% resErrSorted:         estimated fit error from all starting points used for fitting, sorted by increasing fit error
% startingPointsSorted: all starting points used for fitting, sorted by increasing fit error

% notes:



% initial parameter values
nStartPoints_initial=3; % the higher the number the longer it runs, but less chances for local minima, should always be odd
spRangeIntercept_initial=10; % initial range of intercept start points, in std of max var predictor
spRangeSlopes_initial=5; % initial range of slope start points, in std of each predictor
boundsRangeSlopes_initial=30; % intial bounds on slope values, in 1/std of each predictor (smaller std means greater change in y for 1 unit of x)
    
% arrange data - predictor/dependent values must be in columns, participants (or other random unit) in rows
if size(xData,1)==1
    xData=xData';
end
if size(binaryAccData,1)==1
    binaryAccData=binaryAccData';
end

% check illegal values
if sum(~(binaryAccData==1|binaryAccData==0))>0
    error('illegal values in dependent variable (accuracy data). must contain only 1 or 0')
end
if size(xData,1)~=size(binaryAccData,1)
    error('predictor and dependent variable data must have same number of observations')
end

nPred=size(xData,2); % get number of slope predictors (i.e. not counting the intercept)

% initialize adjustable params
nStartPoints=nStartPoints_initial; % the higher the number the longer it runs, but less chances for local minima, should always be odd
spRangeIntercept=spRangeIntercept_initial;
spRangeSlopes=spRangeSlopes_initial;
boundsRangeSlopes=boundsRangeSlopes_initial*ones(1,nPred);

% config fitting algorithm
options=optimoptions(@fmincon,'Display','none','MaxFunctionEvaluations', 10000, 'MaxIterations', 20000, ...
    'ConstraintTolerance', 1.0000e-8, 'OptimalityTolerance', 1.0000e-10);

paramWeightsForConstraint=[zeros(1,nPred+1) 1 -1]; % linear inequality constraint that upper asymptote is always larger than lower asymptote

% define the be optimized (=minimized) function
errFunToMinimize = @(params)err_psychFunFit(params, xData, binaryAccData);

% end preparations


% while loop until estimated params are not suspicious to be local minima/at bounds
conditionsMet=0;
startpointIncreaseCounter=0;
while ~conditionsMet
    
    conditionsCheck=[1 1]; % assume conditions are met, for start points and bounds       
    
    % define multiple starting points for each model param (to avoid local minimum)
    interStartPoints=linspace(-spRangeIntercept*max(std(xData)),spRangeIntercept*max(std(xData)),nStartPoints);
    
    % create matrix of starting points for slopes
    slopeStartPointsMat=[];
    for n=1:nPred
        startPointsCurrPred=linspace(-spRangeSlopes*std(xData(:,n)),spRangeSlopes*std(xData(:,n)),nStartPoints);
        slopeStartPointsMat=[slopeStartPointsMat kron(repmat(startPointsCurrPred',nStartPoints^(n-1),1), ones(nStartPoints^(nPred-n),1))];
    end
    
    if asymptoteBounds(1)==asymptoteBounds(2) % constrained lower asymptote
        asymptoteBounds(2)=asymptoteBounds(2)+1e-10;
        lowAsyStartPoints=asymptoteBounds(1);
    else
        lowAsyStartPoints=linspace(asymptoteBounds(1), asymptoteBounds(2), nStartPoints);
    end
    
    if asymptoteBounds(3)==asymptoteBounds(4) % constrained upper asymptote
        asymptoteBounds(3)=asymptoteBounds(4)-1e-10;
        upperAsyStartPoints=asymptoteBounds(4);
    else
        upperAsyStartPoints=linspace(asymptoteBounds(3), asymptoteBounds(4), nStartPoints);
    end
    
    % define bounds on fitted parameters
    paramLowerBounds=[-Inf -boundsRangeSlopes.*(1./std(xData)) asymptoteBounds(1) asymptoteBounds(3)];
    paramUpperBounds=[Inf boundsRangeSlopes.*(1./std(xData)) asymptoteBounds(2) asymptoteBounds(4)];
    
    % run fitting algorithm with different starting points
    sl_la=[repmat(slopeStartPointsMat, length(lowAsyStartPoints), 1) kron(lowAsyStartPoints', ones(size(slopeStartPointsMat,1),1))];
    sl_la_ua=[repmat(sl_la, length(upperAsyStartPoints), 1) kron(upperAsyStartPoints', ones(size(sl_la,1),1))];
    all_start_points=[kron(interStartPoints', ones(size(sl_la_ua,1),1)) repmat(sl_la_ua, length(interStartPoints), 1)];
    nStartingPoints=size(all_start_points,1);
    
    resFit=zeros(size(all_start_points));
    resErr=zeros(nStartingPoints, 1);
    parfor sp=1:nStartingPoints
        [estParamsCurr, errCurr]=fmincon(errFunToMinimize, all_start_points(sp,:),paramWeightsForConstraint,0,[],[],paramLowerBounds, paramUpperBounds, [], options);
        resFit(sp,:)=estParamsCurr;
        resErr(sp)=errCurr;
    end
    
    % get parameters that give best fit
    [resErrSorted,i]=sort(resErr); % sort according to smallest error
    resFitSorted=resFit(i,:);
    startingPointsSorted=all_start_points(i,:);
    
    bestFitParams=resFitSorted(1,:);
    
    % check if optimization and estimated params are suspicious
    % 1. parameter estimates identical to starting point
    bestResThresh=min(20, round(0.2*size(startingPointsSorted,1)));
    localMinOnStartPoint=mean(abs(startingPointsSorted(1:bestResThresh,1:end-2)-resFitSorted(1:bestResThresh,1:end-2))<0.0001);
    if max(localMinOnStartPoint)>0.8
        disp('WARNING: parameter estimates identical to starting values, increasing range of starting points')
        nStartPoints=nStartPoints+2; % the higher the number the longer it runs, but less chances for local minima, should always be odd
        spRangeIntercept=2*spRangeIntercept;
        spRangeSlopes=2*spRangeSlopes;
        conditionsCheck(1)=0;
        startpointIncreaseCounter=startpointIncreaseCounter+1;
    end
    
    % 2. slopes identical to bounds (not necessary for other params)
    slopeEstimatesAtBound=abs(bestFitParams(2:end-2)-paramLowerBounds(2:end-2))<eps|abs(bestFitParams(2:end-2)-paramUpperBounds(2:end-2))<eps;
    if sum(slopeEstimatesAtBound)~=0
        disp('WARNING: slope estimates identical to bounds, increasing bounds')
        boundsRangeSlopes(slopeEstimatesAtBound)=2*boundsRangeSlopes(slopeEstimatesAtBound);
        conditionsCheck(2)=0;
    end
    
    % its not impossible that the estimates are on the start point, but need to verify its not a local minimum
    if startpointIncreaseCounter>=3
        conditionsCheck(1)=1;
        disp('WARNING: parameter estimates identical to starting values even after increasing range of starting points, check its not a local minimum ')        
    end    
    
    
    if ~ismember(0, conditionsCheck)
        conditionsMet=1;
    end
        
end
    

%%

function fitErr=err_psychFunFit(inParams, xVals, accData)

% calulactes -logLikelihood of data from psychometric function model simulated using a given parameter set

pC=inParams(end-1)+(inParams(end)-inParams(end-1))./(1+exp(-1*([ones(size(xVals,1),1), xVals]*inParams(1:end-2)')));

pC(accData==0)=1-pC(accData==0); % likelihood is reversed for 0 (for a given probability of "correct", the likelihood of an "incorrect" is 1-pC)
pC(pC==0)=eps;

fitErr=-sum(log(pC));


