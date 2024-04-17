function [bestFitParams, resFitSorted]=fitPsychometric(xVals, accData)

% fit psychometric function using fmincon optimization with parameter constraints
% Inputs -
% xVals: vector of x values of psychometric function
% accData: vector of accuracy values for each level of x
%
% Outputs -
% bestFitParams: best estimate (averaged across best fits) of shift, slope, and asymptotes
% resFitSorted - estimated parameters from all starting points used for fitting

% note:
% 1. there is a parameter that controls whether data is from objective or subjective task, for the lower asymptote
% 2. slope must be positive here, can flip descending data 
% 3. the bounds around the asypmtotes can be made tighter or broader



% params for fitting
nStartingPoints=5; % the higher the number the longer it runs, but less chances for local minima
stabilityRangeFitting=3; % defines across how many fits with lowest error the param estimates should be averaged, for stability against noise

objective_subjective=2; %1=objective, 2=subjective


% define boundaries for each model param [lower upper]
shiftBounds=[min(xVals) max(xVals)]; % cannot estimate function whose midpoint is not in range
slopeBounds=[0 Inf]; % slope must be positive
maxAsympBounds=[0.95 1];

if objective_subjective==1
    minAsympBounds=[0.45 0.55]; % for objective
else
    minAsympBounds=[0 0.05]; % for subjective
end

paramLowerBounds=[shiftBounds(1) slopeBounds(1) maxAsympBounds(1) minAsympBounds(1)];
paramUpperBounds=[shiftBounds(2) slopeBounds(2) maxAsympBounds(2) minAsympBounds(2)];


% define range for starting point for each model param (to avoid local minimum)
sh0_range=linspace(shiftBounds(1), shiftBounds(2), nStartingPoints);
sl0_range=2.^(0:nStartingPoints);
mxA0_range=linspace(maxAsympBounds(1), maxAsympBounds(2), nStartingPoints);
mnA0_range=linspace(minAsympBounds(1), minAsympBounds(2), nStartingPoints);


% config fitting algorithm
options=optimoptions(@fmincon,'Display','notify','MaxFunctionEvaluations', 5000, 'MaxIterations', 2000, ...
    'ConstraintTolerance', 1.0000e-10, 'OptimalityTolerance', 1.0000e-10);  


% define the be optimized (=minimized) function. This would be the error 
% (= mean square difference) of the actual data from the fitted model given a set of parameters
% use function handle to send fixed parameters to the error function
errFunToMinimize = @(params)err_psychFunFit(params, xVals, accData); 


% run fitting algorithm with different starting points
resFit=zeros((nStartingPoints+1)*nStartingPoints^3,5);
counter=1;
for sh0=sh0_range
    for sl0=sl0_range
        for mxA0=mxA0_range
            for mnA0=mnA0_range
                [estParamsCurr, errCurr]=fmincon(errFunToMinimize, [sh0 sl0 mxA0 mnA0],[],[],[],[],paramLowerBounds, paramUpperBounds, [], options);
                resFit(counter,:)=[estParamsCurr errCurr];
                counter=counter+1;
            end
        end
    end
end


% sort according to smallest error
[~,i]=sort(resFit(:,end));
resFitSorted=resFit(i,:);

bestFitParams=mean(resFitSorted(1:stabilityRangeFitting, 1:end-1)); % average across best fits for stability


figure; hold all
plot(xVals,accData, 'o'); % plot data

% estimated model with best fit params
estimatedPsycFunction = bestFitParams(4)+(bestFitParams(3)-bestFitParams(4))./(1+exp(-bestFitParams(2)*(xVals-bestFitParams(1))));

plot(xVals, estimatedPsycFunction, 'r')



%%

function fitErr=err_psychFunFit(inParams, xVals, accData)

% calulactes least square error of data from psychometric function model simulated using a given parameter set
shift=inParams(1);
slope=inParams(2);
maxAsymp=inParams(3);
minAsymp=inParams(4);

estimatedPsycFunction = minAsymp+(maxAsymp-minAsymp)./(1+exp(-slope*(xVals-shift)));

fitErr=mean((accData-estimatedPsycFunction).^2);

