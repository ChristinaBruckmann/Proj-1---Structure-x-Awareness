%% FitPsychometricAssaf

filename='ResultsSubjectOfficial10.mat';
load(filename)

%% Logistic Model like in simulation code

newfitmethod=1; % New fit is latest code from Assaf, old fit is like I did in the simulation.


% Model Functions
logisticFunob='(0.5./(1+exp(-a*(x-b))))+0.5'; %Objective
logisticFunsu='1./(1+exp(-a*(x-b)))'; % Subjective

% Fit for each condition

samplepoints=1:nContrasts;

if newfitmethod

    % Fit curved
    for idxcond=1:3
        if idxcond==1
            modelsubjrhy=fitPsychometric(samplepoints,subjmeancondcontrast(idxcond,:),2); %output: midpoint, slope, upper bound, lower bound (optional residuals)
            modelobjrhy=fitPsychometric(samplepoints,objmeancondcontrast(idxcond,:),1);
        elseif idxcond==2
            modelsubjint=fitPsychometric(samplepoints,subjmeancondcontrast(idxcond,:),2);
            modelobjint=fitPsychometric(samplepoints,objmeancondcontrast(idxcond,:),1);
        elseif idxcond==3
            modelsubjirr=fitPsychometric(samplepoints,subjmeancondcontrast(idxcond,:),2);
            modelobjirr=fitPsychometric(samplepoints,objmeancondcontrast(idxcond,:),1);
        end
    end

    % Save Values
    modelfit.subjrhy=modelsubjrhy;
    modelfit.objrhy=modelobjrhy;
    modelfit.subjint=modelsubjint;
    modelfit.objint=modelobjint;
    modelfit.subjirr=modelsubjirr;
    modelfit.objirr=modelobjirr;

    % Estimated Model with Best Fit
    EstPsycFunc_subjrhy = modelfit.subjrhy(4)+(modelfit.subjrhy(3)-modelfit.subjrhy(4))./(1+exp(-modelfit.subjrhy(2)*(samplepoints-modelfit.subjrhy(1))));
    EstPsycFunc_objrhy = modelfit.objrhy(4)+(modelfit.objrhy(3)-modelfit.objrhy(4))./(1+exp(-modelfit.objrhy(2)*(samplepoints-modelfit.objrhy(1))));
    EstPsycFunc_subjint = modelfit.subjint(4)+(modelfit.subjint(3)-modelfit.subjint(4))./(1+exp(-modelfit.subjint(2)*(samplepoints-modelfit.subjint(1))));
    EstPsycFunc_objint = modelfit.objint(4)+(modelfit.objint(3)-modelfit.objint(4))./(1+exp(-modelfit.objint(2)*(samplepoints-modelfit.objint(1))));
    EstPsycFunc_subjirr = modelfit.subjirr(4)+(modelfit.subjirr(3)-modelfit.subjirr(4))./(1+exp(-modelfit.subjirr(2)*(samplepoints-modelfit.subjirr(1))));
    EstPsycFunc_objirr = modelfit.objirr(4)+(modelfit.objirr(3)-modelfit.objirr(4))./(1+exp(-modelfit.objirr(2)*(samplepoints-modelfit.objirr(1))));

    % Plot Curves
    if plots
        figure; % subjective curves (assess fit by comparing to actual values)

        subplot(1,3,1)
        plot(subjmeancondcontrast(1,:),'o')
        hold on
        plot(EstPsycFunc_subjrhy)
        title('Rhythm Subjective')
        hold off

        subplot(1,3,2)
        plot(subjmeancondcontrast(2,:),'o')
        hold on
        plot(EstPsycFunc_subjint)
        title('Interval Subjective')
        hold off

        subplot(1,3,3)
        plot(subjmeancondcontrast(3,:),'o')
        hold on
        plot(EstPsycFunc_subjirr)
        title('Irregular Subjective')
        hold off

        figure; % objective curves (assess fit by comparing to actual values)

        subplot(1,3,1)
        plot(objmeancondcontrast(1,:),'o')
        hold on
        plot(EstPsycFunc_objrhy)
        title('Rhythm Objective')
        hold off

        subplot(1,3,2)
        plot(objmeancondcontrast(2,:),'o')
        hold on
        plot(EstPsycFunc_objint)
        title('Interval Objective')
        hold off

        subplot(1,3,3)
        plot(objmeancondcontrast(3,:),'o')
        hold on
        plot(EstPsycFunc_objirr)
        title('Irregular Objective')
        hold off
    end

    % Compare different conditions to each other:
    figure;
    subplot(1,2,1) %subjective
    plot(EstPsycFunc_subjrhy)
    hold on
    plot(EstPsycFunc_subjint)
    plot(EstPsycFunc_subjirr)
    title('Subjective')
    legend('Rhythm','Interval','Irregular')
    hold off

    subplot(1,2,2)  % objective
    plot(EstPsycFunc_objrhy)
    hold on
    plot(EstPsycFunc_objint)
    plot(EstPsycFunc_objirr)
    title('Objective')
    legend('Rhythm','Interval','Irregular')
    hold off

    % Plot midpoints
    % Extract Midpoints
    mpsubj(1)=modelsubjrhy(1);
    mpobj(1)=modelobjrhy(1);
    mpsubj(2)=modelsubjint(1);
    mpobj(2)=modelobjint(1);
    mpsubj(3)=modelsubjirr(1);
    mpobj(3)=modelobjirr(1);

    if plots
        midpointsplot=figure;% plot midpoints

        subplot(1,2,1)
        plot(1:3,mpsubj,'o','LineWidth',2)
        title('Midpoints Subjective Visibility')
        ylabel('b')
        xlim([0 4])
        xticks(1:3)
        xticklabels({'Rhythm','Interval','Irregular'})

        subplot(1,2,2)
        plot(1:3,mpobj,'o','LineWidth',2)
        title('Midpoints Objective Performance')
        ylabel('b')
        xlim([0 4])
        xticks(1:3)
        xticklabels({'Rhythm','Interval','Irregular'})
    end

else % Old method
    for idxcond=1:3
        if idxcond==1
            modelsubjrhy=fit([0 samplepoints]',[0 subjmeancondcontrast(idxcond,:)]',logisticFunsu);
            modelobjrhy=fit([0 samplepoints]',[0.5 objmeancondcontrast(idxcond,:)]',logisticFunob);
        elseif idxcond==2
            modelsubjint=fit([0 samplepoints]',[0 subjmeancondcontrast(idxcond,:)]',logisticFunsu);
            modelobjint=fit([0 samplepoints]',[0.5 objmeancondcontrast(idxcond,:)]',logisticFunob);
        elseif idxcond==3
            modelsubjirr=fit([0 samplepoints]',[0 subjmeancondcontrast(idxcond,:)]',logisticFunsu);
            modelobjirr=fit([0 samplepoints]',[0.5 objmeancondcontrast(idxcond,:)]',logisticFunob);
        end
    end

    modelfit.subjrhy=modelsubjrhy;
    modelfit.objrhy=modelobjrhy;
    modelfit.subjint=modelsubjint;
    modelfit.objint=modelobjint;
    modelfit.subjirr=modelsubjirr;
    modelfit.objirr=modelobjirr;

    % Plot data points and model fits

    plot(xVals, estimatedPsycFunction, 'r')

    if plots
        subjplot=figure; % subjective curves

        subplot(1,3,1)
        plot(subjmeancondcontrast(1,:),'o')
        hold on
        plot(modelsubjrhy)
        title('Rhythm Subjective')
        hold off

        subplot(1,3,2)
        plot(subjmeancondcontrast(2,:),'o')
        hold on
        plot(modelsubjint)
        title('Interval Subjective')
        hold off

        subplot(1,3,3)
        plot(subjmeancondcontrast(3,:),'o')
        hold on
        plot(modelsubjirr)
        title('Irregular Subjective')
        hold off

        objplot=figure; % objective curves

        subplot(1,3,1)
        plot(objmeancondcontrast(1,:),'o')
        hold on
        plot(modelobjrhy)
        title('Rhythm Objective')
        hold off

        subplot(1,3,2)
        plot(objmeancondcontrast(2,:),'o')
        hold on
        plot(modelobjint)
        title('Interval Objective')
        hold off

        subplot(1,3,3)
        plot(objmeancondcontrast(3,:),'o')
        hold on
        plot(modelobjirr)
        title('Irregular Objective')
        hold off
    end

    % Extract Midpoints
    mpsubj(1)=modelsubjrhy.b;
    mpobj(1)=modelobjrhy.b;
    mpsubj(2)=modelsubjint.b;
    mpobj(2)=modelobjint.b;
    mpsubj(3)=modelsubjirr.b;
    mpobj(3)=modelobjirr.b;

    if plots
        midpointsplot=figure;% plot midpoints

        subplot(1,2,1)
        plot(1:3,mpsubj,'o','LineWidth',2)
        title('Midpoints Subjective Visibility')
        ylabel('b')
        xlim([0 4])
        xticks(1:3)
        xticklabels({'Rhythm','Interval','Irregular'})

        subplot(1,2,2)
        plot(1:3,mpobj,'o','LineWidth',2)
        title('Midpoints Objective Performance')
        ylabel('b')
        xlim([0 4])
        xticks(1:3)
        xticklabels({'Rhythm','Interval','Irregular'})
    end
end

% Save
save(filename, 'mpsubj','mpobj','modelfit',"-append")
%% Assaf's function used above
function [bestFitParams, resFitSorted]=fitPsychometric(xVals, accData, objective_subjective)

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

% objective_subjective=2; %1=objective, 2=subjective (added as input now)


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


% define the to be optimized (=minimized) function. This would be the error 
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


% figure; hold all
% plot(xVals,accData, 'o'); % plot data

% % estimated model with best fit params
% estimatedPsycFunction = bestFitParams(4)+(bestFitParams(3)-bestFitParams(4))./(1+exp(-bestFitParams(2)*(xVals-bestFitParams(1))));
% 
% plot(xVals, estimatedPsycFunction, 'r')

end

%%
function fitErr=err_psychFunFit(inParams, xVals, accData)

% calulactes least square error of data from psychometric function model simulated using a given parameter set
shift=inParams(1);
slope=inParams(2);
maxAsymp=inParams(3);
minAsymp=inParams(4);

estimatedPsycFunction = minAsymp+(maxAsymp-minAsymp)./(1+exp(-slope*(xVals-shift)));

fitErr=mean((accData-estimatedPsycFunction).^2);

end

