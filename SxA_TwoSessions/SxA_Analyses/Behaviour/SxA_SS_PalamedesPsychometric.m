%% SxA Fits with Palamedes
clear
clc
% Load Data
cd 'Y:\el-Christina\SxA\SxA_Data\Behaviour Preprocessed'
load('SxA_ResultsSubject122_Total.mat','dataforpsignifit')

%Stimulus intensities
for c=1:3
ObjData(c,:,:) = dataforpsignifit{1, c}.obj; 
SubjData(c,:,:)= dataforpsignifit{1, c}.subj; 
end

PF = @PAL_Logistic;  % Choose Psychometric Function


%% Fit Objective Curves
paramsFree = [1 1 1 1];  %1: free parameter, 0: fixed parameter
searchGrid.alpha = 0.01:.001:.11;
searchGrid.beta = logspace(0,3,101);
%searchGrid.beta = logspace(log10(0.8), log10(2), 101);
searchGrid.beta=[1];
searchGrid.gamma = [0.4:0.01:0.6];  % guess rate / lower asympotote
searchGrid.lambda = [0.01:0.02];  % lapse rate / deviation from upper asymptote being 1

for c=1:3
    StimLevels=ObjData(c,:,1);
    NumPos=ObjData(c,:,2);
    OutOfNum=ObjData(c,:,3);
[paramsValuesobj(c,:), LLobj(c), ~] = PAL_PFML_Fit(StimLevels,NumPos, ...
    OutOfNum,searchGrid,paramsFree,PF);

disp('done:')
message = sprintf('Threshold estimate: %6.4f',paramsValuesobj(c,1));
disp(message);
message = sprintf('Slope estimate: %6.4f\r',paramsValuesobj(c,2));
disp(message);
end





for c=1:3
%Create simple plot
ProportionCorrectObserved=ObjData(c,:,2)./ObjData(c,:,3); 
StimLevelsFineGrain=[min(StimLevels):max(StimLevels)./1000:max(StimLevels)];
ProportionCorrectModel = PF(paramsValuesobj(c,:),StimLevelsFineGrain);
 
figure('name','Maximum Likelihood Psychometric Function Fitting');
axes
hold on
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','color',[0 .7 0],'linewidth',4);
plot(ObjData(c,:,1),ProportionCorrectObserved,'k.','markersize',40);
set(gca, 'fontsize',16);
set(gca, 'Xtick',StimLevels);
axis([min(StimLevels) max(StimLevels) .4 1]);
xlabel('Stimulus Intensity');
ylabel('proportion correct');
end

% Get goodness of fit and throw warning if it is a bad fit

B=1000; %Number of simulations to perform to determine Goodness-of-Fit

disp('Determining Goodness-of-fit.....');

for c=1:3
    ProportionCorrectObserved=ObjData(c,:,2)./ObjData(c,:,3);
    StimLevelsFineGrain=[min(StimLevels):max(StimLevels)./1000:max(StimLevels)];
    ProportionCorrectModel = PF(paramsValuesobj(c,:),StimLevelsFineGrain);
    [Dev pDev] = PAL_PFML_GoodnessOfFit(StimLevels, NumPos, OutOfNum, ...
        paramsValuesobj(c,:), paramsFree, B, PF, 'searchGrid', searchGrid);

    disp('done:');

    %Put summary of results on screen
    fprintf('Deviance: %6.4f\n',Dev);
    if pDev>0.05
        fprintf('p-value: %6.4f \n',pDev);
    else
        fprintf('p-value: %6.4f\n',pDev);
        warning('Suboptimal fit. Consider excluding.')
    end
end
%% Fit Subjective Curves





%Threshold and Slope are free parameters, guess and lapse rate are fixed
paramsFree = [1 1 1 1];  %1: free parameter, 0: fixed parameter
 
%Parameter grid defining parameter space through which to perform a
%brute-force search for values to be used as initial guesses in iterative
%parameter search.
searchGrid.alpha = 0.01:.001:.11;
searchGrid.beta = logspace(0,3,101);
searchGrid.gamma = 0.5;  %scalar here (since fixed) but may be vector
searchGrid.lambda = 0.02;  %ditto

%Perform fit
disp('Fitting function.....');
[paramsValues LL exitflag] = PAL_PFML_Fit(StimLevels,NumPos, ...
    OutOfNum,searchGrid,paramsFree,PF);

disp('done:')
message = sprintf('Threshold estimate: %6.4f',paramsValues(1));
disp(message);
message = sprintf('Slope estimate: %6.4f\r',paramsValues(2));
disp(message);