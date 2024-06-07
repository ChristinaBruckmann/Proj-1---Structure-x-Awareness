%% Comparing different logistic regression algorithms
% Christina Bruckmann 05.03.24
% This code simulates data based on logistic regression (with either one or multiple predictors)
% and compares how different fitting procedures handle this data

% Needs CompareLogFunc_MultPred.R and CompareLogFunc_OnePred.R

clear
clc
close all

simulate=1; % Simulate or use real data?
multpred=0; % Multiple predictors?

%% Load or Simulate Data
if simulate
    %% Simulate Data

    % Contrast values
    contrastVals=1:10;
    contrastThreshold=5;
    contrastData=sort(repmat(contrastVals',10,1)); 

    % Alpha values
    aDataVals=-5:4;
    aData=repmat(aDataVals', 10,1);

    % Other Parameters 
    simb=[-5,1,0]; % chosen weights (intercept, contrast, alpha) - to get a reasonable curve, put intercept at -5 and contast at 1
    lapserate=0; % lowers the 9th point by this value to check the impact of lapses

    % pC=1./(1+exp(-1*(simb(1)+simb(2)*(contrastVals'-contrastThreshold)))); % Old way of calculating, doesn't work well with some fitting tools
    if multpred
        pC=1./(1+exp(-1*(simb(1)+simb(2)*contrastData+simb(3)*aData)));
    else
        pC=1./(1+exp(-1*(simb(1)+simb(2)*contrastVals')));
    end

    pC=round(pC,2); % Rounding here because some methods use trials as input and some percent correct - like this it can be used for both without rounding differences
    pC(9)=pC(9)-lapserate;

    % Show data
    if multpred
        figure; surf(contrastVals,aDataVals,reshape(pC,10,10)); 
    else
        figure; plot(contrastVals,pC,'.','MarkerSize',20);
    end

    % Format data for Psignifit (3 rows: stim levels, percent correct, total trials)
    if multpred
        curridx=1;
        dataforpsignifit(:,1)=contrastVals';
        for con=contrastVals
            dataforpsignifit(con,2)=mean(pC(curridx:curridx+9)*100);
            curridx=curridx+10;
        end
        dataforpsignifit(:,3)=100;
    else
        dataforpsignifit=[contrastVals',pC.*100];
        dataforpsignifit(:,3)=100;
    end

    % Format Data for R (needs it in trial-format)
    % For now 100 trials per contrast (and alpha level), not realistic
    % Did this, so that the pC is exactly what I feed the other functions without rounding errors.
    if multpred
        rdata1=[contrastData aData pC];
        rdata=[];
        for ntrials=1:length(pC)
            currrow=rdata1(ntrials,:,:); % select one row
            temp=repmat(currrow,100,1); % copy row 10 times (starting is 100 trials, output will be 1000)
            temp(:,3)=0; % set default to 0
            ncorrect=round((currrow(3)*100));
            temp(1:ncorrect,3)=1; % add as many correct responses as required by the pC
            rdata=[rdata; temp]; % merge (contrast, alpha, binary response)
        end
    else
        start=1;
        trialspercon=100;
        rdata=zeros(trialspercon*length(contrastVals),2);
        for con=contrastVals
            rdata(start:(start+trialspercon-1),1)=con;
            trials=zeros(trialspercon,1);
            trials(1:(pC(con)*100))=1;
            rdata(start:(start+trialspercon-1),2)=trials;
            start=start+trialspercon;
        end
    end

else
    %% Load Data (not yet adjusted for multiple predictors)
    s=input('Subject? ');
    loadname=sprintf('SxA_ResultsSubject%i_Total.mat',s);
    load(loadname,'dataforpsignifit','alldataclean')

    % Subj Rhyhm only
    dataforpsignifit=dataforpsignifit{1, 1}.subj;

    % For Assaf's code
    pC= dataforpsignifit(:,2)./dataforpsignifit(:,3);  % percent correct for Assaf's code
    contrastVals = dataforpsignifit(:,1)'; % contrast vals for Assaf's code

    % Data for R (all trials)
    rhyidx=find(alldataclean{:,"Condition"}==1);
    rdata=table2array([alldataclean(rhyidx,"Contrast Level"),alldataclean(rhyidx,"Binary Visibility")]);
end

%% %%%% Run Different Fitting Methods %%%%%%%

if ~multpred % psignifit only works for 1 predictor, still need to check how to adjust assaf's code for 2
    %% Psignifit
    % Prepare Data for Psignifit [x-value, number correct, number of trials]
    options=struct;
    options.expType="YesNo";
    options.sigmoidName='logistic';
    % Run Psignifit
    psignires=psignifit(dataforpsignifit,options);
    plotPsych(psignires);
    psignioutput=[psignires.Fit(1), getSlopePC(psignires, 0.5)];
    plotsModelfit(psignires)
    %% Assaf's Code
    [bestFitParams, resFitSorted]=fitPsychometric(contrastVals', pC, 2);
    assafoutput=bestFitParams(1:2);
end

%% R
save('SimDataforR','rdata')
%Rinit('R.matlab','C:\Program Files\R\R-4.1.3\bin\R.exe',"C:\Program Files\R\R-4.1.3\library")
rdone=0;
while ~rdone
    rdone=input('Executed R Script? ');
end
% Run code in R (ask IT to add the R integration!!)

% After Running code in R, come back here.
Routput= readmatrix('R_output.txt');

% Convert (b0+b1x or s(x-t)? b0=-st, b1=s)
Routput(1)=-Routput(2)*Routput(1);
%% Matlab Built-In

if multpred
    nominalResp=categorical(rdata(:,3)+1); % needed as this function doesnt take 0 as category value
    [b,~,stats]=mnrfit(rdata(:,1:2), nominalResp);
    matlaboutput=[b(1),b(2),b(3)]; %
else
    nominalResp=categorical(rdata(:,2)+1); % needed as this function doesnt take 0 as category value
    [b,~,stats]=mnrfit(rdata(:,1), nominalResp);
    matlaboutput=[b(1),b(2)]; %
end

%% Palamedes
% Prepare data for Palamedes (struct with three fields (x,y,n)
dataforpala=struct;
dataforpala.x=dataforpsignifit(:,1)';
dataforpala.y=dataforpsignifit(:,2)';
dataforpala.n=dataforpsignifit(:,3)';

PF=@PAL_Logistic;
paramsFree = [1 1 0 1]; % All parameters are free, except for guess rate

% Need to define reasonable areas where the function should search
searchGrid.alpha = [0:.01:9];    %threshold between 0 an 10
searchGrid.beta = 10.^[-1:.01:2]; %slope
searchGrid.gamma = [0:.01:.06]; % guess rate (lower boundary I think)
searchGrid.lambda = [0:.01:.06]; % lapse rate (upper boundary I think)

palares=PAL_PFML_Fit(dataforpala.x, dataforpala.y,dataforpala.n,searchGrid,paramsFree,PF); % fit
ProportionCorrectObserved=dataforpala.y./dataforpala.n; 
[Dev, pDev, DevSim, converged] = PAL_PFML_GoodnessOfFit(dataforpala.x, dataforpala.y,dataforpala.n, palares, paramsFree, 1000, PF);
StimLevelsFineGrain=[min(dataforpala.x):max(dataforpala.x)./1000:max(dataforpala.x)];
ProportionCorrectModel = PF(palares,StimLevelsFineGrain);
 
figure('name','Maximum Likelihood Psychometric Function Fitting');
axes
hold on
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',4);
plot(dataforpala.x,ProportionCorrectObserved,'k.','markersize',40);
set(gca, 'fontsize',16);
set(gca, 'Xtick',dataforpala.x);
axis([min(dataforpala.x) max(dataforpala.x) 0 1]);
xlabel('Stimulus Intensity');
ylabel('proportion correct');

%% Plot All
figure;
if ~multpred
    alloutput=[psignioutput;assafoutput;Routput';abs(matlaboutput);palares(1:2)]; 
    plot(1:2,alloutput)
    legend('Psigni','Assaf','R','Matlab','Palamedes')
    xticklabels({'Threshold','Slope'})
    xlim([0,3])
    xticks([1,2])
    if simulate
        hold on
        plot(1:2,[abs(simb(1)),simb(2)],'.','MarkerSize',20)
        legend('Psigni','Assaf','R','Matlab','Palamedes','True Value of Simulation')
    end
else
    alloutput=[Routput';matlaboutput]; 
    plot(1:3,alloutput)
    legend('R','Matlab')
    xticklabels({'Threshold','Contrast','Alpha'})
    xlim([0,4])
    xticks([1,2,3])
    if simulate
        hold on
        plot(1:3,[simb(1),simb(2),simb(3)],'.','MarkerSize',20)
        legend('R','Matlab','True Value of Simulation')
    end
end