%% Contrast-Response Simulation
% Gives overallres.error (step size X samples at each point) with average fit error values across participants
% and overallres.saved (step size X samples at each point) with average
% samples saved by cutting from the top
%% Set-Up
sca;
close all;
clear; clc
w = warning ('off','all');

%% Parameters

params.avmidpointsu=0.5; % average midpoint
params.avslopesu=20; % average slope
params.avmidpointob=0.25; % average midpoint
params.avslopeob=20; % average slope
params.nparticipants=50; % how many participants?
params.noiseStd=0; % standard deviation of noise

% Group Stats Parameters

% Sample Points
params.posssamplesteps=[0.18 0.16 0.14 0.12 0.1 0.08]; % between 0 and 1
params.nsamples=[8 10 12 14 16] ; % how often is each point sampled?

% params.posssamplesteps=[0.08]; % between 0 and 1
% params.nsamples=[1000] ; % how often is each point sampled?

params.stopsamplecount= 4; % repetitions of same results before removal of sample point only from above
params.stopabove=0.7; % contrast from which (upwards)sample points can be eliminated 0-1
params.eliminate=1;

% Preallocate
overallres.errorsu=zeros(length(params.posssamplesteps),length(params.nsamples));
overallres.errorob=zeros(length(params.posssamplesteps),length(params.nsamples));
overallres.saved=zeros(length(params.posssamplesteps),length(params.nsamples));
%% Simulate
% Creates two structs under overallres (error and saved) of the size
% (n different steps x n different sample sizes) as an average of
% nparticipants
params.ncontrasts=zeros(1,length(params.nsamples));
tic
for i=1:length(params.posssamplesteps) %cycle through different step sizes
    groupres.errorsu=zeros(1,length(params.nsamples));
    groupres.errorob=zeros(1,length(params.nsamples));
    groupres.saved=zeros(1,length(params.nsamples));
    params.samplestep=params.posssamplesteps(i); % select step size for this round
    params.ncontrasts(i)=length(0:params.samplestep:1);
    for ii=1:length(params.nsamples) % for each step size, cycle through all sample sizes
        [simulationresults,overallsamplepoints,totalresultssu,totalresultsob,subjmean,objmean]=FullSimulationWIP2(params,ii);
        toc

        % get mean across participants for different sample sizes
        groupres.errorsu(ii) = mean(simulationresults(1,:)); % average fit error subj
        groupres.errorob(ii) = mean(simulationresults(2,:)); % average fit error obj
        groupres.saved(ii) = mean(simulationresults(3,:)); % samples saved through elimination
        groupres.stdsu(ii) = std(simulationresults(1,:));
        groupres.stdob(ii) = std(simulationresults(2,:));
        groupres.stdsaved(ii) = std(simulationresults(3,:));
    end

    % % save group results for different sample steps
    overallres.errorsu(i,:)=groupres.errorsu;
    overallres.errorob(i,:)=groupres.errorob;
    overallres.saved(i,:)=groupres.saved;
    overallres.stdsu(i,:)=groupres.stdsu;
    overallres.stdob(i,:)=groupres.stdob;
    overallres.stdsaved(i,:)=groupres.stdsaved;
end
%% Plot
% plots n sample steps on the x axis (fewer to more steps),
%fit error on y axis and different lines represent different sample sizes at each point.

subplot(1,3,1) % Subjective Error
for i=1:length(params.nsamples)
    varplot(params.ncontrasts,overallres.errorsu(:,i),[overallres.errorsu(:,i)-overallres.stdsu(:,i)],[overallres.errorsu(:,i)+overallres.stdsu(:,i)])
    hold on
    leg=string(params.nsamples);
    legend(leg)
    title('Fitting Error Subjective')
end
hold off

subplot(1,3,2) % Objective Error
for i=1:length(params.nsamples)
    %plot(params.ncontrasts,overallres.errorob)
    varplot(params.ncontrasts,overallres.errorob(:,i),[overallres.errorob(:,i)-overallres.stdob(:,i)],[overallres.errorob(:,i)+overallres.stdob(:,i)])
    hold on
    leg=string(params.nsamples);
    legend(leg)
    title('Fitting Error Objective')
end
hold off

subplot(1,3,3) % Saved Samples
plot(params.ncontrasts,overallres.saved,'--')
legend(leg)
title('Trials saved through dynamic removal')
beep


%% 

%%

%% %%%%%%% Simulation Function %%%%%%%

function [simulationresults,overallsamplepoints,totalresultssu,totalresultsob,subjmean,objmean]=FullSimulationWIP2(params,ii)

% % To-Do
% Shuffle participants' midpoints and slopes in a Gaussian way

% Parameters
maxval=1;
% overallsamplepoints=0:params.samplestep:1; % sample points (contrasts)

% % Assaf added - create the conrast space from 0.01 to 0.99, to later add 0 and 1 as fixed values that dont depend on performance
nSteps=round(1/params.samplestep)+1;
overallsamplepoints=linspace(0.01,0.99,nSteps);

% Preallocate
simulationresults=zeros(3,params.nparticipants); % (1,:) is fitting error subjective  / (2,:) fitting error objective / (3,:) is samples saved
partslopessu=zeros(1,params.nparticipants);
partmidsu=zeros(1,params.nparticipants);
partslopesob=zeros(1,params.nparticipants);
partmidob=zeros(1,params.nparticipants);

% Create True Curves for Participants
for iii=1:params.nparticipants
    pslopesu=params.avslopesu+(randi([-5,5]));
    pmidpointsu=params.avmidpointsu+((0.2-(-0.2)).*rand + (-0.2));
    partslopessu(1,iii)=pslopesu;
    partmidsu(1,iii)=pmidpointsu;
    pslopeob=params.avslopesu+(randi([-5,5]));
    pmidpointob=params.avmidpointob+((0.1-(-0.1)).*rand + (-0.1));
    partslopesob(1,iii)=pslopeob;
    partmidob(1,iii)=pmidpointob;
end

% Model Function
logisticFunsu='1./(1+exp(-a*(x-b)))';
logisticFunob='(0.5./(1+exp(-a*(x-b))))+0.5'; %necessary?
%% Obtain Data
for j=1:params.nparticipants
    % Preallocate / Reset
    sampleresults=zeros(2,length(overallsamplepoints));
    totalresultssu=zeros(params.nsamples(ii),length(overallsamplepoints));
    totalresultsob=zeros(params.nsamples(ii),length(overallsamplepoints));
    currsamplecount=zeros(1,length(overallsamplepoints));
    lastsampleresults=zeros(2,length(overallsamplepoints));
    currsamplepoints=overallsamplepoints;

    % Choose Participant slope/midpoint
    suslope=partslopessu(j);
    sumidpoint=partmidsu(j);
    obslope=partslopesob(j);
    obmidpoint=partmidob(j);

    % Collect Participant Data
    for jj=1:params.nsamples(ii)
        for jjj=1:length(overallsamplepoints)
            currpoint=currsamplepoints(jjj);
            % noise inside logistic function - predicted value between 0-1 but bias
%             percentaccuracysu=maxval./(1+exp(-suslope*(currpoint-sumidpoint+params.noiseStd*randn))); % choose accuracy and add noise
%             percentaccuracyob=(0.5./(1+exp(-obslope*(currpoint-obmidpoint+params.noiseStd*randn))))+0.5; %objective curve between 0.5 and 1

            % noise outside logistic function - beyond 0,1 but no bias
            percentaccuracysu=maxval./(1+exp(-suslope*(currpoint-sumidpoint)))+params.noiseStd*randn; % choose accuracy and add noise
            percentaccuracyob=(0.5./(1+exp(-obslope*(currpoint-obmidpoint))))+0.5+params.noiseStd*randn; %objective curve between 0.5 and 1

            % Calculate Subjective Response
            rs=rand;
            if ~isnan(currpoint) % if current sample point is still included, calculate response
                if rs<=percentaccuracysu
                    responsesu=1;
                else
                    responsesu=0;
                end
            elseif isnan(currpoint)
                responsesu=NaN; % if point has already reached max, do not sample
            end
            sampleresults(1,jjj)=responsesu;

            % Calculate Objective Response
            ro=rand; %(should this be done again? or same value as for sub?)
            if ~isnan(currpoint) % if current sample point is still included, calculate response
                if ro<=percentaccuracyob
                    responseob=1;
                else
                    responseob=0;
                end
            elseif isnan(currpoint)
                responseob=NaN; % if point has already reached max, do not sample
            end
            sampleresults(2,jjj)=responseob;

            % Eliminate repeat sample points
            if params.eliminate==1
                if lastsampleresults(1,jjj)==responsesu && lastsampleresults(2,jjj)==responseob && currpoint>params.stopabove
                    currsamplecount(jjj)=currsamplecount(jjj)+1; % increase counter
                    if currsamplecount(jjj)==params.stopsamplecount % and check if max value is reached
                        currsamplepoints(jjj)=NaN; % if maxvalue is reached, remove sample point
                    else
                    end
                else
                    currsamplecount(jjj)=0; % reset counter
                end
            else
            end
        end
        totalresultssu(jj,:)=sampleresults(1,:);
        totalresultsob(jj,:)=sampleresults(2,:);
        lastsampleresults(1,:)=totalresultssu(jj,:);
        lastsampleresults(2,:)=totalresultsob(jj,:);
    end
    % Mean
    subjmean(j,:)=mean(totalresultssu,1,'omitnan'); % % assaf added - document all individual participants to function output
    objmean(j,:)=mean(totalresultsob,1,'omitnan'); % % same

    x=0:0.01:1; % true curve steps


    if params.eliminate==1
        samplessaved=numel(totalresultssu)-(numel(totalresultssu(~isnan(totalresultssu))));
    else
        samplessaved=0;
    end
    %% Fit Model and Plot
    % % Assaf added - condition to only accept fit if its within range
    done=0;
    while ~done
        modelsu=fit([0 overallsamplepoints 1]',[0 subjmean(j,:) 1]',logisticFunsu); % % assaf added - fixed performance at 0 and 1 to 0 and 1
        if modelsu.b>=0 && modelsu.b<=1
            fittingErrorsu=modelsu.b-sumidpoint;
            done=1;
        end
    end

    done=0;
    while ~done
        modelob=fit([0 overallsamplepoints 1]',[0.5 objmean(j,:) 1]',logisticFunob); % % assaf added - fixed performance at 0 and 1 to 0.5 and 1
        if modelob.b>=0 && modelob.b<=1
            fittingErrorob=modelob.b-obmidpoint;
            done=1;
        else
            objmean(j,:) % % assaf added - display in command window when stuck            
        end
    end
    % %

    simulationresults(1,j)=fittingErrorsu;
    simulationresults(2,j)=fittingErrorob;
    simulationresults(3,j)=samplessaved;
    save('insidefunctionvar')
end
end