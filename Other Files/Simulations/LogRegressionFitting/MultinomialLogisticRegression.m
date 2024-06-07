% Multiple Logistic Regression Simulation
clear 
clc
b=[0,1,0]; % chosen weights (intercept, contrast, a and interaction c*x
% o=1;
% m=5;
contrastData=sort(repmat([1:10]',10,1)); % controlled facor (contrast levels)
aData=10*randn(100,1); % random factor (e.g. spontaneous alpha)
%logodds=1./(1+exp(b(1)+b(2)*contrastData+b(3)*aData)); % Calculated performance based on factors and
%pC=log(logodds);
pC=1./(exp(b(1)+b(2)*contrastData+b(3)*aData)); % Calculated performance based on factors and 
pC=exp(b(1)+b(2)*contrastData)/1+(b(1)+b(2)*contrastData); % Calculated performance based on factors and 
%pC=(exp(-(contrasts-m)/o))./(o*(1+exp(-(contrasts-m)/o)).^2);
% pC=1./(1+exp(-b(1)*(contrasts-m)));
% pC=1./(1+exp(-o*(contrasts-m)));
% test=1./(1+exp(-b(1)-(b(2)*contrastData)));

response=pC>0.5;

plot(log(pC))

% Data Table (Response, contrast, x value)
xdata=[contrastData aData];
ydata=categorical(response);
datmat=[contrastData aData response];

[B,dev,stats] = mnrfit(xdata,ydata);
fitglm(xdata,ydata)

plot(contrastData,pC)

% Non-log
exp(B)
round(stats.p,5)

% FitMNR
[B,dev,stats] = fitmnr(xdata,ydata);