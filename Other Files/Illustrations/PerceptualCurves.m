x=0:0.01:1;
midpoint=0.7;
maxval=1;
slope=20;
%


% generate true curve
truecurve(1,:)=maxval./(1+exp(-slope*(x-midpoint)));
% plot(x,truecurve,'LineWidth',2)

hold on

midpoint=0.6;
truecurve(2,:)=maxval./(1+exp(-slope*(x-midpoint)));
% plot(x,truecurve,'g','LineWidth',2)

midpoint=0.5;
truecurve(3,:)=maxval./(1+exp(-slope*(x-midpoint)));

midpoint=0.65;
truecurve(4,:)=maxval./(1+exp(-slope*(x-midpoint)));

% plot(x,truecurve,'g--','LineWidth',2)
% 
% xlabel('Stimulus Intensity')
% ylabel('Perception')
createfigure(x, truecurve)

function createfigure(X1, YMatrix1)
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',2,'Parent',axes1);
set(plot1(1),'Color',[0 1 0]);
set(plot1(2),'LineStyle','--','Color',[0.301960784313725 0.745098039215686 0.933333333333333]);

% Create ylabel
ylabel('Perception','FontSize',18);

% Create xlabel
xlabel('Stimulus Intensity','FontSize',18);

hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'LineWidth',2,'XTick',[0 1],'YTick',[0 1]);
end