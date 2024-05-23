% Figure for how well people used the temporal structure (0-100)

% Data imported from qualtrics
% interval_survey=[71    71    60    92    19    30    71    60    92    60    29     8    80    81    72    46 72    81    69    65    71    71    60];
% rhythm_survey=[57    62   100    67    40    50    71    70    73    70    52    65    90   100    65    77 91    70    66    81    94    78    81];

interval_survey=table2array(StructurexAwareness22(3:end,2));
rhythm_survey=table2array(StructurexAwareness22(3:end,3));

% Merge
survey_means=[mean(rhythm_survey) mean(interval_survey)];
figure;
%Plot
b=bar(survey_means,'Facecolor', 'flat');
hold on
% plot(1,rhythm_survey,'ob')
% plot(2,interval_survey,'or')
xticklabels({'Rhythm', 'Interval'})
ylabel('Rating (0-100)')
ylim([0 105])
title('How well did you manage to use the...?')
% General idea
% matrix is the matrix of values you want to plot
% rows are subjects, columns conditions
boxplotdata=[rhythm_survey interval_survey];
hold on


% poi i singoli soggetti

for i=1:size(boxplotdata,1)
    if boxplotdata(i,1) > boxplotdata(i,2)
        plot([1,2],[boxplotdata(i,1),boxplotdata(i,2)],'--ro', "LineWidth", 0.7, "Color", 'k')
    else
        plot([1,2],[boxplotdata(i,1),boxplotdata(i,2)],'-ro', "LineWidth", 0.7, "Color", 'k')
    end
end

b.CData(1,:) = [0.00,0.45,0.74];
b.CData(2,:) = [0.85,0.33,0.10];
