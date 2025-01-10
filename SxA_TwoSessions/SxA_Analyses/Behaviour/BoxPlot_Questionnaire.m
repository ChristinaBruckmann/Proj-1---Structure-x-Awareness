% Figure for how well people used the temporal structure (0-100)

% Data imported from qualtrics
interval_survey = [80 81 81 71 72 69, 72, 60, 65, 71, 43, 81, 72, 61, 52, 45, 75, 58, 63, 40, 80, 75, 60, 35, 55, 70, 70, 45, 30, 30, 50, 49, 91, 62, 71, 59, 73, 64, 90, 70, 70, 43, 68, 73, 91, 71, 39, 62, 91, 81, 45, 40];
rhythm_survey = [90, 70, 100, 94, 65, 66, 91, 81, 81, 78, 72, 81, 100, 77, 75, 59, 60, 72, 75, 69, 90, 90, 45, 30, 65, 70, 62, 56, 70, 62, 50, 59, 83, 74, 90, 91, 100, 58, 90, 60, 87, 52, 71, 64, 100, 92, 63, 72, 91, 81, 57, 60];

% Merge
survey_means = [mean(rhythm_survey) mean(interval_survey)];
figure("Position",[763.4000  232.2000  723.2000  541.6000]);
% Plot
b = bar(survey_means, 'Facecolor', 'flat', 'FaceAlpha', 0.7,'EdgeColor', 'none'); 
hold on

xticklabels({'Rhythm', 'Interval'})
ylabel('Rating (0-100)','FontSize', 14)
ylim([0 105])
title('How well did you manage to use the...?')

% Reshape boxplot data into a 2-column format (rows: participants, columns: conditions)
boxplotdata = [rhythm_survey', interval_survey'];

% Plot each participant's data points connected by lines
for i = 1:size(boxplotdata, 1)
    if boxplotdata(i, 1) > boxplotdata(i, 2)
        plot([1, 2], [boxplotdata(i, 1), boxplotdata(i, 2)], '--ro', 'LineWidth', 0.7, 'Color', 'k')
    else
        plot([1, 2], [boxplotdata(i, 1), boxplotdata(i, 2)], '-ro', 'LineWidth', 0.7, 'Color', 'k')
    end
end

b.CData(1, :) = [0.00, 0.45, 0.74];
b.CData(2, :) = [0.85, 0.33, 0.10];

% Make plot pretty
box off
yticks(0:20:100)
ax = gca;
ax.FontSize = 14;