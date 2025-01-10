function raincloudPlot(data, labels, colors)
% data  -  array with measurement vectors
% labels - array of strings for the group labels
% colors - palette each group 

    numGroups = length(data);
    figure; hold on; % hold on is so that the different components can be overlayed in the same figure

    for i = 1:numGroups
        groupData = data{i};

        % density plot approximated with histogram
        [counts, edges] = histcounts(groupData, 30, 'Normalization', 'pdf'); % adjust bins
        binCenters = (edges(1:end-1) + edges(2:end)) / 2;

        % normalize counts width
        densityWidth = counts / max(counts) * 0.5; % 0.5 is scaling you can adjust depending on width wanted

        % plot the vertical density, you can change this to position the distribution where you want
        fill(i - densityWidth, binCenters, colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

        % scatterplot
        scatter(repmat(i, size(groupData)), groupData, 15, colors(i, :), 'filled','MarkerFaceAlpha', 0.6);
        
        % boxplot
        % summary statistics
        medianVal = median(groupData);
        q1 = prctile(groupData, 25);
        q3 = prctile(groupData, 75);

        % IQR plotted as thick horizontal line
        line([i - 0.1, i + 0.1], [q1, q1], 'Color', 'k', 'LineWidth', 3); % bottom of IQR
        line([i - 0.1, i + 0.1], [q3, q3], 'Color', 'k', 'LineWidth', 3); % top of IQR
        line([i, i], [q1, q3], 'Color', 'k', 'LineWidth', 1); % range connecting line

        % median plotted as filled black circle
        scatter(i, medianVal, 60, 'k', 'filled');
    end

    % plot appearance



    % axis limits based on data ranges so that all group data is fully visible and there is consistent
    % space around the data points for each group.
    ylim([min(cellfun(@min, data)) - 1, max(cellfun(@max, data)) + 1]); % adds a small margin
    xlim([0.5, numGroups + 0.5]); % padding around the grouped x-axis

    % x axis appearance
    % 'xticks'so that each group represented at evenly spaced intervals
    set(gca, 'XTick', 1:numGroups, 'XTickLabel', labels);

    % axis labels
    xtickangle(35); % rotates the label of the axis
    xlabel('xLabel', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('yLabel', 'FontSize', 10, 'FontWeight', 'bold');

    % title
    title('How does this look?', 'FontSize', 12, 'FontWeight', 'bold');

    % grid lines only on the y-axis for easier comparison of data points.
    grid on;
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3); % light, dotted grid

    % background color and box style
    set(gca, 'Color', [0.95, 0.95, 0.95]); % grey background
    box off;

    % you can add legend if you want

    legend_entries = arrayfun(@(i) plot(nan, nan, 'o', 'Color', colors(i, :), 'MarkerFaceColor', colors(i, :)), 1:numGroups, 'UniformOutput', false);
    legend([legend_entries{:}], labels, 'Location', 'bestoutside', 'FontSize', 10);


    hold off; % let's see the magic

end


% % just some random data for two groups
% group1 = randn(100, 1) * 5 + 175; 
% group2 = randn(100, 1) * 6 + 165;
% 
% % combine data into array
% data = {group1, group2};
% labels = {'male', 'female'};
% colors = [0.2 0.6 0.9; 1 0.3 0.2]; % let's be stereotypical: blue-male, red-female
% 
% % see if plot works 
% christina_raincloudPlot_m(data, labels, colors);
