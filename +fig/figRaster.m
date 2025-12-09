function figRaster(stg, h, d, subjInfo, ds, dp, clust)
    arguments (Input)
        stg % structure, global settings
        h % structure, graphics handles
        d % structure, plot description
        subjInfo % structure, subject info
        ds % structure, data to stem
        dp % structure, data to plot (in "continuous" time sampled in time bins)
        clust % struct with clusters
    end
    arguments (Output)
    end

    nm = d.Name;
    figure(h.f.(nm))

    if ~isfield(h.a, nm)
        h.a.(nm) = axes('Units', stg.units, 'Position', [3, 1, d.PositionCm(3) - 3.5, d.PositionCm(4) - 1.5]);
    end
    
    % Signal OK marker
    [x, y] = fig.getValidXY(subjInfo, d, dp);
    y = y + (stg.numSubj - subjInfo.ksubj)*2 + 0.5;
    h.p.(nm)(subjInfo.ksubj, 1) = plot(x, y, 'Marker', 'none', 'LineWidth', 1.5, 'Color', stg.subjColor(subjInfo.ksubj, :));
    clear x y
    hold on
    
    % Events
    numev = numel(ds.(d.EventName).OnsDt);
    x = days(ds.(d.EventName).OnsDt - subjInfo.dob); % X data common for polynomial fitting and plotting
    x = repelem(x, 3);
    y1 = zeros(numev, 1);
    y = NaN(3*numev, 1);
    y(1 : 3 : 3*numev) = y1 + (stg.numSubj - subjInfo.ksubj)*2 + 0.5;
    y(2 : 3 : 3*numev) = y1 + (stg.numSubj - subjInfo.ksubj)*2 + 1.5;
    y(3 : 3 : 3*numev) = NaN(numev, 1);
    h.p.(nm)(subjInfo.ksubj, 2) = plot(x, y, 'Marker', 'none', 'LineWidth', 0.5, 'Color', stg.subjColor(subjInfo.ksubj, :));
    clear x y
    
    % Clusters
    for kcl = 1 : numel(clust)
        x(1) = days(clust(kcl).OnsDt(1) - subjInfo.dob);
        x(2) = days(clust(kcl).OnsDt(end) - subjInfo.dob);
        yOffset = (clust(kcl).nested)*0.2;
        y = [1 1]*((stg.numSubj - subjInfo.ksubj)*2 + 0.5 + yOffset);
        h.p.(nm)(subjInfo.ksubj, 3, kcl) = plot(x, y, 'Marker', '.', 'MarkerSize', 4, 'LineWidth', 1.5, 'Color', 'k');
    end
    
    xlabel('Age (days)')
    % xlim([45 155]); % Shows all mice complete
    h.a.(nm).XLim(1) = 0;
    h.a.(nm).YAxis.Visible = 'off';
    % text(-0.1*range(h.a.(nm).XLim), (stg.numSubj - subjInfo.ksubj)*2 + 0.5, ['Mouse ', num2str(subjInfo.ksubj)], 'Interpreter', 'none', 'FontWeight', 'bold')
    % text(-0.02*range(h.a.(nm).XLim), (stg.numSubj - subjInfo.ksubj)*2 + 0.5, subjInfo.subjNm, 'Interpreter', 'none', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')
    text(-0.02*(max(h.a.(nm).XLim)-min(h.a.(nm).XLim)) + h.a.(nm).XLim(1), (stg.numSubj - subjInfo.ksubj)*2 + 1, subjInfo.subjNm,...
        'Interpreter', 'none', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Color', stg.subjColor(subjInfo.ksubj, :), 'FontWeight', 'bold');
    h.a.(nm).FontSize = stg.axFontSize;
    h.a.(nm).Layer = 'top';
    h.a.(nm).YLim = [0, (stg.numSubj - 1)*2 + 1.5; ];
    h.a.(nm).Box = 'off';
    h.a.(nm).Units = 'normalized';
end