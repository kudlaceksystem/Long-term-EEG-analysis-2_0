function figKaroly(stg, h, d, subjInfo, ds, dp, ~)
    nm = d.Name;
    
    figure(h.f.(nm))
    [spx, spy, spWi, spHe, ~, numc] = fig.getSubplotXYWH(stg, h, d, stg.margGlob, stg.marg);
    h.f.(nm).Units = "centimeters";
    h.a.(nm)(subjInfo.ksubj, 1) = axes("Units", "centimeters", "Position", ...
        [spx(mod(subjInfo.ksubj - 1, numc) + 1), spy(ceil(subjInfo.ksubj/numc)), spWi, spHe], "NextPlot", "add");
    
    % Signal OK marker
    [~, ~, xx] = fig.getValidXY(subjInfo, d, dp);
    dxx = xx(1, 2 : end) - xx(2, 1 : end - 1);
    dropStSub = find(dxx > 1/24);
    drop(1, :) = xx(2, dropStSub);
    drop(2, :) = xx(1, dropStSub + 1);
    clear x y
    
    % Karoly plot proper
    % Prepare sz data
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

    onsD = days(ds.(d.EventName).OnsDt - subjInfo.dob);
    x = rem(onsD, 1)*24;
    y = floor(onsD);
    ymi = min(y);
    yma = max(y);
    % Night time shading
    patch([0 6 6 0], [ymi - 1, ymi - 1, yma + 1, yma + 1], 0.92*[1 1 1], 'EdgeColor', 'none');
    hold on
    patch([18 24 24 18], [ymi - 1, ymi - 1, yma + 1, yma + 1], 0.92*[1 1 1], 'EdgeColor', 'none');
    % % % scatter(x, y, 5*ones(size(x)), 'k', 'Marker', 'o', 'MarkerFaceColor', 'k')
    scatter(x, y, 5*ones(size(x)), 'filled', 'MarkerEdgeColor', stg.subjColor(subjInfo.ksubj, :), 'MarkerFaceColor', stg.subjColor(subjInfo.ksubj, :))
    
    % Dropouts
    for kpa = 1 : size(drop, 2)
        patch([0 24 24 0], [drop(1, kpa), drop(1, kpa), drop(2, kpa), drop(2, kpa)], 'k', 'EdgeColor', 'none')
    end
        h.a.(nm)(subjInfo.ksubj).XLim = [0 24];
        h.a.(nm)(subjInfo.ksubj).YLim = [ymi - 1, yma + 1];
        h.a.(nm)(subjInfo.ksubj).XTick = [0 12 24];
        h.a.(nm)(subjInfo.ksubj).Box = stg.box;
        if subjInfo.ksubj > stg.numSubj - stg.sbNCol
            xlabel('Time of day (hours)')
        end
        if mod(subjInfo.ksubj - 1, stg.sbNCol) == 0
            ylabel('Age (days)')
        end
        title(subjInfo.subjNm, 'Interpreter', 'none', 'Color', 'k', 'FontWeight', 'bold');
        % title(subjInfo.subjNm, 'Interpreter', 'none', 'Color', stg.subjColor(ksubj, :), 'FontWeight', 'bold');
        % title(['Mouse ', num2str(ksubj)], 'Interpreter', 'none', 'Color', 'k', 'FontWeight', 'bold');
        h.a.(nm)(subjInfo.ksubj, 1).FontSize = stg.axFontSize;
        h.a.(nm)(subjInfo.ksubj, 1).Layer = 'top';
end