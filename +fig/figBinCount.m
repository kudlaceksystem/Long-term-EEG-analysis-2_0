function figBinCount(stg, h, d, subjInfo, ds, dp, ~)
    nm = d.Name;
    evnm = d.EventName; % Name of the event to analyze
    vanm = d.EventValidSrc; % Name of the dp field to get the data on validity of the source data
    binlenDu = days(3); % Length of couning bin
    binshiftDu = days(1);
    first = ds.(evnm).OnsDt(1);
    last = ds.(evnm).OnsDt(end);
    numbin = floor((last - first - binlenDu)/binshiftDu) + 1;
    % Initialize the bin edges for the histogram
    binDt = first : binshiftDu : last;
    binCounts = NaN(numbin, 1);
    for kb = 1 : numbin
        dpTax = dp.(vanm).tax; % dpTax contains ends of the dp bins
        dpBinlenDu = mode(diff(dpTax)); % Maybe in future we may allow overlapping windows in dp and this line will need to be adapted.
        dpBinlenS = seconds(dpBinlenDu);
        minRequiredValidS = 0.8*dpBinlenS; % 20% missing is tolerated
        validSt = find(dpTax > binDt(kb), 1, "first"); % Find the first dp bin which has the end later than counting bin starts
        validEn = find(dpTax - dpBinlenDu < binDt(kb+1), 1, "last"); % Find the last dp bin which has the start before the conting bin ends
        valid = dp.(vanm).ValidS(validSt : validEn, :);
        if any(valid > dpBinlenS, "all") % Just check that the validS is shorter than the bin length (otherwise there is a bug somewhere in getData)
            disp(valid)
            disp(dpBinlenS)
            pause
        end
        valid = ~all(isnan(valid) | valid <= minRequiredValidS, 2);
        if all(valid)
            binCounts(kb, :) = sum(ds.(evnm).OnsDt >= binDt(kb) & ds.(evnm).OnsDt < binDt(kb+1));
        else
            binCounts(kb, :) = NaN;
        end
    end

    % Stats
    numbins = numel(binCounts);
    numev = sum(binCounts, 1, "omitmissing");
    
    % Compute histogram
    edges = 0 : 20;
    hc = histcounts(binCounts, edges);
    hcNorm = hc/numev;
    
    % Plot figure
    figure(h.f.(nm))
    [spx, spy, spWi, spHe, ~, numc] = fig.getSubplotXYWH(stg, h, d, stg.margGlob, stg.marg);
    h.f.(nm).Units = "centimeters";
    h.a.(nm)(subjInfo.ksubj, 1) = axes("Units", "centimeters", "Position", ...
        [spx(mod(subjInfo.ksubj - 1, numc) + 1), spy(ceil(subjInfo.ksubj/numc)), spWi, spHe], "NextPlot", "add");
    
    x = edges; % X data common for polynomial fitting and plotting
    x = repelem(x, 3);
    x = x(2 : end-1);
    y = NaN(size(x));
    y(1 : 3 : end) = 0;
    y(2 : 3 : end) = hcNorm;
    y(3 : 3 : end) = hcNorm;
    facecolor = 1 - 0.2*(1 - stg.subjColor(subjInfo.ksubj, :));
    h.p.(nm)(subjInfo.ksubj, 1) = patch(x, y, facecolor, 'LineWidth', 0.5, 'EdgeColor', 'k');
    clear x y

    % % % onsD = days(ds.(evnm).OnsDt - subjInfo.dob);
    % % % x = rem(onsD, 1)*24;
    % % % y = floor(onsD);
    % % % ymi = min(y);
    % % % yma = max(y);
    % % % % Night time shading
    % % % patch([0 6 6 0], [ymi - 1, ymi - 1, yma + 1, yma + 1], 0.92*[1 1 1], 'EdgeColor', 'none');
    % % % hold on
    % % % patch([18 24 24 18], [ymi - 1, ymi - 1, yma + 1, yma + 1], 0.92*[1 1 1], 'EdgeColor', 'none');
    % % % % % % scatter(x, y, 5*ones(size(x)), 'k', 'Marker', 'o', 'MarkerFaceColor', 'k')
    % % % scatter(x, y, 5*ones(size(x)), 'filled', 'MarkerEdgeColor', stg.subjColor(subjInfo.ksubj, :), 'MarkerFaceColor', stg.subjColor(subjInfo.ksubj, :))
    
    % % % % % Dropouts
    % % % % for kpa = 1 : size(drop, 2)
    % % % %     patch([0 24 24 0], [drop(1, kpa), drop(1, kpa), drop(2, kpa), drop(2, kpa)], 'k', 'EdgeColor', 'none')
    % % % % end
    % % % %     h.a.(nm)(subjInfo.ksubj).XLim = [0 24];
    % % % %     h.a.(nm)(subjInfo.ksubj).YLim = [ymi - 1, yma + 1];
    % % % %     h.a.(nm)(subjInfo.ksubj).XTick = [0 12 24];
    % % % %     h.a.(nm)(subjInfo.ksubj).Box = stg.box;
    % % % %     if subjInfo.ksubj > stg.numSubj - stg.sbNCol
    % % % %         xlabel('Time of day (hours)')
    % % % %     end
    % % % %     if mod(subjInfo.ksubj - 1, stg.sbNCol) == 0
    % % % %         ylabel('Age (days)')
    % % % %     end
    % % % %     title(subjInfo.subjNm, 'Interpreter', 'none', 'Color', 'k', 'FontWeight', 'bold');
    % % % %     % title(subjInfo.subjNm, 'Interpreter', 'none', 'Color', stg.subjColor(ksubj, :), 'FontWeight', 'bold');
    % % % %     % title(['Mouse ', num2str(ksubj)], 'Interpreter', 'none', 'Color', 'k', 'FontWeight', 'bold');
    % % % %     h.a.(nm)(subjInfo.ksubj, 1).FontSize = stg.axFontSize;
    % % % %     h.a.(nm)(subjInfo.ksubj, 1).Layer = 'top';
end