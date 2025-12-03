function [clust, eventBelongsToClust, stats] = getClusters(subjInfo, ds, dp, cld, ksubj)
    arguments (Input)
        subjInfo
        ds % Data to stem, typically, this will contain the phenomena to cluster
        dp % Data to plot, I do not know now, why we need it here :-)
        cld % Cluster description - one line of the clDesc structure
        ksubj % Which subject we are analyzing now
    end
    arguments (Output)
        clust % Structure with the cluster data (TODO003 POSSIBLY CHANGE IT TO A TABLE)
        eventBelongsToClust % Column vector indicating for each event how many clusters it belongs to.
        stats % One-row table with statistics, can be appended to a table containing all subjects
    end

    % global stg
    onsDt = ds.(cld.EventName).OnsDt; % We always cluster by onset datetime here
    clust = [];
    stats = table;
    if numel(onsDt) < cld.MinNumInClus
        eventBelongsToClust = zeros(size(onsDt));
        stats.clNumClust = 0;
        stats.clNumClustNonNested = 0;
        stats.clFracIntraClustSz = 0;
        stats.clNumSzPerClust = NaN;
        stats.clClusterDuration = NaN;
        stats.clInterclusPeriod = NaN;
        % % % % % stats.clInterclusDiffDur = NaN;
        % % % % % stats.clInterclusDiffRac = NaN;
        % % % % % stats.clInterclusDiffPow = NaN;
        % % % % % stats.clInterclusDiffPos = NaN;
        return
    end
    kc = 0; % Cluster number
    if cld.ExclClAtEdges % Exclude clusters at edges of the recording?
        onsDt = [subjInfo.anStartDt; onsDt; subjInfo.anEndDt]; % Adding dummy events at the analysis start and end
    else
        onsDt = [0; onsDt; years(3000)]; % Adding dummy events at the year 0 and year 3000
    end
    eventBelongsToClust = zeros(size(onsDt)); % Does the event belong to any cluster?
    intMult = cld.InterclusterMultiplier; % Intercluster period must be intMult times longer than the longest ISI within the cluster
    minNumEv = cld.MinNumInClus; % Minimum number of events in the cluster
    numev = length(onsDt); % Number of events (including the artificially added events at the beginning and end of the onsDt vector
    evGrSt = 1; % Event group start index.
    while evGrSt <= numev - minNumEv
        evGrDt = onsDt(evGrSt : evGrSt + minNumEv); % Event group. At this moment, the shortest which could be considered a cluster. Due to the added dummy event at the beginning of OnsDt, the group index needs to be higher by 1 than it would have to be in ds.(*).OnsDt to point at the same event.
        ieiGrDu = diff(evGrDt); % Inter-event intervals (IEIs) within the group in duration class
        % If the first IEI in the group is too short to be intercluster interval (ICI) or any of the within-cluster IEI is > stg.MaxWithinClusIeiD,
        % increase the evGrSt and try again with another group. IMPORTANT: the dummy events at the start and end can make a big difference.
        if ieiGrDu(1) < intMult*max(ieiGrDu(2 : end)) || any(ieiGrDu(2 : end) > cld.MaxWithinClusIeiD)
            evGrSt = evGrSt + 1; % Try another group starting on the next event
        else % The first IEI in the group is long enough to be considered as ICI. Let's try if there is an ICI also after the event group
            addEvTF = true; % Continue adding more events?
            numAddEv = 1; % Number of events to be added to the original group. It increases later.
            while addEvTF
                evGrDt = onsDt(evGrSt : evGrSt + minNumEv + numAddEv); % Add an event
                ieiGrDu = diff(evGrDt); % Get IEI of the, now bigger, event group
                % If the within-cluster IEIs are not too long (the first IEI is still considered ICI) AND last IEI is not too long compared to within-cluster IEI
                % so it does not terminate the current cluster, add more events.
                if intMult*max(ieiGrDu(2 : end)) <= ieiGrDu(1) && ieiGrDu(end) < intMult*max(ieiGrDu(2 : end-1))
                    numAddEv = numAddEv + 1; % Add one more event
                    if evGrSt + minNumEv + numAddEv > numev % But if there is no more event, stop adding events and check if we have a cluster
                        evGrDt7linesAbove = evGrDt;
                        evGrDt = onsDt(evGrSt : evGrSt + minNumEv + numAddEv - 1); % Add an event
                        if evGrDt7linesAbove ~= evGrDt; error('_jk Logical error in the code. evGrDt7linesAbove should be the same as evGrDt. Solve it.'); end
                        ieiGrDu = diff(evGrDt); % Get IEI
                        if ieiGrDu(1) > intMult*max(ieiGrDu(2 : end - 1)) && ieiGrDu(end) > intMult*max(ieiGrDu(2 : end - 1))
                            kc = kc + 1;
                            eventBelongsToClust(evGrSt + 1 : evGrSt + minNumEv + numAddEv - 1) = eventBelongsToClust(evGrSt + 1 : evGrSt + minNumEv + numAddEv - 1) + 1; ...
                                % Increase the number by 1. The number will indicate level of nestedness of the cluster to which the event belongs.
                            clust(kc).Subscript = (evGrSt + 1 : evGrSt + minNumEv + numAddEv - 1) - 1; %#ok<AGROW>
                            clust(kc).OnsDt = onsDt(evGrSt + 1 : evGrSt + minNumEv + numAddEv - 1); ...
                                %#ok<AGROW> % Cluster begins after the ICI period and ends by the beginning of the ICI
                            clust(kc).ds = ds.(cld.EventName)(evGrSt : evGrSt + minNumEv + numAddEv - 2, :); ...
                                %#ok<AGROW> % Note, that here we are subscripting into the original table of events with no dummy event at the beginning.
                            clust(kc).ksubj = ksubj; %#ok<AGROW>
                            clust(kc).subjclustn = kc; %#ok<AGROW>
                            clust(kc).subjNm = subjInfo.subject; %#ok<AGROW>
                            clust(kc).anStartDt = subjInfo.anStartDt; %#ok<AGROW>
                            clust(kc).anEndDt = subjInfo.anEndDt; %#ok<AGROW>
                        else % Remove this when debugging is finished
                            if cld.ExclClAtEdges
                                disp('Last cluster not separated from the anEndDt enough. The cluster not added.')
                            end
                        end
                        evGrSt = numev; % Here we could maybe plug in e.g. numev to stop the outer loop immediatelly
                        addEvTF = false; % Terminate the inner while loop, i.e. stop adding events to this group
                    end
                else % The added event is after too long IEI so it either just disqualifies the first IEI to define the start of the cluster or it is long enough to terminate the cluster.
                    if ieiGrDu(end) < intMult*max(ieiGrDu(2 : end-1)) % _EITHER_ the last IEI is not long enough to be considered terminating ICI. So this group will never be a cluster.
                        addEvTF = false; % Terminate the inner while loop, i.e. stop adding events to this group.
                        evGrSt = evGrSt + 1; % Let's try with a new event group.
                    else % _OR_ it is long enough. So events (2 : end-1) of the group form a cluster separated by at least intMult*max(intraClusterISI) on both sides
                        kc = kc + 1;
                        eventBelongsToClust(evGrSt + 1 : evGrSt + minNumEv + numAddEv - 1) = eventBelongsToClust(evGrSt + 1 : evGrSt + minNumEv + numAddEv - 1) + 1;
                        clust(kc).Subscript = (evGrSt + 1 : evGrSt + minNumEv + numAddEv - 1) - 1; %#ok<AGROW>
                        clust(kc).OnsDt = onsDt(evGrSt + 1 : evGrSt + minNumEv + numAddEv - 1); ...
                            %#ok<AGROW> % Cluster begins after the ICI period and ends by the beginning of the ICI. onsN has a dummy at the beginning.
                        clust(kc).ds = ds.(cld.EventName)(evGrSt : evGrSt + minNumEv + numAddEv - 2, :); ...
                            %#ok<AGROW> % Does not have the dummy at the beginning zero, hence the different subscripts
                        clust(kc).ksubj = ksubj; %#ok<AGROW>
                        clust(kc).subjclustn = kc; %#ok<AGROW>
                        clust(kc).subjNm = subjInfo.subjNm; %#ok<AGROW>
                        clust(kc).anStartDt = subjInfo.anStartDt; %#ok<AGROW>
                        clust(kc).anEndDt = subjInfo.anEndDt; %#ok<AGROW>
                        if evGrSt + minNumEv + numAddEv < numev % If more events exist
                            numAddEv = numAddEv + 1; % Add one more events to the group so if there are multiple clusters starting with the same event, they get detected.
                        else
                            evGrSt = evGrSt + 1;
                            addEvTF = false; % Terminate the inner while loop, i.e. stop adding events to this group.
                        end
                    end
                end
            end
        end
    end
    eventBelongsToClust = eventBelongsToClust(2 : end-1);
    onsDtCorrespondTF = all(clust(kc).ds.OnsDt == clust(kc).OnsDt);
    if ~onsDtCorrespondTF
        error('_jk OnsDt directly in clust(kc) and in clust(kc).ds do not correspond.')
    end
    
    % Remove too long clusters
    clusterTooLongTF = false(numel(clust), 1);
    for kc = 1 : length(clust)
        clusterTooLongTF(kc) = clust(kc).OnsDt(end) - clust(kc).OnsDt(1) > cld.MaxClusterDur;
        if clusterTooLongTF
            eventBelongsToClust(onsDt == clust(kc).OnsDt) = eventBelongsToClust(onsDt == clust(kc).OnsDt) - 1;
        end
    end
    clust = clust(~clusterTooLongTF);
    
    % Find significant dropouts
    tax = dp.tax;
    binLenDu = mode(diff(tax));
    significanceThresholdDu = max([seconds(1/24), 1.01*binLenDu, min(diff(ds.(cld.EventName).OnsDt))]); % The 1.01 constant is to accommodate tiny inaccuracies in the analysis block durations due to rounding errors.
    taxN = (1 : numel(tax))'; % Number the time points
    taxTbl = table(tax, taxN); % tax and taxN are different class so put them in a table
    notNanInd = any(~isnan(dp.(cld.EventValidSrc).ValidS), 2); % At least in one channel, there has to be non-NaN value - then the time point is considered to have valid data
    taxNotNanTbl = taxTbl(notNanInd, :); % Remove the rows corresponding to dropouts (marked by siCharTbl.sz == NaN)
    taxNotNanDiff = diff(taxNotNanTbl.tax); % Where there were NaNs in the siCharTbl.sz, the tax values were removed so there is a long gap, thus high diff.
    signDropOnsSubInTaxNotNan = find(taxNotNanDiff > significanceThresholdDu); % Subscripts into the taxNotNanTbl
    signDropOnsSubInTax =  taxNotNanTbl.taxN(signDropOnsSubInTaxNotNan); % Extract the number added to tax 4 rows above
    signDropOffSubInTax =  taxNotNanTbl.taxN(signDropOnsSubInTaxNotNan + 1); % Same for the first row after the dropout
    signDropOnsDt = tax(signDropOnsSubInTax); % The extracted number is the index into siCharTbl.tax. siCharTbl.tax contains the ends of analysis blocks so this is the end of the last correct block.
    signDropOffDt = tax(signDropOffSubInTax - 1); % - 1 so that it is the end of the last corrupted (dropout) block
    
    % Remove clusters less than the required multiple of within-cluster ISI from the recording limits or signal dropouts
    clusterRemoveTF = false(size(clust));
    for kc = 1 : numel(clust)
        clOnsDt = clust(kc).OnsDt; % Onsets of event within this cluster
        clSep = intMult*max(diff(clOnsDt)); % Minimum required separation of the cluster from other events
        if (clust(kc).OnsDt(1) - subjInfo.anStartDt) < clSep || ... % If the cluster begins too soon after the recording start (which should not happen thanks to the dummies above) OR
                (subjInfo.anEndDt - clust(kc).OnsDt(end)) < clSep % the recording ends too soon after the cluster end
            clusterRemoveTF(kc) = true;
            warning('_jk Cluster at the beginning or end of recording detected and removed.')
            pause
        end
        for kdrop = 1 : numel(signDropOnsDt)
            if (clOnsDt(1) - signDropOffDt(kdrop)) < clSep   &&   signDropOnsDt(kdrop) - clOnsDt(end) < clSep % Too soon after dropout end OR dropout onset too soon after cluster end (or one or both differences are even negative which indicates the dropout even reaches or is contained inside the cluster)
                clusterRemoveTF(kc) = true;
            end
        end
    end
    clusterRemoveSub = find(clusterRemoveTF);
    for kc = 1 : numel(clusterRemoveSub)
        eventBelongsToClust(clust(clusterRemoveSub(kc)).Subscript) = eventBelongsToClust(clust(clusterRemoveSub(kc)).Subscript) - 1;
    end
    clust = clust(~clusterRemoveTF);
    
    %% Nestedness
    clNesTbl =  ...
        table('Size', [0, 5], 'VariableNames', ["Edge", "OnOff", "ClSub", "Dur", "Nes"], ... % Edge (i.e. onset or offset of a cluster), 1 for onset and -1 for offset, subscript, duration of the cluster, nestedness
        'VariableTypes', ["datetime", "double", "double", "duration", "double"]); % Using doubles is not computationally optimal but it makes the rest of the code easy to write the difference in the performance is negligible
    for k = 1 : length(clust)
        clNesTbl.Edge(2*(k-1) + 1) = clust(k).OnsDt(1); % Onset of the first seizure in the cluster
        clNesTbl.OnOff(2*(k-1) + 1) = 1; % Polarity of the edge (1 for onset, -1 for offset)
        clNesTbl.ClSub(2*(k-1) + 1) = k; % Original subscript of the cluster
        clNesTbl.Dur(2*(k-1) + 1, 4) = clust(k).OnsDt(end) - clust(k).OnsDt(1); % Duration of the cluster
        clNesTbl.Edge(2*(k-1) + 2) = clust(k).OnsDt(end); % Offset of the cluster
        clNesTbl.OnOff(2*(k-1) + 2) = -1; % Polarity of the edge (1 for onset, -1 for offset)
        clNesTbl.ClSub(2*(k-1) + 2) = k; % Original subscript of the cluster
        clNesTbl.Dur(2*(k-1) + 2, 4) = clust(k).OnsDt(end) - clust(k).OnsDt(1);
    end
    clNesTbl = sortrows(clNesTbl, ["Edge", "Dur"], {'ascend', 'descend'});
    clNesTbl.Nes = cumsum(clNesTbl.OnOff);
    nes = [clNesTbl.ClSub(clNesTbl.OnOff == 1), clNesTbl.Nes(clNesTbl.OnOff == 1)]; % Take only onsets, and only the subscript into cluster structure and the nestedness value
    for kc = 1 : length(clust)
        clust(nes(kc, 1)).nested = nes(kc, 2);
    end

    % Check if the nestedness is calculated correctly
    eventBelongsToClust2 = zeros(size(eventBelongsToClust));
    for kc = 1 : length(clust)
        eventBelongsToClust2(clust(kc).Subscript) = max(eventBelongsToClust2(clust(kc).Subscript), clust(kc).nested);
    end
    if ~all(eventBelongsToClust == eventBelongsToClust2)
        error('_jk Nestedness calculation mismatch detected.');
    end

    %% Statistics
    stats = table;
    if isempty(clust)
        stats.clNumClust = 0; % Total number of cluster of this subject
        stats.clNumClustNonNested = 0; % Number of non-nested clusters (with nestedness == 1)
        stats.clFracWithinClustEv = 0; % Fraction of events occurring within a cluster (nestedness >= 1)
        stats.clNumEvPerClust = NaN; % Mean number of events per cluster
        stats.clClusterDurDu = NaN; % Mean cluster duration in duration
        stats.clBetwClusPerDu = NaN; % Mean between-cluster period in duration
        % % % stats.clBetwClusDiffDur = NaN; % Mean of differences of the duration of the last event in the previous cluster and first in the next IMPLEMENT IT IN ANOTHER FUNCTION
        % % % stats.clBetwClusDiffRac = NaN;
        % % % stats.clBetwClusDiffPow = NaN;
        % % % stats.clBetwClusDiffPos = NaN;
        return
    end
    clEvOnsDt = [];
    betwClusPeriod = duration(NaN, NaN, NaN); % NaN duration
    % % % betwClusDiffDur = NaN;
    % % % betwClusDiffRac = NaN;
    % % % betwClusDiffPow = NaN;
    % % % betwClusDiffPos = NaN;
    
    clustNonNested = clust([clust.nested] == 1); % Non-nested clusters have nestedness 1 (not 0 like they used to have)
    for kc = 1 : length(clustNonNested)
        clEvOnsDt = [clEvOnsDt; clustNonNested(kc).OnsDt]; %#ok<AGROW>
        if kc >= 2
            betwClusPeriod(kc-1) = clustNonNested(kc).OnsDt(1) - clustNonNested(kc-1).OnsDt(end);
            % % % betwClusDiffDur(kc-1) = clustNonNested(kc).szCharTbl.szDurS(1) - clustNonNested(kc-1).szCharTbl.szDurS(end);
            % % % betwClusDiffRac(kc-1) = clustNonNested(kc).szCharTbl.szRac(1) - clustNonNested(kc-1).szCharTbl.szRac(end);
            % % % betwClusDiffPow(kc-1) = clustNonNested(kc).szCharTbl.szPow(1) - clustNonNested(kc-1).szCharTbl.szPow(end);
            % % % betwClusDiffPos(kc-1) = clustNonNested(kc).szCharTbl.postIctPow(1) - clustNonNested(kc-1).szCharTbl.postIctPow(end);
        end
    end
    % Make sure that each seizure is counted maximum of once (which should be ensured by taking only the non-nested clusters but let's double check).
    clEvOnsDt2 = clEvOnsDt;
    clEvOnsDt = unique(clEvOnsDt); % Remove duplicates (some seizures may be part of more clusters)
    if ~all(clEvOnsDt2 == clEvOnsDt)
        error('_jk Even when only non-nested clusters are considered, there were probably still some seizures counted multiple times.')
    end
    numWithinClEv = length(clEvOnsDt); % Number of within-cluster events
    numAllEv = numev - 2; % Subtract the dummy events
    stats.clNumClust = length(clust);
    stats.clNumClustNonNested = length(clustNonNested);
    stats.clFracWithinClustEv = numWithinClEv/numAllEv; % Proportion of seizures that occurred within a cluster. Easier to determine here than later in subjectStats
    stats.clNumEvPerClust = mean(arrayfun(@(x) numel(x.OnsDt), clust));
    stats.clClusterDur = mean(arrayfun(@(x) x.OnsDt(end) - x.OnsDt(1), clust));
    stats.clInterclusPeriod = mean(betwClusPeriod);
    % % % stats.clBetwClusDiffDur = mean(BetwClusDiffDur);
    % % % stats.clBetwClusDiffRac = mean(BetwClusDiffRac);
    % % % stats.clBetwClusDiffPow = mean(BetwClusDiffPow);
    % % % stats.clBetwClusDiffPos = mean(interclusDiffPos);
end