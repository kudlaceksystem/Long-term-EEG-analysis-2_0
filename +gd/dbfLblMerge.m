function lblSet = dbfLblMerge(lblSet, minSeparationS, pointTF)
    % Merges labels separated less than minSeparationS. Does not care about in which channel the label is.
    % lblSet ........... label set in which labels closer than minSeparationS will be merged
    % minSeparationS ... minimum separation to keep them separate, if closer, merge
    % instant .......... if the label is of type "point", it has no duration, so we will set End = Start
    lblSet = sortrows(lblSet, "Start");
    nummrk = height(lblSet); % Number of markers
    mergeWithPreviousTF = false(nummrk, 1);
    for k = 2 : nummrk
        d = lblSet.Start(k) - lblSet.End(1 : k - 1); % Current marker separation from all the previous
        mergeWithPreviousTF(k) = any(d < seconds(minSeparationS)); % Plugs in true the current marker was insufficiently separated from any of the previous
    end
    mergeWithPreviousSub = find(mergeWithPreviousTF);
    nummrg = numel(mergeWithPreviousSub); % Number of markers to merge with previous
    for k = 1 : nummrg
        toBeMerged = mergeWithPreviousSub(nummrg - k + 1);
        mergedWith = mergeWithPreviousSub(nummrg - k + 1) - 1;
        lblSet.End(mergedWith) = max(lblSet.End(mergedWith), lblSet.End(toBeMerged)); % We go backwards
        lblSet(toBeMerged, :) = [];
    end
    lblSet.ClassName = repelem("Merged", height(lblSet.ClassName), 1);
    if pointTF
        lblSet.End = lblSet.Start;
    end
    a = lblSet.Start;
    b = unique(lblSet.Start);
    if length(a) ~= length(b)
        '====dbfLblMerge===='
        a_ = a
        b_ = b
        pause
        '==================='
    end
end
