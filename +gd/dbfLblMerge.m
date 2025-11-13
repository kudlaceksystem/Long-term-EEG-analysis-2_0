function lblSet = dbfLblMerge(lblSet, minSeparationS, pointTF)
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
        lblSet.End(nummrg - k) = max(lblSet.End(nummrg - k), lblSet.End(nummrg - k + 1)); % We go backwards
        lblSet.End(nummrg - k + 1) = [];
    end
    lblSet.ClassName = "Merged";
    if pointTF
        lblSet.End = lblSet.Start;
    end
end
