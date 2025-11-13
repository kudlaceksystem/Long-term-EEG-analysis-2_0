function lblSet = dbLblMerge(lblSet, minSeparationS, pointTF)
    % lblSet ........... label set in which labels closer than minSeparationS will be merged
    % minSeparationS ... minimum separation to keep them separate, if closer, merge
    % instant .......... if the label is of type "point", it has no duration, so we will set End = Start
    lblSet = sortrows(lblSet, "Start");
    d = lblSet.Start(2 : end) - lblSet.End(1 : end - 1);
    whichMergeWithNext = find(d < seconds(minSeparationS));
    lblSet.End(whichMergeWithNext) = lblSet.End(whichMergeWithNext + 1);
    lblSet(whichMergeWithNext + 1, :) = [];
    if pointTF
        lblSet.End = lblSet.Start;
    end
end