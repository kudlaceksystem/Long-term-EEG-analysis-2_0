function lblSet = dbLblMerge(lblSet, minSeparationS, pointTF)
    % lblSet ........... label set in which labels closer than minSeparationS will be merged
    % minSeparationS ... minimum separation to keep them separate, if closer, merge
    % instant .......... if the label is of type "point", it has no duration, so we will set End = Start
    lblSet = sortrows(lblSet, "Start");
    d = lblSet.Start(2 : end) - lblSet.End(1 : end - 1);
    whichMergeWithPrevious = find(d < seconds(minSeparationS));
    lblSet.End(whichMergeWithPrevious - 1) = lblSet.End(whichMergeWithPrevious);
    lblSet(whichMergeWithPrevious, :) = [];
    if pointTF
        lblSet.End = lblSet.Start;
    end
end