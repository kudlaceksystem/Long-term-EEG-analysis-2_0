function lblSet = dbfLblMergeCh(lblSet, minSeparationS, pointTF)
    % Merges labels separated less than minSeparationS. Works in each channel separately.
    % lblSet ........... label set in which labels closer than minSeparationS will be merged
    % minSeparationS ... in seconds, minimum separation to keep them separate, if closer, merge
    % instant .......... if the label is of type "point", it has no duration, so we will set End = Start
    lblSetMerged = lblSet;
    uch = unique(lblSet.Channel);
    numMrkAfterMrg = zeros(size(uch)); % Number of resulting markers after merging of given channel
    for kch = 1 : numel(uch)
        lblSetCh = lblSet(lblSet.Channel == uch(kch));
        lblSetCh = sortrows(lblSetCh, "Start");
        nummrk = height(lblSetCh); % Number of markers
        mergeWithPreviousTF = false(nummrk, 1);
        for k = 2 : nummrk
            d = lblSetCh.Start(k) - lblSetCh.End(1 : k - 1); % Current marker separation from all the previous
            mergeWithPreviousTF(k) = any(d < seconds(minSeparationS)); % Plugs in true the current marker was insufficiently separated from any of the previous
        end
        mergeWithPreviousSub = find(mergeWithPreviousTF); % Convert logical index to subscript
        nummrg = numel(mergeWithPreviousSub); % Number of markers to merge with previous
        for k = 1 : nummrg
            lblSetCh.End(nummrg - k) = max(lblSetCh.End(nummrg - k), lblSetCh.End(nummrg - k + 1)); % We go backwards
            lblSetCh.End(nummrg - k + 1) = [];
        end
        lblSetCh.ClassName = repelem("Merged", height(lblSetCh.ClassName), 1);
        numMrkAfterMrg(kch) = height(lblSetCh);
        sum(numMrkAfterMrg(1 : (kch - 1)))
        st = sum(numMrkAfterMrg(1 : (kch - 1))) + 1;
        en = sum(numMrkAfterMrg(1 : kch));
        lblSetMerged(st : en, :) = lblSetCh;
    end
    lblSetMerged = lblSetMerged(1 : sum(numMrkAfterMrg));
    lblSetMerged = sortrows(lblSetMerged, "Start");
    if pointTF
        lblSetMerged.End = lblSetMerged.Start;
    end
end
