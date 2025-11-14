function validS = dpfGetValidAmountCh(ll, d, binLimDt)
    % ll ........... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % d ............ relevant line of the dpDesc structure
    % minSepS ...... markers closer than minSepS will be merged
    minSepS = d.MinSepS;
    sigDurS = seconds(ll.sigInfo.SigEnd(1) - ll.sigInfo.SigStart(1));
    % Check if sigDur is the same in all channels
    if any(seconds(ll.sigInfo.SigEnd - ll.sigInfo.SigStart) - sigDurS > 1)
        error('_jk getValidAmountCh sigDur inconsistent across channels')
    end
    numch = height(ll.sigInfo); % Number of channels
    invalidS = zeros(1, numch);
    for kch = 1 : numch
        lblSetInvalCh = ll.lblSet((ismember(ll.lblSet.ClassName, d.ExLblCh) & ll.lblSet.Channel == kch) | ...
            ismember(ll.lblSet.ClassName, d.ExLblAn), :); % lblSet containing only the label classes which invalidate the signal for the intended analysis
        lblSetInvalCh = gd.dbfLblMerge(lblSetInvalCh, minSepS, false);
        lblSetInvalCh = lblSetInvalCh(lblSetInvalCh.End > binLimDt(1) & lblSetInvalCh.Start < binLimDt(2), :);
        if height(lblSetInvalCh) > 0
            lblSetInvalCh.Start(1) = max(lblSetInvalCh.Start(1), binLimDt(1));
            lblSetInvalCh.End(end) = min(lblSetInvalCh.End(end), binLimDt(2));
            invalidS(1, kch) = sum(seconds(lblSetInvalCh.End - lblSetInvalCh.Start));
        else
            invalidS(1, kch) = 0;
        end
    end
    validS = sigDurS*ones(1, numch) - invalidS; % Calculate the valid signal duration
    if any(validS < 0)
        disp(validS)
        warning('_jk dpfGetValidAmountCh validS negative')
        pause
    end
end
