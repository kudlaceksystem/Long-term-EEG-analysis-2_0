function validS = dpGetValidAmountCh(ll, ~, invalidity, binLimDt, minSepS)
    % ll ........... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % ~ ... clnm ... class names - which label classes from the lblSet should be counted (often just one of them) - not used in this function
    % invalidity ... class names - which label classes from the lblSet should be used as a marker of invalid (contaminated signal)
    % binLimDt ..... limits of the bin
    % minSepS ...... markers closer than minSepS will be merged
    sigDurS = seconds(sigInfo.SigEnd(1) - sigInfo.SigStart(1));
    % Check if sigDur is the same in all channels
    if any(seconds(sigInfo.SigEnd - sigInfo.SigStart) - sigDurS > 1)
        error('_jk getValidAmountCh sigDur inconsistent across channels')
    end
    numch = height(sigInfo); % Number of channels
    invalidS = zeros(1, numch);
    for kch = 1 : numch
        lblSetInvalCh = ll.lblSet(ismember(ll.lblSet.ClassName, invalidity) & ll.lblSet.Channel == ch, :); % lblSet containing only the label classes which invalidate the signal for the intended analysis
        lblSetInvalCh = getData.dbLblMerge(lblSetInvalCh, minSepS, false);
        lblSetInvalCh = lblSetInvalCh(lblSetInvalCh.End > binLimDt(1) & lblSetInvalCh.Start < binLimDt(2));
        lblSetInvalCh.Start(1) = binLimDt(1);
        lblSetInvalCh.End(end) = binLimDt(2);
        invalidS(1, kch) = sum(seconds(lblSetInvalCh.End - lblSetInvalCh.Start));
    end
    validS = sigDurS*ones(1, numch) - invalidS; % Calculate the valid signal duration
end
