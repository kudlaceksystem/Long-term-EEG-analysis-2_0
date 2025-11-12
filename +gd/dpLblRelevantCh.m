function [lblSetRelevantCh, lblSetInvalCh] = dpLblRelevantCh(ll, clnm, invalidity, binLimDt, ch)
    % ll ........... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % clnm ......... class names - which label classes from the lblSet should be counted (often just one of them)
    % invalidity ... class names - which label classes from the lblSet should be used as a marker of invalid (contaminated signal)
    % binLimDt ..... time limits of the current bin in datetime
    % ch ........... channels that should be kept only
    lblSetRelevantCh = ll.lblSet(ismember(ll.lblSet.ClassName, clnm) & ll.lblSet.Channel == ch, :); % lblSet containing only the required classes and channels
    lblSetInvalCh = ll.lblSet(ismember(ll.lblSet.ClassName, invalidity) & ll.lblSet.Channel == ch, :); % lblSet containing only the label classes which invalidate the signal for the intended analysis
    for kco = 1 : height(lblSetInvalCh)
        lblSetRelevantCh(lblSetRelevantCh.End > lblSetInvalCh.Start(kco) & lblSetRelevantCh.Start < lblSetInvalCh.End(kco), :) = [];
    end
    lblSetRelevantCh = lblSetRelevantCh(lblSetRelevantCh.Start > binLimDt(1) & lblSetRelevantCh.Start < binLimDt(2), :);
    lblSetInvalCh = lblSetInvalCh(lblSetInvalCh.Start > binLimDt(1) & lblSetInvalCh.Start < binLimDt(2), :);
end
