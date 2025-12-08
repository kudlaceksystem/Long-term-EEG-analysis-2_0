function [lblSetRelevantCh, lblSetInvalCh] = dpfLblRelevant(ll, d, binLimDt, ch)
    % ll ........... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % d ............ relevant line of the dpDesc structure
    % binLimDt ..... time limits of the current bin in datetime
    % ch ........... channels that should be kept only
    lblSetRelevantCh = ll.lblSet(ismember(ll.lblSet.ClassName, d.MainLbl) & ll.lblSet.Channel == ch, :); % lblSet containing only the required classes and channels
    lblSetInvalCh = ll.lblSet((ismember(ll.lblSet.ClassName, d.ExLblCh) & ll.lblSet.Channel == ch) | ...
        ismember(ll.lblSet.ClassName, d.ExLblAn), :); % lblSet containing only the label classes which invalidate the signal for the intended analysis
    for kco = 1 : height(lblSetInvalCh)
        lblSetRelevantCh(lblSetRelevantCh.End > lblSetInvalCh.Start(kco) & lblSetRelevantCh.Start < lblSetInvalCh.End(kco), :) = [];
    end
    lblSetRelevantCh = lblSetRelevantCh(lblSetRelevantCh.Start > binLimDt(1) & lblSetRelevantCh.Start < binLimDt(2), :);
    lblSetInvalCh = lblSetInvalCh(lblSetInvalCh.Start > binLimDt(1) & lblSetInvalCh.Start < binLimDt(2), :);
end
