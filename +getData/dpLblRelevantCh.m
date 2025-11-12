function lblSetRelevantCh = dpGetRelevantCh(ll, clnm, invalidity, binLimDt, ch)
    lblSetRelevant = ll.lblSet(ismember(ll.lblSet.ClassName, clnm) & ll.lblSet.Channel == ch, :); % lblSet containing only the required classes and channels
    lblSetInv = ll.lblSet(ismember(ll.lblSet.ClassName, invalidity), :); % lblSet containing only the label classes which invalidate the signal for the intended analysis
    for kco = 1 : height(lblSetInv)
        lblSetRelevant(lblSetRelevant.End > lblSetInvalCh.Start(kco) & lblSetClCh.Start < lblSetInvalCh.End(kco), :) = [];
    end
    lblSetRelevant = lblSetRelevant(lblSetRelevant.Start > binLimDt(1) & lblSetRelevant.Start < binLimDt(2));
end
