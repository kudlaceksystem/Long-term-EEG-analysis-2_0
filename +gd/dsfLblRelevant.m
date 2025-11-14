function lblSetRelevant = dsfLblRelevant(ll, d)
    % ll ........... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % d ............ relevant line of the dpDesc structure
    lblSetRelevant = ll.lblSet(ismember(ll.lblSet.ClassName, d.MainLbl), :); % lblSet containing only the required classes
    if isempty(d.ExLblAn)
        d.ExLblAn = "";
    end
    lblSetInv = ll.lblSet(ismember(ll.lblSet.ClassName, d.ExLblAn), :); % lblSet containing only the label classes which invalidate the signal for the intended analysis
    for kco = 1 : height(lblSetInv)
        lblSetRelevant(lblSetClCh.Start > lblSetInvalCh.Start(kco) & lblSetClCh.Start < lblSetInvalCh.End(kco), :) = [];
    end
end
