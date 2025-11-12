function lblSetRelevant = dsLblRelevant(ll, clnm, invalidity)
    % ll ........... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % clnm ......... class names - which label classes from the lblSet should be counted (often just one of them)
    % invalidity ... class names - which label classes from the lblSet should be used as a marker of invalid (contaminated signal)
    lblSetRelevant = ll.lblSet(ismember(ll.lblSet.ClassName, clnm), :); % lblSet containing only the required classes
    lblSetInv = ll.lblSet(ismember(ll.lblSet.ClassName, invalidity), :); % lblSet containing only the label classes which invalidate the signal for the intended analysis
    for kco = 1 : height(lblSetInv)
        lblSetRelevant(lblSetClCh.Start > lblSetInvalCh.Start(kco) & lblSetClCh.Start < lblSetInvalCh.End(kco), :) = [];
    end
end
