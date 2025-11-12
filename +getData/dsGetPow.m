function pow = dsGetPow(ll, clnm, invalidity)
    % ll ........... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % clnm ......... class names - which label classes from the lblSet should be counted (often just one of them)
    % invalidity ... class names - which label classes from the lblSet should be used as a marker of invalid (contaminated signal)
    lblSetRelevant = getData.dsGetRelevant(ll, clnm, invalidity);
    pow = ones(size(lblSetRelevant.Start));
end
