function onsDt = dsGetOnsDt(ll, clnm, invalidity, minSepS)
    % ll ........... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % clnm ......... class names - which label classes from the lblSet should be counted (often just one of them)
    % invalidity ... class names - which label classes from the lblSet should be used as a marker of invalid (contaminated signal)
    lblSetRelevant = getData.dsLblRelevant(ll, clnm, invalidity);
    pointTF = ~(any(ll.lblDef.LabelType(ismember(ll.lblDef.ClassName, clnm)) == "roi"));
    lblSetRelevant = getData.dbLblMerge(lblSetRelevant, minSepS, pointTF);
    onsDt = lblSetRelevant.Start;
end