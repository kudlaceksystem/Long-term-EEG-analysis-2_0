function c = dpGetCountCh(ll, clnm, invalidity, binLimDt, minSepS)
    % ll ........... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % clnm ......... class names - which label classes from the lblSet should be counted (often just one of them)
    % invalidity ... class names - which label classes from the lblSet should be used as a marker of invalid (contaminated signal)
    % binLimDt ...... limits of the bin
    % This function also implements removal of invalid portions of EEG, typically contaminated by artifacts which make
    % it impossible to analyze the count of given phenomenon. It is hardcoded here since it is this functions
    % responsibility to do it. It should use a function which will make sure that corresponding values will be written
    % also in the dp.Valid table.
    numch = height(ll.sigInfo); % Number of channels
    pointTF = ~(any(ll.lblDef.LabelType(ismember(ll.lblDef.ClassName, clnm)) == "roi"));
    c = NaN(1, numch);
    for kch = 1 : numch
        lblSetRelevantCh = getData.dpLblRelevantCh(ll, clnm, invalidity, binLimDt);
        lblSetRelevantCh = getData.dbLblMerge(lblSetRelevantCh, minSepS, pointTF);
        c(1, kch) = height(lblSetRelevantCh);
    end
end
