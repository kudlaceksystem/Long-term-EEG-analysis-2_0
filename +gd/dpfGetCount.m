function c = dpfGetCount(ll, d, binLimDt)
    % ll ........... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % d ............ relevant line of the dpDesc structure
    % binLimDt ..... limits of the bin
    % This function uses gd.dpfLblRelevant which also implements removal of events detected during invalid portions of EEG,
    % typically contaminated by artifacts which make it impossible to analyze the count of given phenomenon.
    pointTF = ~(any(ll.lblDef.LabelType(ismember(ll.lblDef.ClassName, d.MainLbl)) == "roi"));
    lblSetRelevant = gd.dpfLblRelevant(ll, d, binLimDt);
    lblSetRelevant = gd.dbfLblMerge(lblSetRelevant, d.MinSepS, pointTF);
    c = height(lblSetRelevant);
end
