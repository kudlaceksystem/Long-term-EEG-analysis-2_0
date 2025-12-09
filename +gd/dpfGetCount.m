function c = dpfGetCount(ll, d, binLimDt)
    % ll ........... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % d ............ relevant line of the dpDesc structure
    % binLimDt ..... limits of the bin
    % This function uses gd.dpfLblRelevantCh which also implements removal of invalid portions of EEG, typically contaminated by artifacts which make
    % it impossible to analyze the count of given phenomenon.
    % % % numch = height(ll.sigInfo); % Number of channels
    pointTF = ~(any(ll.lblDef.LabelType(ismember(ll.lblDef.ClassName, d.MainLbl)) == "roi"));
    c = NaN(1, numch);
    % % % % for kch = 1 : numch
    lblSetRelevant = gd.dpfLblRelevant(ll, d, binLimDt);
        lblSetRelevantCh = gd.dbfLblMerge(lblSetRelevantCh, d.MinSepS, pointTF);
        c(1, kch) = height(lblSetRelevantCh);
    end
end
