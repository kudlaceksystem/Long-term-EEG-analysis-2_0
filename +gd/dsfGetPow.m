function pow = dsfGetPow(ll, d)
    % ll ........... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % d ............ relevant line of the dsDesc structure
    lblSetRelevant = gd.dsfLblRelevant(ll, d);
    pointTF = ~(any(ll.lblDef.LabelType(ismember(ll.lblDef.ClassName, d.MainLbl)) == "roi"));
    lblSetRelevant = gd.dbfLblMerge(lblSetRelevant, d.MinSepS, pointTF);
    pow = ones(size(lblSetRelevant.Start));
end
