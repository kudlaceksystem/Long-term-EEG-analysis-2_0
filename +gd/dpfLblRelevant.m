function [lblSetRelevant, lblSetInval] = dpfLblRelevant(ll, d, binLimDt)
    % ll ........... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % d ............ relevant line of the dpDesc structure
    % binLimDt ..... time limits of the current bin in datetime
    lblSetRelevant = ll.lblSet([], :);
    lblSetInval = ll.lblSet([], :);
    uch = unique(ll.lblSet.Channel);
    for kch = 1 : numel(uch)
        lblSetRelevantCh = ll.lblSet(ismember(ll.lblSet.ClassName, d.MainLbl) & ll.lblSet.Channel == uch(kch), :); % lblSet containing only the required classes and channels
        lblSetInvalCh = ll.lblSet((ismember(ll.lblSet.ClassName, d.ExLblCh) & ll.lblSet.Channel == uch(kch)) | ...
            ismember(ll.lblSet.ClassName, d.ExLblAn), :); % lblSet containing only the label classes which invalidate the signal for the intended analysis
        for kco = 1 : height(lblSetInvalCh)
            lblSetRelevantCh(lblSetRelevantCh.End > lblSetInvalCh.Start(kco) & lblSetRelevantCh.Start < lblSetInvalCh.End(kco), :) = [];
        end
        lblSetRelevantCh = lblSetRelevantCh(lblSetRelevantCh.Start > binLimDt(1) & lblSetRelevantCh.Start < binLimDt(2), :);
        lblSetInvalCh = lblSetInvalCh(lblSetInvalCh.Start > binLimDt(1) & lblSetInvalCh.Start < binLimDt(2), :);
        lblSetRelevantCh.Channel = ones(height(lblSetRelevantCh), 1);
        lblSetInvalCh.Channel = ones(height(lblSetInvalCh), 1);
        lblSetRelevant = [lblSetRelevant; lblSetRelevantCh]; %#ok<AGROW>
        lblSetInval = [lblSetInval; lblSetInvalCh]; %#ok<AGROW>
    end
end
