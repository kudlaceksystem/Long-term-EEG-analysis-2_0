function c = dpGetCountCh(ll, clnm, invalidity, binLimDt)
    % ll ........... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % clnm ......... class names - which label classes from the lblSet should be counted (often just one of them)
    % invalidity ... class names - which label classes from the lblSet should be used as a marker of invalid (contaminated signal)
    % binLimDt ...... limits of the bin
    % This function also implements removal of invalid portions of EEG, typically contaminated by artifacts which make
    % it impossible to analyze the count of given phenomenon. It is hardcoded here since it is this functions
    % responsibility to do it. It should use a function which will make sure that corresponding values will be written
    % also in the dp.Valid table.
    lblSetRelevant = dpGetRelevantCh(ll, clnm, invalidity, binLimDt)
    lblSetCl = ll.lblSet(ismember(ll.lblSet.ClassName, clnm), :); % lblSet containing only the required classes
    lblSetInv = ll.lblSet(ismember(ll.lblSet.ClassName, invalidity), :); % lblSet containing only the label classes which invalidate the signal for the intended analysis
    numch = height(ll.sigInfo); % Number of channels
    c = NaN(1, numch); % Count
    for kch = 1 : numch
        lblSetClCh = lblSetCl(lblSetCl.Channel == kch); % lblSet containing only the required class and channel
        % Remove IEDs occurring when the signal is invalid (contaminated)
        lblSetInvalCh = lblSetInv(lblSetInv.Channel == kch, :);
        for kco = 1 : height(lblSetInvalCh)
            lblSetClCh(lblSetClCh.Start > lblSetInvalCh.Start(kco) & lblSetClCh.Start < lblSetInvalCh.End(kco), :) = [];
        end
        c(1, kch) = sum(lblSetClCh.Start > binLimDt(1) & lblSetClCh.Start < binLimDt(2));
    end
end
    lblSetRelevant = lblSetRelevant(lblSetRelevant.Start > binLimDt(1) & lblSetRelevant.End < binLimDt(2));
