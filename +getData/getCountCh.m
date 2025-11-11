function c = getCountCh(ll, clnm, invalidity, binLimN)
    % ll ........... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % clnm ......... class name - which label class from the lblSet should be used
    % invalidity ... class name - which label class from the lblSet should be used
    % This function also implements removal of invalid portions of EEG, typically contaminated by artifacts which make
    % it impossible to analyze the count of given phenomenon. It is hardcoded here since it is this functions
    % responsibility to do it. It should use a function which will make sure that corresponding values will be written
    % also in the dp.Valid table.
    lblSetCl = ll.lblSet(ll.lblSet.ClassName == clnm); % lblSet containing only the required class
    lblSetContam = getEegInvalidLblSet(ll, invalidity); % Label set containing the labels of contaminating phenomena (e.g. artifacts or seizures)
    numch = height(ll.sigInfo); % Number of channels
    c = NaN(1, numch); % Count
    for kch = 1 : numch
        lblSetClCh = lblSetCl(lblSetCl.Channel == kch); % lblSet containing only the required class and channel
        % Remove IEDs occurring when the signal is contaminated
        lblSetContamCh = lblSetContam(lblSetContam.Channel == kch, :);
        for kco = 1 : height(lblSetContamCh)
            lblSetClCh(lblSetClCh.Start > lblSetContamCh.Start(kco) & lblSetClCh.Start < lblSetContamCh.End(kco), :) = [];
        end
        c(1, kch) = sum(lblSetClCh.Start > binLimN(1) & lblSetClCh.Start < binLimN(2));
    end
end
