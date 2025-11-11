function c = getCountCh(sigInfo, lblSet, clnm, invalidity)
    % sigInfo .. table defined by OSEL system, contains the start and end of the signal file that was labelled
    % lblSet ... table defined by OSEL system, contains onset, duration, etc.
    % clnm ..... class name - which label class from the lblSet should be used
    % This function also implements removal of invalid portions of EEG, typically contaminated by artifacts which make
    % it impossible to analyze the count of given phenomenon. It is hardcoded here since it is this functions
    % responsibility to do it. It should use a function which will make sure that corresponding values will be written
    % also in the dp.Valid table.

    numch = height(sigInfo); % Number of channels
    c = NaN(1, numch);
    for kch = 1 : numch
        c(1, kch) = sum(lblSet.ClassName == clnm & lblSet.Channel == kch);
    end
end
