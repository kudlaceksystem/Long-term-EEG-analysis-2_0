function c = getCountCh(sigInfo, lblSet, clnm)
    % lblSet ... table defined by OSEL system, contains onset, duration, etc.
    % clnm ..... class name - which label class from the lblSet should be used
    numch = height(sigInfo); % Number of channels
    c = NaN(1, numch);
    for kch = 1 : numch
        c(1, kch) = sum(lblSet.ClassName == clnm & lblSet.Channel == kch);
    end
end
