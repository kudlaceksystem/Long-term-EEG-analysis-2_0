function lblSetInv = getEegInvalidLblSet(ll, clnm, ~)
    % ll ....... contents of OSEL label file, i.e. sigInfo, lblDef, lblSet
    % lblSet ... table defined by OSEL system, contains onset, duration, etc.
    % clnm ..... class name - which label class from the lblSet should be used
    if clnm == ""
        clnmAll = unique(lblSet.ClassName);
        whichInvalidate = contains(clnmAll, ...
            ["Seizure", "seizure", "SEIZURE", "S", "art", "Art", "EMG", "emg", "Emg"]); % Edit this according to what invalidates the signal for your analysis
        clnm = clnmAll(whichInvalidate);
    end
    lblSetInv = lblSet(ismember(lblSet.ClassName, clnm));
end
