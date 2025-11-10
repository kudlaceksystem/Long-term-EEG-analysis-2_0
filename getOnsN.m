function onsN = getOnsN(~, lblSet, clnm)
    % lblSet ... table defined by OSEL system, contains onset, duration, etc.
    % clnm ..... class name - which label class from the lblSet should be used
    onsN = datenum(lblSet.Start(lblSet.ClassName == clnm, :)); %#ok<DATNM>
end