function onsN = getOnsN(~, lblSet, clnm)
    % The first argument could be sigInfo, but it is not used. It is needed for the consistency with other functions
    % lblSet ... table defined by OSEL system, contains onset, duration, etc.
    % clnm ..... class name - which label class from the lblSet should be used
    onsN = datenum(lblSet.Start(lblSet.ClassName == clnm, :)); %#ok<DATNM>
end