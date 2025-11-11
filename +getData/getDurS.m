function durS = getDurS(~, lblSet, clnm)
    % The first argument could be sigInfo, but it is not used. It is needed for the consistency with other functions
    % lblSet ... table defined by OSEL system, contains onset, duration, etc.
    % clnm ..... class name - which label class from the lblSet should be used
    st = lblSet.Start(lblSet.ClassName == clnm, :);
    en = lblSet.End(lblSet.ClassName == clnm, :);
    durS = seconds(en - st);
end