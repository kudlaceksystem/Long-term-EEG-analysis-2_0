function durS = getDurS(~, lblSet, clnm)
    % lblSet ... table defined by OSEL system, contains onset, duration, etc.
    % clnm ..... class name - which label class from the lblSet should be used
    st = lblSet.Start(lblSet.ClassName == clnm, :);
    en = lblSet.End(lblSet.ClassName == clnm, :);
    durS = seconds(en - st);
end