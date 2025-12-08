function r = dpbGetRatePhCh(binTbl, nm)
binTbl
    % nm ........... name of the dp field to work on
    if any(binTbl.(nm).Count < 0, "all")
        disp(binTbl.(nm).Count)
        warning('_jk dpbGetRatePhCh Count negative.')
        pause
    end
    if any(binTbl.(nm).ValidS < 0, "all")
        disp(binTbl.(nm).ValidS)
        warning('_jk dpbGetRatePhCh ValidS negative.')
        pause
    end
    r = sum(binTbl.(nm).Count, 1)./sum(binTbl.(nm).ValidS, 1);
    r = r*3600;
end
