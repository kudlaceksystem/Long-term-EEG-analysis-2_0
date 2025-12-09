function r = dpbGetRatePh(binTbl, nm)
    % Get event rate in the unit of events per hour.
    % binTbl ....... dp structure containing only data belonging to given bin
    % nm ........... name of the dp field to work on
    %% Settings
    % The rate is computed as count divided by time in seconds and then multiplied by 3600 to convert to events per hour
    countColName = "Count"; % Name of the column containing data on count
    timeSColName = "ValidS"; % Name of the column containing data on time (period in which events were counted)
    if any(binTbl.(nm).(countColName) < 0, "all")
        disp(binTbl.(nm).(countColName))
        warning('_jk dpbGetRatePhCh Count negative.')
        pause
    end
    if any(binTbl.(nm).(timeSColName) < 0, "all")
        disp(binTbl.(nm).(timeSColName))
        warning('_jk dpbGetRatePhCh ValidS negative.')
        pause
    end
    r = sum(binTbl.(nm).(countColName), 1)./sum(binTbl.(nm).(timeSColName), 1);
    r = r*3600;
end
