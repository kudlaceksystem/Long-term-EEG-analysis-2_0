function dobTable = getSubjectList(stg, dobpn)
    arguments (Input)
        stg
        dobpn
    end
    arguments (Output)
        dobTable
    end
    % datetime below often throws warnings. I want to turn them off for this function.
    origWarningState = warning('query', 'all'); % Save the current warning state
    warning('off', 'all'); % Turn off all warnings

    % Now the function proper
    videoEEGdata = readtable(dobpn);
    numSubj = height(videoEEGdata);
    Subject = strings(numSubj, 1);
    Birth = datetime.empty(numSubj, 0);
    Birth.TimeZone = "UTC";
    Sex = false(numSubj, 1);
    for k = 1 : numSubj
        Subject(k, 1) = string(videoEEGdata.Mouse{k});
        r = regexp(videoEEGdata.Birth(k), '\d\d\d\d-\d\d-\d\d', 'match');
        if ~isempty(r{1})
            dt = datetime(r{1}, 'InputFormat', 'uuuu-MM-dd', 'Format', 'uuuu-MM-dd', 'TimeZone', stg.timeZoneStr);
            dt.TimeZone = "UTC";
            if year(dt) < 1000
                dt.Year = dt.Year + 2000;
            end
        end
        
        r = regexp(videoEEGdata.Birth(k), '\d+-...-\d+', 'match');
        if ~isempty(r{1})
            dt = datetime(r{1}, 'InputFormat', 'dd-MMM-yyyy', 'Format', 'uuuu-MM-dd', 'TimeZone', stg.timeZoneStr);
            dt.TimeZone = "UTC";
            if year(dt) < 1000
                dt.Year = dt.Year + 2000;
            end
        end
        
        r = regexp(videoEEGdata.Birth(k), '\d+\.\d+\.\d+', 'match');
        if ~isempty(r{1})
            dt = datetime(r{1}, 'InputFormat', 'dd.MM.uuuu', 'Format', 'uuuu-MM-dd', 'TimeZone', stg.timeZoneStr);
            dt.TimeZone = "UTC";
            if year(dt) < 1000
                dt.Year = dt.Year + 2000;
            end
        end
        Birth(k, 1) = dt;
        Sex(k, 1) = true;
    end
    % Restore the original warning state
    warning(origWarningState); 
    dobTable = table(Subject, Birth, Sex);
end