function [ll, firstAftChngSub, firstAftChngStr, prevSigInfo, stdFileIntOut, afterTimeChngTFOut] = convertLblToUTC(stg, lblpn, klbl, ll, prevSigInfo, stdFileIntIn, afterTimeChngTFIn)
    % Assumes it is called with ll in chronologic order.
    % To keep times correctly even within the hour of transition from daylight saving time to normal, lblSet is first
    % converted to relative times to SigStart, then the SigStart is converted to TimeZone where data was rcorded, then
    % to UTC and one hour is subtracted if needed (if the recording was before the time change). Then, the lblSet is
    % compute again from the SigStart and relative times. Similarly also the SigEnd is computed.
    % Function outputs also variables firstAftChngSub, firstAftChngStr, prevSigInfo, stdFileIntOut, afterTimeChngTFOut
    % which need to be remembered between the calls to the function or will be used later in getData.
    afterTimeChange = afterTimeChngTFIn;
    firstAftChngSub = [];
    firstAftChngStr = [];
    stdFileIntOut = stdFileIntIn;
    % Get times relative to signal file start
    sigDur = ll.sigInfo.SigEnd(1) - ll.sigInfo.SigStart(1);
    markerTime = ll.lblSet.Start - ll.sigInfo.SigStart(1);
    markerDur = ll.lblSet.End - ll.lblSet.Start;
    
     % Set the time zone where the data was recorded
    ll.sigInfo = dt.tblSetTimeZone(ll.sigInfo, stg.recTimeZoneStr); % First set the time zone where the data was recorded
    
    % Check if it is during the transition from the daylight saving time to the standard time
    if dt.dtIsAmbiguous(ll.sigInfo.SigStart(1))
        fileInterval = ll.sigInfo.SigStart(1) - prevSigInfo.SigStart(1);
        if fileInterval < 0.5*stdFileIntIn
            afterTimeChange = true;
            firstAftChngSub = klbl;
            firstAftChngStr = lblpn{klbl};
            % % % % % disp('_jk First file after change from DST to standard found.')
            % % % % % pause
        end
        if afterTimeChange
            addToUTC = hours(0); % By default, the conversion to UTC assumes the time after the change (i.e. during the standard time already)
        else
            addToUTC = hours(-1); % If we are not yet after the change, we have to subtract one hour, see above
        end
    else
        addToUTC = hours(0);
        if ~isempty(prevSigInfo)
            stdFileIntOut = ll.sigInfo.SigStart(1) - prevSigInfo.SigStart(1);
        end
    end

    prevSigInfo = ll.sigInfo;
    
    % Change TimeZone to UTC
    ll.sigInfo = dt.tblSetTimeZone(ll.sigInfo, "UTC"); % Then, convert it to the UTC
    ll.sigInfo.SigStart = ll.sigInfo.SigStart + addToUTC;
    ll.sigInfo.SigEnd = ll.sigInfo.SigStart + sigDur;

    % Set lblSet
    ll.lblSet = dt.tblSetTimeZone(ll.lblSet, "UTC");
    ll.lblSet.Start = ll.sigInfo.SigStart(1) + markerTime;
    ll.lblSet.End = ll.lblSet.Start + markerDur;

    afterTimeChngTFOut = afterTimeChange;
end