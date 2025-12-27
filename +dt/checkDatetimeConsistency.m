function isConsistentTF = checkDatetimeConsistency(lblpn, klbl, ll)
    % Check if sigInfo.FileName corresponds to file name
    fndattimStr = regexp(lblpn{klbl}, '\d\d\d\d\d\d_\d\d\d\d\d\d', 'match');
    filenameDt = datetime(fndattimStr{1}, 'InputFormat', 'yyMMdd_HHmmss');
    fcdattimStr = regexp(ll.sigInfo.FileName, '\d\d\d\d\d\d_\d\d\d\d\d\d', 'match');
    filecontDt = datetime(fcdattimStr{1}, 'InputFormat', 'yyMMdd_HHmmss');
    filenameConsWithFilen = filenameDt == filecontDt;
    if ~filenameConsWithFilen
        disp(klbl)
        disp(lblpn{klbl})
        disp(ll.sigInfo)
        warning('_jk File name and file contents do not correspond in terms of the date and time.')
        pause
    end
    
    % Check if SigStart corresponds to file name
    fndattimStr = regexp(lblpn{klbl}, '\d\d\d\d\d\d_\d\d\d\d\d\d', 'match');
    filenameDt = datetime(fndattimStr{1}, 'InputFormat', 'yyMMdd_HHmmss');
    filecontDt = ll.sigInfo.SigStart;
    sigStartConsWithFilen = all(filenameDt == filecontDt);
    if ~sigStartConsWithFilen
        disp(klbl)
        disp(lblpn{klbl})
        disp(ll.sigInfo)
        warning('_jk SigStart does not correspond to file name.')
        pause
    end
    isConsistentTF = filenameConsWithFilen && sigStartConsWithFilen;
end