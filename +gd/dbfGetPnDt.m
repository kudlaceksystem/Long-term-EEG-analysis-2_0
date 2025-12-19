function [pn, Dt, N] = dbfGetPnDt(stg, p)
    % Get file path, name, start date in datetime and datenum
    d = dir([p, '\*.mat']);
    n = {d.name}';
    pn = fullfile(p, n);
    N = cellfun(@(x) datenum(regexp(x, '\d\d\d\d\d\d_\d\d\d\d\d\d', 'match'), 'yymmdd_HHMMSS'), n, 'UniformOutput', true); %#ok<DATNM>
    Dt = cellfun(@(x) datetime(regexp(x, '\d\d\d\d\d\d_\d\d\d\d\d\d', 'match'), 'InputFormat', 'yyMMdd_HHmmss'), n, 'UniformOutput', true);
    Dt.TimeZone = stg.timeZoneStr;
    Dt.TimeZone = "UTC";
end
