function [x, y, xx] = getValidXY(subjInfo, d, dp)
    % x .... [blockStart, blockEnd, NaN, blockStart, blockEnd, NaN, ...], time from the date of birth (dob)
    % y .... [0, 0, NaN, 0, 0, NaN, NaN, NaN, NaN, 0, 0, NaN, ...] if there are three NaNs, it indicates no data for given block (dropout)
    % xx ... first row beginnings, second row ends of the valid blocks, time from the date of birth (dob)
    nm = d.EventValidSrc;
    t = dp.(nm).tax;
    ta = [t(1) - (t(2) - t(1)); t]; % Add the beginning of the first block
    x1 = days(ta - subjInfo.dob);
    x = NaN(3*(numel(x1) - 1), 1);
    x(1 : 3 : end - 2) = x1(1 : end - 1);
    x(2 : 3 : end - 1) = x1(2 : end);
    
    yTF = false(size(x));
    yTF(1 : 3 : end - 2) = all(isnan(dp.(nm).ValidS), 2); % If it is NaN in all channels, make y true
    yTF(2 : 3 : end - 1) = all(isnan(dp.(nm).ValidS), 2);
    yTF(1 : 3 : end - 2) = all(dp.(nm).ValidS == 0, 2); % If there is 0 seconds of valid signal, make y true
    yTF(2 : 3 : end - 1) = all(dp.(nm).ValidS == 0, 2);
    y = zeros(size(yTF));
    y(yTF) = NaN; % Where y is true, make it NaN
    y(3 : 3 : end) = NaN; % Separation of bins
    
    xx(1, :) = x(1 : 3 : end - 2);
    xx(2, :) = x(2 : 3 : end - 1);
    yy(1, :) = y(1 : 3 : end - 2);
    xx = xx(:, ~isnan(yy(1, :))); % First row beginnings, second row ends of the valid blocks
end
