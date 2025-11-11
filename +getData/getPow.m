function pow = getPow(~, lblSet, clnm)
    % The first argument could be sigInfo, but it is not used. It is needed for the consistency with other functions
    % lblSet ... table defined by OSEL system, contains onset, duration, etc.
    % clnm ..... class name - which label class from the lblSet should be used
%% This function is not finished. It just fills in ones.
% We should decide if it should be this function which will load the signal file or if it is going to receive it in
% input. I think the second option will be better but it will require increasing the number of input arguments in all
% functions (possibly just ~ in some of them). We need consistent input arguments across all get functions so that the
% whole system works.


pow = ones(sum(lblSet.ClassName == clnm), 1);

end