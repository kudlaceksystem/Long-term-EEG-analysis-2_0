function lblSetInv = getEegInvalidLblSet(~, lblSet, clnm)
    % The first argument could be sigInfo, but it is not used. It is needed for the consistency with other functions
    % lblSet ... table defined by OSEL system, contains onset, duration, etc.
    % clnm ..... class name - which label class from the lblSet should be used
    if clnm == ""
        clnmAll = unique(lblSet.ClassName);
        whichInvalidate = contains(clnmAll, ...
            ["Seizure", "seizure", "SEIZURE", "S", "art", "Art", "EMG", "emg", "Emg"]); % Edit this according to what invalidates the signal for your analysis
        clnm = clnmAll(whichInvalidate);
    end
    % % % % % % % % % % % % % % % % % % numch = height(sigInfo); % Number of channels
    lblSetInv = lblSet(ismember(lblSet.ClassName, clnm));
    % % % % % % % % % % invalidS = zeros(1, numch);
    % % % % % % % % % % for kch = 1 : numch
    % % % % % % % % % %     ls = lblSetInv(lblSetInv.Channel == kch, :); % Keep only the labels that we need
    % % % % % % % % % %     for kl = 1 : height(ls)
    % % % % % % % % % %         invalidS(1, kch) = invalidS(1, kch) + seconds(ls.End(kl) - ls.Start(kl));
    % % % % % % % % % %     end
    % % % % % % % % % % end
    % % % % % % % % % % sigDur = seconds(sigInfo.SigEnd(1) - sigInfo.SigStart(1));
    % % % % % % % % % % % Check if sigDur is the same in all channels
    % % % % % % % % % % if any(seconds(sigInfo.SigEnd - sigInfo.SigStart) - sigDur > 1)
    % % % % % % % % % %     error('_jk getValidAmountCh sigDur inconsistent across channels')
    % % % % % % % % % % end
    % % % % % % % % % % validS = sigDur*ones(1, numch) - invalidS; % Calculate the valid signal duration
end
