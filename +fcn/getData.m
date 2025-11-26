function [subjInfo, ds, dp] = getData(dsDesc, dpDesc, lblp, snlp, dobTable, ksubj, subjNmOrig)
    % dsDesc ........ data to stem description
    % dpDesc ........ data to plot description
    % lblp .......... path to label files
    % snlp .......... path to signal files
    % dobTable ...... tabel containing dates of birth (and possibly sex and other data)
    % ksubj ......... which subject are we processing now
    % subjNmOrig .... original subject name (I do not know why it is needed here

    % Get sz data and signal characteristics such as IED rate, amount of EMG artifacts and critical slowing markers.
    global stg
    
    %% NEEDS REWRITING If the data for this subject already exist load it
    if exist([stg.dataFolder, 'Data-', num2str(seconds(dpDesc.BinLenDu)), '-', char(subjNmOrig), '.mat'], 'file')
        load([stg.dataFolder, 'Data-', num2str(seconds(dpDesc.BinLenDu)), '-', char(subjNmOrig), '.mat'],...
            'subjInfo', 'szCharTbl', 'szCharLabel', 'siCharTbl', 'siCharLabel')
        stg.szCharYLbl = szCharLabel;
        stg.siCharYLbl = siCharLabel;
        stg.saCharYLbl = [szCharLabel, siCharLabel];
        stg.ssCharYLbl = [siCharLabel, siCharLabel];
        if isa(subjInfo.dob, 'table')
            subjInfo.dob = datenum(subjInfo.dob{1, 1});
        end
        subjInfo.subjNmOrig = subjInfo.subjNm;
        if ~stg.keepOriginalSubjectName
            subjInfo.subjNmOrig = subjInfo.subjNm;
            subjInfo.subjNm = ['Mouse', num2str(stg.subjNumber(ksubj), '%02d')];
        end
        subjInfo.sex = dobTable{ksubj, 3};
        return
    end

    % Find out if we will need signal data (as of now, the function is prepared for labels only)
    for knm = 1 : length(dsDesc.Name)
        lblOnlyTF(knm) = all([dsDesc.(dsDesc.Name(knm)).SrcData] == "Lbl");
    end
    dsLblOnlyTF = all(lblOnlyTF);
    for knm = 1 : length(dpDesc.Name)
        lblOnlyTF(knm) = all([dpDesc.(dpDesc.Name(knm)).SrcData] == "Lbl");
    end
    dpLblOnlyTF = all(lblOnlyTF);
    lblOnlyTF = dsLblOnlyTF && dpLblOnlyTF;

    % Get the file names
    [lblpn, lblDt] = gd.dbfGetPnDt(lblp); % Get path name and datenum
    if ~lblOnlyTF
        [snlpn, snlDt] = gd.dbfGetPnDt(snlp); % Get path name and datenum
    end
    [subjNm, anStartDt, anEndDt] = getSubjInfo(lblpn, subjNmOrig);
    
    %% Data to stem
    % Initialize the table in a field of the ds structure
    for kn = 1 : numel(dsDesc.Name)
        nm = dsDesc.Name(kn); % Name of the phenomenon to analyze
        dd = dsDesc.(nm); % Data description (only for this phenomenon)
        ds.(nm) = table('Size', [0, length(dd)], 'VariableTypes', [dd.VarType], 'VariableNames', [dd.VarName]);
    end
    % Loop over label files
    fprintf(['\nLabel File No. ', num2str(0, '%06d'), '/', num2str(numel(lblpn), '%06d'), '\n'])
    for klbl = 1 : numel(lblpn)
        if rem(klbl, 20) == 0
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b') % This can delete the previous line
            fprintf(['\nLabel File No. ', num2str(klbl, '%06d'), '/', num2str(numel(lblpn), '%06d')])
            fprintf('\n')
        end
        ll = load(lblpn{klbl}); % Loaded label. All three variables from the label file will become fields of the ll structure
        for kn = 1 : numel(dsDesc.Name) % Over the names of the phenomena
            nm = dsDesc.Name(kn); % Name of the phenomenon we are now analyzing
            dd = dsDesc.(nm); % Data description (only for this phenomenon)
            % Initialize a new table which will be filled in and appended to the ds.(dsDesc.Name(kn)).
            numNewRows = sum(ismember(ll.lblSet.ClassName, dd(1).MainLbl)); % Number of rows (e.g. number of seizures in this label file)
            newRows = table('Size', [numNewRows, length(dd)], 'VariableTypes', [dd.VarType], 'VariableNames', [dd.VarName]); % Initialization of the table.
            for kchar = 1 : size(dsDesc.(nm), 2) % Over characteristics. Fill in new rows for each characteristic of the phenomenon
                d = dd(kchar); % Description of the current characteristic.
                funcHandle = str2func(d.CalcFcn); % Get function handle from the name of the function.
                colnm = d.VarName; % Column name
                switch d.SrcData % As of now, we probably do not have the sl ready. The loading of sl needs to be finished (or at least tested)
                    case "Lbl"
                        newRows.(colnm) = funcHandle(ll, d); % The function must accept the loaded label and characteristic description structure.
                    case "Snl"
                        newRows.(colnm) = funcHandle(ls, d);
                    case "LblSnl"
                        newRows.(colnm) = funcHandle(ll, ls, d);
                end
                newRows.(colnm) = funcHandle(ll, d);
            end
            ds.(nm) = [ds.(nm); newRows];
        end
    end

    %% Data to plot
    % Find out if we will need the signal files
    for knm = 1 : length(dpDesc.Name)
        lblOnlyTF(knm) = all([dpDesc.(dpDesc.Name(knm)).SrcData] == "Lbl");
    end
    lblOnlyTF = all(lblOnlyTF);
    
    % Split the time into bins
    binDt = (anStartDt : dpDesc.BinLenDu : anEndDt)'; % Borders of bins in datenum
    numbin = numel(binDt) - 1; % Number of bins
    dp.tax = binDt(2 : end);  % Each bin should be assigned the timestamp of its end because that is the moment we have had acquired (and processed) all data of the block.
    
    % Initialize the table in a field of the dp structure
    for kn = 1 : numel(dpDesc.Name)
        nm = dpDesc.Name(kn); % Name of the phenomenon to analyze
        dd = dpDesc.(nm); % Data description (only for this phenomenon)
        dp.(nm) = table('Size', [0, length(dd)], 'VariableTypes', [dd.VarType], 'VariableNames', [dd.VarName]); % Initialize with zero number of rows
    end
    
    % Loop over time bins. For each bin, it will load all the files it needs and another loop.
    fprintf(['\nBin No. ', num2str(0, '%06d'), '/', num2str(numbin, '%06d'), '\n'])
    loadedLblpn = ""; % Keep track of the currently loaded label file
    loadedSnlpn = ""; % Keep track of the currently loaded signal file
    for kb = 1 : numbin % Loop over time blocks
        if rem(kb, 20) == 0
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
            fprintf(['\nBin No. ', num2str(kb, '%06d'), '/', num2str(numbin, '%06d')])
            fprintf('\n')
        end
        tol = seconds(0.001); % Tolerance in seconds
        lblfSub = max(1, find(lblDt > binDt(kb) - tol, 1, 'first') - 1) : find(lblDt <= binDt(kb + 1), 1, 'last'); % Subscripts of label files relevant for this bin
        % Check if the label and signal files correspond
        if ~dpLblOnlyTF
            snlfSub = find(snlDt > binDt(kb) + tol, 1, 'first') - 1 : find(snlDt <= binDt(kb + 1), 1, 'last'); % Subscripts of signal file relevant for this bin
            % Check if we have the same label files and signal files
            if any(size(lblfSub) ~= size(snlfSub))
                disp(size(lblfSub))
                disp(size(snlfSub))
                warning('Label file and signal file subscripts have different size')
                pause
            elseif any(lblfSub ~= snlfSub)
                disp(lblfSub)
                disp(snlfSub)
                warning('Label file and signal file subscripts are not equal')
                pause
            end
        end

        % Initialize a new table which will be filled in and appended to the dp.(dpDesc.Name(kn)).
        for kn = 1 : numel(dpDesc.Name) % Over the names of the phenomena
            nm = dpDesc.Name(kn); % Name of the current phenomenon
            dd = dpDesc.(nm); % Structure with the description of the computations on the current phenomenon
            numRows = numel(lblfSub); % Number of rows
            binTables.(nm) = table('Size', [numRows, length(dd)], 'VariableTypes', [dd.VarType], 'VariableNames', [dd.VarName]); % Table of data from all files belonging to current time bin
        end
        
        % Loop over files within this block
        for klf = 1 : numel(lblfSub) % k-th label file (out of those relevant for this block)
            if loadedLblpn ~= string(lblpn{lblfSub(klf)}) % Check if the required data file is already loaded. If not, load it.
                ll = load(lblpn{lblfSub(klf)}, 'sigInfo', 'lblDef', 'lblSet');
                loadedLblpn = string(lblpn{lblfSub(klf)}); % Update which file is currently loaded.
            end
            % If signal data are needed, load them and check if te label and signal files correspond
            if ~dpLblOnlyTF
                load(snlpn{snlfSub(klf)}, 'sigTbl')
                snlChToProcessSub = cellfun(@(x) isempty(x), regexp(sigTbl.ChName, 'Rhd.X-\d', 'match')); % Subscripts of signal channels to be processed (e.g. we may want do ignore accelerometer channels)
                % Check if signal and label correspond
                sigTbl = sigTbl(snlChToProcessSub, :);
                sigInfo = sigInfo(snlChToProcessSub, :);
                if ~all(sigInfo.Subject == subjNm)
                    error('_jk Inconsistency of subjects in label file.')
                end
                if ~all(sigTbl.Subject == subjNm)
                    error('_jk Inconsistency of subjects in signal file.')
                end
                % TODO001 Sometimes, there is a mismatch between time extent of the signal and label file. In the future, fix the data files, so that this does
                % not happen. This is probably an issue in the creation of the data files, not of this code. Possibly make some utility to check and
                % fix the files.
                if abs(seconds(min(sigInfo.SigStart) - min(sigTbl.SigStart))) > 10 || abs(seconds(min(sigInfo.SigEnd) - min(sigTbl.SigEnd))) > 10 % If start or end times of signal and label file differ by more than 1 s
                    disp(['lblpn ', lblpn{lblfSub(klf)}, 10])
                    disp(['snlpn ', snlpn{snlfSub(klf)}, 10])
                    disp(sigInfo)
                    disp(sigTbl)
                    error('_jk Label file and signal file have different time extent')
                end
                sigInfo.SigStart = sigTbl.SigStart; % The signal file's SigStart will be used
                sigInfo.SigEnd = sigTbl.SigEnd; % The signal file's SigEnd will be used
                
                % Check if block start is after the end of given file. I believe, this should never happend unless there is a gap in the recording.
                tol = seconds(60); % Gap of 60 seconds will be tolerated
                if binDt(kb) > max(ll.sigInfo.SigEnd) + tol
                    warning(['Data missing at ', datestr(binDt(kb)), '.'])
                    disp(lblfSub)
                    disp(snlfSub)
                    disp(['Block start: ', char(binDt(kb))]); %#ok<*DATST>
                    disp(['sigInfo.SigEnd: ', char(min(ll.sigInfo.SigEnd))]);
                    disp(['Difference: ', char(min(ll.sigInfo.SigEnd) - binDt(kb))])
                    continue
                end
            end

            % Main calculation
            for kn = 1 : numel(dpDesc.Name) % Over the names of the phenomena
                nm = dpDesc.Name(kn); % Name of the current phenomenon
                dd = dpDesc.(nm); % Description of all the calculations on the current phenomenon
                for kchar = 1 : size(dpDesc.(nm), 2) % Fill in new rows for each characteristic of the phenomenon
                    d = dd(kchar); % Description of the current characteristic of the phenomenon
                    if d.CalcLvl == "file"
                        funcHandle = str2func(d.CalcFcn);
                        colnm = d.VarName; % Column name
                        numch = height(ll.sigInfo);
                        switch d.SrcData
                            case "Lbl"
                                % In contrast to ds calculation, here we include also bin limits so that the function can disregard data not belonging to the current bin
                                binTables.(nm).(colnm)(klf, 1 : numch) =...
                                    funcHandle(ll, d, [binDt(kb), binDt(kb+1)]);
                            case "Snl"
                                binTables.(nm).(colnm)(klf, 1 : numch) =...
                                    funcHandle(ls, d, [binDt(kb), binDt(kb+1)]);
                            case "LblSnl"
                                binTables.(nm).(colnm)(klf, 1 : numch) =...
                                    funcHandle(ll, ls, d, [binDt(kb), binDt(kb+1)]);
                        end
                    end
                end
            end
        end % Over files within the block. Here, we will aggregate the data from individual files, to get data for the given bin.
        for kn = 1 : numel(dpDesc.Name) % Over the names of the phenomena
            nm = dpDesc.Name(kn); % Name of the current phenomenon
            dd = dpDesc.(nm); % Description of all the calculations on the current phenomenon
            for kchar = 1 : size(dpDesc.(nm), 2) % Fill in new rows for each characteristic of the phenomenon
                d = dd(kchar); % Description of the current characteristic of the phenomenon
                switch d.CalcLvl
                    case "file"
                        if ~isempty(binTables.(nm)) % If it is not empty, sum over the first dimension (columns).
                            binTableFinal.(nm) = sum(binTables.(nm), 1);
                        else % If it is empty, the function sum would create one NaN in each cell of the table, which would be inconsistent with the dimensions of other cells if we process more than one channel
                            for kcol = 1 : width(dp.(nm))
                                nn{1, kcol} = NaN(1, numel(dp.(nm){1, kcol})); %#ok<AGROW> % Create a cell array of NaNs to plug it in the binTables.(nm)
                            end
                            binTableFinal.(nm)(end, :) = nn; % Plug in the created NaN cells into the table. The NaNs have to be there because the row exists in the time axis dp.tax.
                        end
                    case "bin"
                        funcHandle = str2func(d.CalcFcn);
                        colnm = d.VarName; % Column name
                        numch = height(ll.sigInfo); % Number of channels
                        binTableFinal.(nm).(colnm)(1, 1 : numch) = funcHandle(binTables, nm); % Each cell of the table can contain either a scalar or a row vector (if there are more channels)
                end
            end
            dp.(nm) = [dp.(nm); binTableFinal.(nm)]; % The binTableFinal has one row containg the data for current time bin. Append it to the main table dp.(nm)
        end
    end
    
    %% Get subject info
    subjInfo.subjNm = subjNm;
    subjInfo.subjNmOrig = subjInfo.subjNm;
    subjInfo.anStartDt = anStartDt;
    subjInfo.anEndDt = anEndDt;
    subjNumber = regexp(subjNm, '\D\D\d\d\d+', 'match');
    subjNumber = subjNumber{1}(3 : end);
    whichSubj = find(contains(string(dobTable{:, 1}), subjNumber));
    if ~isscalar(whichSubj)
        error(['_jk Date of birth table has multiple subjects which have ', num2str(subjNumber), ' in their name.'])
    end
    subjInfo.dob = dobTable.Birth(whichSubj);
    subjInfo.sex = dobTable{string(dobTable{:, 1}) == subjNm, 3};
    if ~stg.keepOriginalSubjectName
        subjInfo.subjNm = ['Mouse', num2str(stg.subjNumber(ksubj), '%02d')];
    end
    %% Saving of the ds and dp structures not sorted out yet
% % % %     save([stg.dataFolder, 'Data-', num2str(stg.dpBinLenS), '-', char(subjNmOrig), '.mat'], 'subjInfo', 'szCharTbl', 'szCharLabel', 'siCharTbl', 'siCharLabel')
    
    %% Nested functions - possibly put them also in +gd?
    function [subjNm, anStartDt, anEndDt] = getSubjInfo(lblpn, subjNmOrig, varargin)
        load(lblpn{1}, 'sigInfo', 'lblDef', 'lblSet') %#ok<NASGU>
        % There can be multiple subjects in one lbl3 file. Keep only channels containing the data on the subject.
        ss = strsplit(subjNmOrig, 'ET'); % ET stands for ear tag. Sometimes it is included in the subject name
        whichChannelsLbl = find(contains(sigInfo.Subject, ss{end}));
        sigInfo = sigInfo(whichChannelsLbl, :);
        % % % % % % % % % % lblSet = lblSet(ismember(lblSet.Channel, whichChannelsLbl), :);
        % Check that all rows belong to the same subject
        subjNm = sigInfo.Subject(1);
        if ~all(sigInfo.Subject == subjNm)
            disp(sigInfo)
            error('_jk Multiple subjects in label file.')
        end
        anStartDt = min(sigInfo.SigStart); % Analysis start determined by label files
        lastLbl = load(lblpn{end});
        lastSigInfo = lastLbl.sigInfo(whichChannelsLbl, :);
        anEndDt = max(lastSigInfo.SigEnd); % Analysis end determined by signal files
        if numel(varargin) == 0
            return
        else
            snlpn = varargin{1};
        end
        % Now the same with signal data (e.g. markers of critical slowing)
        load(snlpn{1}, 'sigTbl')
        whichChannelsSnl = contains(sigTbl.Subject, ss{end});
        sigTbl = sigTbl(whichChannelsSnl, :);
        subjNm = sigTbl.Subject(1);
        if ~all(sigTbl.Subject == subjNm)
            error('_jk Multiple subjects in signal file.')
        end
        % % % % % % % % % % % anStartDt = min(sigInfo.SigStart); % Analysis start determined by label files
        anStartDtSig = min(sigTbl.SigStart); % Analysis start determined by signal files
        % % % % % % % % % % % % lastLbl = load(lblpn{end});
        lastSigInfo = lastLbl.sigInfo(whichChannelsLbl, :);
        if any(sigInfo.Subject ~= lastSigInfo.Subject)
            error('_jk Last label file has diffent channels than the first file.')
        end
        lastSnl = load(snlpn{end});
        lastSigTbl = lastSnl.sigTbl(whichChannelsSnl, :);
        if any(sigTbl.Subject ~= lastSigTbl.Subject)
            error('_jk Last signal file has diffent channels than the first file.')
        end
        % % % % % % % % % % anEndDt = max(lastSigInfo.SigEnd); % Analysis end determined by signal files
        anEndDtSig = max(lastSigTbl.SigEnd); % Analysis end determined by label files
        if abs(anEndDt - anEndDtSig) > seconds(1) || abs(anStartDt - anStartDtSig) > seconds(1) % If they differ by more than a second
            disp('Analysis start difference:')
            disp((anStartDt - anStartDtSig)*3600*24)
            disp('Analysis end difference:')
            disp((anEndDt - anEndDtSig)*3600*24)
            error('_jk Label data and signal data have different time extent')
        end
    end
%% OLD NESTED FUNCTION - NEED REWRITING OR DELETING
% % % % % % % %     function szRacine = getRacine(szOnsetN, lblSet)
% % % % % % % %         szRacine = [];
% % % % % % % %         if ~isempty(szOnsetN)
% % % % % % % %             szRacine = NaN(size(szOnsetN));
% % % % % % % %             racInd = lblSet.ClassName == "RacineJK01";
% % % % % % % %             if any(racInd)
% % % % % % % %                 lblSetRac = lblSet(racInd, :);
% % % % % % % %                 for ksz = 1 : size(szOnsetN, 1)
% % % % % % % %                     whichRac = find(abs(datenum(lblSetRac.Start) - szOnsetN(ksz)) < 20/3600/24, 1);
% % % % % % % %                     if isempty(whichRac)
% % % % % % % %                         szRacine(1, ksz) = NaN;
% % % % % % % %                     else
% % % % % % % %                         szRacine(1, ksz) = lblSetRac.Value(whichRac);
% % % % % % % %                     end
% % % % % % % %                 end
% % % % % % % %             end
% % % % % % % %         end
% % % % % % % %     end
% % % % % % % %     function szPower = getPower(szOnsetN, szOffsetN, sigTbl)
% % % % % % % %         szPower = [];
% % % % % % % %         if ~isempty(szOnsetN)
% % % % % % % % %             szPower = NaN(size(sigInfo, 1), size(szOnsetN, 2));
% % % % % % % %             szPower = NaN(stg.numEegCh, size(szOnsetN, 2));
% % % % % % % %             for kc = 1 : stg.numEegCh
% % % % % % % %                 [stg.flt.num, stg.flt.den] = butter(stg.flt.szN, [stg.flt.szF1, stg.flt.szF2]/(sigTbl.Fs(kc)/2)); % Get the filter which will be used to remove slow waves and EMG from the seizure signal
% % % % % % % %                 snlSz = filtfilt(stg.flt.num, stg.flt.den, double(sigTbl.Data{kc})); % Filter the whole signal to (hopefully) avoid artifacts at the beginning and end of the seizure
% % % % % % % %                 for ksz = 1 : size(szOnsetN, 1)
% % % % % % % %                     szStartSub = floor((szOnsetN(ksz) - datenum(sigTbl.SigStart(kc)))*3600*24*sigTbl.Fs(kc)) + 1;
% % % % % % % %                     szEndSub = floor((szOffsetN(ksz) - datenum(sigTbl.SigStart(kc)))*3600*24*sigTbl.Fs(kc));
% % % % % % % %                     if szEndSub > numel(snlSz)
% % % % % % % %                         numSamp2 = (szOffsetN(ksz) - snlN(snlfSub(klf) + 1)) * 3600*24*sigTbl.Fs(kc); % Number of samples that need to be taken from the second signal file
% % % % % % % %                         if numSamp2 > 0
% % % % % % % %                             l2 = load(snlpn{snlfSub(klf) + 1}); % Get next file. Not very efficient but this situation is rare.
% % % % % % % %                             snlSz2 = filtfilt(stg.flt.num, stg.flt.den, double(l2.sigTbl.Data{kc}));
% % % % % % % %                             snlSz2 = snlSz2(1 : numSamp2);
% % % % % % % %                             szSnl = [snlSz, snlSz2];
% % % % % % % %                         else
% % % % % % % %                             szEndSub = numel(snlSz);
% % % % % % % %                             szSnl = snlSz(szStartSub : szEndSub);
% % % % % % % %                         end
% % % % % % % %                     else
% % % % % % % %                         szSnl = snlSz(szStartSub : szEndSub);
% % % % % % % %                     end
% % % % % % % %                     szPower(kc, ksz) = sum(szSnl.*szSnl)/size(szSnl, 2);
% % % % % % % %                 end
% % % % % % % %             end
% % % % % % % %         end
% % % % % % % %     end
% % % % % % % %     function postIctPower = getPostIctPower(szOffsetN, sigTbl)
% % % % % % % %         postIctDurS = 5;
% % % % % % % %         if sigTbl.Subject(1) == "jc20190313_2" && abs(sigTbl.SigStart(1) - datetime('190417_165755', 'InputFormat', 'yyMMdd_HHmmss')) < seconds(10)
% % % % % % % %             postIctDurS = 3;
% % % % % % % %         end
% % % % % % % %         postIctPower = [];
% % % % % % % %         if ~isempty(szOffsetN)
% % % % % % % % %             postIctPower = NaN(size(sigInfo, 1), size(szOffsetN, 2));
% % % % % % % %             postIctPower = NaN(stg.numEegCh, size(szOffsetN, 2));
% % % % % % % %             for kc = 1 : stg.numEegCh
% % % % % % % %                 [stg.flt.num, stg.flt.den] = butter(stg.flt.szN, [stg.flt.szF1, stg.flt.szF2]/(sigTbl.Fs(kc)/2)); % Get the filter which will be used to remove slow waves and EMG from the seizure signal
% % % % % % % %                 snlSz = filtfilt(stg.flt.num, stg.flt.den, double(sigTbl.Data{kc})); % Filter the whole signal to (hopefully) avoid artifacts at the beginning and end of the seizure
% % % % % % % %                 for ksz = 1 : size(szOffsetN, 1)
% % % % % % % %                     postIctStartSub = floor((szOffsetN(ksz) - datenum(sigTbl.SigStart(kc)))*3600*24*sigTbl.Fs(kc)) + 1;
% % % % % % % %                     postIctEndSub = floor((szOffsetN(ksz) - datenum(sigTbl.SigStart(kc)))*3600*24*sigTbl.Fs(kc) + postIctDurS*sigTbl.Fs(kc));
% % % % % % % % %                     if postIctEndSub > numel(snlSz)
% % % % % % % % %                         l2 = load(snlpn{snlfSub(klf) + 1}); % Get next file. Not very efficient but this situation is rare.
% % % % % % % % %                         snlSz2 = filtfilt(stg.flt.num, stg.flt.den, double(l2.sigTbl.Data{kc}));
% % % % % % % % %                         snlSz = [snlSz, snlSz2];
% % % % % % % % %                     end
% % % % % % % % %                     szSnl = snlSz(postIctStartSub : postIctEndSub);
% % % % % % % % %                     postIctPower(kc, ksz) = sum(szSnl.*szSnl)/size(szSnl, 2);
% % % % % % % % 
% % % % % % % % 
% % % % % % % %                     if postIctEndSub > numel(snlSz)
% % % % % % % %                         numSamp2 = (szOffsetN(ksz) + postIctDurS/3600/24 - snlN(snlfSub(klf) + 1)) * 3600*24*sigTbl.Fs(kc); % Number of samples that need to be taken from the second signal file
% % % % % % % %                         if numSamp2 > 0
% % % % % % % %                             l2 = load(snlpn{snlfSub(klf) + 1}); % Get next file. Not very efficient but this situation is rare.
% % % % % % % %                             snlSz2 = filtfilt(stg.flt.num, stg.flt.den, double(l2.sigTbl.Data{kc}));
% % % % % % % %                             snlSz2 = snlSz2(1 : numSamp2);
% % % % % % % %                             szSnl = [snlSz, snlSz2];
% % % % % % % %                         else
% % % % % % % %                             postIctEndSub = numel(snlSz);
% % % % % % % %                             szSnl = snlSz(postIctStartSub : postIctEndSub);
% % % % % % % %                         end
% % % % % % % %                     else
% % % % % % % %                         szSnl = snlSz(postIctStartSub : postIctEndSub);
% % % % % % % %                     end
% % % % % % % %                     postIctPower(kc, ksz) = sum(szSnl.*szSnl)/size(szSnl, 2);
% % % % % % % %                 end
% % % % % % % %             end
% % % % % % % %         end
% % % % % % % %     end
% % % % % % % %     function numMrk = countMrk(lblSetMrk, sigStartN, sigEndN, blStart, blEnd)
% % % % % % % %         lblSetMrk = sortrows(lblSetMrk, "Start"); % Sort them
% % % % % % % %         iedStartN = datenum(lblSetMrk.Start);
% % % % % % % %         remInd = [false; diff(iedStartN) < stg.minIedSepS/3600/24]; % Get indices to...
% % % % % % % %         lblSetMrk(remInd, :) = []; % ...remove IEDs too early after previous one. It deletes possible duplicates or polyspikes.
% % % % % % % %         lblSetMrk = lblSetMrk(iedStartN >= sigStartN & iedStartN <= sigEndN, :); % Remove markers belonging to other signal files (this could happen due to a bug);
% % % % % % % %         lblSetMrk = lblSetMrk(iedStartN > blStart & iedStartN < blEnd, :);
% % % % % % % %         numMrk = size(lblSetMrk, 1);
% % % % % % % %     end
% % % % % % % %     function [mrkOnN, mrkOffN] = tToOON(t, tStartN, tFs)
% % % % % % % %         t = [0; double(t(:)); 0];
% % % % % % % %         mrkOnN = (find(diff(t) == 1) - 1)/tFs/3600/24 + tStartN;
% % % % % % % %         mrkOffN = (find(diff(t) == -1) - 1)/tFs/3600/24 + tStartN;
% % % % % % % %     end
end
