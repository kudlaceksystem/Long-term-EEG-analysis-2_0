function [subjInfo, ds, dp] = getData(stg, dsDesc, dpDesc, lblp, snlp, dobTable, ksubj, subjNmOrig)
    % dsDesc ........ data to stem description
    % dpDesc ........ data to plot description
    % lblp .......... path to label files
    % snlp .......... path to signal files
    % dobTable ...... tabel containing dates of birth (and possibly sex and other data)
    % ksubj ......... which subject are we processing now
    % subjNmOrig .... original subject name (I do not know why it is needed here

    % Get sz data and signal characteristics such as IED rate, amount of EMG artifacts and critical slowing markers.
    % % % % % % % global stg
    
    %% NEEDS REWRITING If the data for this subject already exist load it
    % if exist([stg.dataFolder, 'Data-', num2str(seconds(dpDesc.BinLenDu)), '-', char(subjNmOrig), '.mat'], 'file')
    %     load([stg.dataFolder, 'Data-', num2str(seconds(dpDesc.BinLenDu)), '-', char(subjNmOrig), '.mat'],...
    %         'subjInfo', 'szCharTbl', 'szCharLabel', 'siCharTbl', 'siCharLabel')
    %     stg.szCharYLbl = szCharLabel;
    %     stg.siCharYLbl = siCharLabel;
    %     stg.saCharYLbl = [szCharLabel, siCharLabel];
    %     stg.ssCharYLbl = [siCharLabel, siCharLabel];
    %     if isa(subjInfo.dob, 'table')
    %         subjInfo.dob = datenum(subjInfo.dob{1, 1});
    %     end
    %     subjInfo.subjNmOrig = subjInfo.subjNm;
    %     if ~stg.keepOriginalSubjectName
    %         subjInfo.subjNmOrig = subjInfo.subjNm;
    %         subjInfo.subjNm = ['Mouse', num2str(stg.subjNumber(ksubj), '%02d')];
    %     end
    %     subjInfo.sex = dobTable{ksubj, 3};
    %     return
    % end

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
    [lblpn, lblDt] = gd.dbfGetPnDt(stg, lblp); % Get path name and datenum
    if ~lblOnlyTF
        [snlpn, snlDt] = gd.dbfGetPnDt(stg, snlp); % Get path name and datenum
    end
    [subjNm, anStartDt, anEndDt, chName] = getSubjInfo(stg, lblpn, subjNmOrig);
    
    %% Data to stem
    % Initialize the table in a field of the ds structure
    for kn = 1 : numel(dsDesc.Name)
        nm = dsDesc.Name(kn); % Name of the phenomenon to analyze
        dd = dsDesc.(nm); % Data description (only for this phenomenon)
        ds.(nm) = table('Size', [0, length(dd)], 'VariableTypes', [dd.VarType], 'VariableNames', [dd.VarName]);
        ds.(nm) = helper.tblSetTimeZone(ds.(nm), stg.timeZoneStr);
        ds.(nm) = helper.tblSetTimeZone(ds.(nm), "UTC");
    end
    % Loop over label files
    fprintf(['\nLabel File No. ', num2str(0, '%06d'), '/', num2str(numel(lblpn), '%06d'), '\n'])
    clear prevEnd
    for klbl = 1 : numel(lblpn)
        if rem(klbl, 20) == 0
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b') % This can delete the previous line
            fprintf(['\nLabel File No. ', num2str(klbl, '%06d'), '/', num2str(numel(lblpn), '%06d')])
            fprintf('\n')
        end
lblpn{klbl}
        ll = load(lblpn{klbl}) % Loaded label. All three variables from the label file will become fields of the ll structure
        ll.sigInfo = helper.tblSetTimeZone(ll.sigInfo, stg.timeZoneStr);
        ll.sigInfo = helper.tblSetTimeZone(ll.sigInfo, "UTC");
        ll.lblSet = helper.tblSetTimeZone(ll.lblSet, stg.timeZoneStr);
        ll.lblSet = helper.tblSetTimeZone(ll.lblSet, "UTC");

        % Keep only channels belonging to this animal
        chToKeep = find(ll.sigInfo.Subject == string(subjNmOrig));
        ll.sigInfo = ll.sigInfo(chToKeep, :);
        ll.lblSet = ll.lblSet(ismember(ll.lblSet.Channel, chToKeep), :);

        % Check if new file begins after the end of the previous file
        if exist('prevSigInfo', 'var')
            siginfo_ = ll.sigInfo.SigStart
            prevsiginfo_ = prevSigInfo.SigEnd
            currentStartMinusPrevEnd = ll.sigInfo.SigStart(1) - prevSigInfo.SigEnd(1);
            if currentStartMinusPrevEnd < minutes(-40)
                klbl_ = klbl
                prev_lblpn_ = lblpn{klbl-1}
                prevSigInfo_ = prevSigInfo
                prevTz_ = prevSigInfo.SigStart.TimeZone
                lblpn_klbl_ = lblpn{klbl}
                thisSigInfo_ = ll.sigInfo
                thisTz_ = ll.sigInfo.SigStart.TimeZone
                currentStartMinusPrevEnd_ = currentStartMinusPrevEnd
                'now change the time'
                ll.sigInfo.SigStart = ll.sigInfo.SigStart + hours(1);
                ll.sigInfo.SigEnd = ll.sigInfo.SigEnd + hours(1);
                newTime_ = ll.sigInfo
                pause
            end
        end
        prevSigInfo = ll.sigInfo;

        % Check if sigInfo.FileName corresponds to file name
        fndattimStr = regexp(lblpn{klbl}, '\d\d\d\d\d\d_\d\d\d\d\d\d', 'match');
        filenameDt = datetime(fndattimStr{1}, 'InputFormat', 'yyMMdd_HHmmss', 'TimeZone', stg.timeZoneStr);
        filenameDt.TimeZone = "UTC";
        fcdattimStr = regexp(ll.sigInfo.FileName, '\d\d\d\d\d\d_\d\d\d\d\d\d', 'match');
        filecontDt = datetime(fcdattimStr{1}, 'InputFormat', 'yyMMdd_HHmmss', 'TimeZone', stg.timeZoneStr);
        filecontDt.TimeZone = "UTC";
        if filenameDt ~= filecontDt
                klbl
                lblpn{klbl}
                ll.sigInfo
                pause
        end

        % Check if SigStart corresponds to file name
        fndattimStr = regexp(lblpn{klbl}, '\d\d\d\d\d\d_\d\d\d\d\d\d', 'match');
        filenameDt = datetime(fndattimStr{1}, 'InputFormat', 'yyMMdd_HHmmss', 'TimeZone', stg.timeZoneStr);
        filenameDt.TimeZone = "UTC";
        % % % % % % % % % % fcdattimStr = regexp(ll.sigInfo.FileName, '\d\d\d\d\d\d_\d\d\d\d\d\d', 'match');
        filecontDt = ll.sigInfo.SigStart;
        if filenameDt ~= filecontDt
                klbl
                lblpn{klbl}
                ll.sigInfo
                pause
        end
        
        for kn = 1 : numel(dsDesc.Name) % Over the names of the phenomena
            nm = dsDesc.Name(kn); % Name of the phenomenon we are now analyzing
            dd = dsDesc.(nm); % Data description (only for this phenomenon)
            % Initialize a new table which will be filled in and appended to the ds.(dsDesc.Name(kn)).

            numNewRows = sum(ismember(ll.lblSet.ClassName, dd(1).MainLbl)); % Number of rows (e.g. number of seizures in this label file)
            newRows = table('Size', [numNewRows, length(dd)], 'VariableTypes', [dd.VarType], 'VariableNames', [dd.VarName]); % Initialization of the table.
            newRows = helper.tblSetTimeZone(newRows, stg.timeZoneStr);
            newRows = helper.tblSetTimeZone(newRows, "UTC");
            for kchar = 1 : size(dsDesc.(nm), 2) % Over characteristics. Fill in new rows for each characteristic of the phenomenon
                d = dd(kchar); % Description of the current characteristic.
                funcHandle = str2func(d.CalcFcn); % Get function handle from the name of the function.
                colnm = d.VarName; % Column name
                switch d.SrcData % As of now, we probably do not have the ls ready. The loading of ls needs to be finished (or at least tested)
                    case "Lbl"
                        y = funcHandle(ll, d);
                        newRows.(colnm)(1 : numel(y)) = y; % The function must accept the loaded label and characteristic description structure.
                    case "Snl"
                        newRows.(colnm) = funcHandle(ls, d);
                    case "LblSnl"
                        newRows.(colnm) = funcHandle(ll, ls, d);
                end
            end
            ds.(nm) = [ds.(nm); newRows];
        end
    end
    for kn = 1 : numel(dsDesc.Name) % Over the names of the phenomena
        nm = dsDesc.Name(kn); % Name of the phenomenon we are now analyzing
        ds.(nm) = ds.(nm)(~isnat(ds.(nm){:, 1}), :);
        % Remove duplicates
        [~, sub] = unique(ds.(nm){:, 1});
        ds.(nm) = ds.(nm)(sub, :);
        % TODO006 Merge too close events here again (they may be in different files, yet still too close
        if ~any(strcmp(ds.(nm).Properties.VariableNames, 'DurDu'))
            ds.(nm).DurDu = seconds(zeros(height(ds.(nm)), 1));
        end
        if all(ds.(nm).DurDu == seconds(0))
            pointTF = true;
        else
            pointTF = false;
        end
        ds.(nm) = sortrows(ds.(nm), "OnsDt");
        onsDt = ds.(nm).OnsDt;
        offDt = ds.(nm).OnsDt + ds.(nm).DurDu;
        nummrk = height(ds.(nm)); % Number of markers
        mergeWithPreviousTF = false(nummrk, 1);
        for k = 2 : nummrk
            d = onsDt(k) - offDt(1 : k - 1); % Current marker separation from all the previous
            mergeWithPreviousTF(k) = any(d < seconds(dsDesc.(nm)(1).MinSepS)); % Plugs in true the current marker was insufficiently separated from any of the previous
        end
        mergeWithPreviousSub = find(mergeWithPreviousTF);
        nummrg = numel(mergeWithPreviousSub); % Number of markers to merge with previous
        for k = 1 : nummrg
            toBeMerged = mergeWithPreviousSub(nummrg - k + 1);
            mergedWith = mergeWithPreviousSub(nummrg - k + 1) - 1;
            offDt(mergedWith) = max(offDt(mergedWith), offDt(toBeMerged));
            offDt(toBeMerged) = [];
            onsDt(toBeMerged) = [];
            ds.(nm)(toBeMerged, :) = [];
        end
        ds.(nm).OnsDt = onsDt;
        ds.(nm).DurDu = offDt - onsDt;
        if pointTF
            ds.(nm).DurDu = seconds(zeros(height(ds.(nm)), 1));
        end
    end
    
    %% Data to plot
    % Find out if we will need the signal files
    for knm = 1 : length(dpDesc.Name)
        lblOnlyTF(knm) = all([dpDesc.(dpDesc.Name(knm)).SrcData] == "Lbl");
    end
    lblOnlyTF = all(lblOnlyTF);

    %% Which bin lengths are there in the dpDesc?
    binlenAllDu = seconds(zeros(numel(dpDesc.Name), 1));
    for kn = 1 : numel(dpDesc.Name)
        binlenAllDu(kn) = dpDesc.(dpDesc.Name(kn))(1).BinLenDu;
    end
    binlenUnDu = unique(binlenAllDu);
    
    for kbinlen = 1 : numel(binlenUnDu)
        % Split the time into bins
        binDt = (anStartDt : binlenUnDu(kbinlen) : anEndDt)'; % Edges of bins in datenum
        numbin = numel(binDt) - 1; % Number of bins
        
        % Initialize the table in a field of the dp structure
        for kn = 1 : numel(dpDesc.Name)
            nm = dpDesc.Name(kn); % Name of the phenomenon to analyze
            dd = dpDesc.(nm); % Data description (only for this phenomenon)
            if dd(1).BinLenDu == binlenUnDu(kbinlen)
                dp.(nm) = table('Size', [0, length(dd)], 'VariableTypes', [dd.VarType], 'VariableNames', [dd.VarName]); % Initialize with zero number of rows
                dp.(nm) = helper.tblSetTimeZone(dp.(nm), stg.timeZoneStr);
                dp.(nm) = helper.tblSetTimeZone(dp.(nm), "UTC");
            end
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
                if dd(1).BinLenDu == binlenUnDu(kbinlen)
                    numRows = numel(lblfSub); % Number of rows
                    binTables.(nm) = table('Size', [numRows, length(dd)], 'VariableTypes', [dd.VarType], 'VariableNames', [dd.VarName]); % Table of data from all files belonging to current time bin
                    binTables.(nm) = helper.tblSetTimeZone(binTables.(nm), stg.timeZoneStr);
                    binTables.(nm) = helper.tblSetTimeZone(binTables.(nm), "UTC");
                end
            end
            
            % Get channel names
            ll = load(lblpn{1}, 'sigInfo');
            % Keep only channels belonging to this animal
            chToKeep = find(ll.sigInfo.Subject == string(subjNmOrig));
            ll.sigInfo = ll.sigInfo(chToKeep, :);
            channelNames = ll.sigInfo.ChName;
            clear ll

            % Loop over files within this block
            for klf = 1 : numel(lblfSub) % k-th label file (out of those relevant for this block)
                if loadedLblpn ~= string(lblpn{lblfSub(klf)}) % Check if the required data file is already loaded. If not, load it.
                    ll = load(lblpn{lblfSub(klf)}, 'sigInfo', 'lblDef', 'lblSet');
                    ll.sigInfo = helper.tblSetTimeZone(ll.sigInfo, stg.timeZoneStr);
                    if diff(isdst(ll.sigInfo.SigStart(1), ll.sigInfo.SigStart(1) + hours(1)))
                        'Daylight fuck'
                        ll.sigInfo
                        pause
                    end
                    ll.sigInfo = helper.tblSetTimeZone(ll.sigInfo, "UTC");
                    ll.lblSet = helper.tblSetTimeZone(ll.lblSet, stg.timeZoneStr);
                    ll.lblSet = helper.tblSetTimeZone(ll.lblSet, "UTC");

                    % Keep only channels belonging to this animal
                    chToKeep = find(ll.sigInfo.Subject == string(subjNmOrig));
                    ll.sigInfo = ll.sigInfo(chToKeep, :);
                    ll.lblSet = ll.lblSet(ismember(ll.lblSet, chToKeep), :);
                    % Check channel names
                    if numel(ll.sigInfo.ChName) ~= numel(channelNames)
                        error('_jk getData: Number of channels inconsistent.')
                    end
                    if ~all(ll.sigInfo.ChName == channelNames)
                        error('_jk getData: Channel order inconsistent.')
                    end
                    % Update which file is currently loaded.
                    loadedLblpn = string(lblpn{lblfSub(klf)});
                end
                % If signal data are needed, load them and check if the label and signal files correspond
                if ~dpLblOnlyTF
                    load(snlpn{snlfSub(klf)}, 'sigTbl')
                    % Solve the daylight saving time issues here as well
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
                    if dd(1).BinLenDu == binlenUnDu(kbinlen)
                        for kchar = 1 : size(dpDesc.(nm), 2) % Fill in new rows for each characteristic of the phenomenon
                            d = dd(kchar); % Description of the current characteristic of the phenomenon
                            if d.CalcLvl == "file"
                                funcHandle = str2func(d.CalcFcn(1));
                                colnm = d.VarName; % Column name
                                % % % % % % % numch = height(ll.sigInfo);
                                switch d.SrcData
                                    case "Lbl"
                                        % In contrast to ds calculation, here we include also bin limits so that the function can disregard data not belonging to the current bin
                                        y = funcHandle(ll, d, [binDt(kb), binDt(kb+1)]);
                                    case "Snl"
                                        y = funcHandle(ls, d, [binDt(kb), binDt(kb+1)]);
                                    case "LblSnl"
                                        y = funcHandle(ll, ls, d, [binDt(kb), binDt(kb+1)]);
                                end
                                binTables.(nm).(colnm)(klf, 1 : numel(y)) = y;
                            end
                        end
                    end
                end
            end % Over files within the block. Here, we will aggregate the data from individual files, to get data for the given bin.
            for kn = 1 : numel(dpDesc.Name) % Over the names of the phenomena
                nm = dpDesc.Name(kn); % Name of the current phenomenon
                dd = dpDesc.(nm); % Description of all the calculations on the current phenomenon
                if dd(1).BinLenDu == binlenUnDu(kbinlen)
                    binTableFinal.(nm) = binTables.(nm)([], :); % Initialize
                    for kchar = 1 : size(dpDesc.(nm), 2) % Fill in new rows for each characteristic of the phenomenon
                        d = dd(kchar); % Description of the current characteristic of the phenomenon
                        colnm = d.VarName; % Column name
                        switch d.CalcLvl
                            case "file"
                                if ~isempty(binTables.(nm)) % If it is not empty, sum over the first dimension.
                                    funcHandle = str2func(d.CalcFcn(2));
                                    warning('off', 'MATLAB:table:RowsAddedExistingVars')
                                    y = funcHandle(binTables.(nm).(colnm), 1);
                                    binTableFinal.(nm).(colnm)(1, 1 : numel(y)) = y;
                                    warning('on', 'MATLAB:table:RowsAddedExistingVars')
                                else % If it is empty, the function sum would create one NaN in each cell of the table. If we process >1 channel, a single NaN would be inconsistent with the dimensions of other cells.
                                    warning('off', 'MATLAB:table:RowsAddedExistingVars')
                                    y = NaN(1, numel(dp.(nm).(colnm)(1, :)));
                                    binTableFinal.(nm).(colnm)(1, 1 : numel(y)) = y;
                                    warning('on', 'MATLAB:table:RowsAddedExistingVars')
                                end
                            case "bin"
                                funcHandle = str2func(d.CalcFcn);
                                numch = height(ll.sigInfo); % Number of channels
                                warning('off', 'MATLAB:table:RowsAddedExistingVars')
                                y = funcHandle(binTables, nm); % Each cell of the table can contain either a scalar or a row vector (if there are more channels)
                                binTableFinal.(nm).(colnm)(1, 1 : numel(y)) = y;
                                warning('on', 'MATLAB:table:RowsAddedExistingVars')
                        end
                    end
                    validS = binTableFinal.(nm).ValidS
                    if validS > seconds(dpDesc.(nm)(1).BinLenDu)
                        y
                        validS
                        pause
                    end
                    dp.(nm) = [dp.(nm); binTableFinal.(nm)]; % The binTableFinal has one row containg the data for current time bin. Append it to the main table dp.(nm)
                end
            end
        end
        for kn = 1 : numel(dpDesc.Name) % Over the names of the phenomena
            nm = dpDesc.Name(kn); % Name of the current phenomenon
            dd = dpDesc.(nm); % Description of all the calculations on the current phenomenon
            if dd(1).BinLenDu == binlenUnDu(kbinlen)
                dp.(nm).tax = binDt(2 : end);  % Each bin should be assigned the timestamp of its end because that is the moment we have had acquired (and processed) all data of the block.
                dp.(nm) = movevars(dp.(nm), 'tax', 'Before', 1);
            end
        end
    end

    %% Get subject info
    subjInfo.ksubj = ksubj;
    subjInfo.subjNm = subjNm;
    subjInfo.subjNmOrig = subjInfo.subjNm;
    subjInfo.anStartDt = anStartDt;
    subjInfo.anEndDt = anEndDt;
    subjInfo.chName = chName;
    subjNumber = regexp(subjNm, '\D\D\d\d\d+', 'match');
    subjNumber = subjNumber{1}(3 : end);
    whichSubj = find(contains(string(dobTable{:, 1}), subjNumber));
    if isempty(whichSubj)
        error(['_jk Mouse ', num2str(subjNumber), ' not found in the table of dates or birth.'])
    end
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
    function [subjNm, anStartDt, anEndDt, chName] = getSubjInfo(stg, lblpn, subjNmOrig, varargin)
        lll = load(lblpn{1}, 'sigInfo', 'lblDef', 'lblSet');
        lll.sigInfo = helper.tblSetTimeZone(lll.sigInfo, stg.timeZoneStr);
        lll.sigInfo = helper.tblSetTimeZone(lll.sigInfo, "UTC");
        % There can be multiple subjects in one lbl3 file. Keep only channels containing the data on the subject.
        ss = strsplit(subjNmOrig, 'ET'); % ET stands for ear tag. Sometimes it is included in the subject name
        whichChannelsLbl = find(contains(lll.sigInfo.Subject, ss{end}));
        sigInfo = lll.sigInfo(whichChannelsLbl, :);
        chName = lll.sigInfo.ChName;

        % % % % % % % % % % lblSet = lblSet(ismember(lblSet.Channel, whichChannelsLbl), :);
        % Check that all rows belong to the same subject
        subjNm = lll.sigInfo.Subject(1);
        if ~all(lll.sigInfo.Subject == subjNm)
            disp(lll.sigInfo)
            error('_jk Multiple subjects in label file.')
        end
        anStartDt = min(lll.sigInfo.SigStart); % Analysis start determined by label files
        lastLbl = load(lblpn{end}, 'sigInfo');
        lastLbl.sigInfo = helper.tblSetTimeZone(lastLbl.sigInfo, stg.timeZoneStr);
        lastLbl.sigInfo = helper.tblSetTimeZone(lastLbl.sigInfo, "UTC");
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
        if any(lll.sigInfo.Subject ~= lastSigInfo.Subject)
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
