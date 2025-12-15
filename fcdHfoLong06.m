% This script takes label files and optionally signal files and creates long-term profiles of labeled patterns and transients
analyzeIndividualSubjects = 1;
analyzePopulation = 1;
if analyzeIndividualSubjects
    % close all
    % clear
    analyzeIndividualSubjects = 1;
    analyzePopulation = 1;
end


%% ANALYSIS
% TODO001 Different time extents of label and signal files. Needs to be fixed in the data not in this script.
% TODO002 If file is already loaded do not load it again
% TODO003 Change from struct to table in getClusters
% TODO004 DO SUBJECT STATS LATER
% TODO005 get rid of global variables
% In siCharCl fix the y-axis labels
% Add before-, during- and after-cluster raw IED rate and possibly also sz chars
% Add violin plots of baselines after the exponentials
% P-values and significance in individual subjects
% Evaluate the trends as fold change?
% Special functions for fitting Poisson or power-law distribution to the signal characteristics
% Circadian distribution of lead seizures or cluster onsets and offsets
% Tukey
% HSD
% Add first seizure occurrence time in the descriptiveStats
% Try different cluster definition settings
% How many droupouts, how long, how long in total?

%% PLOT FORMATTING
% ShowStat simulated similarity

%% CODE CLEANLINESS
% Naming of the phase and radius or angle and modulus or theta and R
% Input into function should never be fields of a structure. Input the whole structure and choose fields within the function.
% Call sample "sample" and not "population"?
% Get rid of global variables

%% DATA

%% ANALYSIS IDEAS
% Forecasting
% Add seizure sizes according to Osorio 2010.
% Osorio 2010: Omori law. Is the IED rate decay after seizures exponential or power-law?
% Compute skewness
% Analyze all Isa's signal features
% Does sz rate differ between males and females? Any other differences? Is there estral cycle-related cyclicity in females?
% Add Bohdana's mice, e.g. 867
% Sleep analysis (CVUT student?)
% How to statistically analyze the circadian distribution

%% %%%%%%%%%% %%
%% CODE START %%
%% %%%%%%%%%% %%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%% NEW SETTINGS %%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
global stg % TODO005 get rid of global variables

%% Select subjects
subjList = {'BH002390'};
path0 = 'r:\Kudlacek\FCD HFO\HFO long-term profile'; % Without '\' at the end
path1 = {
    'Testing snl and lbl'
};
subjToPlot = {
        'BH002390';
};
pathEeg3 = {
    'BH002390_smrx converted data full'
};
pathLbl3 = {
    'BH002390_label_TEST with dropout'
};


%% Data to stem
dsDesc.Name = "Seizure"; % Types of tables that will be created. Typically, each table belongs to one type of event (e.g. seizure, sleep epoch, some behavioral event)
exLblAn = [];
minSepSzS = 60;
d(1).VarName    = "OnsDt";
d(1).VarType    = "datetime";
d(1).CalcFcn    = "gd.dsfGetOnsDt";
d(1).SrcData    = "Lbl";
d(1).MainLbl    = ["Seizure", "seizure", "SEIZURE", "S", "s"];
d(1).ExLblAn    = exLblAn; % Labels to exclude in all channels if present in any
d(1).MinSepS    = minSepSzS;
d(1).PlotTitle  = "Seizure occurrence";
d(1).YAxisLabel = "";
d(2).VarName    = "DurDu";
d(2).VarType    = "duration";
d(2).CalcFcn    = "gd.dsfGetDurDu";
d(2).SrcData    = "Lbl";
d(2).MainLbl    = ["Seizure", "seizure", "SEIZURE", "S", "s"];
d(2).ExLblAn    = exLblAn; % Labels to exclude in all channels if present in any
d(2).MinSepS    = minSepSzS;
d(2).PlotTitle  = "Seizure duration";
d(2).YAxisLabel = "Sz dur (s)";
d(3).VarName    = "Pow";
d(3).VarType    = "double";
d(3).CalcFcn    = "gd.dsfGetPow";
d(3).SrcData    = "Lbl";
d(3).MainLbl    = ["Seizure", "seizure", "SEIZURE", "S", "s"];
d(3).ExLblAn    = exLblAn; % Labels to exclude in all channels if present in any
d(3).MinSepS    = minSepSzS;
d(3).PlotTitle  = "Seizure signal power";
d(3).YAxisLabel = "Sz power (a.u.)";
dsDesc.(dsDesc.Name(1)) = d;
clear d

%% Data to plot
% % % dpDesc.BinLenDu = seconds(6*3600);
dpDesc.Name = ["Seizure21600"; "Ied3600"]; % Tables that will be created. Typically, each table belongs to one characteristic of signal (e.g. IED rate, mean IED amplitude, signal power, delta/theta ratio, etc.)

% Seizure
binlenDu = seconds(6*3600);
szLbl = ["Seizure", "seizure", "SEIZURE", "S"];
exLblCh = ""; % Has to be string - use empty string
exLblAn = "";
d(1).VarName    = "ValidS";
d(1).VarType    = "double";
d(1).BinLenDu   = binlenDu;
d(1).CalcLvl    = "file";
d(1).CalcFcn    = ["gd.dpfGetValidAmount", "sum"]; % If ClcLvl is "file", specify function to apply on individual files and function to merge data from multiple files to bin
d(1).SrcData    = "Lbl";
d(1).MainLbl    = szLbl;
d(1).ExLblCh    = exLblCh; % Labels to exclude in individual channels
d(1).ExLblAn    = exLblAn; % Labels to exclude in all channels if present in any
d(1).MinSepS    = minSepSzS;
d(1).PlotTitle  = "Total duration of usable rec";
d(1).YAxisLabel = "Usable rec (s)";
d(2).VarName    = "Count";
d(2).VarType    = "double";
d(2).BinLenDu   = binlenDu;
d(2).CalcLvl    = "file";
d(2).CalcFcn    = ["gd.dpfGetCount", "sum"]; % If ClcLvl is "file", specify function to apply on individual files and function to merge data from multiple files to bin
d(2).SrcData    = "Lbl";
d(2).MainLbl    = szLbl;
d(2).ExLblCh    = exLblCh; % Labels to exclude in individual channels
d(2).ExLblAn    = exLblAn; % Labels to exclude in all channels if present in any
d(2).MinSepS    = minSepSzS;
d(2).PlotTitle  = "Sz count";
d(2).YAxisLabel = "Sz count";
d(3).VarName    = "RatePh";
d(3).VarType    = "double";
d(3).BinLenDu   = binlenDu;
d(3).CalcLvl    = "bin";
d(3).CalcFcn    = "gd.dpbGetRatePh"; % If ClcLvl is "file", specify function to apply on individual files and function to merge data from multiple files to bin
d(3).SrcData    = "Lbl";
d(3).MainLbl    = szLbl;
d(3).ExLblCh    = exLblCh; % Labels to exclude in individual channels
d(3).ExLblAn    = exLblAn; % Labels to exclude in all channels if present in any
d(3).MinSepS    = minSepSzS;
d(3).PlotTitle  = "Sz rate";
d(3).YAxisLabel = "Sz/hour";
dpDesc.(dpDesc.Name(1)) = d;
clear d

% Ied
binlenDu = seconds(3600);
exLblCh = ["art", "Art", "EMG", "emg", "Emg"];
exLblAn = ["Seizure", "seizure", "SEIZURE", "S"];
minSepIedS = 0.1;
d(1).VarName    = "ValidS";
d(1).VarType    = "double";
d(1).BinLenDu   = binlenDu;
d(1).CalcLvl    = "file";
d(1).CalcFcn    = ["gd.dpfGetValidAmountCh", "sum"];
d(1).SrcData    = "Lbl";
d(1).MainLbl    = "IED_Janca";
d(1).ExLblCh    = exLblCh; % Labels to exclude in individual channels
d(1).ExLblAn    = exLblAn; % Labels to exclude in all channels if present in any
d(1).MinSepS    = minSepIedS;
d(1).PlotTitle  = "Total duration of usable rec";
d(1).YAxisLabel = "Usable rec (s)";
d(2).VarName    = "Count";
d(2).VarType    = "double";
d(2).BinLenDu   = binlenDu;
d(2).CalcLvl    = "file";
d(2).CalcFcn    = ["gd.dpfGetCountCh", "sum"];
d(2).SrcData    = "Lbl";
d(2).MainLbl    = "IED_Janca";
d(2).ExLblCh    = exLblCh; % Labels to exclude in individual channels
d(2).ExLblAn    = exLblAn; % Labels to exclude in all channels if present in any
d(2).MinSepS    = minSepIedS;
d(2).PlotTitle  = "IED count";
d(2).YAxisLabel = "IED count";
d(3).VarName    = "RatePh";
d(3).VarType    = "double";
d(3).BinLenDu   = binlenDu;
d(3).CalcLvl    = "bin";
d(3).CalcFcn    = "gd.dpbGetRatePhCh";
d(3).SrcData    = "Lbl";
d(3).MainLbl    = "IED_Janca";
d(3).ExLblCh    = exLblCh; % Labels to exclude in individual channels
d(3).ExLblAn    = exLblAn; % Labels to exclude in all channels if present in any
d(3).MinSepS    = minSepIedS;
d(3).PlotTitle  = "IED rate";
d(3).YAxisLabel = "IEDs/hour";
dpDesc.(dpDesc.Name(2)) = d;
clear d

%% Clusters description
clDesc(1).EventName = "Seizure"; % This has a different structure than dsDesc and dpDesc
clDesc(1).EventValidSrc = "Seizure21600";
clDesc(1).MinNumInClus = 4; % Minimum required number of given phenomena in the cluster
clDesc(1).InterclusterMultiplier = 2; % The intercluster period must be stg.InterclusterMultiplier times longer than the longest intracluster inter-event interval
clDesc(1).MaxClusterDur = 7; % Maximum cluster duration in days
clDesc(1).MaxWithinClusIeiD = 2; % Maximum inter-event interval within the cluster in days, if longer, it is not a cluster
clDesc(1).ExclClAtEdges = true;

%% Figures description
% General settings
stg.numSubj = numel(subjToPlot);
stg.sbNCol = max(1, ceil(sqrt(stg.numSubj)) - 1); % Subplots of subjects: number of columns
stg.sbNRow = ceil(stg.numSubj/stg.sbNCol); % Subplots subjects - number of rows
stg.figWidth1Cm = 8.5;
stg.figWidth2Cm = stg.figWidth1Cm*2 + 0.5;
stg.margGlob = [2 0.5 0 0.2]; % Left, bottom, right, top
stg.marg = [0.6 0.6 0.5 0.5]; % Left, bottom, right, top
stg.margGlobCi = [0 0 0 0]; % Left, bottom, right, top
stg.margCi = [0.1 0.5 0.1 0.5]; % Left, bottom, right, top
stg.margGlobSlopeBox = [0 0 0 0];
stg.margSlopeBox = [0.15 0.1 0.1 0.3];

% List the figures you wish to plot
figDesc.Name = ["SzRaster"; "SzKaroly"];

% SzRaster
kfig = 1;
d.Name          = figDesc.Name(kfig);
d.FigFcn        = "fig.figRaster"; % Function to use
d.EventName     = "Seizure"; % Phenomenon to stem
d.EventChar     = "DurDu"; % Characteristic to display as the height of the stems
d.EventValidSrc = "Seizure21600";
d.PositionCm    = [5, 5, stg.figWidth2Cm, stg.numSubj*0.7 + 1.5]; % Position in centimeters
figDesc.(figDesc.Name(kfig)) = d;
clear d

% SzKaroly
kfig = kfig + 1;
d.Name          = figDesc.Name(kfig);
d.FigFcn        = "fig.figKaroly"; % Function to use
d.EventName     = "Seizure"; % Phenomenon to stem
d.EventValidSrc = "Seizure21600";
% % % d.PositionCm    = [5, 5, stg.figWidth2Cm, stg.numSubj*5 + 1.5]; % Position in centimeters
d.subplotHeCm   = 5;
figDesc.(figDesc.Name(kfig)) = d;
clear d




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%% END OF NEW SETTINGS %%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%







colorfulSubjects = true;
stg.uniformSubjectColor = [0.8 0.1 0.1];
fcn.getSubjAndSubjClr(subjToPlot, subjList, colorfulSubjects);

% General
stg.dataFolder = 'DataEmgNotExcluded/';
stg.removeEmgContaminatedTF = false;
stg.keepOriginalSubjectName = true;
stg.numEegCh = 4; % Number of EEG channels (other channels may be analysis results)

% Seizures
stg.minIsiS = 10; % Minimal required ISI, otherwise consider them one seizure. Used to be 120, changed to 60 seconds on 2023-01-17
stg.minSzVal = 3;
stg.afterSzMarginS = 10;
% Filter for the EEG signal for power computation
stg.flt.szN = 3; % Order
stg.flt.szF1 = 3; % Low-cut
stg.flt.szF2 = 50; % High-cut
stg.isiHistXScale = "linear";
stg.isiHistYScale = "log";
stg.isiPlotExpFitTF = false;
stg.isiPlotPwlFitTF = false;


stg.leadSzTimeS = 4*3600; % Duration of the period before the seizure that must be seizure free
stg.fitSzDurS = 2*3600; % Duration of the fitted region before or after the seizure in seconds

stg.minIedSepS = 0.1; % Minimum IED separation otherwise delete the later one, in seconds
stg.szClNm = ["SEIZURE", "Seizure", "seizure", "s", "sz"];
stg.artClNm = ["jkArtifact01", "highAmpArtifact01"];
% stg.iedClNm = "IED_Janca30Hz5";
% stg.iedClNm = "IED_Janca";
stg.iedClNm = "fast ripple";
stg.emgClNm = "EMG_det01";
stg.fsLbl = 10; % Hz. When manipulating the labels, they are sometimes converted to (binary or m-ary) signal. Here we set the Fs for this signal.
% % % stg.snlDecontaminationCh = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 -1 -1]; % Which channel of the label pertains to the analysis signal. -1 stands for any.
stg.snlDecontaminationCh = [1 2 3 4];
stg.whWinLen = 1; % Whole signal analysis window length in days
stg.clDurMultiple = 1; % How many cluster durations before and after the cluster we should use for pre- and post-cluster analysis
stg.curFitPopNorm = true;

% Simulated signal characteristics
stg.normalizePopCur = true; % Normalize fitted curves (sums of exponentials) in population analysis
stg.simSiCharFltOrder = 3; % 2 exponentials or 3 exponentials
% stg.simMovAveLen = 3600/dpDesc.BinLenDu*4;
stg.simMovAveLen = 1;

% Statistics
stg.withinSubjectStat = @mean;
stg.acrossSubjectStat = @median;
stg.sdTF = true; % If true, SD is reported. If false, SEM is reported.
stg.Q = 0.1; % False discovery rate for Benjamini-Hochberg correction
stg.multCompCorr = "BenjaminiHochberg"; % Method of multiple comparisons correction
stg.numBootstrapSamples = 1000;
stg.andrzejak = false; % Use circ_t which includes a re-normalization according to Andrzejak et al. (2023), Chaos

%% Plot parameters
stg.plotCharsStacked = false;
% Full configuration
% % % stg.szCharToPlot =      ["szOnsN"; "szDurS"; "szRac"; "szPow"; "postIctPow"];
% % % stg.szCharNameShort =   ["Sz/day"; "Dur (s)"; "Racine"; "Pwr (mV^2)"; "Post (mV^2)"];
% % % stg.szCharYLim =        ["nonneg", "nonneg"; "nonneg", "nonneg"; "0", "6"; "nonneg", "nonneg"; "nonneg", "nonneg"];
% % % stg.siCharToPlot =      ["emg";     "ied";    "vrn";   "ac1";  "crc"];
% % % stg.siCharYLim =        ["nonneg", "nonneg"; "nonneg", "nonneg"; "nonneg", "nonneg"; "nonneg", "nonneg"; "nonneg", "nonneg"];
% % % stg.siCharBinWeights =  ["emgValidS"; "iedValidS"; "iedValidS"; "iedValidS"; "crcValidS"];
% % % stg.siCharNameShort =   ["EMG"; "IED"; "Var"; "ACF"; "CRC"];

% Long-term paper configuration
stg.szCharToPlot =      ["szOnsN"; "szDurS"; "szPow"];
stg.szCharNameShort =   ["Sz/day"; "Dur (s)"; "Pwr (mV^2)"];
stg.szCharYLim =        ["nonneg", "nonneg"; "nonneg", "nonneg"; "nonneg", "nonneg"];
stg.siCharToPlot =      "ied";
stg.siCharYLim =        ["nonneg", "nonneg"];
stg.siCharBinWeights =  "iedValidS";
stg.siCharNameShort =   "IED";


stg.saCharToPlot =      [stg.szCharToPlot; stg.siCharToPlot];
stg.ssCharToPlot =      [stg.siCharToPlot; stg.siCharToPlot];
stg.saCharYLim =        [stg.szCharYLim; stg.siCharYLim];
stg.ssCharYLim =        [stg.siCharYLim; stg.siCharYLim];
stg.saCharNameShort =   [stg.szCharNameShort; stg.siCharNameShort];
stg.ssCharNameShort =   [stg.siCharNameShort; stg.siCharNameShort];
stg.statsPlotType = ["violin"; "swarm"]; % box, violin, swarm
stg.szCharHe = 1.5 + 2/numel(stg.szCharToPlot); % Height of the subplot of a single seizure characteristic
stg.siCharHe = 1.5 + 2/numel(stg.siCharToPlot); % Height of the subplot of a single signal characteristic
stg.saCharHe = 1.5 + 2/numel(stg.saCharToPlot); % Height of the subplot of a single signal characteristic
stg.ssCharHe = 1.5 + 2/numel(stg.saCharToPlot); % Height of the subplot of a single signal characteristic
stg.fitColor = [1 0.7 0; 0.0 1.0 0.0; 0.0 0.0 1.0];
stg.curColor = [0.2 0.8 0.4];
stg.simColor = [0.3 1 0.1];
stg.simPopColor = [0.5 0.6 1];
% % stg.fitColor = [1 0.8 0; 0.8 0 0; 0 0 0.9]; % SalyOlomouc
stg.axFontSize = 7;
% stg.statFontSize = 6.66;
stg.statFontSize = 7;


stg.units = 'centimeters';
stg.box = 'on';
stg.lnWiFit = 0.5;
stg.lnWiFitMean = 1.5;
stg.subjColorSubjMeanMult = 0;

global h %#ok<*GVMIS>
h = struct; h.f = []; h.a = [];

%% Prepare for plotting
setFormat; % Set plot colors etc.
h = fig.prepareFigures(stg, h, figDesc);
%% Get data from each subject, analyze them
dobTable = fcn.getSubjectList('Video-EEG data.xlsx'); % Get list of subjects including their date of birth
if analyzeIndividualSubjects % If you have all the subject data in RAM, you may want to skip loading individual subjects
    for ksubj = 1 : stg.numSubj
        lblp = [path0, '\', path1{ksubj}, '\', subjToPlot{ksubj}, '\', pathLbl3{ksubj}]; % Get label path
        snlp = [path0, '\', path1{ksubj}, '\', subjToPlot{ksubj}, '\', pathEeg3{ksubj}]; % Get signal path
        % % % [subjInfo, ds, dp] = fcn.getData(dsDesc, dpDesc, lblp, snlp, dobTable, ksubj, subjToPlot{ksubj}); % Subject info, seizure properties table, signal characteristics table, signal characteristics y-axis labels
        % % % [clust, szBelongsToClust, clustStats] = fcn.getClusters(subjInfo, ds, dp, clDesc(1), ksubj);
        %% TODO004 DO SUBJECT STATS LATER
        subjStats(ksubj, :) = subjectStats(stg, subjInfo, ds, dp, clustStats); %#ok<SAGROW>

        for kfig = 1 : numel(figDesc.Name)
            d = figDesc.(figDesc.Name(kfig));
            funcHandle = str2func(d.FigFcn); % Get function handle from the name of the function.
            funcHandle(stg, h, d, subjInfo, ds, dp, clust)
        end
        % % % % Seizure occurrence analysis
        % % % [szRate, szRateBinlen] = plotSzRate(subjInfo, szCharTbl, siCharTbl, ksubj); % szRate and binlen are used in the plotSzPsd function
        % % % [szPsdPax(ksubj, :), szPsd(ksubj, :)] = plotSzPsd(subjInfo, szRate, szRateBinlen, ksubj); %#ok<SAGROW>
        % % % [isiH{ksubj, 1}, isiStats(ksubj, :), isiHist(ksubj, :)] = plotSzIsiHist(subjInfo, szCharTbl, ksubj); %#ok<SAGROW>
        % % % 
        % % % % Seizure characteristics analyses
        % % % plotSzChar(ksubj, subjInfo, szCharTbl, siCharTbl)
        % % % [szCharFit.whole(ksubj, :), szCharFit.wholeX(ksubj, :), szCharFit.wholeY(ksubj, :)]...
        % % %     = plotSzCharWhFit(ksubj, subjInfo, szCharTbl, siCharTbl);
        % % % [szCharFit.clDu(ksubj, :), szCharFit.clDuX(ksubj, :), szCharFit.clDuY(ksubj, :)]...
        % % %     = plotSzCharCl(ksubj, subjInfo, szCharTbl, siCharTbl, clust);
        % % % [szCharFit.clDu(ksubj, :), szCharFit.clDuX(ksubj, :), szCharFit.clDuY(ksubj, :)]...
        % % %     = plotSzCharClFit(ksubj, subjInfo, szCharTbl, siCharTbl, clust);
        % % % [szCharFit.circ(ksubj, :), szChar.circEd, szChar.circR(ksubj, :)]...
        % % %     = plotSzCharCiFit(ksubj, subjInfo, szCharTbl, siCharTbl);
        % % % clustSzCharTbl(ksubj, :) = clTerminDur(szCharTbl, clust, szBelongsToClust); %#ok<SAGROW>
        % % % 
        % % % % Seizure and signal characteristics analyses in one figure
        % % % plotSaChar(ksubj, subjInfo, szCharTbl, siCharTbl)
        % % % % % % plotSaCharCWT(ksubj, subjInfo, szCharTbl, siCharTbl)
        % % % [saCharFit.whole(ksubj, :), saCharFit.wholeX(ksubj, :), saCharFit.wholeY(ksubj, :)]...
        % % %     = plotSaCharWhFit(ksubj, subjInfo, szCharTbl, siCharTbl);
        % % % [saCharFit.clBe(ksubj, :), saCharFit.clDu(ksubj, :), saCharFit.clAf(ksubj, :)]...
        % % %     = plotSaCharClFit(ksubj, subjInfo, szCharTbl, siCharTbl, clust);
        % % % plotSaCharCl(ksubj, subjInfo, szCharTbl, siCharTbl, clust);
        % % % [saCharFit.circ(ksubj, :), saChar.circP(ksubj, :), saChar.circR(ksubj, :), saChar.ray(ksubj, :), saChar.omn(ksubj, :)]...
        % % %     = plotSaCharCiFit(ksubj, subjInfo, szCharTbl, siCharTbl);
        % % % 
        % % % % Signal characteristics analyses
        % % % plotSiChar(subjInfo, szCharTbl, siCharTbl, ksubj);
        % % % % % % [siPaxTbl(ksubj, :), siPsdTbl(ksubj, :)]...
        % % %     % % % = plotSiCharPsd(subjInfo, siCharTbl, ksubj); %#ok<SAGROW>
        % % % siCharFit.whole(ksubj, :)...
        % % %     = plotSiCharWhFit(subjInfo, szCharTbl, siCharTbl, ksubj);
        % % % [siChar.clBeX, siChar.clBeY(ksubj, :), siChar.clDuX, siChar.clDuY(ksubj, :), siChar.clAfX, siChar.clAfY(ksubj, :)]...
        % % %     = plotSiCharCl(subjInfo, siCharTbl, clust, ksubj);
        % % % [siCharFit.clBe(ksubj, :), siCharFit.clDu(ksubj, :), siCharFit.clAf(ksubj, :)]...
        % % %     = plotSiCharClFit(subjInfo, szCharTbl, siCharTbl, clust, ksubj);
        % % % [siCharFit.circ(ksubj, :), siChar.circP(ksubj, :), siChar.circR(ksubj, :)]...
        % % %     = plotSiCharCiFit(subjInfo, szCharTbl, siCharTbl, ksubj);
        % % % [siCharSzDiffBe(ksubj, :), siCharSzDiffAf(ksubj, :)]...
        % % %     = plotSiCharSzBeAfVsOther(subjInfo, szCharTbl, siCharTbl, ksubj);
        % % % [siChar.szBeX, siChar.szBeY(ksubj, :), siChar.szAfX, siChar.szAfY(ksubj, :)]...
        % % %     = plotSiCharSz(subjInfo, szCharTbl, siCharTbl, ksubj);
        % % % [siCharFit.szBe(ksubj, :), siCharFit.szAf(ksubj, :)]...
        % % %     = plotSiCharSzFit(subjInfo, szCharTbl, siCharTbl, ksubj); % The fitDurS nad leadTimeS should be > 2*dp.BinLenS
        % % % [siCharCur.baseline(ksubj, :), siCharCur.ma(ksubj, :), siCharCur.exp(ksubj, :), siCharCur.expFltB(ksubj, :), siCharCur.expFltA(ksubj, :), siCharCur.tauH(ksubj, :), siCharCur.pwl(ksubj, :), ...
        % % %     siCharCur.fe(ksubj, :), siCharCur.fp(ksubj, :)]...
        % % %     = plotSiCharSzCur(subjInfo, szCharTbl, siCharTbl, isiHist(ksubj, :), ksubj);
        % % % [risingTbl(ksubj, :), risingClTbl(ksubj, :), risingNonClTbl(ksubj, :)] = siCharRisingAroundSz(szCharTbl, siCharTbl, szBelongsToClust); %#ok<SAGROW>
    end
end
if analyzePopulation
    %% Population figures
    % Descriptive statistics
    descriptiveStats(subjStats)
    
    % Seizure occurrence analysis
    % % % plotSzPsdAllPop(szPsdPax, szPsd)
    plotSzIsiHistAll(isiHist)
    plotSzIsiHistPop(cell2mat(isiH))
    
    % Seizure characteristics analyses
    plotSzCharWhFitAllPop(szCharFit.whole, szCharFit.wholeX, szCharFit.wholeY, szCharTbl)
    plotSzCharClFitAllPop(szCharFit.clDu, szCharTbl)
    plotSzCharCiFitAllPop(szCharFit.circ, szChar.circEd, szChar.circR, szCharTbl)
    clTerminDurPop(clustSzCharTbl);

    % Seizure and signal characteristics analyses
    plotSaCharWhFitAllPop(saCharFit.whole, saCharFit.wholeX, saCharFit.wholeY, szCharTbl, siCharTbl);
    plotSaCharClFitAllPop(saCharFit.clBe, saCharFit.clDu, saCharFit.clAf, szCharTbl, siCharTbl)
    plotSaCharCiFitAllPop(saCharFit.circ, saChar.circP, saChar.circR, szCharTbl, siCharTbl)
    
    % Signal characteristics analyses
    % % % plotSiCharPsdAllPop(siPaxTbl, siPsdTbl, siCharTbl)
    plotSiCharWhFitAllPop(siCharFit.whole, siCharTbl);
    plotSiCharClAllPop(siChar.clBeX, siChar.clBeY, siChar.clDuX, siChar.clDuY, siChar.clAfX, siChar.clAfY, siCharTbl)
    plotSiCharClFitAllPop(siCharFit.clBe, siCharFit.clDu, siCharFit.clAf, siCharTbl)
    plotSiCharCiFitAllPop(siCharFit.circ, siChar.circP, siChar.circR, siCharTbl)
    plotSiCharSzAllPop(siChar.szBeX, siChar.szBeY, siChar.szAfX, siChar.szAfY, siCharTbl)
    plotSiCharSzFitAllPop(siCharFit.szBe, siCharFit.szAf, siCharTbl)
    [siCharCurPop.baseline, siCharCurPop.ma, siCharCurPop.expFltB, siCharCurPop.expFltA, siCharCurPop.fe] = plotSiCharSzCurAllPop(siCharCur, siChar.szAfX, siChar.szAfY, siCharTbl);
    siCharRisingAroundSzPop(risingTbl, risingClTbl, risingNonClTbl);
    if stg.numSubj == 14
        save('siCharCurPop.mat', 'siCharCurPop')
    end
end

load('siCharCurPop.mat', 'siCharCurPop')
for ksubj = 1 : stg.numSubj
    lblp = [path0, '\', path1{ksubj}, '\', subjToPlot{ksubj}, '\', pathLbl3{ksubj}];
    snlp = [path0, '\', path1{ksubj}, '\', subjToPlot{ksubj}, '\', pathEeg3{ksubj}]; % Taking "EMG not removed" because the contaminated data are removed later
    % [subjInfo, szCharTbl, siCharTbl] = getData(lblp, snlp, dobTable, ksubj, subjToPlot{ksubj}); % Subject info, seizure properties table, signal characteristics table, signal characteristics y-axis labels
    % [clust, szBelongsToClust, clustStats] = extractClusters(subjInfo, szCharTbl, siCharTbl, ksubj);
    % % % subjStats(ksubj, :) = subjectStats(subjInfo, szCharTbl, siCharTbl, clustStats); %#ok<SAGROW>

    % Seizure and signal characteristics analyses in one figure
    [siCharTbl] = addSimulatedData(siCharTbl, siCharCur, siCharCurPop, ksubj);
    [siCharSimSim.rmsnrmse(ksubj, :), siCharSimSim.rho(ksubj, :), siCharSimSim.pval(ksubj, :), ...
        siCharSimSim.rmsnrmsePop(ksubj, :), siCharSimSim.rhoPop(ksubj, :), siCharSimSim.pvalPop(ksubj, :)] ...
        = simulatedDataSimilarity(siCharTbl);
    plotSsChar(ksubj, subjInfo, szCharTbl, siCharTbl)
    plotSsExplainConvolution(subjInfo, siCharTbl, siCharCurPop)
    % [sfCharFit.whole(ksubj, :), sfCharFit.wholeX(ksubj, :), sfCharFit.wholeY(ksubj, :)]...
    %     = plotSfCharWhFit(ksubj, subjInfo, szCharTbl, siCharTbl, siCharCur);
    % [saCharFit.clBe(ksubj, :), saCharFit.clDu(ksubj, :), saCharFit.clAf(ksubj, :)]...
    %     = plotSaCharClFit(ksubj, subjInfo, szCharTbl, siCharTbl, clust);
    % [saCharFit.circ(ksubj, :), saChar.circP(ksubj, :), saChar.circR(ksubj, :)]...
    %     = plotSaCharCiFit(ksubj, subjInfo, szCharTbl, siCharTbl);
end
plotSsExplainSumOfExp(siCharCurPop)
plotSimSim(siCharSimSim)

% Print figures
printFigures


%% %%%%%%%%%%%%% %%
%%   FUNCTIONS   %%
%% %%%%%%%%%%%%5 %%


%% %%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%   HERE I FINISHED SO FAR   %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% %%
function [szRate, binlen] = plotSzRate(subjInfo, szCharTbl, siCharTbl, ksubj)
    global stg
    global h
    binlen = 0.25; % In days
    typeOfChart = "plot"; % "patch" or "plot"
    [szRate, ed, xOk, yOk] = calcSzRate(subjInfo, szCharTbl, siCharTbl, binlen);
    [plotTF, plotName] = createFigInd(3);
    if plotTF
        % Calculate axes positions
        % % % % % % % % margGlob = [1.8 0.6 0 0]; % Left, bottom, right, top
        % % % % % % % % marg = [0.7 0.7 0.5 0.5]; % Left, bottom, right, top
        [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlob, stg.marg);
        h.f.(plotName).Units = stg.units;
        h.a.(plotName)(ksubj, 1) = axes('Units', stg.units, 'Position', ...
            [spx(mod(ksubj - 1, numc) + 1), spy(ceil(ksubj/numc)), spWi, spHe], 'NextPlot', 'add');
        
        % Signal OK marker
        % [x, y, xx] = getSnlOkXY(subjInfo, siCharTbl);
        plot(xOk, yOk, 'Marker', 'none', 'LineWidth', 1.5, 'Color', stg.subjColor(ksubj, :));
        hold on
        
        switch typeOfChart
            case "patch"
                x = zeros(1, 3*(numel(ed)-1) + 1);
                x(1 : 3 : end-3) = ed(1 : end-1);
                x(2 : 3 : end-2) = ed(1 : end-1);
                x(3 : 3 : end-1) = ed(2 : end);
                x(end) = ed(end);
                y = zeros(1, 3*(numel(ed)-1) + 1);
                y(2 : 3 : end-2) = szRate(1 : end);
                y(3 : 3 : end-1) = szRate(1 : end);
                xp = x(~isnan(y));
                yp = y(~isnan(y));
                h.pa.(plotName)(ksubj, 1) = patch(xp, yp, stg.subjColor(ksubj, :));
            case "plot"
                x = (ed(1 : end-1) + ed(2 : end))/2; % Timestamp data points at bin centers
                y = szRate;
                % x = ed(2 : end); % Timestamp data points at bin ends (makes sense for causal forecasting)
                % x = x - x(1);
                h.p.(plotName)(ksubj, 1) = plot(x, y, 'Marker', 'none', 'MarkerSize', 4, 'LineWidth', 0.5, 'Color', stg.subjColor(ksubj, :));
        end
        h.a.(plotName)(ksubj, 1).Box = stg.box;
        if ksubj > stg.numSubj - stg.sbNCol
            xlabel('Age (days)')
        end
        if mod(ksubj - 1, stg.sbNCol) == 0
            ylabel('Sz rate (day^{-1})')
        end
        title(subjInfo.subjNm, 'Interpreter', 'none', 'Color', 'k', 'FontWeight', 'bold');
        % title(subjInfo.subjNm, 'Interpreter', 'none', 'Color', stg.subjColor(ksubj, :), 'FontWeight', 'bold');
        h.a.(plotName)(ksubj, 1).FontSize = stg.axFontSize;
        h.a.(plotName)(ksubj, 1).Layer = 'top';
    end
end
function [pax, psd] = plotSzPsd(subjInfo, szRate, binlen, ksubj)
    % psd ....... power spectral density
    % pax ....... period axis, btw it is sorted from the high periods (low freq) to the shorter (high freq)
    % subjInfo .. to get subject name
    % szRate .... seizure rate calculated in plotSzRate
    % binlen .... length of the bins in which seizure rate was calculated, defined in plotSzRate
    % ksubj ..... subject index
    global stg
    global h
    szRate = fillmissing(szRate, 'linear');
    maxPer = 7; % In days
    szRate = szRate - mean(szRate);
    [psd, fax, psdci] = pwelch(szRate, fix(maxPer/binlen), fix(maxPer/2/binlen), fix(maxPer/binlen), 1/binlen); % The Fs is in days^-1
    pax = 1./fax(fax >= 1/maxPer); % Period axis
    psd = psd(fax >= 1/maxPer);
    psdci = psdci(fax >= 1/maxPer);
    [plotTF, plotName] = createFigInd(3);
    if plotTF
        % Calculate axes positions
        % % % % % % % % margGlob = [1.8 0.6 0 0]; % Left, bottom, right, top
        % % % % % % % % marg = [0.7 0.7 0.5 0.5]; % Left, bottom, right, top
        [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlob, stg.marg);
        h.f.(plotName).Units = stg.units;
        h.a.(plotName)(ksubj, 1) = axes('Units', stg.units, 'Position', ...
            [spx(mod(ksubj - 1, numc) + 1), spy(ceil(ksubj/numc)), spWi, spHe], 'NextPlot', 'add');
        x = pax;
        y = psd;
        e = psdci;
        h.pa.(plotName)(ksubj, 1) = patch([x; flipud(x)], 10*log10([y - e; flipud(y) + flipud(e)]), 1 - (1 - stg.subjColor(ksubj, :))/8, 'EdgeColor', 'none');
        % h.pa.(plotName)(ksubj, 1) = patch([x; flipud(x)], [y - e; flipud(y) + flipud(e)], 1 - (1 - stg.subjColor(ksubj, :))/8, 'EdgeColor', 'none');
        hold on
        h.p.(plotName)(ksubj, 1) = plot(x, 10*log10(y), 'Marker', 'x', 'MarkerSize', 3, 'LineWidth', 0.5, 'Color', stg.subjColor(ksubj, :));
        % h.p.(plotName)(ksubj, 1) = plot(x, (y), 'Marker', 'x', 'MarkerSize', 3, 'LineWidth', 0.5, 'Color', stg.subjColor(ksubj, :));
        h.p.(plotName)(ksubj, 2) = plot([1, 1], h.a.(plotName)(ksubj).YLim, ':k', 'LineWidth', 1);
        h.a.(plotName)(ksubj, 1).XGrid = 'on';
        h.a.(plotName)(ksubj, 1).Box = stg.box;
        if ksubj > stg.numSubj - stg.sbNCol
            xlabel('Period (days)')
        end
        if mod(ksubj - 1, stg.sbNCol) == 0
            ylabel(['Sz rate', 10, 'PSD (dB)'])
        end
        title(subjInfo.subjNm, 'Interpreter', 'none', 'Color', 'k', 'FontWeight', 'bold');
        % title(subjInfo.subjNm, 'Interpreter', 'none', 'Color', stg.subjColor(ksubj, :), 'FontWeight', 'bold');
        h.a.(plotName)(ksubj, 1).FontSize = stg.axFontSize;
        h.a.(plotName)(ksubj, 1).Layer = 'top';
    end
end
function plotSzPsdAllPop(szPsdPax, szPsd)
    global stg
    global h
    positionCm = [10, 25, stg.singleColumnWidth, 7];
    [plotTF, plotName] = createFigPos(positionCm);
    if plotTF
        h.a.(plotName) = axes('Units', 'centimeters', 'Position', [2.2, 1.3, positionCm(3) - (2.2+0.5), positionCm(4) - (1.3+0.5)]);
        x = szPsdPax(1, :)';
        yNormalizationMatrix = (szPsd(szPsdPax == 1)*ones(1, size(szPsd, 2)));
        y = 10*log10(szPsd./yNormalizationMatrix)';
        plot(x, y, 'LineWidth', 0.5);
        colororder(stg.subjColor);
        hold on;
        y = 10*log10(mean(szPsd./yNormalizationMatrix, 1))';
        plot(x, y, 'k', 'LineWidth', 2)
        xlabel('Period (days)')
        ylabel('Sz rate PSD (dB)')
        h.a.(plotName).FontSize = stg.axFontSize;
        h.a.(plotName).Layer = 'top';
        h.a.(plotName).Box = stg.box;
        h.a.(plotName).Units = 'normalized';
    end
end
function [isiH, stats, isiHist] = plotSzIsiHist(subjInfo, szCharTbl, ksubj)
    global stg
    global h
    isiHist = table('Size', [1 2], 'VariableNames', {'BinEdges', 'Probability'}, 'VariableTypes', {'double', 'double'});
    % Get histogram counts for real data
    isiH = diff(szCharTbl.szOnsN)*24; % Inter-seizure interval in hours
    switch stg.isiHistXScale
        case "log"
            ed = logspace(-2, 3, 17); % Bin edges
        case "linear"
            ed = logspace(-2, log10(48), 17); % Bin edges
    end
    isiHist.BinEdges = ed;
    hc = histcounts(isiH, ed, 'Normalization', 'pdf');
    isiHist.Probability = hc;
    
    % Fit exponential distribution
    xe = (ed(1 : end-1) + ed(2 : end))./2;
    pd = fitdist(isiH, 'Exponential'); % Probability distribution
    ye = pdf(pd, xe);
    ye = ye/sum(ye);
    
    stats = table;
    [stats.chiIsiExpH, stats.chiIsiExpP] = chi2gof(isiH, 'CDF', pd, 'EMin', 1, 'NBins', 20);
    [stats.corIsiR, stats.corIsiP] = corr(isiH(1 : end-1), isiH(2 : end), 'type', 'Spearman');

    [plotTF, plotName] = createFigInd(3);
    if plotTF
        % Calculate axes positions
        % % % % % % % % margGlob = [1.8 0.6 0 0]; % Left, bottom, right, top
        % % % % % % % % marg = [0.7 0.7 0.5 0.5]; % Left, bottom, right, top
        [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlob, stg.marg);
        h.f.(plotName).Units = stg.units;
        h.a.(plotName)(ksubj, 1) = axes('Units', stg.units, 'Position', ...
            [spx(mod(ksubj - 1, numc) + 1), spy(ceil(ksubj/numc)), spWi, spHe], 'NextPlot', 'add');
        
        % Plot the real data
        switch stg.isiHistXScale
            case "log"
                xd = log10(ed); % Data
                xe = log10(xe); % Exponential fit
            case "linear"
                xd = ed;
        end
        y = hc;
        xd = repelem(xd, 3);
        xd = xd(1 : end-1);
        yd(1 : 3 : length(y)*3 - 2) = 0;
        yd(2 : 3 : length(y)*3 - 1) = y;
        yd(3 : 3 : length(y)*3) = y;
        yd = [0, yd, 0];
        switch stg.isiHistYScale
            case "log"
                flor = 1e-4;
                yd(yd <= flor) = flor;
                yd = log10(yd);
                ye = log10(ye);
            case "linear"
                flor = 0;
        end
        patch(xd, yd, stg.subjColor(ksubj, :));
        hold on
    
        % Plot the exponential fit
        if stg.isiPlotExpFitTF
            plot(xe, ye, 'Color', 0*[0.8 0.1 1], 'Marker', 'none', 'LineStyle', '--', 'LineWidth', 1)
        end
        ylim([log10(flor), 1])
        hax = gca;
        switch stg.isiHistXScale
            case "log"
                hax.XTickLabel = cellfun(@(x) ['10^{', x, '}'], hax.XTickLabel, 'UniformOutput', false);
            case "linear"
                hax.XLim = [0 48];
        end
        switch stg.isiHistYScale
            case "log"
                hax.YTickLabel = cellfun(@(x) ['10^{', x, '}'], hax.YTickLabel, 'UniformOutput', false);
        end
        hax.Box = stg.box;
        if ksubj > stg.numSubj - stg.sbNCol
            xlabel('ISI (hours)')
        end
        if mod(ksubj - 1, stg.sbNCol) == 0
            ylabel('Probability')
        end
        title([subjInfo.subjNm], 'Interpreter', 'none', 'Color', 'k')
        % title([subjInfo.subjNm], 'Interpreter', 'none', 'Color', stg.subjColor(ksubj, :))
        % title(['Mouse ', num2str(ksubj)], 'Interpreter', 'none', 'Color', stg.subjColor(ksubj, :))
        h.a.(plotName)(ksubj, 1).FontSize = stg.axFontSize;
        h.a.(plotName)(ksubj, 1).Layer = 'top';
    end
end
function plotSzIsiHistAll(isiHist)
    global stg
    global h
    % Plotting
    positionCm = [10, 25, 8.5, 7];
    [plotTF, plotName] = createFigPos(positionCm);
    if plotTF
        h.a.(plotName) = axes('Units', 'centimeters', 'Position', [2.2, 1.3, positionCm(3) - (2.2+0.5), positionCm(4) - (1.3+0.5)]);
        % Plot the real data
        ed1 = isiHist.BinEdges(1, :);
        for ksubj = 1 : height(isiHist)
            ed = isiHist.BinEdges(ksubj, :);
            if ~all(ed == ed1, 'all')
                error('_jk Different bin edges.')
            end
            prob(ksubj, :) = isiHist.Probability(ksubj, :); %#ok<AGROW>
        end
        switch stg.isiHistXScale
            case "log"
                xd = log10(ed1); % Data
            case "linear"
                xd = ed1;
        end
        xd = repelem(xd, 3);
        xd = xd(1 : end-1);
        y = mean(prob, 1);
        yd(1 : 3 : length(y)*3 - 2) = 0;
        yd(2 : 3 : length(y)*3 - 1) = y;
        yd(3 : 3 : length(y)*3) = y;
        yd = [0, yd, 0];
        e = std(prob, 1)/sqrt(size(prob, 1));
        switch stg.isiHistYScale
            case "log"
                flor = 1e-4;
                yd(yd <= flor) = flor;
                yd = log10(yd);
                y = log10(y);
                e = log10(y);
            case "linear"
                flor = 0;
        end
        patch(xd, yd, [0.8 0.8 0.8]);
        hold on
        errorbar((ed1(1:end-1) + ed1(2:end))/2, y, e)
        hold on
       
        % Format the graph
        hax = h.a.(plotName);
        switch stg.isiHistXScale
            case "log"
                hax.XTickLabel = cellfun(@(x) ['10^{', x, '}'], hax.XTickLabel, 'UniformOutput', false);
            case "linear"
                hax.XLim = [0 48];
        end
        switch stg.isiHistYScale
            case "log"
                hax.YLim = ([log10(flor), 1]);
                hax.YTickLabel = cellfun(@(x) ['10^{', x, '}'], hax.YTickLabel, 'UniformOutput', false);
        end
        xlabel('ISI (hours)')
        ylabel('Probability')
        hax.FontSize = stg.axFontSize;
        hax.Layer = 'top';
        hax.Box = stg.box;
        hax.Units = 'normalized';
    end
end
function plotSzIsiHistPop(isiH)
    % Source of useful functions:
    % https://github.com/wyz996/power-law
    % NOTE! xmin influences the alpha, see below.
    % Cite them in the paper!
    global stg
    global h
    % % For degugging
    % close all; h.f = [];
    % 
    % % Artificial data generation (for debugging)
    % % Poisson process (has exponential distribution of ISIs)
    % N = 1000; % Number of artificial ISIs
    % muu = 9.5; % Mean ISI (for the Poisson process)
    % poissonEvents = sort(rand(N, 1)*N*muu); % Poisson process
    % isiHe = diff(poissonEvents); % Inter-seizure intervals in hours (exponentially distributed)
    % % Power-law distributed ISI
    % isiHp = ((1-rand(N,1)).^(-1/(2.5-1)))*0.1 + 0.0*rand(N,1); % Power-law distributed data according to the link above
    % isiH = isiHp

    
    % Data analysis
    switch stg.isiHistXScale
        case "log"
            ed = logspace(-2, 3, 17); % Bin edges
        case "linear"
            ed = logspace(-2, log10(48), 17); % Bin edges
    end
    
    % % % ed = logspace(-2, 3, 9); % Bin edges
    hc = histcounts(isiH, ed, 'Normalization', 'pdf');
    x = (ed(1 : end-1) + ed(2 : end))/2;

    % Fit exponential distribution
    pd = fitdist(isiH, 'Exponential'); % Probability distribution
    yfe = pdf(pd, x);
    [exponChiH, exponChiP] = chi2gof(isiH, 'CDF', pd, 'EMin', 1, 'NBins', 20); %#ok<ASGLU>
    [exponKsH, exponKsP] = kstest(isiH, 'CDF', pd); %#ok<ASGLU>
    % [corrR, corrP] = corr(isiH(1 : end-1), isiH(2 : end), 'type', 'Spearman'); % Correlation of subsequent ISI. This implementation is wrong since multiple subjects are within the isiH vector.
    
    % Fit power-law distribution
    [powerLawAlpha, powerLawXmin] = plfit(isiH, 'xmin', 0.5, 'range', (1.1 : 0.1 : 2.5)'); % The xmin value influences the alpha a lot
    % [powerLawAlpha, powerLawXmin] = plfit(isiH, 'range', (1.1 : 0.1 : 2.5)')
    [powerKsP, powerGof] = plpva(isiH, powerLawXmin, 'range', (1.1 : 0.1 : 2.5)', 'reps', 100, 'silent'); %#ok<ASGLU>
    yfp = (powerLawAlpha-1)/powerLawXmin * (x/powerLawXmin).^-powerLawAlpha; % For plotting and chi2gof test
    % [pwrlwChiH, pwrlwChiP] = chi2gof(xf, 'Ctrs', xf, 'Frequency', hc, 'Expected', yfp, 'NParams', 2) % I do not know how to do chi2 test. This is probably wrong.

    % Plotting
    positionCm = [10, 25, 8.5, 7];
    [plotTF, plotName] = createFigPos(positionCm);
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp('ISI histogram:')
        disp(['(', plotName, ')'])
        disp(['    Exponential', 10, ...
            '        Chi-squared   p = ', num2str(exponChiP, '%.2g'), 10, ...
            '        KS test       p = ', num2str(exponKsP, '%.2g'), 10, ...
            '    Power-law', 10, ...
            '        Parameter alpha = ', num2str(powerLawAlpha, '%.2g'), 10, ...
            '        Parameter  xmin = ', num2str(powerLawXmin, '%.2g'), 10, ...
            '        KS test       p = ', num2str(powerKsP, '%.2g')]);
    end
    if plotTF
        h.a.(plotName) = axes('Units', 'centimeters', 'Position', [2.2, 1.3, positionCm(3) - (2.2+0.5), positionCm(4) - (1.3+0.5)]);
        % Plot the real data
        switch stg.isiHistXScale
            case "log"
                xd = log10(ed); % Data
                x = log10(x); % Exponential fit
            case "linear"
                xd = ed;
        end

        y = hc;
        xd = repelem(xd, 3);
        xd = xd(1 : end-1);
        yd(1 : 3 : length(y)*3 - 2) = 0;
        yd(2 : 3 : length(y)*3 - 1) = y;
        yd(3 : 3 : length(y)*3) = y;
        yd = [0, yd, 0];
        switch stg.isiHistYScale
            case "log"
                flor = 1e-4;
                yd(yd <= flor) = flor;
                yd = log10(yd);
                ye = log10(yfe);
                yp = log10(yfp);
            case "linear"
                flor = 0;
                ye = yfe;
                yp = yfp;
        end

        patch(xd, yd, [0.8 0.8 0.8]);
        hold on
        % Plot the exponential fit
        % % % % % % ye = log10(yfe);
        if stg.isiPlotExpFitTF
            plot(x, ye, 'Color', 0*[0.8 0.1 1], 'Marker', 'none', 'LineStyle', '--', 'LineWidth', 1)
        end
        % Plot the power-law fit
        if stg.isiPlotPwlFitTF
            plot(x, yp, 'Color', 0*[0.8 0.1 1], 'Marker', 'none', 'LineStyle', ':', 'LineWidth', 1)
        end
        % Format the graph
        hax = h.a.(plotName);
        switch stg.isiHistXScale
            case "log"
                hax.XTickLabel = cellfun(@(x) ['10^{', x, '}'], hax.XTickLabel, 'UniformOutput', false);
            case "linear"
                hax.XLim = [0 48];
        end
        switch stg.isiHistYScale
            case "log"
                hax.YLim = ([log10(flor), 1]);
                hax.YTickLabel = cellfun(@(x) ['10^{', x, '}'], hax.YTickLabel, 'UniformOutput', false);
        end
        % % % % % % hax.YTickLabel = cellfun(@(x) ['10^{', x, '}'], hax.YTickLabel, 'UniformOutput', false);
        % % % % % % hax.XTickLabel = cellfun(@(x) ['10^{', x, '}'], hax.XTickLabel, 'UniformOutput', false);
        xlabel('ISI (hours)')
        ylabel('Probability')
        hax.FontSize = stg.axFontSize;
        hax.Layer = 'top';
        hax.Box = stg.box;
        hax.Units = 'normalized';
        % Text
        if stg.isiPlotPwlFitTF
            text(hax.XLim(1) + 0.75*(diff(hax.XLim)), hax.YLim(1) + 0.75*(diff(hax.YLim)), ['$\alpha=', num2str(powerLawAlpha, 3), '$'], 'Interpreter', 'latex')
        end
    end
end
% Plot seizure characteristics
function plotSzChar(ksubj, subjInfo, szCharTbl, siCharTbl)
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSzCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl); % Refer to the created axes using h.f.(thisFcnName)
        formatAxesInd(plotName, ksubj, subjInfo, szCharTbl.Properties.VariableNames, 'Age (days)', "analysisPeriod")
    end
end
function [fitTbl, xxFitTbl, yyFitTbl] = plotSzCharWhFit(ksubj, subjInfo, szCharTbl, siCharTbl)
    global stg
    global h
    numChar = numel(stg.szCharToPlot);
    [fitTbl, ed, ~, yyTbl, xxFitTbl, yyFitTbl] = fitSzCharWh(subjInfo, szCharTbl, siCharTbl); % Compute the trend
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSzCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl); % Refer to the created axes using h.f.(thisFcnName)
        
        % Plot the trend
        for kchar = 1 : numChar
            y = yyTbl{1, kchar};
            % Plotting patch (bar graph)
            hax = h.a.(plotName)(ksubj, kchar);
            pax = [ed(1), ed(1), repelem(ed(2:end-1), 3), ed(end), ed(end), ed(1)];
            pay = repelem(y, 3);
            pay(1 : 3 : end) = 0;
            pay = [pay, 0, 0]; %#ok<AGROW>
            pay(isnan(pay)) = 0;
            patch(hax, pax, pay, 1 - 0.2*(1 - stg.subjColor(ksubj, :)))
            if kchar == 1
                h.p.(plotName)(ksubj, kchar, 2).YData(2 : 3 : end - 1) = ceilToEven1sig(max(pay)); % Change the height of seizures in the plot of sz rate
            end
            hold on
            % Plotting fitted line
            x = xxFitTbl{1, kchar};
            y = yyFitTbl{1, kchar};
            plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', 1, 'Color', stg.subjColor(ksubj, :)*stg.subjColorSubjMeanMult);
            hax.Children([1 2 3 4]) = hax.Children([1 3 4 2]);
            slopeText = num2str(fitTbl{1, kchar + 0*numChar}, '%.2g');
            if all(sign(fitTbl{1, kchar + [1, 2]*numChar})) == sign(fitTbl{1, kchar + 0*numChar})
                fontWeight = 'bold';
            else
                fontWeight = 'normal';
            end
            textbp(slopeText, 'Parent', hax, 'FontSize', stg.axFontSize, 'FontWeight', fontWeight)
        end
        formatAxesInd(plotName, ksubj, subjInfo, szCharTbl.Properties.VariableNames, 'Age (days)', "analysisPeriod")
    end
end
function plotSzCharWhFitAllPop(fitTbl, xxFitTbl, yyFitTbl, szCharTbl) %#ok<INUSD>
    global stg
    global h
    normalizedTF = 1;
    numChar = numel(stg.szCharToPlot);
    str = strings(1, numChar); medFit = NaN(1, numChar); sigRanP = NaN(1, numChar);
    for kchar = 1 : numChar
        [str(kchar), ~, ~, medFit(kchar), ~, sigRanP(kchar)] = getBasicStats(fitTbl.(stg.szCharToPlot(kchar) + "fitSlope"), 'Slope');
    end
    [sigRanP, sigRanH] = getSignificanceFromP(sigRanP, stg.Q, stg.multCompCorr);
    positionCm = [10, 25, stg.singleColumnWidth, min(25, stg.szCharHe*numChar)];
    [plotTF, plotName] = createFigPos(positionCm);
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp('Whole recording:')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, szCharTbl.Properties.VariableNames, kchar);
            disp(['    ', ylbl])
            disp(str(kchar))
        end
    end
    if plotTF
        if ~isfield(h.a, plotName)
            createAxesAllPopStack(plotName)
        end
        strP = strings(1, numChar);
        for kchar = 1 : numChar
            for ksubj = 1 : stg.numSubj
                hax = h.a.(plotName)(1, kchar);
                if normalizedTF
                    x = [0 1]; %#ok<*UNRCH>
                    fitSlope = fitTbl{ksubj, kchar};
                    y = x*fitSlope;
                else
                    x = xxFitTbl{ksubj, kchar};
                    y = yyFitTbl{ksubj, kchar};
                end
                plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', stg.lnWiFit, 'Color', stg.subjColor(ksubj, :)); % Not normalized
            end
            if normalizedTF
                y = x*medFit(1, kchar);
                plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            else
                x = stg.acrossSubjectStat(xxFitTbl{:, kchar}, 1, 'omitmissing');
                y = stg.acrossSubjectStat(yyFitTbl{:, kchar}, 1, 'omitmissing');
                plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            end
            switch sign(medFit(kchar))
                case -1
                    signMed = '-';
                case 0
                    signMed = '0';
                case 1
                    signMed = '+';
            end
            strP(kchar) = string(['  p=', num2str(sigRanP(1, kchar), 2), '(', signMed, ')']);
        end
        formatAxesAllPop(plotName, szCharTbl.Properties.VariableNames, 'Normalized time', "none", ["any", "any"])
        for kchar = 1 : numChar
            hax = h.a.(plotName)(kchar);
            xt = hax.XLim(1);
            % yt = mean(hax.YLim);
            yt = hax.YLim(2);
            ht = text(hax, xt, yt, strP(kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','left', 'VerticalAlignment','top', 'BackgroundColor','none');
            if sigRanH(1, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
        end
    end
end
function [mFitTbl, mxxFitTbl, myyFitTbl] = plotSzCharCl(ksubj, subjInfo, szCharTbl, siCharTbl, cl) %#ok<DEFNU,INUSD>
% NOT FINISHED
    % ksubj ...... subject index
    % subjInfo ... table with subject info
    % szCharTbl .. table with seizure characteristics
    % siCharTbl .. table with signal characteristics
    % cl ......... structure with clusters of given mouse
    % mFitTbl .... table with mean slopes, offsets and confidence intervlas (the CI mean has no real meaning though)
    % mxFit0 ..... mean cluster duration (probably not further used)
    % myFit0 ..... mean difference between value of the fit at the offset - onset of the cluster for each characteristic (probably not further used)
    global stg
    global h %#ok<NUSED>
    numChar = numel(stg.szCharToPlot); %#ok<NASGU>
    [fitTbl, eded, ~, yyTbl, xxFitTbl, yyFitTbl] = fitSzCharCl(subjInfo, siCharTbl, cl);
    mFitTbl = stg.withinSubjectStat(fitTbl, "omitmissing");
    vrnm = xxFitTbl.Properties.VariableNames;
    if isempty(xxFitTbl)
        mxxFitTbl = varfun(@(x) [x; NaN, NaN], xxFitTbl, 'OutputFormat', 'table'); % Essentially mean cluster duration
        myyFitTbl = varfun(@(x) [x; NaN, NaN], yyFitTbl, 'OutputFormat', 'table'); % Essentially end - start of the fitted line
        myyTbl = varfun(@(x) [x; NaN(1, size(yyTbl{1, 1}, 2))], yyTbl, 'OutputFormat', 'table'); % Average of the data
    else
        % mxxFitTbl = stg.withinSubjectStat(varfun(@(x) diff(x, [], 2), xxFitTbl, 'OutputFormat', 'table'), 1); % Returns scalar x
        % myyFitTbl = stg.withinSubjectStat(varfun(@(x) diff(x, [], 2), yyFitTbl, 'OutputFormat', 'table'), 1); % Returns scalar y
        mxxFitTbl = stg.withinSubjectStat(varfun(@(x) x - [x(:, 1), x(:, 1)], xxFitTbl, 'OutputFormat', 'table'), 1); % Returns vector [0, x]
        myyFitTbl = stg.withinSubjectStat(varfun(@(x) x - [x(:, 1), x(:, 1)], yyFitTbl, 'OutputFormat', 'table'), 1); % Returns vector [0, y]
        myyTbl = stg.withinSubjectStat(varfun(@(x) x, yyTbl, 'OutputFormat', 'table'), 1); % Average of the data
    end
    mxxFitTbl.Properties.VariableNames = vrnm;
    myyFitTbl.Properties.VariableNames = vrnm;
    myyTbl.Properties.VariableNames = vrnm;
    [plotTF, plotName] = createFigInd(1);
    if plotTF

                x = 0 : size(eded, 2) - 2;
                x = 0.33 + 0.33*(x/max(x));
plotColor = ones(size(yyTbl, 1), 1)*stg.subjColor(ksubj, :);
numRealOut = plotPeriEventData(plotName, ksubj, x, yyTbl, myyTbl, plotColor); %#ok<NASGU>

        % % formatAxesInd(plotName, ksubj, subjInfo, szCharTbl.Properties.VariableNames, 'Age (days)', "analysisPeriod")
    end
end
function [mFitTbl, mxxFitTbl, myyFitTbl] = plotSzCharClFit(ksubj, subjInfo, szCharTbl, siCharTbl, cl)
    % ksubj ...... subject index
    % subjInfo ... table with subject info
    % szCharTbl .. table with seizure characteristics
    % siCharTbl .. table with signal characteristics
    % cl ......... structure with clusters of given mouse
    % mFitTbl .... table with mean slopes, offsets and confidence intervlas (the CI mean has no real meaning though)
    % mxFit0 ..... mean cluster duration (probably not further used)
    % myFit0 ..... mean difference between value of the fit at the offset - onset of the cluster for each characteristic (probably not further used)
    global stg
    global h
    numChar = numel(stg.szCharToPlot);
    [fitTbl, eded, ~, yyTbl, xxFitTbl, yyFitTbl] = fitSzCharCl(subjInfo, siCharTbl, cl);
    mFitTbl = stg.withinSubjectStat(fitTbl, "omitmissing");
    vrnm = xxFitTbl.Properties.VariableNames;
    if isempty(xxFitTbl)
        mxxFitTbl = varfun(@(x) [x; NaN, NaN], xxFitTbl, 'OutputFormat', 'table'); % Essentially mean cluster duration
        myyFitTbl = varfun(@(x) [x; NaN, NaN], yyFitTbl, 'OutputFormat', 'table'); % Essentially end - start of the fitted line
    else
        % mxxFitTbl = stg.withinSubjectStat(varfun(@(x) diff(x, [], 2), xxFitTbl, 'OutputFormat', 'table'), 1); % Returns scalar x
        % myyFitTbl = stg.withinSubjectStat(varfun(@(x) diff(x, [], 2), yyFitTbl, 'OutputFormat', 'table'), 1); % Returns scalar y
        mxxFitTbl = stg.withinSubjectStat(varfun(@(x) x - [x(:, 1), x(:, 1)], xxFitTbl, 'OutputFormat', 'table'), 1); % Returns vector [0, x]
        myyFitTbl = stg.withinSubjectStat(varfun(@(x) x - [x(:, 1), x(:, 1)], yyFitTbl, 'OutputFormat', 'table'), 1); % Returns vector [0, y]
    end
    mxxFitTbl.Properties.VariableNames = vrnm;
    myyFitTbl.Properties.VariableNames = vrnm;
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSzCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl); % Refer to the created axes using h.f.(thisFcnName)
        % Plot the fit
        numCl = numel(cl);
        for kchar = 1 : numChar
            maxHeight = 0;
            for kcl = 1 : numCl
                % Bar graph (it is actually a patch but looks like a bar graph)
                y = yyTbl{kcl, kchar};
                ed = eded(kcl, :);
                % Plotting
                hax = h.a.(plotName)(ksubj, kchar);
                pax = [ed(1), ed(1), repelem(ed(2:end-1), 3), ed(end), ed(end), ed(1)];
                pay = repelem(y, 3);
                pay(1 : 3 : end) = 0;
                pay = [pay, 0, 0]; %#ok<AGROW>
                pay(isnan(pay)) = 0;
                maxHeight = max([maxHeight, pay]);
                h.pa.(plotName)(ksubj, kchar, 1) = patch(hax, pax, pay, 1 - 0.2*(1 - stg.subjColor(ksubj, :)), 'Tag', 'szCharClDuBar');
                hold on
                % Fit plot
                x = xxFitTbl{kcl, kchar};
                y = yyFitTbl{kcl, kchar};
                maxHeight = max([maxHeight, y]);
                h.p.(plotName)(ksubj, kchar, 3) = plot(hax, x, y, ...
                    'Marker', 'none', 'LineStyle', '-', 'LineWidth', stg.lnWiFitMean, 'Color', stg.subjColor(ksubj, :)*stg.subjColorSubjMeanMult, 'Tag', 'szCharClDuFit');
                uistack(h.pa.(plotName)(ksubj, kchar), 'bottom')
                uistack(h.p.(plotName)(ksubj, kchar, 3), 'top')
            end
            if kchar == 1
                h.p.(plotName)(ksubj, kchar, 2).YData(2 : 3 : end-1) = ceilToEven1sig(maxHeight); % Change the height of seizures in the plot of sz rate
            end
        end
        formatAxesInd(plotName, ksubj, subjInfo, szCharTbl.Properties.VariableNames, 'Age (days)', "analysisPeriod")
    end
end
function plotSzCharClFitAllPop(fitTbl, szCharTbl)
    global stg
    global h
    numChar = numel(stg.szCharToPlot);
    str = strings(1, numChar); medFit = NaN(1, numChar); sigRanP = NaN(1, numChar);
    for kchar = 1 : numChar
        [str(kchar), ~, ~, medFit(kchar), ~, sigRanP(kchar)] = getBasicStats(fitTbl.(stg.szCharToPlot(kchar) + "fitSlope"), 'Slope');
    end
    [sigRanP, sigRanH] = getSignificanceFromP(sigRanP, stg.Q, stg.multCompCorr);
    positionCm = [10, 25, stg.singleColumnWidth, min(25, stg.szCharHe*numChar)];
    [plotTF, plotName] = createFigPos(positionCm);
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp('Cluster means (animal-wise):')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, szCharTbl.Properties.VariableNames, kchar);
            disp(['    ', ylbl])
            disp(str(kchar))
        end
    end
    if plotTF
        if ~isfield(h.a, plotName)
            createAxesAllPopStack(plotName)
        end
        x = [0 1];
        strP = strings(1, numChar);
        for kchar = 1 : numChar
            hax = h.a.(plotName)(kchar);
            for ksubj = 1 : stg.numSubj
                fitSlope = fitTbl{ksubj, kchar};
                y = x*fitSlope;
                plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', stg.lnWiFit, 'Color', stg.subjColor(ksubj, :));
            end
            y = x*medFit(1, kchar);
            plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            switch sign(medFit(kchar))
                case -1
                    signMed = '-';
                case 0
                    signMed = '0';
                case 1
                    signMed = '+';
            end
            strP(kchar) = string(['  p=', num2str(sigRanP(1, kchar), 2), '(', signMed, ')']);
        end
        formatAxesAllPop(plotName, szCharTbl.Properties.VariableNames, 'Normalized cluster time', "none", ["any", "any"])
        for kchar = 1 : numChar
            hax = h.a.(plotName)(kchar);
            % xt = mean(hax.XLim);
            xt = 0;
            % yt = mean(hax.YLim);
            yt = hax.YLim(2);
            ht = text(hax, xt, yt, strP(kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','left', 'VerticalAlignment','top', 'BackgroundColor','none');
            if sigRanH(1, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
        end
    end
end
function [cirTbl, ed, rrTbl] = plotSzCharCiFit(ksubj, subjInfo, szCharTbl, siCharTbl)
    global stg
    global h
    numBetw = 10; % Only for smoothness of the circular sections
    numChar = numel(stg.szCharToPlot);
    [cirTbl, ed, ~, rrTbl, ppCirTbl, rrCirTbl] = fitSzCharCi(szCharTbl, siCharTbl); % Compute the circular statistics
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSzCharDataCirc(plotName, ksubj, szCharTbl); % Refer to the created axes using h.f.(thisFcnName)
        % Plot the trend
        for kchar = 1 : numChar
            % Plotting patch
            r = rrTbl{1, kchar};
            if kchar == 1 % Here, the size of the individual seizure vectors has no meaning, so we normalize histogram to max of the histogram
                r = r/max(r);
            else
                y1 = szCharTbl{:, stg.szCharToPlot(kchar)};
                r = r/max(y1);
            end
            hax = h.a.(plotName)(ksubj, kchar);
            circleCenter = max(r)/1000; % Pure zero minimum radius makes trouble probably due to rounding errors
            pap = NaN(1, (numel(ed)-1)*(numBetw+3) + 1); par = pap;
            pap(1) = ed(1);
            par(1) = 0;
            for ke = 1 : numel(ed) - 1
                pap((ke-1)*(numBetw+3) + (2 : 2+numBetw+1)) = linspace(ed(ke), ed(ke+1), numBetw+2);
                pap((ke-1)*(numBetw+3) + 2+numBetw+2) = ed(ke+1);
                par((ke-1)*(numBetw+3) + (2 : 2+numBetw+1)) = r(ke);
                par((ke-1)*(numBetw+3) + 2+numBetw+2) = circleCenter;
            end
            pap = pap*2*pi;
            par(isnan(par)) = circleCenter;
            par(end) = circleCenter;
            [x, y] = pol2cart(-pap - pi/2, par);
            z = -5*ones(size(x));
            h.pa.(plotName)(ksubj, kchar, 2) = patch(hax, x, y, z, 1 - 0.2*(1 - stg.subjColor(ksubj, :)), 'EdgeColor', stg.subjColor(ksubj, :)*stg.subjColorSubjMeanMult);
            % Plotting resultant vector
            resVecWi = 1.5;
            prv = ppCirTbl{1, kchar};
            rrv = rrCirTbl{1, kchar}*hax.XLim(2);
            [x, y] = pol2cart(-prv - pi/2, rrv);
            h.p.(plotName)(ksubj, kchar, 3) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', stg.subjColor(ksubj, :)*0);
            % Arrowhead
            [x(1), y(1)] = pol2cart(-(prv(2)-5/6*pi) - pi/2, rrv(2)/4);
            x(1) = sum(x);
            y(1) = sum(y);
            h.p.(plotName)(ksubj, kchar, 4) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', stg.subjColor(ksubj, :)*0);
            x(1) = 0; y(1) = 0;
            [x(1), y(1)] = pol2cart(-(prv(2)+5/6*pi) - pi/2, rrv(2)/4);
            x(1) = sum(x);
            y(1) = sum(y);
            h.p.(plotName)(ksubj, kchar, 5) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', stg.subjColor(ksubj, :)*0);
            % % % % R scaling
            % % % if kchar == 1
            % % %     factor = ceilToEven1sig(max([r, rrCirTbl{1, kchar}]));
            % % %     for kplot = 1 : 2
            % % %         h.p.(plotName)(ksubj, kchar, kplot).XData = factor*h.p.(plotName)(ksubj, kchar, kplot).XData;
            % % %         h.p.(plotName)(ksubj, kchar, kplot).YData = factor*h.p.(plotName)(ksubj, kchar, kplot).YData;
            % % %     end
            % % %     h.pa.(plotName)(ksubj, kchar, 1).XData = factor*h.pa.(plotName)(ksubj, kchar, 1).XData;
            % % %     h.pa.(plotName)(ksubj, kchar, 1).YData = factor*h.pa.(plotName)(ksubj, kchar, 1).YData;
            % % %     h.a.(plotName)(ksubj, kchar).XLim = factor*h.a.(plotName)(ksubj, kchar).XLim;
            % % %     h.a.(plotName)(ksubj, kchar).YLim = factor*h.a.(plotName)(ksubj, kchar).YLim;
            % % % end
            % Text
            xt = 0;
            yt = hax.YLim(1);
            text(hax, xt, yt, ['r=', num2str(rrv(2), 1)], 'Interpreter', 'tex', 'FontSize', stg.axFontSize, 'FontWeight', 'normal',...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
            % ylbl = getYlbl(plotName, szCharTbl.Properties.VariableNames, kchar);
            ylbl = stg.szCharNameShort(kchar);
            yt = hax.YLim(2);
            text(hax, xt, yt, ylbl, 'Interpreter', 'tex', 'FontSize', stg.axFontSize, 'FontWeight', 'normal',...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        end
        % % % % % % % % margGlob = [0 0 0 0]; marg = [0.1 0.5 0.1 0.5];
        [~, ~, spWi, ~, ~, ~] = getSubplotXYWH(plotName, stg.margGlobCi, stg.margCi);
        hax = h.a.(plotName)(ksubj, 1);
        x = hax.XLim(1) + diff(hax.XLim)*(0.9*spWi/hax.Position(3))/2;
        y = hax.YLim(2) + 0.3*diff(hax.YLim);
        text(hax, x, y, subjInfo.subjNm, 'Interpreter', 'none', 'FontSize', stg.axFontSize,...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        for kchar = 1 : numChar
            set(h.a.(plotName)(ksubj, kchar), 'Units', 'normalized');
        end
    end
end
function plotSzCharCiFitAllPop(cirTbl, ed, rrTbl, szCharTbl)
    global stg
    global h
    numBetw = 10; % Only for smoothness of the circular sections
    numChar = numel(stg.szCharToPlot);
    str = strings(1, numChar); resM = NaN(1, numChar); resR = NaN(1, numChar); rayP = NaN(1, numChar); omnP = NaN(1, numChar);
    for kchar = 1 : numChar
        p = cirTbl{:, kchar};
        p = p(~isnan(p));
        resM(kchar) = circ_mean(p);
        if ~stg.andrzejak
            resR(kchar) = circ_r(p);
        else
            resR(kchar) = circ_t(p);
        end
        rayP(kchar) = circ_rtest(p);
        omnP(kchar) = circ_otest(p);
        clear p
        xNm = 'Time of day';
        str(kchar) = ...
            [repelem(' ', 1, 13-numel(xNm)), xNm, ' = ', num2str(mod(resM(kchar)/2/pi*24, 24), '%.2g'), ', L=', num2str(resR(kchar), '%.2g'), 10, ...
             '        Ray p = ', num2str(rayP(kchar), '%.2g'), 10, ...
             '       Omni p = ', num2str(omnP(kchar), '%.2g'), 10];
    end
    [~, rayH] = getSignificanceFromP(rayP, stg.Q, stg.multCompCorr);
    % [~, omnH] = getSignificanceFromP(omnP, stg.Q, stg.multCompCorr);
    positionCm = [10, 25, stg.singleColumnWidth, min(25, stg.szCharHe*numChar)];
    [plotTF, plotName] = createFigPos(positionCm);
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp('Circadian statistics (subject-wise):')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, szCharTbl.Properties.VariableNames, kchar);
            disp(['    ', ylbl])
            disp(str(kchar))
        end
    end
    if plotTF
        % Calculate axes positions
        % % % % % % % % margGlob = [0 0 0 0]; % Left, bottom, right, top
        % % % % % % % % marg = [0.1 0.5 0.1 0.5];
        [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlobCi, stg.margCi);
        h.f.(plotName).Units = stg.units;
        
        % Plot the data
        numCol = ceil((numChar+1)^0.8);
        numRow = ceil((numChar+1)/numCol);
        axWiHe = 0.8*min(spWi/numCol, spHe/numRow);
        for kchar = 1 : numChar
            h.a.(plotName)(1, kchar) = axes('Units', stg.units, 'Position', [...
                spx(mod(1 - 1, numc) + 1) + mod(kchar-1, numCol)*spWi/numCol + 0.1*spWi/numCol,...
                spy(ceil(1/numc)) + (numRow-ceil(kchar/numCol))*spHe/numRow + 0.0*spHe/numRow,...
                axWiHe, ...
                axWiHe],...
                'NextPlot', 'add', 'Visible', 'off', 'XLimMode', 'manual', 'YLimMode', 'manual');
            hax = h.a.(plotName)(1, kchar);
            
            r1 = NaN(stg.numSubj, numel(rrTbl{1, 1}));
            for ksubj = 1 : height(cirTbl)
                % Subject's circular bar graph
                r1(ksubj, :) = rrTbl{ksubj, kchar};
                r1(ksubj, :) = r1(ksubj, :)/max(r1(ksubj, :));
                p = NaN(1, (numel(ed)-1)*(numBetw+2) + 1);
                r = NaN(1, (numel(ed)-1)*(numBetw+2) + 1);
                for ke = 1 : numel(ed) - 1
                    p((ke-1)*(numBetw+2) + 1  :  ke*(numBetw+2)) = linspace(ed(ke), ed(ke+1), numBetw+2);
                    r((ke-1)*(numBetw+2) + 1  :  ke*(numBetw+2)) = r1(ksubj, ke);
                end
                p(end) = ed(1);
                p = p*2*pi;
                r(end) = r1(ksubj, 1);

                [x, y] = pol2cart(-p - pi/2, r);
                h.p.(plotName)(ksubj, kchar, 1) = plot(hax, x, y, 'Color', 1 - 0.5*(1 -stg.subjColor(ksubj, :)));
                % Subject's mean resultant vector
                p = [0, cirTbl{ksubj, kchar}];
                % r = rrCirTbl{ksubj, kchar};
                r = [0 1];
                [x, y] = pol2cart(-p - pi/2, r);
                h.p.(plotName)(ksubj, kchar, 1) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-',...
                    'LineWidth', stg.lnWiFit, 'Color', stg.subjColor(ksubj, :)); % Not normalized
            end
            % Circular bar graph mean
            r1m = mean(r1, 1, 'omitmissing');
            pm = NaN(1, (numel(ed)-1)*(numBetw+2) + 1);
            rm = NaN(1, (numel(ed)-1)*(numBetw+2) + 1);
            for ke = 1 : numel(ed) - 1
                pm((ke-1)*(numBetw+2) + 1  :  ke*(numBetw+2)) = linspace(ed(ke), ed(ke+1), numBetw+2);
                rm((ke-1)*(numBetw+2) + 1  :  ke*(numBetw+2)) = r1m(1, ke);
            end
            pm(end) = ed(1);
            pm = pm*2*pi;
            rm(end) = rm(1);
            [x, y] = pol2cart(-pm - pi/2, rm);
            h.p.(plotName)(1, kchar, 2) = plot(hax, x, y, 'Color', 'k', 'LineWidth', 1.5);
            % Circle
            circleP = (0 : 1 : 360)/360*2*pi;
            circleR = 1;
            % if ~all(isnan(rrCirTbl{:, kchar}))
            %     circleR = ceilToEven1sig(max(rrCirTbl{:, kchar}, [], 'all'));
            % else
            %     circleR = 1;
            % end
            [x, y] = pol2cart(-circleP - pi/2, circleR);
            h.p.(plotName)(1, kchar, 3) = plot(x, y, 'k');
            hax.XLim = max([x, y])*[-1 1];
            hax.YLim = hax.XLim;
            % Night shading
            circleP = (-90 : 1 : 90)/360*2*pi;
            % colo = [linspace(1, 0, 90), 1, linspace(0, 1, 90)]'*[1 1 1];
            colo = [linspace(1, 0.7, 90), 1, linspace(0.7, 1, 90)]'*[1 1 1];
            colo = permute(colo, [1 3 2]);
            [x, y] = pol2cart(-circleP - pi/2, circleR);
            z = -10*ones(size(x));
            h.pa.(plotName)(1, kchar, 1) = patch(x, y, z, colo, 'EdgeColor', 'none');
            % Mean resultant vector
            resVecWi = 1.5;
            prv = [0, resM(kchar)];
            rrv = [0, resR(kchar)];
            [x, y] = pol2cart(-prv - pi/2, rrv);
            h.p.(plotName)(1, kchar, 4) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', 'k');
            % Arrowhead
            [x(1), y(1)] = pol2cart(-(prv(2)-5/6*pi) - pi/2, rrv(2)/4);
            x(1) = sum(x);
            y(1) = sum(y);
            h.p.(plotName)(1, kchar, 5) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', 'k');
            x(1) = 0; y(1) = 0;
            [x(1), y(1)] = pol2cart(-(prv(2)+5/6*pi) - pi/2, rrv(2)/4);
            x(1) = sum(x);
            y(1) = sum(y);
            h.p.(plotName)(1, kchar, 6) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', 'k');
            x = 0;
            y = hax.YLim(2);
            % text(hax, x, y, num2str(y), 'Interpreter', 'tex', 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
            % ylbl = getYlbl(plotName, szCharTbl.Properties.VariableNames, kchar);
            ylbl = stg.szCharNameShort(kchar);
            text(hax, x, y, ylbl, 'Interpreter', 'tex', 'FontSize', stg.axFontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
            strP = string(['  p=', num2str(rayP(1, kchar), 2)]);
            % ht = text(hax, xrv(2)/2, yrv(2)/2, strP, 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment', 'middle', 'BackgroundColor','w');
            ht = text(hax, 0, -1, strP, 'FontSize', stg.axFontSize, 'HorizontalAlignment','center', 'VerticalAlignment', 'top', 'BackgroundColor','none');
            if rayH(1, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
        end
    end
end
function clustSzCharTbl = clTerminDur(szCharTbl, cl, szBelongsToClust)
    global stg
    numChar = numel(stg.szCharToPlot);
    numCl = size(cl, 2);
    clustSzCharTbl = table('Size', [1, 4*numChar], 'VariableTypes', repelem("double", 4*numChar));
    clustSzCharTbl{1, :} = deal(NaN);
    for kchar = 1 : numChar
        clustSzCharTbl.Properties.VariableNames{1, (kchar-1)*4 + 1} = char(stg.szCharToPlot(kchar) + "WithinCl");
        clustSzCharTbl.Properties.VariableNames{1, (kchar-1)*4 + 2} = char(stg.szCharToPlot(kchar) + "TerminCl");
        clustSzCharTbl.Properties.VariableNames{1, (kchar-1)*4 + 3} = char(stg.szCharToPlot(kchar) + "BetweenCl");
        clustSzCharTbl.Properties.VariableNames{1, (kchar-1)*4 + 4} = char(stg.szCharToPlot(kchar) + "TermWithCl");
    end
    for kchar = 1 : numChar
        withinClust = NaN(numCl, 1);
        terminClust = NaN(numCl, 1);
        terminWithinClust = NaN(numCl, 1);
        for kcl = 1 : numCl
            withinClust(kcl) = mean(cl(kcl).szCharTbl.(stg.szCharToPlot(kchar))(1 : end - 1), "omitnan");
            terminClust(kcl) = cl(kcl).szCharTbl.(stg.szCharToPlot(kchar))(end);
            terminWithinClust(kcl) = terminClust(kcl) - withinClust(kcl);
        end
        clustSzCharTbl{1, (kchar-1)*4 + 1} = mean(withinClust, "omitnan");
        clustSzCharTbl{1, (kchar-1)*4 + 2} = mean(terminClust, "omitnan");
        clustSzCharTbl{1, (kchar-1)*4 + 3} = mean(szCharTbl.(stg.szCharToPlot(kchar))(~szBelongsToClust), "omitnan");
        clustSzCharTbl{1, (kchar-1)*4 + 4} = mean(terminWithinClust, "omitnan");
    end
end
function clTerminDurPop(clustSzCharTbl)
    disp('Terminating seizure of cluster - duration')
    withinCl = clustSzCharTbl{:, 5};
    getBasicStats(withinCl, 'withinCl')
    terminCl = clustSzCharTbl{:, 6};
    getBasicStats(terminCl, 'terminCl')
    betweenCl = clustSzCharTbl{:, 7};
    getBasicStats(betweenCl, 'betweenCl')
    termWithCl = clustSzCharTbl{:, 8};
    getBasicStats(termWithCl, 'termWithCl')
    termin_within_sigP = signrank(terminCl - withinCl)
    termiWith_sigP = signrank(termWithCl)

    disp('Terminating seizure of cluster - power')
    withinCl = clustSzCharTbl{:, 9};
    getBasicStats(withinCl, 'withinCl')
    terminCl = clustSzCharTbl{:, 10};
    getBasicStats(terminCl, 'terminCl')
    betweenCl = clustSzCharTbl{:, 11};
    getBasicStats(betweenCl, 'betweenCl')
    termWithCl = clustSzCharTbl{:, 12};
    getBasicStats(termWithCl, 'termWithCl')
    termin_within_sigP = signrank(terminCl - withinCl)
    termWith_sigP = signrank(termWithCl)
end
% Plot seizure and signal characteristics analyses in one figure
function plotSaChar(ksubj, subjInfo, szCharTbl, siCharTbl)
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSaCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl); % Refer to the created axes using h.f.(thisFcnName)
        szCharTblVarNames = szCharTbl.Properties.VariableNames;
        siCharTblVarNames = siCharTbl.Properties.VariableNames;
        charTblVarNames = [szCharTblVarNames, siCharTblVarNames];
        formatAxesInd(plotName, ksubj, subjInfo, charTblVarNames, 'Age (days)', "analysisPeriod")
    end
end
function plotSaPsd(ksubj, subjInfo, szCharTbl, siCharTbl)
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSaCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl); % Refer to the created axes using h.f.(thisFcnName)
        szCharTblVarNames = szCharTbl.Properties.VariableNames;
        siCharTblVarNames = siCharTbl.Properties.VariableNames;
        charTblVarNames = [szCharTblVarNames, siCharTblVarNames];
        formatAxesInd(plotName, ksubj, subjInfo, charTblVarNames, 'Age (days)', "analysisPeriod")
    end
end
function plotSaCharCWT(ksubj, subjInfo, szCharTbl, siCharTbl)
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSaCharDataCWT(plotName, ksubj, subjInfo, szCharTbl, siCharTbl); % Refer to the created axes using h.f.(thisFcnName)

        % % % szCharTblVarNames = {'szOnsN'};
        % % % siCharTblVarNames = [siCharTbl.Properties.VariableNames, siCharTbl.Properties.VariableNames];
        % % % charTblVarNames = [szCharTblVarNames, siCharTblVarNames];
        % % % formatAxesInd(plotName, ksubj, subjInfo, charTblVarNames, 'Age (days)', "analysisPeriod")
    end
end
function [fitTbl, xxFitTbl, yyFitTbl] = plotSaCharWhFit(ksubj, subjInfo, szCharTbl, siCharTbl)
    global stg
    global h
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    [fitTbl, ed, ~, yySzTbl, ~, ~, xxFitTbl, yyFitTbl] = fitSaCharWh(subjInfo, szCharTbl, siCharTbl); % Compute the trend
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSaCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl); % Refer to the created axes using h.f.(thisFcnName)
        
        % Plot the trend
        for kchar = 1 : numSzChar
            y = yySzTbl{1, kchar};
            % Plotting patch (bar graph)
            hax = h.a.(plotName)(ksubj, kchar);
            pax = [ed(1), ed(1), repelem(ed(2:end-1), 3), ed(end), ed(end), ed(1)];
            pay = repelem(y, 3);
            pay(1 : 3 : end) = 0;
            pay = [pay, 0, 0]; %#ok<AGROW>
            pay(isnan(pay)) = 0;
            patch(hax, pax, pay, 1 - 0.2*(1 - stg.subjColor(ksubj, :)))
            if kchar == 1
                h.p.(plotName)(ksubj, kchar, 2).YData(2 : 3 : end - 1) = ceilToEven1sig(max(pay)); % Change the height of seizures in the plot of sz rate
            end
        end
        for kchar = 1 : numChar
            % Plotting fitted line
            x = xxFitTbl{1, kchar};
            y = yyFitTbl{1, kchar};
            cil = fitTbl{1, kchar + 1*numChar};
            cih = fitTbl{1, kchar + 2*numChar};
            if cil*cih > 0 % Both have the same sign and none is zero
                lnstyle = '-';
            else
                lnstyle = '--';
            end
            hax = h.a.(plotName)(ksubj, kchar);
            plot(hax, x, y, 'Marker', 'none', 'LineStyle', lnstyle, 'LineWidth', stg.lnWiFitMean, 'Color', stg.subjColor(ksubj, :)*stg.subjColorSubjMeanMult);
            if kchar ~= numChar
                hax.Children([1 2 3 4]) = hax.Children([1 3 4 2]);
            end
        end
        charTblVarNames = [szCharTbl.Properties.VariableNames, siCharTbl.Properties.VariableNames];
        formatAxesInd(plotName, ksubj, subjInfo, charTblVarNames, 'Age (days)', "analysisPeriod")
    end
end
function plotSaCharWhFitAllPop(fitTbl, xxFitTbl, yyFitTbl, szCharTbl, siCharTbl) %#ok<INUSD>
    global stg
    global h
    normalizedTF = true;
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    charTblVarNames = [szCharTbl.Properties.VariableNames, siCharTbl.Properties.VariableNames];
    str = strings(1, numChar); medFit = NaN(1, numChar); sigRanP = NaN(1, numChar);
    for kchar = 1 : numChar
        [str(kchar), ~, ~, medFit(kchar), ~, sigRanP(kchar)] = getBasicStats(fitTbl.(stg.saCharToPlot(kchar) + "fitSlope"), 'Slope');
    end
    [sigRanP, sigRanH] = getSignificanceFromP(sigRanP, stg.Q, stg.multCompCorr);
    if stg.plotCharsStacked
        figHeight = min(25, stg.saCharHe*numChar);
    else
        figHeight = 5; % Centimeters
    end
    positionCm = [10, 10, stg.singleColumnWidth, figHeight];
    [plotTF, plotName] = createFigPos(positionCm);
    
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp('Whole recording:')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, charTblVarNames, kchar);
            disp(['    ', ylbl])
            disp(str(kchar))
            cil = fitTbl{:, kchar + 1*numChar};
            cih = fitTbl{:, kchar + 2*numChar};
            numPos = sum(cil > 0 & cih > 0);
            numNeg = sum(cil < 0 & cih < 0);
            numNS = sum(cil < 0 & cih > 0);
            disp(['        ', 'Positive: ', num2str(numPos), ', Negative: ', num2str(numNeg), ', Non-sig: ', num2str(numNS)])
        end
    end
    if plotTF
        if ~isfield(h.a, plotName)
            if stg.plotCharsStacked
                createAxesAllPopStack(plotName)
            else
                createAxesAllPopSlopeBox(plotName)
            end
        end
        strP = strings(1, numChar);
        for kchar = 1 : numChar
            for ksubj = 1 : stg.numSubj
                hax = h.a.(plotName)(1, kchar);
                if normalizedTF
                    x = [0 1];
                    fitSlope = fitTbl{ksubj, kchar};
                    y = x*fitSlope;
                else
                    x = xxFitTbl{ksubj, kchar};
                    y = yyFitTbl{ksubj, kchar};
                end
                plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', stg.lnWiFit, 'Color', stg.subjColor(ksubj, :)); % Not normalized
            end
            if normalizedTF
                y = x*medFit(1, kchar);
                plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            else
                x = stg.acrossSubjectStat(xxFitTbl{:, kchar}, 1, 'omitmissing');
                y = stg.acrossSubjectStat(yyFitTbl{:, kchar}, 1, 'omitmissing');
                plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            end
            signMed = getSignChar(medFit(kchar));
            strP(kchar) = string(['  p=', num2str(sigRanP(1, kchar), 2), '(', signMed, ')']);
        end
        % Box plot, violin plot, swarm chart
        if ~stg.plotCharsStacked
            for kchar = 1 : numChar
                hax = h.a.(plotName)(2, kchar);
                y = fitTbl{:, kchar};
                for kp = 1 : numel(stg.statsPlotType)
                    switch stg.statsPlotType(kp)
                        case "box"
                            boxplot(hax, y, 'Colors', 'k', 'Widths', 8);
                        case "violin"
                            violinMedianConf(hax, y, stg.uniformSubjectColor, stg.numBootstrapSamples, 0.05)
                        case "swarm"
                            swarmchart(hax, ones(size(y)), y, 3, stg.uniformSubjectColor, "filled")
                    end
                    [yl, ~] = getYLimYTick([min(y), max(y)]);
                    hax.YLim = yl; hax.YTickLabel = []; hax.XTick = []; hax.TickLength = 1.5*hax.TickLength;
                    hax.Box = 'on';
                end
            end
        end
        % Formatting
        if stg.plotCharsStacked
            formatAxesAllPop(plotName, charTblVarNames, 'Normalized time', "none", ["any", "any"])
        else
            for kchar = 1 : numChar
                hax = h.a.(plotName)(1, kchar);
                if kchar > numChar - 2
                    xlabel(hax, 'Normalized time')
                end
                y = fitTbl{:, kchar};
                [yl, ~] = getYLimYTick([min(y), max(y)]);
                hax.YLim = yl;
                ylabel(hax, getYlbl(plotName, charTblVarNames, kchar))
                hax.Box = 'on';
            end
        end
        for kchar = 1 : numChar
            hax = h.a.(plotName)(1, kchar);
            % xt = mean(hax.XLim);
            xt = 0;
            % yt = mean(hax.YLim);
            yt = hax.YLim(2);
            ht = text(hax, xt, yt, strP(kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','left', 'VerticalAlignment','top', 'BackgroundColor','none');
            if sigRanH(1, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
        end
    end
end
function plotSaCharCl(ksubj, subjInfo, szCharTbl, siCharTbl, cl)
    % ksubj ...... subject index
    % subjInfo ... table with subject info
    % szCharTbl .. table with seizure characteristics
    % siCharTbl .. table with signal characteristics
    % cl ......... structure with clusters of given mouse
    % mFitTbl .... table with mean slopes, offsets and confidence intervlas (the CI mean has no real meaning though)
    % mxFit0 ..... mean cluster duration (probably not further used)
    % myFit0 ..... mean difference between value of the fit at the offset - onset of the cluster for each characteristic (probably not further used)
    global stg
    global h
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    fitPosition = ["before", "during", "after"];
    [plotTF, plotName] = createFigInd(1);
    numPoints = 10; % Each cluster have different length. Number of points to which we will interpolate the characteristic.
    myyInterpTblCell = cell(1, numel(fitPosition)); xCell = cell(1, numel(fitPosition));
    
    for kfp = 1 : numel(fitPosition) % k fit position
        xxInterpTbl = table; yyInterpTbl = table;

        % Seizure char
        [~, ~, xxTbl, yyTbl, xxFitTbl, ~] = fitSzCharCl(subjInfo, siCharTbl, cl);
        vrnmSz = xxFitTbl.Properties.VariableNames;
        if fitPosition(kfp) == "during" && ~isempty(xxTbl)
            if isempty(xxFitTbl)
                myyTbl = varfun(@(x) [x; NaN(1, size(yyTbl{1, 1}, 2))], yyTbl, 'OutputFormat', 'table'); % Average of the data
            else
                myyTbl = stg.withinSubjectStat(varfun(@(x) x, yyTbl, 'OutputFormat', 'table'), 1); % Average of the data
            end
            myyTbl.Properties.VariableNames = vrnmSz;
            xxInterpTbl = table; yyInterpTbl = table;
            for kcl = 1 : height(xxTbl)
                for kchar = 1 : numSzChar
                    % % % xxInterpTbl{kcl, kchar} = linspace(eded(kcl, 1), eded(kcl, end), numPoints);
                    xxInterpTbl{kcl, kchar} = linspace(xxTbl{kcl, kchar}(1), xxTbl{kcl, kchar}(end), numPoints);
                    yyInterpTbl{kcl, kchar} = interp1(xxTbl{kcl, kchar}, yyTbl{kcl, kchar}, xxInterpTbl{kcl, kchar}, "linear");
                end
            end
        else
            heiXxTbl = height(xxTbl);
            if heiXxTbl > 1
                numClust = heiXxTbl;
            else
                numClust = 1;
            end
            for kcl = 1 : numClust
                for kchar = 1 : numSzChar
                    xxInterpTbl{kcl, kchar} = NaN(1, numPoints);
                    yyInterpTbl{kcl, kchar} = NaN(1, numPoints);
                end
            end
        end
        
        % Signal char
        [~, xxTbl, yyTbl, ~, ~] = fitSiCharCl(subjInfo, siCharTbl, cl, fitPosition(kfp));
        vrnmSi = xxTbl.Properties.VariableNames;
        if isempty(xxTbl)
            for kchar = 1 : numSiChar
                xxInterpTbl{1, numSzChar + kchar} = NaN(1, numPoints);
                yyInterpTbl{1, numSzChar + kchar} = NaN(1, numPoints);
            end
            plotColor = stg.subjColor(ksubj, :);
        else
            for kcl = 1 : height(xxTbl)
                for kchar = 1 : numSiChar
                    if numel(xxTbl{kcl, kchar}{1}) < 2
                        xxInterpTbl{kcl, numSzChar + kchar} = NaN(1, numPoints);
                        yyInterpTbl{kcl, numSzChar + kchar} = NaN(1, numPoints);
                    else
                        xxInterpTbl{kcl, numSzChar + kchar} = linspace(xxTbl{kcl, kchar}{1}(1), xxTbl{kcl, kchar}{1}(end), numPoints);
                        yyInterpTbl{kcl, numSzChar + kchar} = interp1(xxTbl{kcl, kchar}{1}, yyTbl{kcl, kchar}{1}, xxInterpTbl{kcl, numSzChar + kchar}, "linear");
                    end
                end
            end
            plotColor = ones(size(yyTbl, 1), 1)*stg.subjColor(ksubj, :);
        end

        % Common
        vrnm = [vrnmSz, vrnmSi];
        xxInterpTbl.Properties.VariableNames = vrnm;
        yyInterpTbl.Properties.VariableNames = vrnm;
        myyInterpTblCell{kfp} = stg.withinSubjectStat(yyInterpTbl, 'omitmissing');
        switch fitPosition(kfp)
            case "before"
                xCell{kfp} = linspace(0, 0.33, numPoints);
            case "during"
                xCell{kfp} = linspace(0.34, 0.66, numPoints);
            case "after"
                xCell{kfp} = linspace(0.67, 1, numPoints);
        end

        % Plot
        if plotTF
            numCl = plotPeriEventData(plotName, ksubj, xCell{kfp}, yyInterpTbl, myyInterpTblCell{kfp}, plotColor); % Refer to the created axes using h.a.(plotName)
'Some variables'
numCl
numChar
thisAxes_ = h.a.(plotName)
ksubj
kchar

            for kchar = 1 : numChar
                hax = h.a.(plotName)(ksubj, kchar);
                h.t.(plotName)(ksubj, kchar, kfp, 1) = text(hax, mean(xCell{kfp}), h.a.(plotName)(ksubj, kchar).YLim(2), ['n=', num2str(numCl(kchar))],...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', stg.statFontSize); % A four dimensional matrix. Cool, isn't it.
            end
            h.t.(plotName)(ksubj, kchar, kfp, 2) = text(mean(xCell{kfp}), 0, fitPosition(kfp),...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', stg.statFontSize);
        end
    end
    % Formatting
    if plotTF
        set(h.a.(plotName), 'XLim', [0 1])
        charTblVarNames = [szCharTbl.Properties.VariableNames, siCharTbl.Properties.VariableNames];
        formatAxesInd(plotName, ksubj, subjInfo, charTblVarNames, 'Normalized time', "none")
        set(h.a.(plotName), 'XTick', [])
        for kchar = 1 : numChar
            yl = h.a.(plotName)(ksubj, kchar).YLim;
            h.pa.(plotName)(ksubj, kchar, 1) = patch(h.a.(plotName)(ksubj, kchar), [0 0.33 0.33 0], [yl(1) yl(1) yl(2) yl(2)],...
                1 - 0.15*(1-stg.fitColor(1, :)), 'ZData', [-1 -1 -1 -1], 'EdgeColor', 'none');
            h.pa.(plotName)(ksubj, kchar, 2) = patch(h.a.(plotName)(ksubj, kchar), [0.34 0.66 0.66 0.34], [yl(1) yl(1) yl(2) yl(2)],...
                1 - 0.15*(1-stg.fitColor(2, :)), 'ZData', [-1 -1 -1 -1], 'EdgeColor', 'none');
            h.pa.(plotName)(ksubj, kchar, 3) = patch(h.a.(plotName)(ksubj, kchar), [0.67 1 1 0.67], [yl(1) yl(1) yl(2) yl(2)],...
                1 - 0.15*(1-stg.fitColor(3, :)), 'ZData', [-1 -1 -1 -1], 'EdgeColor', 'none');
            for kfp = 1 : numel(fitPosition)
                h.t.(plotName)(ksubj, kchar, kfp, 1).Position(2) = h.a.(plotName)(ksubj, kchar).YLim(2);
            end
        end
    end
end
function [mFitTblBe, mFitTblDu, mFitTblAf] = plotSaCharClFit(ksubj, subjInfo, szCharTbl, siCharTbl, cl)
    % ksubj ...... subject index
    % subjInfo ... table with subject info
    % szCharTbl .. table with seizure characteristics
    % siCharTbl .. table with signal characteristics
    % cl ......... structure with clusters of given mouse
    % mFitTbl .... table with mean slopes, offsets and confidence intervlas (the CI mean has no real meaning though)
    % mxFit0 ..... mean cluster duration (probably not further used)
    % myFit0 ..... mean difference between value of the fit at the offset - onset of the cluster for each characteristic (probably not further used)
    global stg
    global h
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    fitPosition = ["before", "during", "after"];
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSaCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl); % Refer to the created axes using h.a.(thisFcnName)
    end
    numFitPos = numel(fitPosition);
    mFitTblCell = cell(1, numFitPos);
    maxHeight = zeros(1, numChar);
    for kfp = 1 : numFitPos % Fit position is before, during or after the cluster
        [fitTbl, eded, ~, yySzTbl, ~, ~, xxFitTbl, yyFitTbl] = fitSaCharCl(subjInfo, siCharTbl, cl, fitPosition(kfp));
        mFitTblCell{kfp} = stg.withinSubjectStat(fitTbl, "omitmissing");
        if plotTF
            for kchar = 1 : numChar
                if isempty(xxFitTbl)
                    continue
                end
                for kcl = 1 : size(xxFitTbl, 1)
                    % if clust(kcl).nested
                    %     continue
                    % end
                    hax = h.a.(plotName)(ksubj, kchar);
                    if kchar <= numSzChar && kfp == 2
                        % Bar graph (it is actually a patch but looks like a bar graph)
                        y = yySzTbl{kcl, kchar};
                        ed = eded(kcl, :);
                        % Plotting
                        hax = h.a.(plotName)(ksubj, kchar);
                        pax = [ed(1), ed(1), repelem(ed(2:end-1), 3), ed(end), ed(end), ed(1)];
                        pay = repelem(y, 3);
                        pay(1 : 3 : end) = 0;
                        pay = [pay, 0, 0]; %#ok<AGROW>
                        pay(isnan(pay)) = 0;
                        maxHeight(kchar) = max([maxHeight(kchar), pay]);
                        h.pa.(plotName)(ksubj, kchar, kcl) = patch(hax, pax, pay, 1 - 0.2*(1 - stg.subjColor(ksubj, :)), 'Tag', 'szCharClDuBar');
                        % Fit plot
                        x = xxFitTbl{kcl, kchar};
                        y = yyFitTbl{kcl, kchar};
                        maxHeight(kchar) = max([maxHeight(kchar), y]);
                        h.p.(plotName)(ksubj, kchar, 4) = plot(hax, x, y, ...
                           'Marker', 'none', 'LineWidth', 2, 'LineStyle', '-', 'Color', stg.fitColor(kfp, :), 'Marker', 'none');
                        % % % h.p.(plotName)(ksubj, kchar, 3) = plot(hax, x, y, ...
                        % % %     'Marker', 'none', 'LineStyle', '--', 'LineWidth', 1, 'Color', stg.subjColor(ksubj, :)*stg.subjColorSubjMeanMult, 'Tag', 'szCharClDuFit');
                        if kchar <= numSzChar && kfp == 2
                            uistack(h.pa.(plotName)(ksubj, kchar, kcl), 'bottom')
                        end
                        uistack(h.p.(plotName)(ksubj, kchar, 4), 'top')
                    else
                        h.p.(plotName)(ksubj, kchar, 4) = plot(hax, xxFitTbl{kcl, kchar}, yyFitTbl{kcl, kchar}, ...
                           'Marker', 'none', 'LineWidth', 2, 'LineStyle', '-', 'Color', stg.fitColor(kfp, :), 'Marker', 'none');
                    end
                end
                if kchar == 1
                    h.p.(plotName)(ksubj, kchar, 2).YData(2 : 3 : end-1) = ceilToEven1sig(maxHeight(kchar)); % Change the height of seizures in the plot of sz rate
                end
            end
        end
    end
    if plotTF
        charTblVarNames = [szCharTbl.Properties.VariableNames, siCharTbl.Properties.VariableNames];
        formatAxesInd(plotName, ksubj, subjInfo, charTblVarNames, 'Age (days)', "analysisPeriod")
        if stg.clusterExampleMouseJc20190509_2 && (subjInfo.subjNmOrig == "jc20190509_2" || subjInfo.subjNm == "Mouse13")
            timeLimits = [82, 86];
            ymax = NaN(1, numChar); % Initialize variable for YLim(2) of each axes
            % The first characteristic is the seizure rate. We will use the patch (bar graph) data to get y-limits
            numPatch = size(h.pa.(plotName)(ksubj, 1, :), 3); % How many clusters
            allPatchXData = NaN(1000, 1); % Initialize
            allPatchYData = NaN(1000, 1);
            counter = 0;
            for kp = 1 : numPatch
                numData = numel(h.pa.(plotName)(ksubj, 1, kp).XData);
                allPatchXData(counter + 1 : counter + numData) = h.pa.(plotName)(ksubj, 1, kp).XData;
                allPatchYData(counter + 1 : counter + numData) = h.pa.(plotName)(ksubj, 1, kp).YData;
                counter = counter + numData;
            end
            allPatchXData = allPatchXData(~isnan(allPatchXData));
            allPatchYData = allPatchYData(~isnan(allPatchYData));
            ind = allPatchXData > timeLimits(1) & allPatchXData < timeLimits(2); % Choose only data within the time limits
            ymax(1) = max(allPatchYData(ind));
            % Seizure characteristics will use individual seizure data which are in h.p.(plotName)(ksubj, kchar, 2)
            for kchar = 2 : numSzChar
                szInd = h.p.(plotName)(ksubj, kchar, 2).XData > timeLimits(1) & h.p.(plotName)(ksubj, kchar, 2).XData < timeLimits(2);
                ymax(kchar) = max(h.p.(plotName)(ksubj, kchar, 2).YData(szInd));
            end
            % Signal characteristics will use data which are in h.p.(plotName)(ksubj, kchar, 3)
            for kchar = numSzChar + (1 : numSiChar)
                szInd = h.p.(plotName)(ksubj, kchar, 3).XData > timeLimits(1) & h.p.(plotName)(ksubj, kchar, 3).XData < timeLimits(2);
                ymax(kchar) = max(h.p.(plotName)(ksubj, kchar, 3).YData(szInd));
            end
            % Set the axes limits and ticks
            for kchar = 1 : numChar
                h.a.(plotName)(ksubj, kchar).XLim = timeLimits;
                h.a.(plotName)(ksubj, kchar).YLim(2) = ceilToEven1sig(ymax(kchar));
                h.a.(plotName)(ksubj, kchar).YTick = [h.a.(plotName)(ksubj, kchar).YLim(1), mean(h.a.(plotName)(ksubj, kchar).YLim), h.a.(plotName)(ksubj, kchar).YLim(2)];
            end
        
        end
    end
    mFitTblBe = mFitTblCell{1};
    mFitTblDu = mFitTblCell{2};
    mFitTblAf = mFitTblCell{3};
end
function plotSaCharClFitAllPop(fitTblBe, fitTblDu, fitTblAf, szCharTbl, siCharTbl)
    global stg
    global h
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    charTblVarNames = [szCharTbl.Properties.VariableNames, siCharTbl.Properties.VariableNames];
    numSubj = stg.numSubj;
    str = strings(3, numChar); medFit = NaN(3, numChar); sigRanP = NaN(3, numChar);
    for kchar = 1 : numChar
        [str(1, kchar), ~, ~, medFit(1, kchar), ~, sigRanP(1, kchar)] = getBasicStats(fitTblBe.(stg.saCharToPlot(kchar) + "fitSlope"), 'Slope');
        [str(2, kchar), ~, ~, medFit(2, kchar), ~, sigRanP(2, kchar)] = getBasicStats(fitTblDu.(stg.saCharToPlot(kchar) + "fitSlope"), 'Slope');
        [str(3, kchar), ~, ~, medFit(3, kchar), ~, sigRanP(3, kchar)] = getBasicStats(fitTblAf.(stg.saCharToPlot(kchar) + "fitSlope"), 'Slope');
    end
    [sigRanP, sigRanH] = getSignificanceFromP(sigRanP, stg.Q, stg.multCompCorr);
    positionCm = [20, 10, stg.singleColumnWidth, min(25, stg.saCharHe*numChar)];
    [plotTF, plotName] = createFigPos(positionCm);
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp('Before cluster:')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, charTblVarNames, kchar);
            disp(['    ', ylbl])
            disp(str(1, kchar))
        end
        disp('During cluster:')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, charTblVarNames, kchar);
            disp(['    ', ylbl])
            disp(str(2, kchar))
        end
        disp('After cluster:')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, charTblVarNames, kchar);
            disp(['    ', ylbl])
            disp(str(3, kchar))
        end
    end
    if plotTF
        if ~isfield(h.a, plotName)
            if stg.plotCharsStacked
                createAxesAllPopStack(plotName)
            else
                createAxesAllPopSlopeBox(plotName)
            end
        end
        strP = strings(3, numChar); xt = NaN(3, 1);
        for kchar = 1 : numChar
            if stg.plotCharsStacked
                hax(1) = h.a.(plotName)(1, kchar);
                hax(2) = h.a.(plotName)(1, kchar);
                hax(3) = h.a.(plotName)(1, kchar);
                x = [0 0.33; 0.34 0.66; 0.37 1];
            else
                if kchar <= numSzChar
                    hax(1) = h.a.(plotName)(1, kchar);
                    hax(2) = h.a.(plotName)(1, kchar);
                    hax(3) = h.a.(plotName)(1, kchar);
                else
                    hax(1) = h.a.(plotName)(1, numChar + 2*(kchar-numSzChar-1) + 1);
                    hax(2) = h.a.(plotName)(1, kchar);
                    hax(3) = h.a.(plotName)(1, numChar + 2*(kchar-numSzChar-1) + 2);
                end
                x = [0 1; 0 1; 0 1];
            end
            for ksubj = 1 : numSubj
                % Before cluster
                fitSlope = fitTblBe{ksubj, kchar + 0*numChar};
                y = [-1 0]*fitSlope;
                plot(hax(1), x(1,:), y, 'Marker', 'none', 'LineWidth', stg.lnWiFit, 'LineStyle', '-', 'Color', stg.fitColor(1, :), 'Marker', 'none');
                % During cluster
                fitSlope = fitTblDu{ksubj, kchar + 0*numChar};
                y = [0 1]*fitSlope;
                plot(hax(2), x(2,:), y, 'Marker', 'none', 'LineWidth', stg.lnWiFit, 'LineStyle', '-', 'Color', stg.fitColor(2, :), 'Marker', 'none');
                % After cluster
                fitSlope = fitTblAf{ksubj, kchar + 0*numChar};
                y = [0 1]*fitSlope;
                plot(hax(3), x(3,:), y, 'Marker', 'none', 'LineWidth', stg.lnWiFit, 'LineStyle', '-', 'Color', stg.fitColor(3, :), 'Marker', 'none');
            end

            %  Plot medians
            % Before cluster
            xt(1) = mean(x(1,:)); % x-coordinate for text with the p-value
            y = [-1 0]*medFit(1, kchar);
            plot(hax(1), x(1,:), y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            signMed = getSignChar(medFit(1, kchar));
            strP(1, kchar) = string(['p=', num2str(sigRanP(1, kchar), 2), '(', signMed, ')']);
            % During cluster
            xt(2) = mean(x(2,:));
            y = [0 1]*medFit(2, kchar);
            plot(hax(2), x(2,:), y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            signMed = getSignChar(medFit(2, kchar));
            strP(2, kchar) = string(['p=', num2str(sigRanP(2, kchar), 2), '(', signMed, ')']);
            % After cluster
            xt(3) = mean(x(2,:));
            y = [0 1]*medFit(3, kchar);
            signMed = getSignChar(medFit(3, kchar));
            plot(hax(3), x(3,:), y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            strP(3, kchar) = string(['p=', num2str(sigRanP(3, kchar), 2), '(', signMed, ')']);
            % hax.XLim = [0 1];
        end
        strP(contains(strP, 'n/a')) = "";

        % Box plot, violin plot, swarm chart
        if ~stg.plotCharsStacked
            % During
            for kchar = 1 : numChar
                hax2 = h.a.(plotName)(2, kchar);
                y = fitTblDu{:, kchar};
                for kp = 1 : numel(stg.statsPlotType)
                    switch stg.statsPlotType(kp)
                        case "box"
                            boxplot(hax2, y, 'Colors', 'k', 'Widths', 8);
                        case "violin"
                            violinMedianConf(hax2, y, stg.fitColor(2, :), stg.numBootstrapSamples, 0.05)
                        case "swarm"
                            swarmchart(hax2, ones(size(y)), y, 3, stg.fitColor(2, :), "filled")
                    end
                    [yl, ~] = getYLimYTick([min(y), max(y)]);
                    hax2.YLim = yl; hax2.YTickLabel = []; hax2.XTick = []; hax2.TickLength = 1.5*hax2.TickLength;
                    hax2.Box = 'on';
                end
            end
            % Before
            for kchar = 1 : numSiChar
                hax2 = h.a.(plotName)(2, numChar + 2*(kchar-1) + 1);
                y = fitTblBe{:, numSzChar + kchar};
                for kp = 1 : numel(stg.statsPlotType)
                    switch stg.statsPlotType(kp)
                        case "box"
                            boxplot(hax2, y, 'Colors', 'k', 'Widths', 8);
                        case "violin"
                            violinMedianConf(hax2, y, stg.fitColor(1, :), stg.numBootstrapSamples, 0.05)
                        case "swarm"
                            swarmchart(hax2, ones(size(y)), y, 3, stg.fitColor(1, :), "filled")
                    end
                    [yl, ~] = getYLimYTick([min(y), max(y)]);
                    hax2.YLim = yl; hax2.YTickLabel = []; hax2.XTick = []; hax2.TickLength = 1.5*hax2.TickLength;
                    hax2.Box = 'on';
                end
            end
            % After
            for kchar = 1 : numSiChar
                hax2 = h.a.(plotName)(2, numChar + 2*(kchar-1) + 2);
                y = fitTblAf{:, numSzChar + kchar};
                for kp = 1 : numel(stg.statsPlotType)
                    switch stg.statsPlotType(kp)
                        case "box"
                            boxplot(hax2, y, 'Colors', 'k', 'Widths', 8);
                        case "violin"
                            violinMedianConf(hax2, y, stg.fitColor(3, :), stg.numBootstrapSamples, 0.05)
                        case "swarm"
                            swarmchart(hax2, ones(size(y)), y, 3, stg.fitColor(3, :), "filled")
                    end
                end
            end
        end
        % Formatting
        if stg.plotCharsStacked
            formatAxesAllPop(plotName, charTblVarNames, 'Normalized time', "none", ["any", "any"])
        else
            for kchar = 1 : numSzChar
                hax = h.a.(plotName)(1, kchar);
                xlabel(hax, 'During')
                y = fitTblDu{:, kchar};
                [yl, ~] = getYLimYTick([min(y), max(y)]);
                hax.YLim = yl;
                hax.YAxis.Exponent = floor(log10(max(abs(yl))));
                ylabel(hax, getYlbl(plotName, charTblVarNames, kchar))
            end
            for kchar = 1 : numSiChar
                y = [fitTblDu{:, numSzChar + kchar}; fitTblBe{:, numSzChar + kchar}; fitTblAf{:, numSzChar + kchar}];
                min(y)
                [yl, ~] = getYLimYTick([min(y), max(y)]);
                hax = h.a.(plotName)(1, numSzChar + kchar);
                hax.YLim = yl; xlabel(hax, 'During'); ylabel(hax, getYlbl(plotName, charTblVarNames, numSzChar + kchar)); hax.YAxis.Exponent = floor(log10(yl(2)));
                hax = h.a.(plotName)(2, numSzChar + kchar);
                hax.YLim = yl; hax.YTickLabel = []; hax.TickLength = 1.5*hax2.TickLength; hax.YAxis.Exponent = floor(log10(yl(2)));
                hax = h.a.(plotName)(1, numChar + 2*(kchar-1) + 1);
                hax.YLim = yl; xlabel(hax, 'Before'); ylabel(hax, getYlbl(plotName, charTblVarNames, numSzChar + kchar)); hax.YAxis.Exponent = floor(log10(yl(2)));
                hax = h.a.(plotName)(2, numChar + 2*(kchar-1) + 1);
                hax.YLim = yl; hax.YTickLabel = []; hax.TickLength = 1.5*hax2.TickLength; hax.YAxis.Exponent = floor(log10(yl(2)));
                hax = h.a.(plotName)(1, numChar + 2*(kchar-1) + 2);
                hax.YLim = yl; xlabel(hax, 'After'); ylabel(hax, getYlbl(plotName, charTblVarNames, numSzChar + kchar)); hax.YAxis.Exponent = floor(log10(yl(2)));
                hax = h.a.(plotName)(2, numChar + 2*(kchar-1) + 2);
                hax.YLim = yl; hax.YTickLabel = []; hax.TickLength = 1.5*hax2.TickLength; hax.YAxis.Exponent = floor(log10(yl(2)));
            end
            set(h.a.(plotName), 'Box', 'on')
            set(h.a.(plotName), 'XTick', [])
            get(h.a.(plotName)(2, :), 'TickLength')
        end
        for kchar = 1 : numChar
            hax = h.a.(plotName)(1, kchar);
            if abs(hax.YLim(1)) > hax.YLim(2)
                yt = hax.YLim(1);
                ht = text(hax, xt(2), yt, strP(2, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'BackgroundColor','none');
            else                
                yt = hax.YLim(2);
                ht = text(hax, xt(2), yt, strP(2, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            end
            if sigRanH(2, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
        end
        for kchar = 1 : numSiChar
            hax = h.a.(plotName)(1, numChar + 2*(kchar-1) + 1);
            if abs(hax.YLim(1)) > hax.YLim(2)
                yt = hax.YLim(1);
                ht = text(hax, xt(1), yt, strP(1, numSzChar + kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'BackgroundColor','none');
            else
                yt = hax.YLim(2);
                ht = text(hax, xt(1), yt, strP(1, numSzChar + kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            end
            if sigRanH(1, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
            hax = h.a.(plotName)(1, numChar + 2*(kchar-1) + 2);
            if abs(hax.YLim(1)) > hax.YLim(2)
                yt = hax.YLim(1);
                ht = text(hax, xt(3), yt, strP(3, numSzChar + kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'BackgroundColor','none');
            else                
                yt = hax.YLim(2);
                ht = text(hax, xt(3), yt, strP(3, numSzChar + kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            end
            if sigRanH(3, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
        end
        h.f.(plotName).Units = 'normalized';
        % % % % % % % % % % % hax.XTick = [];
        % % % % % % % % % % % hax.XLabel = [];
        % % % % % % % % % % % text(0 + 0.33/2, hax.YLim(1) - 0.2*diff(hax.YLim), 'Pre-cluster', 'FontSize',  8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
        % % % % % % % % % % % text(0.34 + 0.33/2, hax.YLim(1) - 0.2*diff(hax.YLim), 'Cluster', 'FontSize',  8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
        % % % % % % % % % % % text(0.67 + 0.33/2, hax.YLim(1) - 0.2*diff(hax.YLim), 'Post-cluster', 'FontSize',  8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
    end
end
function [cirTbl, ppTbl, rrTbl, rayCirTbl, omnCirTbl] = plotSaCharCiFit(ksubj, subjInfo, szCharTbl, siCharTbl)
    global stg
    global h
    numPointsAroundCircle = 97; % Only for smoothness of the circular sections
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    [cirTbl, edSz, ~, rrSzTbl, ~, ppSiTbl, rrSiTbl, ppCirTbl, rrCirTbl, rayCirTbl, omnCirTbl] = fitSaCharCi(szCharTbl, siCharTbl); % Calculate the circular statistics
    % citCir = rayCirTbl{1, :}'; % cit stands for circular test - either Rayleigh or Omnibus - choose here
    citCir = omnCirTbl{1, :}'; % cit stands for circular test - either Rayleigh or Omnibus - choose here
    whichAreNotNan = ~isnan(citCir);
    [citPS, citHS] = getSignificanceFromP(citCir(whichAreNotNan), stg.Q, stg.multCompCorr);
    citP = citCir'; citP(whichAreNotNan) = citPS';
    citH = false(1, 4); citH(whichAreNotNan) = citHS';
    fontWt(~citH) = "normal";
    fontWt(citH) = "bold";
    % Prepare the output (interpolate to same length)
    ppTbl = table('Size', [0, numChar], 'VariableTypes', repelem("cell", 1, numChar), 'VariableNames', stg.saCharToPlot);
    rrTbl = ppTbl;
    for kchar = 1 : numChar
        if kchar <= numSzChar
            % Plotting circular bar graph
            r = rrSzTbl{1, kchar};
            if kchar == 1 % Here, the size of the individual seizure vectors has no meaning, so we normalize histogram to max of the histogram
                r = r/max(r);
            else
                y1 = szCharTbl{:, stg.szCharToPlot(kchar)};
                r = r/max(y1);
            end
            circleCenter = max(r)/1000; % Pure zero minimum radius makes trouble probably due to rounding errors
            numBetw = (numPointsAroundCircle-1)/(numel(edSz)-1) - 3;
            pap = NaN(1, (numel(edSz)-1)*(numBetw+3) + 1); par = pap;
            pap(1) = edSz(1);
            par(1) = 0;
            for ke = 1 : numel(edSz) - 1
                pap((ke-1)*(numBetw+3) + (2 : 2+numBetw+1)) = linspace(edSz(ke), edSz(ke+1), numBetw+2);
                pap((ke-1)*(numBetw+3) + 2+numBetw+2) = edSz(ke+1);
                par((ke-1)*(numBetw+3) + (2 : 2+numBetw+1)) = r(ke);
                par((ke-1)*(numBetw+3) + 2+numBetw+2) = circleCenter;
            end
            pap = pap*2*pi;
            par(isnan(par)) = circleCenter;
            par(end) = circleCenter;
            ppTbl{1, kchar} = {pap};
            rrTbl{1, kchar} = {par};
        else
            ppTbl{1, kchar} = {ppSiTbl{1, kchar - numSzChar}}; %#ok<CCAT1>
            rrTbl{1, kchar} = {rrSiTbl{1, kchar - numSzChar}}; %#ok<CCAT1>
        end
    end
    ppTbl = varfun(@(x) x{1}, ppTbl);
    rrTbl = varfun(@(x) x{1}, rrTbl);

    % Plot
    [plotTF, plotName] = createFigInd(1.1);
    if plotTF
        numCol = plotSaCharDataCirc(plotName, ksubj, szCharTbl, ppSiTbl, rrSiTbl); % Refer to the created axes using h.f.(thisFcnName)
        % Plot the fit
        for kchar = 1 : numChar
            hax = h.a.(plotName)(ksubj, kchar);
            if kchar <= numSzChar
                pap = ppTbl{1, kchar};
                par = rrTbl{1, kchar};
                [x, y] = pol2cart(-pap - pi/2, par);
                z = -5*ones(size(x));
                h.pa.(plotName)(ksubj, kchar, 2) = patch(hax, x, y, z, 1 - 0.2*(1 - stg.subjColor(ksubj, :)), 'EdgeColor', stg.subjColor(ksubj, :)*stg.subjColorSubjMeanMult);
            end
            % Plotting resultant vector
            resVecWi = 1.5;
            prv = ppCirTbl{1, kchar};
            rrv = rrCirTbl{1, kchar}*hax.XLim(2);
            [x, y] = pol2cart(-prv - pi/2, rrv);
            h.p.(plotName)(ksubj, kchar, 3) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', stg.subjColor(ksubj, :)*0);
            % Arrowhead
            [x(1), y(1)] = pol2cart(-(prv(2)-5/6*pi) - pi/2, rrv(2)/4);
            x(1) = sum(x);
            y(1) = sum(y);
            h.p.(plotName)(ksubj, kchar, 4) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', stg.subjColor(ksubj, :)*0);
            x(1) = 0; y(1) = 0;
            [x(1), y(1)] = pol2cart(-(prv(2)+5/6*pi) - pi/2, rrv(2)/4);
            x(1) = sum(x);
            y(1) = sum(y);
            h.p.(plotName)(ksubj, kchar, 5) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', stg.subjColor(ksubj, :)*0);
            
            % % % % R scaling
            % % % if kchar == 1
            % % %     factor = ceilToEven1sig(max([r, rrCirTbl{1, kchar}]));
            % % %     for kplot = 1 : 2
            % % %         h.p.(plotName)(ksubj, kchar, kplot).XData = factor*h.p.(plotName)(ksubj, kchar, kplot).XData;
            % % %         h.p.(plotName)(ksubj, kchar, kplot).YData = factor*h.p.(plotName)(ksubj, kchar, kplot).YData;
            % % %     end
            % % %     h.pa.(plotName)(ksubj, kchar, 1).XData = factor*h.pa.(plotName)(ksubj, kchar, 1).XData;
            % % %     h.pa.(plotName)(ksubj, kchar, 1).YData = factor*h.pa.(plotName)(ksubj, kchar, 1).YData;
            % % %     h.a.(plotName)(ksubj, kchar).XLim = factor*h.a.(plotName)(ksubj, kchar).XLim;
            % % %     h.a.(plotName)(ksubj, kchar).YLim = factor*h.a.(plotName)(ksubj, kchar).YLim;
            % % % end

            % Text
            xt = 0;
            yt = hax.YLim(1);
            if ~isnan(citP(1, kchar))
                rpString = ['r=', num2str(rrv(2), 1), ',p=', num2str(citP(1, kchar), 1)];
            else
                rpString = ['r=', num2str(rrv(2), 1)];
            end
            text(hax, xt, yt, rpString, 'Interpreter', 'tex', 'FontSize', stg.axFontSize, 'FontWeight', fontWt{kchar},...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
            % ylbl = getYlbl(plotName, szCharTbl.Properties.VariableNames, kchar);
            ylbl = stg.saCharNameShort(kchar);
            yt = hax.YLim(2);
            text(hax, xt, yt, ylbl, 'Interpreter', 'tex', 'FontSize', stg.axFontSize, 'FontWeight', 'normal',...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        end
        if rem(numCol, 2) == 0
            hax = h.a.(plotName)(ksubj, numCol/2);
            x = hax.XLim(2);
            y = hax.YLim(2) + 0.3*diff(hax.YLim);
        else
            hax = h.a.(plotName)(ksubj, ceil(numCol/2));
            x = hax.XLim(1) + diff(hax.XLim)/2;
            y = hax.YLim(2) + 0.3*diff(hax.YLim);
        end
        text(hax, x, y, subjInfo.subjNm, 'Interpreter', 'none', 'FontSize', stg.axFontSize + 1,...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        for kchar = 1 : numChar + 1
            set(h.a.(plotName)(ksubj, kchar), 'Units', 'normalized');
        end
    end
end
function plotSaCharCiFitAllPop(cirTbl, ppTbl, rrTbl, szCharTbl, siCharTbl)
    global stg
    global h
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    str = strings(1, numChar); resM = NaN(1, numChar); resR = NaN(1, numChar); rayP = NaN(1, numChar); omnP = NaN(1, numChar);
    for kchar = 1 : numChar
        p = cirTbl{:, kchar};
        p = p(~isnan(p));
        resM(kchar) = circ_mean(p);
        if ~stg.andrzejak
            resR(kchar) = circ_r(p);
        else
            resR(kchar) = circ_t(p);
        end
        rayP(kchar) = circ_rtest(p);
        omnP(kchar) = circ_otest(p);
        clear p
        xNm = 'Time of day';
        str(kchar) = ...
            [repelem(' ', 1, 13-numel(xNm)), xNm, ' = ', num2str(mod(resM(kchar)/2/pi*24, 24), '%.2g'), ', L=', num2str(resR(kchar), '%.2g'), 10, ...
             '        Ray p = ', num2str(rayP(kchar), '%.2g'), 10, ...
             '       Omni p = ', num2str(omnP(kchar), '%.2g'), 10];
    end
    [~, rayH] = getSignificanceFromP(rayP, stg.Q, stg.multCompCorr);
    % omnH = getSignificanceFromP(omnP, stg.Q, stg.multCompCorr);
    positionCm = [10, 25, stg.singleColumnWidth, min(25, stg.saCharHe*numChar*1.1)];
    [plotTF, plotName] = createFigPos(positionCm);
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp('Circadian statistics (animal-wise):')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, [szCharTbl.Properties.VariableNames, siCharTbl.Properties.VariableNames], kchar);
            disp(['    ', ylbl])
            disp(str(kchar))
        end
    end
    if plotTF
        % Calculate axes positions
        [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlobCi, stg.margCi);
        h.f.(plotName).Units = stg.units;
        
        % Plot the data
        numCol = ceil((numChar+1)^0.6);
        numRow = ceil((numChar+1)/numCol);
        axWiHe = 0.7*min(spWi/numCol, spHe/numRow);
        for kchar = 1 : numChar+1
            h.a.(plotName)(1, kchar) = axes('Units', stg.units, 'Position', [...
                spx(mod(1 - 1, numc) + 1) + mod(kchar-1, numCol)*spWi/numCol + max(spWi/numCol - axWiHe)/2,...
                spy(ceil(1/numc)) + (numRow-ceil(kchar/numCol))*spHe/numRow + 0.0*spHe/numRow,...
                axWiHe, ...
                axWiHe],...
                'NextPlot', 'add', 'Visible', 'off', 'XLimMode', 'manual', 'YLimMode', 'manual');
            hax = h.a.(plotName)(1, kchar);
            
            % % % r1 = NaN(stg.numSubj, numel(rrTbl{1, 1}));
            if kchar ~= numChar+1
                r = NaN(stg.numSubj, size(rrTbl{1, kchar}, 2)); p = r;
                for ksubj = 1 : height(cirTbl)
                    % Subject's circular bar graph
                    r(ksubj, :) = rrTbl{ksubj, kchar};
                    p(ksubj, :) = ppTbl{ksubj, kchar};
                    % [x, y] = pol2cart(-p(ksubj, :) - pi/2, r(ksubj, :));
                    % h.p.(plotName)(ksubj, kchar, 1) = plot(hax, x, y, 'Color', 1 - 0.5*(1 -stg.subjColor(ksubj, :)));
                    % Subject's mean resultant vector
                    pv = [0, cirTbl{ksubj, kchar}];
                    rv = [0, cirTbl{ksubj, kchar + 3*numChar}];
                    % % % rv = [0 1];
                    [x, y] = pol2cart(-pv - pi/2, rv);
                    h.p.(plotName)(ksubj, kchar, 1) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-',...
                        'LineWidth', stg.lnWiFit, 'Color', stg.subjColor(ksubj, :)/2);
                    [x(1), y(1)] = pol2cart(-(pv(2)-5/6*pi) - pi/2, rv(2)/4);
                    x(1) = sum(x);
                    y(1) = sum(y);
                    h.p.(plotName)(1, kchar, 5) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', stg.lnWiFit, 'Color', stg.subjColor(ksubj, :)/2);
                    x(1) = 0; y(1) = 0;
                    [x(1), y(1)] = pol2cart(-(pv(2)+5/6*pi) - pi/2, rv(2)/4);
                    x(1) = sum(x);
                    y(1) = sum(y);
                    h.p.(plotName)(1, kchar, 6) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', stg.lnWiFit, 'Color', stg.subjColor(ksubj, :)/2);
               end
                % Circular bar graph mean
                rm = mean(r, 1, 'omitmissing');
                pm = p(1, :);
                pm(end+1) = pm(1); %#ok<AGROW>
                rm(end+1) = rm(1); %#ok<AGROW>
                [x, y] = pol2cart(-pm - pi/2, rm);
                h.p.(plotName)(1, kchar, 2) = plot(hax, x, y, 'Color', 'k', 'LineWidth', 1.5);
            end
            % Circle
            circleP = (0 : 1 : 360)/360*2*pi;
            circleR = 1;
            % if ~all(isnan(rrCirTbl{:, kchar}))
            %     circleR = ceilToEven1sig(max(rrCirTbl{:, kchar}, [], 'all'));
            % else
            %     circleR = 1;
            % end
            [x, y] = pol2cart(-circleP - pi/2, circleR);
            h.p.(plotName)(1, kchar, 3) = plot(x, y, 'k');
            hax.XLim = max([x, y])*[-1 1];
            hax.YLim = hax.XLim;
            % Night shading
            circleP = (-90 : 1 : 90)/360*2*pi;
            % colo = [linspace(1, 0, 90), 1, linspace(0, 1, 90)]'*[1 1 1];
            colo = [linspace(1, 0.7, 90), 1, linspace(0.7, 1, 90)]'*[1 1 1];
            colo = permute(colo, [1 3 2]);
            [x, y] = pol2cart(-circleP - pi/2, circleR);
            z = -10*ones(size(x));
            h.pa.(plotName)(1, kchar, 1) = patch(x, y, z, colo, 'EdgeColor', 'none');
            if kchar ~= numChar+1
                % Mean resultant vector
                resVecWi = 1.5;
                prv = [0, resM(kchar)];
                rrv = [0, resR(kchar)];
                [x, y] = pol2cart(-prv - pi/2, rrv);
                h.p.(plotName)(1, kchar, 4) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', 'k');
                % Arrowhead
                [x(1), y(1)] = pol2cart(-(prv(2)-5/6*pi) - pi/2, rrv(2)/4);
                x(1) = sum(x);
                y(1) = sum(y);
                h.p.(plotName)(1, kchar, 5) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', 'k');
                x(1) = 0; y(1) = 0;
                [x(1), y(1)] = pol2cart(-(prv(2)+5/6*pi) - pi/2, rrv(2)/4);
                x(1) = sum(x);
                y(1) = sum(y);
                h.p.(plotName)(1, kchar, 6) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', 'k');
                x = 0;
                y = hax.YLim(2);
                % text(hax, x, y, num2str(y), 'Interpreter', 'tex', 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
                % ylbl = getYlbl(plotName, szCharTbl.Properties.VariableNames, kchar);
                ylbl = stg.saCharNameShort(kchar);
                text(hax, x, y, ylbl, 'Interpreter', 'tex', 'FontSize', stg.axFontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
                strP = string(['  p=', num2str(rayP(1, kchar), 2)]);
                % ht = text(hax, xrv(2)/2, yrv(2)/2, strP, 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment', 'middle', 'BackgroundColor','w');
                ht = text(hax, 0, -1, strP, 'FontSize', stg.axFontSize, 'HorizontalAlignment','center', 'VerticalAlignment', 'top', 'BackgroundColor','none');
                if rayH(1, kchar)
                    ht.FontWeight = "bold";
                else
                    ht.FontWeight = "normal";
                end
            end
            % Plot explanatory circle
            if kchar == numChar+1
                plotExplanatoryCircle(hax)
            end
        end
        % Title
        if rem(numCol, 2) == 0
            hax = h.a.(plotName)(1, numCol/2);
            x = hax.XLim(2);
            y = hax.YLim(2) + 0.3*diff(hax.YLim);
        else
            hax = h.a.(plotName)(1, ceil(numCol/2));
            x = hax.XLim(1) + diff(hax.XLim)/2;
            y = hax.YLim(2) + 0.3*diff(hax.YLim);
        end
        text(hax, x, y, 'All mice', 'Interpreter', 'none', 'FontSize', stg.axFontSize + 1,...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        for kchar = 1 : numChar + 1
            set(h.a.(plotName)(1, kchar), 'Units', 'normalized');
        end
    end
    % disp(['resR = ', resR])
end
% Plot signal characteristics analyses
function plotSiChar(subjInfo, szCharTbl, siCharTbl, ksubj)
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSiCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl); % Refer to the created axes using h.f.(thisFcnName)
        formatAxesInd(plotName, ksubj, subjInfo, siCharTbl.Properties.VariableNames, 'Age (days)', "analysisPeriod")
    end
end
function [paxTbl, psdTbl] = plotSiCharPsd(subjInfo, siCharTbl, ksubj)
    global stg
    global h
    charToPlot = ["sz", stg.siCharToPlot];
    numChar = numel(charToPlot);
    [paxTbl, psdTbl, psdciTbl] = psdSiChar(siCharTbl);
    [plotTF, plotName] = createFigInd(1.0);
    if plotTF
        % Calculate axes positions
        % % % % % % % % margGlob = [2.2 0.6 0 0]; % Left, bottom, right, top
        % % % % % % % % marg = [0.7 0.7 0.5 0.5]; % Left, bottom, right, top
        [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlob, stg.marg);
        
        % Plot the data
        for kchar = 1 : numChar
            h.a.(plotName)(ksubj, kchar) = axes('Units', stg.units, 'Position', ...
                [spx(mod(ksubj - 1, numc) + 1), spy(ceil(ksubj/numc)) - (kchar - numChar)*spHe/numChar, spWi, spHe/numChar]);
            x = paxTbl{1, kchar}';
            y = psdTbl{1, kchar}';
            e = psdciTbl{1, kchar}';
            h.pa.(plotName)(ksubj, kchar) = patch([x; flipud(x)], 10*log10([y - e; flipud(y) + flipud(e)]), ...
               1 - (1 - stg.subjColor(ksubj, :))/8, 'EdgeColor', 'none');
            hold on
            h.p.(plotName)(ksubj, kchar, 1) = plot(x, 10*log10(y), 'Marker', 'x', 'MarkerSize', 3, 'LineWidth', 0.5, 'Color', stg.subjColor(ksubj, :));
            ymi = [min(10*log10(y)), min(10*log10([y - e; flipud(y) + flipud(e)]))];
            yma = [max(10*log10(y)), max(10*log10([y - e; flipud(y) + flipud(e)]))];
            yl = getYLimYTick([ymi, yma])';
            h.p.(plotName)(ksubj, kchar, 2) = plot([1; 1], yl, ':k', 'LineWidth', 0.5);
            % h.a.(plotName)(ksubj, kchar).XScale = 'log';
            clear x y e
        end
% % % % % asdf = siCharTbl.Properties.VariableNames
        formatAxesInd(plotName, ksubj, subjInfo, {'sz', 'ied'}, 'Period (days)', "none", char([10, 'PSD (dB)']))
    end
end
function plotSiCharPsdAllPop(paxTbl, psdTbl, siCharTbl)
    global stg
    global h
    [numSubj, numChar] = size(psdTbl);
    positionCm = [20, 10, stg.singleColumnWidth, min(25, stg.siCharHe*numChar*1.0 + 1.3)];
    [plotTF, plotName] = createFigPos(positionCm);
    if plotTF
        if ~isfield(h.a, plotName)
            createAxesAllPopStack(plotName)
        end
        for kchar = 1 : numChar
            for ksubj = 1 : numSubj
                hax = h.a.(plotName)(kchar); % Axes handle (I know there are no longer handles to the graphic objects but I still use the name hax)
                pax = paxTbl{ksubj, kchar};
                x = pax;
                psd = psdTbl{ksubj, kchar}; % Power spectral density
                [~, pax1sub] = min(abs(pax - 1));
                yNormalizationMatrix = (psd(pax1sub)*ones(1, size(psd, 2)));
                y = 10*log10(psd./yNormalizationMatrix)';
                plot(hax, x, y, 'LineWidth', 0.5, 'Color', stg.subjColor(ksubj, :));
            end
            psd = mean(psdTbl{:, kchar}, 1, 'omitmissing');
            yNormalizationMatrix = (psd(pax1sub)*ones(1, size(psd, 2)));
            y = 10*log10(psd./yNormalizationMatrix)';
            plot(hax, x, y, 'k', 'LineWidth', 2)
        end
        
        % Axes formatting
        xmi = NaN(1, numChar); xma = NaN(1, numChar); % Initialization
        for kchar = 1 : numChar
            xmi(kchar) = h.a.(plotName)(1, kchar).XLim(1);
            xma(kchar) = h.a.(plotName)(1, kchar).XLim(2);
        end
        xmi = min(xmi);
        xma = max(xma);
        for kchar = 1 : numChar
            hax = h.a.(plotName)(1, kchar);
            hax.XLim = [xmi xma];
            if kchar ~= numChar
                hax.XTickLabel = [];
            end
            hax.TickLength = [0.01 0.01];
            if kchar == numChar
                xlabel('Period (days)')
            end
            % y-axis
            numChildren = numel(hax.Children);
            ymi = NaN(1, numChildren); yma = NaN(1, numChildren); % Initialization
            for kchild = 1 : numChildren
                ymi(kchild) = min(hax.Children(kchild).YData(:));
                yma(kchild) = max(hax.Children(kchild).YData(:));
            end
            ymi = min(ymi);
            yma = max(yma);
            % [yl, yt] = getYLimYTick([ymi yma], stg.siCharYLim(kchar, :)); % Input could be a signal but getYLimYTick takes the min and max of it anyway
            [yl, yt] = getYLimYTick([ymi yma]); % Input could be a signal but getYLimYTick takes the min and max of it anyway
            hax.YLim = yl;
            hax.YTick = yt;
            if mod(ksubj - 1, stg.sbNCol) == 0
                ylbl = [stg.siCharYLbl{strcmp(siCharTbl.Properties.VariableNames, stg.siCharToPlot(kchar))}, 10, ' PSD (dB)'];
                hax.YLabel.String = ylbl;
                if numChar > 1
                    hax.YLabel.Rotation = -35;
                    hax.YLabel.Position(1) = hax.YLabel.Position(1) - (hax.XLim(1) - hax.YLabel.Position(1))/3;
                    hax.YLabel.HorizontalAlignment = 'right';
                    hax.YLabel.VerticalAlignment = 'middle';
                end
            end
            hax.Layer = 'top';
            hax.Box = stg.box;
            hax.FontSize = stg.axFontSize;
            hax.Units = 'normalized';
        end
    end
end
function fitTbl = plotSiCharWhFit(subjInfo, szCharTbl, siCharTbl, ksubj)
    global stg
    global h
    numChar = numel(stg.siCharToPlot);
    [fitTbl, ~, ~, xxFitTbl, yyFitTbl] = fitSiCharWh(subjInfo, siCharTbl);
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSiCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl); % Refer to the created axes using h.f.(thisFcnName)
        % Plot the trends
        for kchar = 1 : numChar
            hax = h.a.(plotName)(ksubj, kchar);
            x = xxFitTbl{1, kchar};
            y = yyFitTbl{1, kchar};
            plot(hax, x, y, 'Marker', 'none', 'LineStyle', '--', 'LineWidth', 2, 'Color', 1 - 0.5*(1 - stg.subjColor(ksubj, :)));
            % % plot(hax, x, y, 'Marker', 'none', 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k'); % SalyOlomouc
        end
        formatAxesInd(plotName, ksubj, subjInfo, siCharTbl.Properties.VariableNames, 'Age (days)', "analysisPeriod")
    end
end
function plotSiCharWhFitAllPop(fitTbl, siCharTbl)
    global stg
    global h
    numChar = numel(stg.siCharToPlot);
    numSubj = stg.numSubj;
    positionCm = [20, 25, stg.singleColumnWidth, min(25, stg.siCharHe*numChar)];
    [plotTF, plotName] = createFigPos(positionCm);
    str = strings(1, numChar); medFit = NaN(1, numChar); sigRanP = NaN(1, numChar);
    for kchar = 1 : numChar
        [str(kchar), ~, ~, medFit(kchar), ~, sigRanP(kchar)] = getBasicStats(fitTbl.(stg.siCharToPlot(kchar) + "fitSlope"), 'Slope');
    end
    [sigRanP, sigRanH] = getSignificanceFromP(sigRanP, stg.Q, stg.multCompCorr);
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp('Whole recording:')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, siCharTbl.Properties.VariableNames, kchar);
            disp(['    ', ylbl])
            disp(str(kchar))
        end
    end
    if plotTF
        if ~isfield(h.a, plotName)
            createAxesAllPopStack(plotName)
        end
        strP = strings(1, numChar);
        for kchar = 1 : numChar
            hax = h.a.(plotName)(kchar);
            for ksubj = 1 : numSubj
                fitSlope = fitTbl{ksubj, kchar + 0*numChar};
                x = [0 1];
                y = x*fitSlope;
                h.p.(plotName)(kchar) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', stg.lnWiFit, 'Color', stg.subjColor(ksubj, :));
            end
            % Medians
            % % % % % % % xt = mean(x);
            y = x*medFit(kchar);
            plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '--', 'Color', 'k', 'Marker', 'none');
            switch sign(medFit(kchar))
                case -1
                    signMed = '-';
                case 0
                    signMed = '0';
                case 1
                    signMed = '+';
            end
            strP(kchar) = string([' p=', num2str(sigRanP(1, kchar), 2), '(', signMed, ')']);
        end
        formatAxesAllPop(plotName, siCharTbl.Properties.VariableNames, 'Normalized time', "none", ["any", "any"])
        for kchar = 1 : numChar
            hax = h.a.(plotName)(kchar);
            % xt = mean(hax.XLim);
            xt = 0;
            % yt = mean(hax.YLim);
            if abs(hax.YLim(1)) > hax.YLim(2)
                yt = hax.YLim(1);
                ht = text(hax, xt, yt, strP(kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','left', 'VerticalAlignment','bottom', 'BackgroundColor','none');
            else                
                yt = hax.YLim(2);
                ht = text(hax, xt, yt, strP(kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','left', 'VerticalAlignment','top', 'BackgroundColor','none');
            end
            if sigRanH(1, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
        end
    end
end
function [xBe, myyBeTbl, xDu, myyDuTbl, xAf, myyAfTbl] = plotSiCharCl(subjInfo, siCharTbl, cl, ksubj)
    % cl ... structure with all clusters of given subject
    global stg
    global h
    numPoints = 4; % Each cluster have different length. Number of points to which we will interpolate the characteristic.
    numChar = numel(stg.siCharToPlot);
    fitPosition = ["before", "during", "after"];
    [plotTF, plotName] = createFigInd(1);
    myyInterpTblCell = cell(1, numel(fitPosition)); xCell = cell(1, numel(fitPosition));
    for kfp = 1 : numel(fitPosition) % k fit position
        [~, xxTbl, yyTbl, ~, ~] = fitSiCharCl(subjInfo, siCharTbl, cl, fitPosition(kfp));
        xxInterpTbl = table; yyInterpTbl = table;
        varNames = xxTbl.Properties.VariableNames;
        if isempty(xxTbl)
            for kchar = 1 : numChar
                xxInterpTbl{1, kchar} = NaN(1, numPoints);
                yyInterpTbl{1, kchar} = NaN(1, numPoints);
            end
            plotColor = stg.subjColor(ksubj, :);
        else
            for kcl = 1 : height(xxTbl)
                for kchar = 1 : numChar
                    if numel(xxTbl{kcl, kchar}{1}) < 2
                        xxInterpTbl{kcl, kchar} = NaN(1, numPoints);
                        yyInterpTbl{kcl, kchar} = NaN(1, numPoints);
                    else
                        xxInterpTbl{kcl, kchar} = linspace(xxTbl{kcl, kchar}{1}(1), xxTbl{kcl, kchar}{1}(end), numPoints);
                        yyInterpTbl{kcl, kchar} = interp1(xxTbl{kcl, kchar}{1}, yyTbl{kcl, kchar}{1}, xxInterpTbl{kcl, kchar}, "linear");
                    end
                end
            end
            plotColor = ones(size(yyTbl, 1), 1)*stg.subjColor(ksubj, :);
        end
        xxInterpTbl.Properties.VariableNames = varNames;
        yyInterpTbl.Properties.VariableNames = varNames;
        myyInterpTblCell{kfp} = stg.withinSubjectStat(yyInterpTbl, 'omitmissing');
        switch fitPosition(kfp)
            case "before"
                xCell{kfp} = linspace(0, 0.33, numPoints);
            case "during"
                xCell{kfp} = linspace(0.34, 0.66, numPoints);
            case "after"
                xCell{kfp} = linspace(0.67, 1, numPoints);
        end
        if plotTF
            numCl = plotPeriEventData(plotName, ksubj, xCell{kfp}, yyInterpTbl, myyInterpTblCell{kfp}, plotColor); % Refer to the created axes using h.a.(plotName)
            h.t.(plotName)(ksubj, kchar, kfp, 1) = text(mean(xCell{kfp}), h.a.(plotName)(ksubj, kchar).YLim(2), ['n=', num2str(numCl(kchar))],...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', stg.statFontSize); % A four dimensional matrix. Cool, isn't it.
            h.t.(plotName)(ksubj, kchar, kfp, 2) = text(mean(xCell{kfp}), 0, fitPosition(kfp),...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', stg.statFontSize);
            h.a.(plotName)(ksubj, kchar).XTick = [];
        end
    end
    if plotTF
        formatAxesInd(plotName, ksubj, subjInfo, siCharTbl.Properties.VariableNames, 'Normalized time', "none")
        for kchar = 1 : numChar
            yl = h.a.(plotName)(ksubj, kchar).YLim;
            h.pa.(plotName)(ksubj, kchar, 1) = patch(h.a.(plotName)(ksubj, kchar), [0 0.33 0.33 0], [yl(1) yl(1) yl(2) yl(2)],...
                1 - 0.15*(1-stg.fitColor(1, :)), 'ZData', [-1 -1 -1 -1], 'EdgeColor', 'none');
            h.pa.(plotName)(ksubj, kchar, 2) = patch(h.a.(plotName)(ksubj, kchar), [0.34 0.66 0.66 0.34], [yl(1) yl(1) yl(2) yl(2)],...
                1 - 0.15*(1-stg.fitColor(2, :)), 'ZData', [-1 -1 -1 -1], 'EdgeColor', 'none');
            h.pa.(plotName)(ksubj, kchar, 3) = patch(h.a.(plotName)(ksubj, kchar), [0.67 1 1 0.67], [yl(1) yl(1) yl(2) yl(2)],...
                1 - 0.15*(1-stg.fitColor(3, :)), 'ZData', [-1 -1 -1 -1], 'EdgeColor', 'none');
            for kfp = 1 : numel(fitPosition)
                h.t.(plotName)(ksubj, kchar, kfp, 1).Position(2) = h.a.(plotName)(ksubj, kchar).YLim(2);
            end
        end
    end
    xBe = xCell{1};
    xDu = xCell{2};
    xAf = xCell{3};
    myyBeTbl = myyInterpTblCell{1};
    myyDuTbl = myyInterpTblCell{2};
    myyAfTbl = myyInterpTblCell{3};
end
function plotSiCharClAllPop(xBe, myyBeTbl, xDu, myyDuTbl, xAf, myyAfTbl, siCharTbl)
    global stg
    global h
    numChar = numel(stg.siCharToPlot);
    positionCm = [20, 10, stg.singleColumnWidth, min(25, stg.siCharHe*numChar)];
    [plotTF, plotName] = createFigPos(positionCm);
    if contains(plotName, 'Pop')
        statFunction = stg.acrossSubjectStat;
    else
        statFunction = stg.withinSubjectStat;
    end
    if plotTF
        if ~isfield(h.a, plotName)
            createAxesAllPopStack(plotName)
        end
        mmyyBeTbl = statFunction(myyBeTbl, [], 'omitmissing');
        plotPeriEventData(plotName, 1, xBe, myyBeTbl, mmyyBeTbl, stg.subjColor);
        mmyyDuTbl = statFunction(myyDuTbl, [], 'omitmissing');
        plotPeriEventData(plotName, 1, xDu, myyDuTbl, mmyyDuTbl, stg.subjColor);
        mmyyAfTbl = statFunction(myyAfTbl, [], 'omitmissing');
        plotPeriEventData(plotName, 1, xAf, myyAfTbl, mmyyAfTbl, stg.subjColor);
        formatAxesAllPop(plotName, siCharTbl.Properties.VariableNames, 'Normalized time', "none", ["any", "any"])
        for kchar = 1 : numChar
            yl = h.a.(plotName)(1, kchar).YLim;
            h.pa.(plotName)(1, kchar, 1) = patch(h.a.(plotName)(1, kchar), [0 0.33 0.33 0], [yl(1) yl(1) yl(2) yl(2)],...
                1 - 0.15*(1-stg.fitColor(1, :)), 'ZData', [-1 -1 -1 -1], 'EdgeColor', 'none');
            h.pa.(plotName)(1, kchar, 2) = patch(h.a.(plotName)(1, kchar), [0.34 0.66 0.66 0.34], [yl(1) yl(1) yl(2) yl(2)],...
                1 - 0.15*(1-stg.fitColor(2, :)), 'ZData', [-1 -1 -1 -1], 'EdgeColor', 'none');
            h.pa.(plotName)(1, kchar, 3) = patch(h.a.(plotName)(1, kchar), [0.67 1 1 0.67], [yl(1) yl(1) yl(2) yl(2)],...
                1 - 0.15*(1-stg.fitColor(3, :)), 'ZData', [-1 -1 -1 -1], 'EdgeColor', 'none');
        end
    end
end
function [mFitTblBe, mFitTblDu, mFitTblAf] = plotSiCharClFit(subjInfo, szCharTbl, siCharTbl, cl, ksubj)
    % cl ... structure with all clusters of given subject
    global stg
    global h
    numChar = numel(stg.siCharToPlot);
    fitPosition = ["before", "during", "after"];
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSiCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl); % Refer to the created axes using h.a.(thisFcnName)
    end
    mFitTblCell = cell(1, numel(fitPosition));
    for kfp = 1 : numel(fitPosition)
        [fitTbl, ~, ~, xxFitTbl, yyFitTbl] = fitSiCharCl(subjInfo, siCharTbl, cl, fitPosition(kfp));
        mFitTblCell{kfp} = stg.withinSubjectStat(fitTbl, "omitmissing");
        if plotTF
            for kchar = 1 : numChar
                if isempty(xxFitTbl)
                    continue
                end
                for kcl = 1 : size(xxFitTbl, 1)
                    % if clust(kcl).nested
                    %     continue
                    % end
                    hax = h.a.(plotName)(ksubj, kchar);
                    plot(hax, xxFitTbl{kcl, kchar}, yyFitTbl{kcl, kchar}, ...
                       'Marker', 'none', 'LineWidth', 2, 'LineStyle', '-', 'Color', stg.fitColor(kfp, :), 'Marker', 'none');
                end
            end
        end
    end
    if plotTF
        formatAxesInd(plotName, ksubj, subjInfo, siCharTbl.Properties.VariableNames, 'Age (days)', "analysisPeriod")
    end
    mFitTblBe = mFitTblCell{1};
    mFitTblDu = mFitTblCell{2};
    mFitTblAf = mFitTblCell{3};
end
function plotSiCharClFitAllPop(fitTblBe, fitTblDu, fitTblAf, siCharTbl)
    global stg
    global h
    numChar = numel(stg.siCharToPlot);
    numSubj = stg.numSubj;
    str = strings(3, numChar); medFit = NaN(3, numChar); sigRanP = NaN(3, numChar);
    for kchar = 1 : numChar
        [str(1, kchar), ~, ~, medFit(1, kchar), ~, sigRanP(1, kchar)] = getBasicStats(fitTblBe.(stg.siCharToPlot(kchar) + "fitSlope"), 'Slope');
        [str(2, kchar), ~, ~, medFit(2, kchar), ~, sigRanP(2, kchar)] = getBasicStats(fitTblDu.(stg.siCharToPlot(kchar) + "fitSlope"), 'Slope');
        [str(3, kchar), ~, ~, medFit(3, kchar), ~, sigRanP(3, kchar)] = getBasicStats(fitTblAf.(stg.siCharToPlot(kchar) + "fitSlope"), 'Slope');
    end
    [sigRanP, sigRanH] = getSignificanceFromP(sigRanP, stg.Q, stg.multCompCorr);
    positionCm = [20, 10, stg.singleColumnWidth, min(25, stg.siCharHe*numChar)];
    [plotTF, plotName] = createFigPos(positionCm);
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp('Before cluster:')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, siCharTbl.Properties.VariableNames, kchar);
            disp(['    ', ylbl])
            disp(str(1, kchar))
        end
        disp('During cluster:')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, siCharTbl.Properties.VariableNames, kchar);
            disp(['    ', ylbl])
            disp(str(2, kchar))
        end
        disp('After cluster:')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, siCharTbl.Properties.VariableNames, kchar);
            disp(['    ', ylbl])
            disp(str(3, kchar))
        end
    end
    if plotTF
        if ~isfield(h.a, plotName)
            createAxesAllPopStack(plotName)
        end
        strP = strings(3, numChar); xt = NaN(3, 1);
        for kchar = 1 : numChar
            hax = h.a.(plotName)(kchar);
            for ksubj = 1 : numSubj
                % Before cluster
                fitSlope = fitTblBe{ksubj, kchar + 0*numChar};
                x = [0 0.33];
                y = [-1 0]*fitSlope;
                plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFit, 'LineStyle', '-', 'Color', stg.fitColor(1, :), 'Marker', 'none');
                % During cluster
                fitSlope = fitTblDu{ksubj, kchar + 0*numChar};
                x = [0.34 0.66];
                y = [0 1]*fitSlope;
                plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFit, 'LineStyle', '-', 'Color', stg.fitColor(2, :), 'Marker', 'none');
                % After cluster
                fitSlope = fitTblAf{ksubj, kchar + 0*numChar};
                x = [0.67 1];
                y = [0 1]*fitSlope;
                plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFit, 'LineStyle', '-', 'Color', stg.fitColor(3, :), 'Marker', 'none');
            end
            % Medians
            % Before cluster
            x = [0 0.33];
            xt(1) = mean(x);
            y = [-1 0]*medFit(1, kchar);
            plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            switch sign(medFit(1, kchar))
                case -1
                    signMed = '-';
                case 0
                    signMed = '0';
                case 1
                    signMed = '+';
            end
            strP(1, kchar) = string(['p=', num2str(sigRanP(1, kchar), 2), '(', signMed, ')']);
            % During cluster
            x = [0.34 0.66];
            xt(2) = mean(x);
            y = [0 1]*medFit(2, kchar);
            plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            switch sign(medFit(2, kchar))
                case -1
                    signMed = '-';
                case 0
                    signMed = '0';
                case 1
                    signMed = '+';
            end
            strP(2, kchar) = string(['p=', num2str(sigRanP(2, kchar), 2), '(', signMed, ')']);
            % After cluster
            x = [0.67 1];
            xt(3) = mean(x);
            y = [0 1]*medFit(3, kchar);
            switch sign(medFit(3, kchar))
                case -1
                    signMed = '-';
                case 0
                    signMed = '0';
                case 1
                    signMed = '+';
            end
            plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            strP(3, kchar) = string(['p=', num2str(sigRanP(3, kchar), 2), '(', signMed, ')']);
        end
        formatAxesAllPop(plotName, siCharTbl.Properties.VariableNames, 'Normalized time', "none", ["any", "any"])
        for kchar = 1 : numChar
            hax = h.a.(plotName)(kchar);
            % % % % yt = mean(hax.YLim);
            % % % % ht = text(hax, xt(1), yt, strP(1, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            if abs(hax.YLim(1)) > hax.YLim(2)
                yt = hax.YLim(1);
                ht = text(hax, xt(1), yt, strP(1, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'BackgroundColor','none');
            else                
                yt = hax.YLim(2);
                ht = text(hax, xt(1), yt, strP(1, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            end
            if sigRanH(1, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
            % % % % ht = text(hax, xt(2), yt, strP(2, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            if abs(hax.YLim(1)) > hax.YLim(2)
                yt = hax.YLim(1);
                ht = text(hax, xt(2), yt, strP(2, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'BackgroundColor','none');
            else                
                yt = hax.YLim(2);
                ht = text(hax, xt(2), yt, strP(2, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            end
            if sigRanH(2, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
            % % % % ht = text(hax, xt(3), yt, strP(3, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            if abs(hax.YLim(1)) > hax.YLim(2)
                yt = hax.YLim(1);
                ht = text(hax, xt(3), yt, strP(3, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'BackgroundColor','none');
            else                
                yt = hax.YLim(2);
                ht = text(hax, xt(3), yt, strP(3, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            end
            if sigRanH(3, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
        end
        hax.XTick = [];
        hax.XLabel = [];
        text(0 + 0.33/2, hax.YLim(1) - 0.2*diff(hax.YLim), 'Pre-cluster', 'FontSize',  8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
        text(0.34 + 0.33/2, hax.YLim(1) - 0.2*diff(hax.YLim), 'Cluster', 'FontSize',  8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
        text(0.67 + 0.33/2, hax.YLim(1) - 0.2*diff(hax.YLim), 'Post-cluster', 'FontSize',  8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
    end
end
function [cirTbl, ppTbl, rrTbl] = plotSiCharCiFit(subjInfo, szCharTbl, siCharTbl, ksubj)
    global stg
    global h
    numChar = numel(stg.siCharToPlot);
    [cirTbl, ~, ppTbl, rrTbl, ppCirTbl, rrCirTbl] = fitSiCharCi(siCharTbl); % Compute the circular statistics
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSiCharDataCirc(plotName, ksubj, szCharTbl, ppTbl, rrTbl); % Refer to the created axes using h.f.(thisFcnName)
        % Plot the trend
        for kchar = 1 : numChar
            hax = h.a.(plotName)(ksubj, kchar);
            % Plotting resultant vector
            resVecWi = 1.5;
            prv = ppCirTbl{1, kchar};
            rrv = rrCirTbl{1, kchar};
            % rrv = [0, hax.XLim(2)];
            [x, y] = pol2cart(-prv - pi/2, rrv);
            h.p.(plotName)(ksubj, kchar, 3) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', stg.subjColor(ksubj, :)*0);
            % Arrowhead
            [x(1), y(1)] = pol2cart(-(prv(2)-5/6*pi) - pi/2, rrv(2)/4);
            x(1) = sum(x);
            y(1) = sum(y);
            h.p.(plotName)(ksubj, kchar, 4) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', stg.subjColor(ksubj, :)*0);
            x(1) = 0; y(1) = 0;
            [x(1), y(1)] = pol2cart(-(prv(2)+5/6*pi) - pi/2, rrv(2)/4);
            x(1) = sum(x);
            y(1) = sum(y);
            h.p.(plotName)(ksubj, kchar, 5) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', stg.subjColor(ksubj, :)*0);
            % Text
            x = 0;
            y = hax.YLim(1);
            text(hax, x, y, ['r=', num2str(rrv(2), 1)], 'Interpreter', 'tex', 'FontWeight', 'normal', 'FontSize', stg.axFontSize,...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
            % ylbl = getYlbl(plotName, szCharTbl.Properties.VariableNames, kchar);
            ylbl = stg.siCharNameShort(kchar);
            y = hax.YLim(2);
            text(hax, x, y, ylbl, 'Interpreter', 'tex', 'FontWeight', 'normal', 'FontSize', stg.axFontSize,...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        end
        % % % % % % % % margGlob = [0 0 0 0]; marg = [0.1 0.5 0.1 0.5];
        [~, ~, spWi, ~, ~, ~] = getSubplotXYWH(plotName, stg.margGlobCi, stg.margCi);
        hax = h.a.(plotName)(ksubj, 1);
        x = hax.XLim(1) + diff(hax.XLim)*(0.9*spWi/hax.Position(3))/2;
        y = hax.YLim(2) + 0.2*diff(hax.YLim);
        text(hax, x, y, subjInfo.subjNm, 'Interpreter', 'none', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        for kchar = 1 : numChar
            set(h.a.(plotName)(ksubj, kchar), 'Units', 'normalized');
        end
    end
end
function plotSiCharCiFitAllPop(cirTbl, ppTbl, rrTbl, siCharTbl)
    global stg
    global h
    numChar = numel(stg.siCharToPlot);
    str = strings(1, numChar); resM = NaN(1, numChar); resR = NaN(1, numChar); rayP = NaN(1, numChar); omnP = NaN(1, numChar);
    for kchar = 1 : numChar
        resM(kchar) = circ_mean(cirTbl{:, kchar});
        % % % % % % % resM(kchar) = circ_mean(cirTbl{:, kchar});
        if ~stg.andrzejak
            resR(kchar) = circ_r(cirTbl{:, kchar});
        else
            resR(kchar) = circ_t(cirTbl{:, kchar});
        end
        rayP(kchar) = circ_rtest(cirTbl{:, kchar});
        omnP(kchar) = circ_otest(cirTbl{:, kchar});
        xNm = 'Hour of day';
        str(kchar) = ...
            [repelem('       ', 1, 13-numel(xNm)), xNm, ' = ', num2str(mod(resM(kchar)/2/pi*24, 24), '%.2g'), ', L=', num2str(resR(kchar), '%.2g'), 10, ...
             '              Ray p = ', num2str(rayP(kchar), '%.2g'), 10, ...
             '             Omni p = ', num2str(omnP(kchar), '%.2g'), 10];
    end
    [~, rayH] = getSignificanceFromP(rayP, stg.Q, stg.multCompCorr);
    % [~, omnH] = getSignificanceFromP(omnP, stg.Q, stg.multCompCorr);
    positionCm = [10, 25, stg.singleColumnWidth, min(25, stg.siCharHe*numChar)];
    [plotTF, plotName] = createFigPos(positionCm);
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp('Circadian statistics (animal-wise):')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, siCharTbl.Properties.VariableNames, kchar);
            disp(['    ', ylbl])
            disp(str(kchar))
        end
    end
    if plotTF
        % Calculate axes positions
        % % % % % % % % margGlob = [0 0 0 0]; % Left, bottom, right, top
        % % % % % % % % marg = [0.1 0.5 0.1 0.5];
        [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlobCi, stg.margCi);
        h.f.(plotName).Units = stg.units;
        % Plot the data
        numCol = ceil(numChar^0.7);
        numRow = ceil(numChar/numCol);
        axWiHe = 0.8*min(spWi/numCol, spHe/numRow);
        for kchar = 1 : numChar
            h.a.(plotName)(1, kchar) = axes('Units', stg.units, 'Position', [...
                spx(mod(1 - 1, numc) + 1) + mod(kchar-1, numCol)*spWi/numCol + 0.1*spWi/numCol,...
                spy(ceil(1/numc)) + (numRow-ceil(kchar/numCol))*spHe/numRow + 0.0*spHe/numRow,...
                axWiHe, ...
                axWiHe],...
                'NextPlot', 'add', 'Visible', 'off', 'XLimMode', 'manual', 'YLimMode', 'manual');
            hax = h.a.(plotName)(1, kchar);
            for ksubj = 1 : height(cirTbl)
                % Plot subject's mean circadian siChar
                p = [ppTbl{ksubj, kchar}, ppTbl{ksubj, kchar}(1)];
                r = [rrTbl{ksubj, kchar}, rrTbl{ksubj, kchar}(1)];
                [x, y] = pol2cart(-p - pi/2, r);
                h.p.(plotName)(ksubj, kchar, 1) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-',...
                    'LineWidth', 0.5*stg.lnWiFit, 'Color', 1 - 0.5*(1 - stg.subjColor(ksubj, :)));
                % Plot subject's resultant vector
                p = [0, cirTbl{ksubj, kchar}];
                r = [0 1];
                [x, y] = pol2cart(-p - pi/2, r);
                h.p.(plotName)(ksubj, kchar, 2) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-',...
                    'LineWidth', 1*stg.lnWiFit, 'Color', stg.subjColor(ksubj, :));
            end
            % Circle
            circleP = (0 : 1 : 360)/360*2*pi;
            circleR = 1;
            % if ~all(isnan(rrCirTbl{:, kchar}))
            %     circleR = ceilToEven1sig(max(rrCirTbl{:, kchar}, [], 'all'));
            % else
            %     circleR = 1;
            % end
            [x, y] = pol2cart(-circleP - pi/2, circleR);
            h.p.(plotName)(1, kchar, 3) = plot(x, y, 'k');
            hax.XLim = max([x, y])*[-1 1];
            hax.YLim = hax.XLim;
            % Night shading
            circleP = (-90 : 1 : 90)/360*2*pi;
            % colo = [linspace(1, 0, 90), 1, linspace(0, 1, 90)]'*[1 1 1];
            colo = [linspace(1, 0.7, 90), 1, linspace(0.7, 1, 90)]'*[1 1 1];
            colo = permute(colo, [1 3 2]);
            [x, y] = pol2cart(-circleP - pi/2, circleR);
            z = -10*ones(size(x));
            h.pa.(plotName)(1, kchar, 1) = patch(x, y, z, colo, 'EdgeColor', 'none');
            % Mean circadian siChar over subjects
            p = mean([ppTbl{1, kchar}, ppTbl{1, kchar}(1)], 1);
            r = mean([rrTbl{:, kchar}, rrTbl{:, kchar}(:, 1)], 1);
            [x, y] = pol2cart(-p - pi/2, r);
            h.p.(plotName)(ksubj, kchar, 1) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-',...
                'LineWidth', 1*stg.lnWiFitMean, 'Color', 'k');
            % Mean resultant vector
            resVecWi = 1.5;
            prv = [0, resM(kchar)];
            rrv = [0, resR(kchar)];
            [x, y] = pol2cart(-prv - pi/2, rrv);
            h.p.(plotName)(1, kchar, 4) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', 'k');
            % Arrowhead
            [x(1), y(1)] = pol2cart(-(prv(2)-5/6*pi) - pi/2, rrv(2)/4);
            x(1) = sum(x);
            y(1) = sum(y);
            h.p.(plotName)(1, kchar, 6) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', 'k');
            x(1) = 0; y(1) = 0;
            [x(1), y(1)] = pol2cart(-(prv(2)+5/6*pi) - pi/2, rrv(2)/4);
            x(1) = sum(x);
            y(1) = sum(y);
            h.p.(plotName)(1, kchar, 6) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', 'k');
            % Text
            x = 0;
            y = hax.YLim(2);
            % text(hax, x, y, num2str(y), 'Interpreter', 'tex', 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
            % ylbl = getYlbl(plotName, siCharTbl.Properties.VariableNames, kchar);
            ylbl = stg.siCharNameShort(kchar);
            text(hax, x, y, ylbl, 'Interpreter', 'tex', 'FontSize', stg.axFontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
            strP = string(['p=', num2str(rayP(1, kchar), 2)]);
            % ht = text(hax, xrv(2)/2, yrv(2)/2, strP, 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment', 'middle', 'BackgroundColor','w');
            ht = text(hax, 0, -1, strP, 'FontSize', stg.axFontSize, 'HorizontalAlignment','center', 'VerticalAlignment', 'top', 'BackgroundColor','none');
            if rayH(1, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
        end
    end
end
function [diffBe, diffAf] = plotSiCharSzBeAfVsOther(subjInfo, szCharTbl, siCharTbl, ksubj)
    global stg
    global h
    fitPosition = ["before", "after"];
    [plotTF, plotName] = createFigInd(1);
    myyTblCell = cell(1, numel(fitPosition)); xpCell = cell(1, numel(fitPosition)); numSz = cell(1, numel(fitPosition));
    for kfp = 1 : numel(fitPosition) % k fit position
        % % % % % % % % % % % [~, xxTbl, yyTbl, ~, ~] = befSzVsOther(subjInfo, szCharTbl, siCharTbl, fitPosition(kfp));
        [fitTbl, xxTbl, yyTbl, xxFitTbl, yyFitTbl, xxOtherTbl, yyOtherTbl] = aroundSzVsOther(subjInfo, szCharTbl, siCharTbl, fitPosition(kfp));
        % Get mean value of the other (not before or not after)
        yOther = varfun(@(x) cellfun(@(y) stg.withinSubjectStat(y, 'omitmissing'), x), yyOtherTbl);
        yOther.Properties.VariableNames = cellfun(@(x) x(5:end), yOther.Properties.VariableNames, 'UniformOutput', false); % Remove the 'Fun_' from variable names
        % Get mean value of the before or after
        yThis = varfun(@(x) mean(x,"all", "omitmissing"), yyTbl);
        yThis.Properties.VariableNames = cellfun(@(x) x(5:end), yThis.Properties.VariableNames, 'UniformOutput', false); % Remove the 'Fun_' from variable names
        % Get average waveform of IED rate
        myyTblCell{kfp} = stg.withinSubjectStat(yyTbl, 'omitmissing');
        xLen = numel(xxTbl{1, 1});
        daysTF = (xLen-1)*dpDesc.BinLenDu/3600/24 > 1.5; % Will the time units in the plot be days or hours?
        normFactor = daysTF*(24-1) + 1;
        if daysTF
            timeUnit = '(days)';
        else
            timeUnit = '(hours)';
        end
        switch fitPosition(kfp)
            case "before"
                xpCell{kfp} = (-(xLen-1)*dpDesc.BinLenDu/3600 : dpDesc.BinLenDu/3600 : 0)/normFactor;
                diffBe = yThis - yOther;
            case "after"
                xpCell{kfp} = (0 : dpDesc.BinLenDu/3600 : (xLen-1)*dpDesc.BinLenDu/3600)/normFactor;
                diffAf = yThis - yOther;
        end
        if plotTF
            plotColor = ones(size(yyTbl, 1), 1)*stg.subjColor(ksubj, :);
            numSz{kfp} = plotPeriEventDataWithOther(plotName, ksubj, xpCell{kfp}, yyTbl, myyTblCell{kfp}, yOther, plotColor); % Refer to the created axes using h.a.(plotName)
        end
    end
    if plotTF
        formatAxesInd(plotName, ksubj, subjInfo, siCharTbl.Properties.VariableNames, timeUnit, "xdata")
        for kfp = 1 : numel(fitPosition) % k fit position
            hax = h.a.(plotName)(ksubj);
            text(mean(xpCell{kfp}), hax.YLim(2), ['n=', num2str(numSz{kfp}(1))], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', stg.statFontSize)
        end
    end
    % % % xBeH = xpCell{1};
    % % % xAfH = xpCell{2};
    myyBeTbl = myyTblCell{1};
    myyAfTbl = myyTblCell{2};
end
function [xBeH, myyBeTbl, xAfH, myyAfTbl] = plotSiCharSz(subjInfo, szCharTbl, siCharTbl, ksubj)
    global stg
    global h
    fitPosition = ["before", "after"];
    [plotTF, plotName] = createFigInd(1);
    myyTblCell = cell(1, numel(fitPosition)); xpCell = cell(1, numel(fitPosition));
    for kfp = 1 : numel(fitPosition) % k fit position
        [~, xxTbl, yyTbl, ~, ~] = fitSiCharSz(subjInfo, szCharTbl, siCharTbl, fitPosition(kfp));
        myyTblCell{kfp} = stg.withinSubjectStat(yyTbl, 'omitmissing');
        xLen = numel(xxTbl{1, 1});
        daysTF = (xLen-1)*stg.dpBinLenS/3600/24 > 1.5;
        normFactor = daysTF*(24-1) + 1;
        if daysTF
            timeUnit = '(days)';
        else
            timeUnit = '(hours)';
        end
        switch fitPosition(kfp)
            case "before"
                xpCell{kfp} = (-(xLen-1)*stg.dpBinLenS/3600 : stg.dpBinLenS/3600 : 0)/normFactor;
            case "after"
                xpCell{kfp} = (0 : stg.dpBinLenS/3600 : (xLen-1)*stg.dpBinLenS/3600)/normFactor;
        end
        if plotTF
            plotColor = ones(size(yyTbl, 1), 1)*stg.subjColor(ksubj, :);
            numSz{kfp} = plotPeriEventData(plotName, ksubj, xpCell{kfp}, yyTbl, myyTblCell{kfp}, plotColor); % Refer to the created axes using h.a.(plotName)
        end
    end
    if plotTF
        formatAxesInd(plotName, ksubj, subjInfo, siCharTbl.Properties.VariableNames, timeUnit, "xdata")
        for kfp = 1 : numel(fitPosition) % k fit position
            hax = h.a.(plotName)(ksubj);
            text(mean(xpCell{kfp}), hax.YLim(2), ['n=', num2str(numSz{kfp}(1))], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', stg.statFontSize)
        end
    end
    xBeH = xpCell{1};
    xAfH = xpCell{2};
    myyBeTbl = myyTblCell{1};
    myyAfTbl = myyTblCell{2};
end
function plotSiCharSzAllPop(xBeH, myyBeTbl, xAfH, myyAfTbl, siCharTbl)
    global stg
    global h
    numChar = numel(stg.siCharToPlot);
    positionCm = [20, 10, stg.singleColumnWidth, min(25, stg.siCharHe*numChar)];
    [plotTF, plotName] = createFigPos(positionCm);
    if contains(plotName, 'Pop')
        statFunction = stg.acrossSubjectStat;
    else
        statFunction = stg.withinSubjectStat;
    end
    if plotTF
        if ~isfield(h.a, plotName)
            createAxesAllPopStack(plotName)
        end
        % Normalize the data in each cell of the table
        for kc = 1 : size(myyBeTbl, 2)
            for kr = 1 : size(myyBeTbl, 1)
                if stg.curFitPopNorm
                    normalizingFactor = max([max(myyBeTbl{kr, kc}), max(myyAfTbl{kr, kc})]);
                else
                    normalizingFactor = 1;
                end
                myyBeTbl{kr, kc} = myyBeTbl{kr, kc}/normalizingFactor;
                myyAfTbl{kr, kc} = myyAfTbl{kr, kc}/normalizingFactor;
            end
        end
        % Average and plot
        mmyyBeTbl = statFunction(myyBeTbl, 1, 'omitmissing');
        plotPeriEventData(plotName, 1, xBeH, myyBeTbl, mmyyBeTbl, stg.subjColor);
        mmyyAfTbl = statFunction(myyAfTbl, 1, 'omitmissing');
        plotPeriEventData(plotName, 1, xAfH, myyAfTbl, mmyyAfTbl, stg.subjColor);
        formatAxesAllPop(plotName, siCharTbl.Properties.VariableNames, 'Time (hours)', "xdata", ["any", "any"])
    end
end
function [mFitTblBe, mFitTblAf] = plotSiCharSzFit(subjInfo, szCharTbl, siCharTbl, ksubj)
    global stg
    global h
    numChar = numel(stg.siCharToPlot);
    fitPosition = ["before", "after"];
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSiCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl); % Refer to the created axes using h.a.(thisFcnName)
    end
    mFitTblCell = cell(1, numel(fitPosition));
    for kfp = 1 : numel(fitPosition)
        [fitTbl, ~, ~, xxFitTbl, yyFitTbl] = fitSiCharSz(subjInfo, szCharTbl, siCharTbl, fitPosition(kfp));
        mFitTblCell{kfp} = stg.withinSubjectStat(fitTbl, "omitmissing");
        if plotTF
            for kchar = 1 : numChar
                if isempty(xxFitTbl)
                    continue
                end
                for ksz = 1 : size(xxFitTbl, 1)
                    % if clust(kcl).nested
                    %     continue
                    % end
                    hax = h.a.(plotName)(ksubj, kchar);
                    plot(hax, xxFitTbl{ksz, kchar}, yyFitTbl{ksz, kchar}, ...
                       'Marker', 'none', 'LineWidth', 2, 'LineStyle', '-', 'Color', stg.fitColor(kfp*2-1, :), 'Marker', 'none');
                end
            end
        end
    end
    if plotTF
        formatAxesInd(plotName, ksubj, subjInfo, siCharTbl.Properties.VariableNames, 'Age (days)', "xdata")
    end
    mFitTblBe = mFitTblCell{1};
    mFitTblAf = mFitTblCell{2};
end
function plotSiCharSzFitAllPop(fitTblBe, fitTblAf, siCharTbl)
    % fitTblBe ... table containing the slopes, offsets and their confidence intervals before seizures, each row is mean of one subject
    % fitTblAf ... same for periods after seizures
    % siCharTbl .. table with signal characteristics, used to get characteristics' names
    global stg
    global h
    numSubj = stg.numSubj;
    numChar = numel(stg.siCharToPlot);
    str = strings(2, numChar); medFit = NaN(2, numChar); sigRanP = NaN(2, numChar);
    for kchar = 1 : numChar
        [str(1, kchar), ~, ~, medFit(1, kchar), ~, sigRanP(1, kchar)] = getBasicStats(fitTblBe.(stg.siCharToPlot(kchar) + "fitSlope"), 'Slope');
        [str(2, kchar), ~, ~, medFit(2, kchar), ~, sigRanP(2, kchar)] = getBasicStats(fitTblAf.(stg.siCharToPlot(kchar) + "fitSlope"), 'Slope');
    end
    [sigRanP, sigRanH] = getSignificanceFromP(sigRanP, stg.Q, stg.multCompCorr);

    if stg.plotCharsStacked
        figHeight = min(25, stg.saCharHe*numChar);
    else
        figHeight = 3 ; % Centimeters
    end
    positionCm = [10, 10, stg.singleColumnWidth, figHeight];
    [plotTF, plotName] = createFigPos(positionCm);
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp(['Before seizure (', num2str(stg.fitSzDurS), ' seconds):'])
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, siCharTbl.Properties.VariableNames, kchar);
            disp(['    ', ylbl])
            disp(str(1, kchar))
        end
        disp(['After seizure (', num2str(stg.fitSzDurS), ' seconds):'])
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, siCharTbl.Properties.VariableNames, kchar);
            disp(['    ', ylbl])
            disp(str(2, kchar))
        end
    end
    medFit = medFit/24; % Convert from change per day to change per hour which is more natural especially since the IED rate is in IEDs/hour
    if plotTF
        if ~isfield(h.a, plotName)
            if stg.plotCharsStacked
                createAxesAllPopStack(plotName)
            else
                createAxesAllPopSlopeBox(plotName)
            end
        end
        strP = strings(2, numChar);
        for kchar = 1 : numChar
            if stg.plotCharsStacked
                hax(1, 1) = h.a.(plotName)(1, kchar);
                hax(1, 2) = h.a.(plotName)(1, kchar);
            else
                hax(1, 1) = h.a.(plotName)(1, kchar);
                hax(1, 2) = h.a.(plotName)(1, kchar + 1);
            end
            for ksubj = 1 : numSubj
                % Before seizure
                fitSlope = fitTblBe{ksubj, kchar + 0*numChar}; % Here in the units of (IEDs/hour)/day
                fitSlope = fitSlope/24; % Converting to (IEDs/hour)/hour
                % x(1, :) = [-fitDurS/3600/24, 0];
                x(1, :) = [-1, 0];
                y{1}(ksubj, :) = x(1, :)*fitSlope;
                y{1}(ksubj, :) = y{1}(ksubj, :) - x(1, 1)*fitSlope; % Normalize to the start of the pre-seizure period (will look more natural)
                % % % % % % x(1, :) = x(1, :)*24;
                plot(hax(1), x(1, :), y{1}(ksubj, :), 'Marker', 'none', 'LineWidth', stg.lnWiFit, 'LineStyle', '-', 'Color', stg.fitColor(1, :), 'Marker', 'none');
                % After seizure
                fitSlope = fitTblAf{ksubj, kchar + 0*numChar}; % Here in the units of (IEDs/hour)/day
                fitSlope = fitSlope/24; % Converting to (IEDs/hour)/hour
                % x(2, :) = [0, fitDurS/3600/24];
                x(2, :) = [0, 1];
                y{2}(ksubj, :) = x(2, :)*fitSlope;
                y{2}(ksubj, :) = y{2}(ksubj, :) - x(2, 2)*fitSlope; % Normalize to the end of the post-seizure period (will look more natural)
                % % % % % % % % % % % x(2, :) = x(2, :)*24;
                plot(hax(2), x(2, :), y{2}(ksubj, :), 'Marker', 'none', 'LineWidth', stg.lnWiFit, 'LineStyle', '-', 'Color', stg.fitColor(3, :), 'Marker', 'none');
            end
            % Median before seizure
            % % % % % % % % x = [-fitDurS/3600/24, 0];
            ym(1, :) = x(1, :)*medFit(1, kchar);
            ym(1, :) = ym(1, :) - x(1, 1)*medFit(1, kchar); % Normalize to the start of the pre-seizure period (will look more natural)
            % % % % x = x*24;
            plot(hax(1), x(1, :), ym(1, :), 'k', 'LineWidth', 2, 'LineStyle','-');
            xt(1) = mean(x(1, :));
            signMed = getSignChar(medFit(1, kchar));
            strP(1, kchar) = string(['p=', num2str(sigRanP(1, kchar), 2), '(', signMed, ')']);
            % Median after seizure
            % % % % % % % % x = [0, fitDurS/3600/24];
            ym(2, :) = x(2, :)*medFit(2, kchar);
            ym(2, :) = ym(2, :) - x(2, 2)*medFit(2, kchar); % Normalize to the end of the post-seizure period (will look more natural)
            % % % % % % % % x = x*24;
            plot(hax(2), x(2, :), ym(2, :), 'k', 'LineWidth', 2, 'LineStyle','-');
            xt(2) = mean(x(2, :));
            signMed = getSignChar(medFit(2, kchar));
            strP(2, kchar) = string(['p=', num2str(sigRanP(2, kchar), 2), '(', signMed, ')']);
        end
        % Box plot, violin plot, swarm chart
        if ~stg.plotCharsStacked
            % Before
            for kchar = 1 : numChar
                hax(2, 1) = h.a.(plotName)(2, kchar);
                % % % % y = fitTblBe{:, kchar}/24;
                for kp = 1 : numel(stg.statsPlotType)
                    switch stg.statsPlotType(kp)
                        case "box"
                            boxplot(hax(2, 1), y{1}(:, 2), 'Colors', 'k', 'Widths', 8);
                        case "violin"
                            violinMedianConf(hax(2, 1), y{1}(:, 2), stg.fitColor(1, :), stg.numBootstrapSamples, 0.05)
                        case "swarm"
                            swarmchart(hax(2, 1), ones(size(y{1}(:, 2))), y{1}(:, 2), 3, stg.fitColor(1, :), "filled")
                    end
                    % % % % [yl, ~] = getYLimYTick([min(y{1}(:, 2)), max(y{1}(:, 2))]);
                    % % % % hax(2, 1).YLim = yl; hax(2, 1).YTickLabel = []; hax(2, 1).XTick = []; hax(2, 1).TickLength = 1.5*hax(2, 1).TickLength;
                    % % % % hax(2, 1).Box = 'on';
                end
            end
            % After
            for kchar = 1 : numChar
                hax(2, 2) = h.a.(plotName)(2, kchar + 1);
                % % % % % y = fitTblAf{:, kchar}/24;
                for kp = 1 : numel(stg.statsPlotType)
                    switch stg.statsPlotType(kp)
                        case "box"
                            boxplot(hax(2, 2), y{2}(:, 1), 'Colors', 'k', 'Widths', 8);
                        case "violin"
                            violinMedianConf(hax(2, 2), y{2}(:, 1), stg.fitColor(3, :), stg.numBootstrapSamples, 0.05)
                        case "swarm"
                            swarmchart(hax(2, 2), ones(size(y{2}(:, 1))), y{2}(:, 1), 3, stg.fitColor(3, :), "filled")
                    end
                    % % % % [yl, ~] = getYLimYTick([min(y{1}(:, 2)), max(y{1}(:, 2))]);
                    % % % % hax(2, 2).YLim = yl; hax(2, 2).YTickLabel = []; hax(2, 2).XTick = []; hax(2, 2).TickLength = 1.5*hax(2, 2).TickLength;
                    % % % % hax(2, 2).Box = 'on';
                end
            end
        end
        if stg.plotCharsStacked
            formatAxesAllPop(plotName, siCharTbl.Properties.VariableNames, 'Time (hours)', "none", ["any", "any"])
        else
            for kchar = 1 : numChar
                yForLim = [y{1}(:, 2); y{2}(:, 1)];
                [yl, ~] = getYLimYTick([min(yForLim), max(yForLim)]);
                xlabel(hax(1, 1), 'Before'); ylabel(hax(1, 1), getYlbl(plotName, siCharTbl.Properties.VariableNames, kchar))
                xlabel(hax(1, 2), 'After'); ylabel(hax(1, 2), getYlbl(plotName, siCharTbl.Properties.VariableNames, kchar));
                set(hax, 'YLim', yl);
                set(hax, 'Box', 'on')
                set(hax(2, :), 'XTick', [])
                set(hax(2, :), 'YTickLabel', [])
                set(hax, 'TickLength', [0.03 0.075])
                for khax = 1 : numel(hax)
                    hax(khax).YAxis.Exponent = floor(log10(yl(2)));
                end
            end
        end
        for kchar = 1 : numChar
            % % % % hax = h.a.(plotName)(kchar);
            % % yt = mean(hax.YLim);
            % % ht = text(hax, xt(1), yt, strP(1, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            if abs(hax(1, 1).YLim(1)) > hax(1, 1).YLim(2)
                yt = hax(1, 1).YLim(1);
                ht = text(hax(1, 1), xt(1), yt, strP(1, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'BackgroundColor','none');
            else                
                yt = hax(1, 1).YLim(2);
                ht = text(hax(1, 1), xt(1), yt, strP(1, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            end
            if sigRanH(1, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
            % % ht = text(hax, xt(2), yt, strP(2, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            if abs(hax(1, 2).YLim(1)) > hax(1, 2).YLim(2)
                yt = hax(1, 2).YLim(1);
                ht = text(hax(1, 2), xt(2), yt, strP(2, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'BackgroundColor','none');
            else                
                yt = hax(1, 2).YLim(2);
                ht = text(hax(1, 2), xt(2), yt, strP(2, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            end

            if sigRanH(2, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
        end
        % text(-fitDurS/3600/2, hax.YLim(1) - 0.8*diff(hax.YLim), 'Pre-seizure', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
        % text(fitDurS/3600/2, hax.YLim(1) - 0.8*diff(hax.YLim), 'Post-seizure', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
    end
end
function [baseline, ma, curExp, curExpFltB, curExpFltA, curExpTauH, curPwl, fe, fp] = plotSiCharSzCur(subjInfo, szCharTbl, siCharTbl, isiHist, ksubj)
    global stg
    global h
    numChar = numel(stg.siCharToPlot);
    [plotTF, plotName] = createFigInd(1);
    [~, xxTbl, yyTbl, ~, ~] = fitSiCharSz(subjInfo, szCharTbl, siCharTbl, "after");
    myyTbl = stg.withinSubjectStat(yyTbl, 'omitmissing');
    xLen = numel(xxTbl{1, 1});
    xH = 0 : stg.dpBinLenS/3600 : (xLen-1)*stg.dpBinLenS/3600;
    xD = xH/24;
    if plotTF
        if stg.fitSzDurS/3600/24 <= 1.5
        % % % if xD(end) <= 1.5
            xp = xH;
            timeUnit = '(hours)';
        else
            xp = xD;
            timeUnit = '(days)';
        end
        plotPeriEventData(plotName, ksubj, xp, [], myyTbl, []);
    end
    baseline = NaN(1, numChar); ma = NaN(1, numChar);
    fe = cell(1, numChar); fp = cell(1, numChar); % Fit objects
    fecoeff = cell(1, numChar); fpcoeff = cell(1, numChar); % Coefficients extracted from the fit objects (just vectors of numbers)
    tauH = cell(1, numChar); % Time constants of the exponential in hours
    fltB = cell(1, numChar); fltA = cell(1, numChar); % Filter coefficients. Filter will serve for simulated IED rate (or other siChar) generation.
    u1 = zeros(size(xH)); % Vector for unit impulse for drawing the impulse response
    u1(1) = 1; % Unit impulse
    for kchar = 1 : numChar
        y = myyTbl{1, kchar};
        toKeep = ~isnan(y);
        if sum(toKeep) < 3
            fecoeff{kchar} = NaN(1, 2*stg.simSiCharFltOrder + 2);
            fltB{kchar} = NaN(1, stg.simSiCharFltOrder);
            fltA{kchar} = NaN(1, stg.simSiCharFltOrder + 1);
            tauH{kchar} = NaN(1, stg.simSiCharFltOrder);
            fpcoeff{kchar} = NaN(1, 4);
        else
            toKeep = ~isnan(y);
            x = xH(toKeep);
            y = y(toKeep);
            baseline(kchar) = getBaseline(y);
            ma(kchar) = max(y);
            % Fit sum of multiple exponentials
            [fecoeff{kchar}, amplitudes, decayFactors, tauH{kchar}, fe{kchar}] = fitExponentials(x, y, baseline(kchar));
            % IIR filter design
            [fltB{kchar}, fltA{kchar}] = designFiltFromExponentials(decayFactors, amplitudes);
            % Power-law fit (c.f. Osorio 2010 - Quakes of the Brain). NOT USED NOW, PREPARED FOR FUTURE.
            [fpcoeff{kchar}, fp{kchar}] = fitPowerLaw(x, y, baseline(kchar));
            if plotTF
                hax = h.a.(plotName)(ksubj, kchar);
                % Prepare the impulse response
                yfl = filter(fltB{kchar}, fltA{kchar}, u1'); % The impulse response of the designed filter
                ypf = yfl + baseline(kchar);
% Plot the real data
ed = isiHist.BinEdges;
switch stg.isiHistXScale
    case "log"
        xd = log10(ed); % Data
    case "linear"
        xd = ed;
end
y = isiHist.Probability;
xd = repelem(xd, 3);
xd = xd(1 : end-1);
yd(1 : 3 : length(y)*3 - 2) = 0;
yd(2 : 3 : length(y)*3 - 1) = y;
yd(3 : 3 : length(y)*3) = y;
yd = [0, yd, 0];
switch stg.isiHistYScale
    case "log"
        flor = 1e-4;
        yd(yd <= flor) = flor;
        yd = log10(yd);
    case "linear"
        flor = 0;
end
xd = xd/24;
yd = yd*100;
yd = yd - min(yd);
yd = yd/max(yd)*max(ypf);
patch(xd, yd, stg.subjColor(ksubj, :), 'EdgeColor', stg.subjColor(ksubj, :), 'LineWidth', 2, 'FaceAlpha', 0.5);
% % % set(gca, 'Children', flipud(get(gca, 'Children')));
hold on

                % Plot the filter impulse response
                plot(hax, xp, ypf, 'LineWidth', 1, 'Color', stg.curColor, 'LineStyle', '-')
            end
        end
    end
    if plotTF
        formatAxesInd(plotName, ksubj, subjInfo, siCharTbl.Properties.VariableNames, ['Time ', timeUnit], "xdata")
        for kchar = 1 : numChar
            if exist('amplitudes', 'var')
                tA = sprintf('$A_{%d}$=%.0f,' , [1 : numel(amplitudes); amplitudes]);
                tA = [tA(1 : end - 1), ' IEDs/hour'];
                tT = sprintf('$\\tau_{%d}$=%.1f,', [1 : numel(tauH{kchar}); tauH{kchar}]);
                tT = [tT(1 : end - 1), ' hours'];
                t = [tA, 10, tT];
                text(hax.XLim(2), hax.YLim(2), t, 'Interpreter', 'latex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 6);
            end
        end
    end

    curExp = cell2table(fecoeff);
    curExp.Properties.VariableNames = stg.siCharToPlot;
    curExpFltB = cell2table(fltB);
    curExpFltB.Properties.VariableNames = stg.siCharToPlot;
    curExpFltA = cell2table(fltA);
    curExpFltA.Properties.VariableNames = stg.siCharToPlot;
    curExpTauH = cell2table(tauH);
    curExpTauH.Properties.VariableNames = stg.siCharToPlot;
    curPwl = cell2table(fpcoeff);
    curPwl.Properties.VariableNames = stg.siCharToPlot;
end
function [baselinePop, maPop, fltBPop, fltAPop, fePop] = plotSiCharSzCurAllPop(siCharCur, tax, myyAfTbl, siCharTbl)
    global stg
    global h
    numChar = numel(stg.siCharToPlot);
    positionCm = [20, 10, 0.8*stg.singleColumnWidth, min(25, stg.siCharHe*numChar)];
    [plotTF, plotName] = createFigPos(positionCm);
    xH = 0 : stg.dpBinLenS/3600 : (numel(tax)-1)*stg.dpBinLenS/3600;
    u1 = zeros(size(xH)); % Vector for unit impulse for drawing the impulse response
    u1(1) = 1; % Unit impulse
    if plotTF
        % % % % % % % [plotTF, plotName] = createFigPos(positionCm);
        if ~isfield(h.a, plotName)
            createAxesAllPopStack(plotName)
        end
    end
    % Prepare some other stuff useful for plotting before we enter the for loop over chars
    xD = xH/24;
    if stg.fitSzDurS/3600/24 <= 1.5
    % % % if xD(end) <= 1.5
        xp = xH;
        timeUnit = '(hours)';
        tauMult = 60;
    else
        xp = xD;
        timeUnit = '(days)';
        tauMult = 1;
    end

    % Initialize variables
    fePop = cell(1, numChar);
    fecoeffPop = cell(1, numChar); fpcoeffPop = cell(1, numChar);
    fltBPop = cell(1, numChar); fltAPop = cell(1, numChar);
    baselinePop = NaN(1, numChar); maPop = NaN(1, numChar);
    str = strings(2, numChar, 2*stg.simSiCharFltOrder + 2);
    tauHPop = cell(1, numChar);
    for kchar = 1 : numChar
        % All individual subjects
        % Statistics of exponential fit
        fecoeffAll = siCharCur.exp{:, kchar};
        fecoeffAllmea = mean(fecoeffAll, 1, 'omitmissing');
        fecoeffAllsem = std(fecoeffAll, [], 1, 'omitmissing')/sqrt(size(fecoeffAll, 1));
        fecoeffAllsdv = std(fecoeffAll, [], 1, 'omitmissing');
        fecoeffAllmed = median(fecoeffAll, 1, 'omitmissing');
        fecoeffAlliqr = iqr(fecoeffAll, 1);
        tauHAll = tauMult*siCharCur.tauH{:, kchar};
        tauHAllmea = mean(tauHAll, 1, 'omitmissing');
        tauHAllsem = std(tauHAll, [], 1, 'omitmissing')/sqrt(size(tauHAll, 1));
        tauHAllsdv = std(tauHAll, [], 1, 'omitmissing');
        tauHAllmed = median(tauHAll, 1, 'omitmissing');
        tauHAlliqr = iqr(tauHAll, 1);
        % Create names of the variables to be reported
        numcoeff = numel(fecoeffAllmea) - 2; % Number of exponentials' coefficients
        coeffname = cell(2, numcoeff + 2);
        for kc = 1 : numcoeff/2
            coeffname{1, 2*(kc-1) + 1} = ['    A', num2str(kc)];
            coeffname{1, 2*(kc-1) + 2} = ['  tau', num2str(kc)];
        end
        coeffname{1, end-1} = '   R^2';
        coeffname{1, end} = '  RMSE';
        coeffname(2, 1 : 4) = {'     K', ' alpha', '   R^2', '  RMSE'};
        coeffname = cellfun(@(x) ['            ', x], coeffname, 'UniformOutput', false);
        % Fill in values of the variables
        for kc = 1 : numcoeff/2 % + 2 to include R^2 and RMSE
            if stg.sdTF
                str(1, kchar, 2*(kc-1) + 1) = string([coeffname{1, 2*(kc-1) + 1}, ' = ', ...
                    num2str(fecoeffAllmea(kc), 2), '±', num2str(fecoeffAllsdv(kc), 2), ...
                    ' (', num2str(fecoeffAllmed(kc), 2), '±', num2str(fecoeffAlliqr(kc), 2), ')']);
                str(1, kchar, 2*(kc-1) + 2) = string([coeffname{1, 2*(kc-1) + 2}, ' = ', ...
                    num2str(tauHAllmea(kc), 2), '±', num2str(tauHAllsdv(kc), 2), ...
                    ' (', num2str(tauHAllmed(kc), 2), '±', num2str(tauHAlliqr(kc), 2), ')']); % In every second position instead of the rate constant we provide tau
            else
                str(1, kchar, 2*(kc-1) + 1) = string([coeffname{1, 2*(kc-1) + 1}, ' = ', ...
                    num2str(fecoeffAllmea(kc), 2), '±', num2str(fecoeffAllsem(kc), 2), ...
                    ' (', num2str(fecoeffAllmed(kc), 2), '±', num2str(fecoeffAlliqr(kc), 2), ')']);
                str(1, kchar, 2*(kc-1) + 2) = string([coeffname{1, 2*(kc-1) + 2}, ' = ', ...
                    num2str(tauHAllmea(kc), 2), '±', num2str(tauHAllsem(kc), 2), ...
                    ' (', num2str(tauHAllmed(kc), 2), '±', num2str(tauHAlliqr(kc), 2), ')']); % In every second position instead of the rate constant we provide tau
            end
        end
        if stg.sdTF
            str(1, kchar, end - 1) = string([coeffname{1, end - 1}, ' = ', ...
                num2str(fecoeffAllmea(end - 1), 2), '±', num2str(fecoeffAllsdv(end - 1), 2), ...
                    ' (', num2str(fecoeffAllmed(end - 1), 2), '±', num2str(fecoeffAlliqr(end - 1), 2), ')']);
            str(1, kchar, end) = string([coeffname{1, end}, ' = ', ...
                num2str(fecoeffAllmea(end), 2), '±', num2str(fecoeffAllsdv(end), 2), ...
                    ' (', num2str(fecoeffAllmed(end), 2), '±', num2str(fecoeffAlliqr(end), 2), ')']);
        else
            str(1, kchar, end - 1) = string([coeffname{1, end - 1}, ' = ', ...
                num2str(fecoeffAllmea(end - 1), 2), '±', num2str(fecoeffAllsem(end - 1), 2), ...
                    ' (', num2str(fecoeffAllmed(end - 1), 2), '±', num2str(fecoeffAlliqr(end - 1), 2), ')']);
            str(1, kchar, end) = string([coeffname{1, end}, ' = ', ...
                num2str(fecoeffAllmea(end), 2), '±', num2str(fecoeffAllsem(end), 2), ...
                    ' (', num2str(fecoeffAllmed(end), 2), '±', num2str(fecoeffAlliqr(end), 2), ')']);
        end
        % Statistics of power-law fit
        fpcoeffAll = siCharCur.pwl{:, kchar};
        fpcoeffAllmea = mean(fpcoeffAll, 1, 'omitmissing');
        fpcoeffAllsem = std(fpcoeffAll, [], 1, 'omitmissing')/sqrt(size(fpcoeffAll, 1));
        fpcoeffAllsdv = std(fpcoeffAll, [], 1, 'omitmissing');
        fpcoeffAllmed = median(fpcoeffAll, 1, 'omitmissing');
        fpcoeffAlliqr = iqr(fpcoeffAll, 1);
        for kc = 1 : numel(fpcoeffAllmea)
            if stg.sdTF
                str(2, kchar, kc) = string([coeffname{2, kc}, ' = ', ...
                    num2str(fpcoeffAllmea(kc), 2), '±', num2str(fpcoeffAllsdv(kc), 2), ...
                    ' (', num2str(fpcoeffAllmed(kc), 2), '±', num2str(fpcoeffAlliqr(kc), 2), ')']);
            else
                str(2, kchar, kc) = string([coeffname{2, kc}, ' = ', ...
                    num2str(fpcoeffAllmea(kc), 2), '±', num2str(fpcoeffAllsem(kc), 2), ...
                    ' (', num2str(fpcoeffAllmed(kc), 2), '±', num2str(fpcoeffAlliqr(kc), 2), ')']);
            end
        end
        % Fit of the averaged data over all subjects, i.e. population statistics (although the term population might be mathematically inaccurate, it is actually sample)
        % Normalize the data in each cell of the table
        myyAfNorm = NaN(size(myyAfTbl, 1), numel(myyAfTbl{1, 1}));
        for ksubj = 1 : size(myyAfTbl, 1)
            if stg.curFitPopNorm
                normalizingFactor = max(myyAfTbl{ksubj, kchar}, [], 2, 'omitmissing');
            else
                normalizingFactor = 1;
            end
            myyAfNorm(ksubj, :) = myyAfTbl{ksubj, kchar}/normalizingFactor;
        end
        % Get the mean or median
        yPop = stg.acrossSubjectStat(myyAfNorm, 1, 'omitmissing'); % Get the sample mean or median
        toKeep = ~isnan(yPop);
        if sum(toKeep) < 3
            fecoeffPop{kchar} = NaN(1, numcoeff); fflcoeff{kchar} = NaN(1, 4); fpcoeff{kchar} = NaN(1, 4);
        else
            x = xH(toKeep); % xH is in hours
            yPop = yPop(toKeep);
            baselinePop(kchar) = getBaseline(yPop);
            maPop(kchar) = max(yPop);
            % Fit sum of multiple exponentials
            [fecoeffPop{kchar}, amplitudes, decayFactors, tauHPop{kchar}, fePop{kchar}] = fitExponentials(x, yPop, baselinePop(kchar), 1);
            % % % % % % % % % % % % % tauHPop{kchar} = tauMult*tauHPop{kchar};
            % IIR filter design
            [fltBPop{kchar}, fltAPop{kchar}] = designFiltFromExponentials(decayFactors, amplitudes);
            % Power-law fit (c.f. Osorio 2010 - Quakes of the Brain)
            [fpcoeffPop{kchar}] = fitPowerLaw(x, yPop, baselinePop(kchar));
            if plotTF
                hax = h.a.(plotName)(1, kchar);
                for ksubj = 1 : size(myyAfTbl, 1)
                    y = myyAfTbl{ksubj, kchar};
                    baseline = getBaseline(y);
                    % % % % Plot the original fitted exponentials. Just to check that the filter works correctly.
                    % % % plot(hax, xp, y, 'LineWidth', 0.1, 'Color', stg.subjColor(ksubj, :))
                    % % % yfe = siCharCur.fe{ksubj, kchar}(xH);
                    % % % yfe = yfe + baseline;
                    % % % hold on
                    % % % plot(hax, xp, yfe, 'LineWidth', 3, 'Color', [0 0 0.9], 'LineStyle', '--')
                    if any(isnan(siCharCur.expFltB{ksubj, kchar}))
                        yfl = NaN(size(u1'));
                    else
                        yfl = filter(siCharCur.expFltB{ksubj, kchar}, siCharCur.expFltA{ksubj, kchar}, u1'); % The impulse response of the designed filter
                        yfl = yfl + baseline;
                        yfl = yfl/max(yfl);
                    end
                    plot(hax, xp, yfl, 'LineWidth', 0.5, 'Color', stg.curColor, 'LineStyle', '-')
                    % % % hold off
                end
                hax = h.a.(plotName)(1, kchar);
                yfl = filter(fltBPop{kchar}, fltAPop{kchar}, u1'); % The impulse response of the designed filter
                yfl = yfl + baselinePop(kchar);
                yfl = yfl/max(yfl);
                plot(hax, xp, yfl, 'LineWidth', 1.5, 'Color', [0 0 0], 'LineStyle', '-')
                % Violin
                positionCmViolin = [20, 10, 6, 2.2];
                createFigPos(positionCmViolin, [plotName, 'Violin']);
                ymax = [0.25, 2, 2];
                for kt = 1 : numcoeff/2
                    h.a.([plotName, 'Violin'])(kchar, kt) = subplot(numChar, numcoeff/2, (kchar-1)*numChar + kt);
                    hax = h.a.([plotName, 'Violin'])(kchar, kt);
                    y = siCharCur.tauH{:, kchar}(:, kt)/24;
                    v = violinMedianConf(hax, y, stg.curColor, stg.numBootstrapSamples, 0.05);
                    hold on
                    swarmchart(hax, ones(size(y)), y, 3, stg.curColor, "filled")
                    line(hax, v.DensityWidth*[-1, 1]/2 + 1, tauHPop{kchar}(1, kt)*[1, 1]/24, 'LineWidth', 1, 'Color', 'k')
                    ylabel(['tau', num2str(kt), ' ', timeUnit])
                    ylim([0 ymax(kt)])
                    % ylim([0 2])
                end
                % yAll = siCharCur.tauH{:, kchar}(:, :);
                % [yl, ~] = getYLimYTick([min(yAll, [], 'all'), max(yAll, [], 'all')]);
                % yl = [0 48];
                % set(h.a.([plotName, 'Violin'])(kchar, :), 'YLim', yl)
                % yt = [0 24 48];
                % set(h.a.([plotName, 'Violin'])(kchar, :), 'YTick', yt)
                set(h.a.([plotName, 'Violin'])(kchar, :), 'Box', 'on')
                if stg.printFigures
                    print('Tau violin.eps', '-depsc', '-vector')
                end
            end
        end
    end
    if plotTF
        formatAxesAllPop(plotName, siCharTbl.Properties.VariableNames, ['Time ', timeUnit], "xdata")
    end
    
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp(['Curve fit after seizure (', num2str(stg.fitSzDurS), ' seconds):'])
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, siCharTbl.Properties.VariableNames, kchar);
            disp(['    ', ylbl])
            disp('        Individual subjects');
            disp(['tau in hours/', num2str(tauMult)])
            if stg.sdTF
                disp('mean±SD (median±IQR)')
            else
                disp('mean±SEM (median±IQR)')
            end
            for km = 1 : 2 % Over models (sum of exponentials, power-law)
                for k = 1 : size(str, 3) % Over coefficients
                    disp(str(km, kchar, k))
                end
            end
            disp('        Population data');
            disp(['            Fit coeff exp:          ', num2str(fecoeffPop{kchar}(1 : end - 2))])
            disp(['                  Fit tau:          ', num2str(tauMult*tauHPop{kchar})])
            disp(['            Fit coeff pwr-law:       ', num2str(fpcoeffPop{kchar}(1 : end - 2))])
            disp(['            Goodness of fit exp:    ', num2str(fecoeffPop{kchar}(end - 1 : end))])
            disp(['            Goodness of fit pwr-law: ', num2str(fpcoeffPop{kchar}(end - 1 : end))])
        end
    end
end
function [risingTbl, risingClTbl, risingNonClTbl] = siCharRisingAroundSz(szCharTbl, siCharTbl, szBelongsToClust)
    global stg
    szBelongsToClust = logical(szBelongsToClust)';
    numChar = numel(stg.siCharToPlot);
    numSz = size(szCharTbl, 1);
    risingTbl = table('Size', [1, 3*numChar], 'VariableTypes', repelem("double", 3*numChar));
    risingTbl{1, :} = deal(NaN);
    for kchar = 1 : numChar
        risingTbl.Properties.VariableNames{1, (kchar-1)*3 + 1} = char(stg.siCharToPlot(kchar) + "NumPos"); % Number of seizures during which the siChar difference was positive
        risingTbl.Properties.VariableNames{1, (kchar-1)*3 + 2} = char(stg.siCharToPlot(kchar) + "NumTotal"); % Number of seizures
        risingTbl.Properties.VariableNames{1, (kchar-1)*3 + 3} = char(stg.siCharToPlot(kchar) + "FracPos"); % Fraction of seizures during which the siChar difference was positive
    end
    risingClTbl = risingTbl;
    risingNonClTbl = risingTbl;
    for kchar = 1 : numChar
        y = siCharTbl{:, stg.siCharToPlot(kchar)};
        y = fillmissing(y, 'linear', 1, 'MaxGap', ceil(4*3600/stg.dpBinLenS));
        y = filter(1/stg.simMovAveLen*ones(1, stg.simMovAveLen), 1, y);
        posTF = NaN(numSz, 1);
        for ksz = 1 : numSz
            be = find(siCharTbl.tax < szCharTbl.szOnsN(ksz), 1, 'last'); % Row of siCharTbl of the value before the seizure
            af = find(siCharTbl.tax > szCharTbl.szOnsN(ksz), 1, 'first'); % Row of siCharTbl of the value after the seizure
            posTF(ksz) = y(af) > y(be);
        end
        risingTbl{1, (kchar-1)*3 + 1} = sum(posTF, "omitnan"); % Number of seizures during which the siChar difference was positive
        risingTbl{1, (kchar-1)*3 + 2} = sum(~isnan(posTF)); % Number of seizures
        risingTbl{1, (kchar-1)*3 + 3} = risingTbl{1, (kchar-1)*3 + 1}/risingTbl{1, (kchar-1)*3 + 2}; % Fraction of seizures during which the siChar difference was positive
        risingClTbl{1, (kchar-1)*3 + 1} = sum(posTF(logical(szBelongsToClust)), "omitnan"); % Number of seizures during which the siChar difference was positive
        risingClTbl{1, (kchar-1)*3 + 2} = sum(~isnan(posTF(logical(szBelongsToClust)))); % Number of seizures
        risingClTbl{1, (kchar-1)*3 + 3} = risingClTbl{1, (kchar-1)*3 + 1}/risingClTbl{1, (kchar-1)*3 + 2}; % Fraction of seizures during which the siChar difference was positive
        risingNonClTbl{1, (kchar-1)*3 + 1} = sum(posTF(~logical(szBelongsToClust)), "omitnan"); % Number of seizures during which the siChar difference was positive
        risingNonClTbl{1, (kchar-1)*3 + 2} = sum(~isnan(posTF(~logical(szBelongsToClust)))); % Number of seizures
        risingNonClTbl{1, (kchar-1)*3 + 3} = risingNonClTbl{1, (kchar-1)*3 + 1}/risingNonClTbl{1, (kchar-1)*3 + 2}; % Fraction of seizures during which the siChar difference was positive
    end
end 
function siCharRisingAroundSzPop(risingTbl, risingClTbl, risingNonClTbl)
    numChar = size(risingTbl, 2)/3;
    for kchar = 1 : numChar
        rising = risingTbl{:, kchar*3} - 0.5
        signrank(rising)
        getBasicStats(risingTbl{:, kchar*3}, 'rising')
        risingCl = risingClTbl{:, kchar*3} - 0.5
        signrank(risingCl)
        getBasicStats(risingClTbl{:, kchar*3}, 'risingCl')
        risingNonCl = risingNonClTbl{:, kchar*3} - 0.5
        signrank(risingNonCl)
        getBasicStats(risingNonClTbl{:, kchar*3}, 'risingNonCl')
    end
end
% Plot seizure and signal characteristics and simulated signal characteristics analyses in one figure
function plotSsChar(ksubj, subjInfo, szCharTbl, siCharTbl)
    global stg
    if stg.numSubj == 1
        positionCm = [10, 25, 2.8*stg.singleColumnWidth, 6];
        [plotTF, plotName] = createFigPos(positionCm);
    else
        [plotTF, plotName] = createFigInd(1);
    end
    if plotTF
        plotSsCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl); % Refer to the created axes using h.f.(thisFcnName)
        % szCharTblVarNames = szCharTbl.Properties.VariableNames;
        siCharTblVarNames = siCharTbl.Properties.VariableNames;
        % charTblVarNames = [szCharTblVarNames, siCharTblVarNames];
        formatAxesInd(plotName, ksubj, subjInfo, siCharTblVarNames, 'Age (days)', "xdata")
    end
end
function [fitTbl, xxFitTbl, yyFitTbl] = plotSfCharWhFit(ksubj, subjInfo, szCharTbl, siCharTbl)
    global stg
    global h
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    [fitTbl, ed, ~, yySzTbl, ~, ~, xxFitTbl, yyFitTbl] = fitSaCharWh(subjInfo, szCharTbl, siCharTbl); % Compute the trend
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSaCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl); % Refer to the created axes using h.f.(thisFcnName)
        
        % Plot the trend
        for kchar = 1 : numSzChar
            y = yySzTbl{1, kchar};
            % Plotting patch (bar graph)
            hax = h.a.(plotName)(ksubj, kchar);
            pax = [ed(1), ed(1), repelem(ed(2:end-1), 3), ed(end), ed(end), ed(1)];
            pay = repelem(y, 3);
            pay(1 : 3 : end) = 0;
            pay = [pay, 0, 0]; %#ok<AGROW>
            pay(isnan(pay)) = 0;
            patch(hax, pax, pay, 1 - 0.2*(1 - stg.subjColor(ksubj, :)))
            if kchar == 1
                h.p.(plotName)(ksubj, kchar, 2).YData(2 : 3 : end - 1) = ceilToEven1sig(max(pay)); % Change the height of seizures in the plot of sz rate
            end
        end
        for kchar = 1 : numChar
            % Plotting fitted line
            x = xxFitTbl{1, kchar};
            y = yyFitTbl{1, kchar};
            hax = h.a.(plotName)(ksubj, kchar);
            plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', stg.lnWiFitMean, 'Color', stg.subjColor(ksubj, :)*stg.subjColorSubjMeanMult);
            if kchar ~= numChar
                hax.Children([1 2 3 4]) = hax.Children([1 3 4 2]);
            end
        end
        charTblVarNames = [szCharTbl.Properties.VariableNames, siCharTbl.Properties.VariableNames];
        formatAxesInd(plotName, ksubj, subjInfo, charTblVarNames, 'Age (days)', "analysisPeriod")
    end
end
function plotSfCharWhFitAllPop(fitTbl, xxFitTbl, yyFitTbl, szCharTbl, siCharTbl) %#ok<INUSD>
    global stg
    global h
    normalizedTF = true;
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    charTblVarNames = [szCharTbl.Properties.VariableNames, siCharTbl.Properties.VariableNames];
    str = strings(1, numChar); medFit = NaN(1, numChar); sigRanP = NaN(1, numChar);
    for kchar = 1 : numChar
        [str(kchar), ~, ~, medFit(kchar), ~, sigRanP(kchar)] = getBasicStats(fitTbl.(stg.saCharToPlot(kchar) + "fitSlope"), 'Slope');
    end
    [sigRanP, sigRanH] = getSignificanceFromP(sigRanP, stg.Q, stg.multCompCorr);
    positionCm = [10, 10, stg.singleColumnWidth, min(25, stg.saCharHe*numChar)];
    [plotTF, plotName] = createFigPos(positionCm);
    
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp('Whole recording:')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, charTblVarNames, kchar);
            disp(['    ', ylbl])
            disp(str(kchar))
        end
    end
    if plotTF
        if ~isfield(h.a, plotName)
            createAxesAllPopStack(plotName)
        end
        strP = strings(1, numChar);
        for kchar = 1 : numChar
            for ksubj = 1 : stg.numSubj
                hax = h.a.(plotName)(1, kchar);
                if normalizedTF
                    x = [0 1];
                    fitSlope = fitTbl{ksubj, kchar};
                    y = x*fitSlope;
                else
                    x = xxFitTbl{ksubj, kchar}; %#ok<UNRCH>
                    y = yyFitTbl{ksubj, kchar};
                end
                plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', stg.lnWiFit, 'Color', stg.subjColor(ksubj, :)); % Not normalized
            end
            if normalizedTF
                y = x*medFit(1, kchar);
                plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            else
                x = stg.acrossSubjectStat(xxFitTbl{:, kchar}, 1, 'omitmissing'); %#ok<UNRCH>
                y = stg.acrossSubjectStat(yyFitTbl{:, kchar}, 1, 'omitmissing');
                plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            end
            switch sign(medFit(kchar))
                case -1
                    signMed = '-';
                case 0
                    signMed = '0';
                case 1
                    signMed = '+';
            end
            strP(kchar) = string(['  p=', num2str(sigRanP(1, kchar), 2), '(', signMed, ')']);
        end
        formatAxesAllPop(plotName, charTblVarNames, 'Normalized time', "none", ["any", "any"])
        for kchar = 1 : numChar
            hax = h.a.(plotName)(kchar);
            % xt = mean(hax.XLim);
            xt = 0;
            % yt = mean(hax.YLim);
            yt = hax.YLim(2);
            ht = text(hax, xt, yt, strP(kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','left', 'VerticalAlignment','top', 'BackgroundColor','none');
            if sigRanH(1, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
        end
    end
end
function [mFitTblBe, mFitTblDu, mFitTblAf] = plotSfCharClFit(ksubj, subjInfo, szCharTbl, siCharTbl, cl)
    % ksubj ...... subject index
    % subjInfo ... table with subject info
    % szCharTbl .. table with seizure characteristics
    % siCharTbl .. table with signal characteristics
    % cl ......... structure with clusters of given mouse
    % mFitTbl .... table with mean slopes, offsets and confidence intervlas (the CI mean has no real meaning though)
    % mxFit0 ..... mean cluster duration (probably not further used)
    % myFit0 ..... mean difference between value of the fit at the offset - onset of the cluster for each characteristic (probably not further used)
    global stg
    global h
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    fitPosition = ["before", "during", "after"];
    [plotTF, plotName] = createFigInd(1);
    if plotTF
        plotSaCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl); % Refer to the created axes using h.a.(thisFcnName)
    end
    mFitTblCell = cell(1, numel(fitPosition));
    for kfp = 1 : numel(fitPosition)
        [fitTbl, eded, ~, yySzTbl, ~, ~, xxFitTbl, yyFitTbl] = fitSaCharCl(subjInfo, siCharTbl, cl, fitPosition(kfp));
        mFitTblCell{kfp} = stg.withinSubjectStat(fitTbl, "omitmissing");
        if plotTF
            for kchar = 1 : numChar
                if isempty(xxFitTbl)
                    continue
                end
                maxHeight = 0;
                for kcl = 1 : size(xxFitTbl, 1)
                    % if clust(kcl).nested
                    %     continue
                    % end
                    hax = h.a.(plotName)(ksubj, kchar);
                    if kchar <= numSzChar && kfp == 2
                        % Bar graph (it is actually a patch but looks like a bar graph)
                        y = yySzTbl{kcl, kchar};
                        ed = eded(kcl, :);
                        % Plotting
                        hax = h.a.(plotName)(ksubj, kchar);
                        pax = [ed(1), ed(1), repelem(ed(2:end-1), 3), ed(end), ed(end), ed(1)];
                        pay = repelem(y, 3);
                        pay(1 : 3 : end) = 0;
                        pay = [pay, 0, 0]; %#ok<AGROW>
                        pay(isnan(pay)) = 0;
                        maxHeight = max([maxHeight, pay]);
                        h.pa.(plotName)(ksubj, kchar, 1) = patch(hax, pax, pay, 1 - 0.2*(1 - stg.subjColor(ksubj, :)), 'Tag', 'szCharClDuBar');
                        % Fit plot
                        x = xxFitTbl{kcl, kchar};
                        y = yyFitTbl{kcl, kchar};
                        maxHeight = max([maxHeight, y]);
                        h.p.(plotName)(ksubj, kchar, 3) = plot(hax, x, y, ...
                           'Marker', 'none', 'LineWidth', 2, 'LineStyle', '-', 'Color', stg.fitColor(kfp, :), 'Marker', 'none');

                        % % % h.p.(plotName)(ksubj, kchar, 3) = plot(hax, x, y, ...
                        % % %     'Marker', 'none', 'LineStyle', '--', 'LineWidth', 1, 'Color', stg.subjColor(ksubj, :)*stg.subjColorSubjMeanMult, 'Tag', 'szCharClDuFit');
                        if kchar <= numSzChar && kfp == 2
                            uistack(h.pa.(plotName)(ksubj, kchar), 'bottom')
                        end
                        uistack(h.p.(plotName)(ksubj, kchar, 3), 'top')
                    else
                        h.p.(plotName)(ksubj, kchar, 3) = plot(hax, xxFitTbl{kcl, kchar}, yyFitTbl{kcl, kchar}, ...
                           'Marker', 'none', 'LineWidth', 2, 'LineStyle', '-', 'Color', stg.fitColor(kfp, :), 'Marker', 'none');
                    end
                end
                if kchar == 1
                    h.p.(plotName)(ksubj, kchar, 2).YData(2 : 3 : end-1) = ceilToEven1sig(maxHeight); % Change the height of seizures in the plot of sz rate
                end
            end
        end
    end
    if plotTF
        charTblVarNames = [szCharTbl.Properties.VariableNames, siCharTbl.Properties.VariableNames];
        formatAxesInd(plotName, ksubj, subjInfo, charTblVarNames, 'Age (days)', "analysisPeriod")
    end
    mFitTblBe = mFitTblCell{1};
    mFitTblDu = mFitTblCell{2};
    mFitTblAf = mFitTblCell{3};
end
function plotSfCharClFitAllPop(fitTblBe, fitTblDu, fitTblAf, szCharTbl, siCharTbl)
    global stg
    global h
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    charTblVarNames = [szCharTbl.Properties.VariableNames, siCharTbl.Properties.VariableNames];
    numSubj = stg.numSubj;
    str = strings(3, numChar); medFit = NaN(3, numChar); sigRanP = NaN(3, numChar);
    for kchar = 1 : numChar
        [str(1, kchar), ~, ~, medFit(1, kchar), ~, sigRanP(1, kchar)] = getBasicStats(fitTblBe.(stg.saCharToPlot(kchar) + "fitSlope"), 'Slope');
        [str(2, kchar), ~, ~, medFit(2, kchar), ~, sigRanP(2, kchar)] = getBasicStats(fitTblDu.(stg.saCharToPlot(kchar) + "fitSlope"), 'Slope');
        [str(3, kchar), ~, ~, medFit(3, kchar), ~, sigRanP(3, kchar)] = getBasicStats(fitTblAf.(stg.saCharToPlot(kchar) + "fitSlope"), 'Slope');
    end
    [sigRanP, sigRanH] = getSignificanceFromP(sigRanP, stg.Q, stg.multCompCorr);
    positionCm = [20, 10, stg.singleColumnWidth, min(25, stg.saCharHe*numChar)];
    [plotTF, plotName] = createFigPos(positionCm);
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp('Before cluster:')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, charTblVarNames, kchar);
            disp(['    ', ylbl])
            disp(str(1, kchar))
        end
        disp('During cluster:')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, charTblVarNames, kchar);
            disp(['    ', ylbl])
            disp(str(2, kchar))
        end
        disp('After cluster:')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, charTblVarNames, kchar);
            disp(['    ', ylbl])
            disp(str(3, kchar))
        end
    end
    if plotTF
        if ~isfield(h.a, plotName)
            createAxesAllPopStack(plotName)
        end
        strP = strings(3, numChar); xt = NaN(3, 1);
        for kchar = 1 : numChar
            hax = h.a.(plotName)(kchar);
            for ksubj = 1 : numSubj
                % Before cluster
                fitSlope = fitTblBe{ksubj, kchar + 0*numChar};
                x = [0 0.33];
                y = [-1 0]*fitSlope;
                plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFit, 'LineStyle', '-', 'Color', stg.fitColor(1, :), 'Marker', 'none');
                % During cluster
                fitSlope = fitTblDu{ksubj, kchar + 0*numChar};
                x = [0.34 0.66];
                y = [0 1]*fitSlope;
                plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFit, 'LineStyle', '-', 'Color', stg.fitColor(2, :), 'Marker', 'none');
                % After cluster
                fitSlope = fitTblAf{ksubj, kchar + 0*numChar};
                x = [0.67 1];
                y = [0 1]*fitSlope;
                plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFit, 'LineStyle', '-', 'Color', stg.fitColor(3, :), 'Marker', 'none');
            end
            % Medians
            % Before cluster
            x = [0 0.33];
            xt(1) = mean(x);
            y = [-1 0]*medFit(1, kchar);
            plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            if isnan(sign(medFit(1, kchar)))
                signMed = 'n/a';
            else
                switch sign(medFit(1, kchar))
                    case -1
                        signMed = '-';
                    case 0
                        signMed = '0';
                    case 1
                        signMed = '+';
                end
            end
            strP(1, kchar) = string(['p=', num2str(sigRanP(1, kchar), 2), '(', signMed, ')']);
            % During cluster
            x = [0.34 0.66];
            xt(2) = mean(x);
            y = [0 1]*medFit(2, kchar);
            plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            if isnan(sign(medFit(2, kchar)))
                signMed = 'n/a';
            else
                switch sign(medFit(2, kchar))
                    case -1
                        signMed = '-';
                    case 0
                        signMed = '0';
                    case 1
                        signMed = '+';
                end
            end
            strP(2, kchar) = string(['p=', num2str(sigRanP(2, kchar), 2), '(', signMed, ')']);
            % After cluster
            x = [0.67 1];
            xt(3) = mean(x);
            y = [0 1]*medFit(3, kchar);
            if isnan(sign(medFit(3, kchar)))
                signMed = 'n/a';
            else
                switch sign(medFit(3, kchar))
                    case -1
                        signMed = '-';
                    case 0
                        signMed = '0';
                    case 1
                        signMed = '+';
                end
            end
            plot(hax, x, y, 'Marker', 'none', 'LineWidth', stg.lnWiFitMean, 'LineStyle', '-', 'Color', 'k', 'Marker', 'none');
            strP(3, kchar) = string(['p=', num2str(sigRanP(3, kchar), 2), '(', signMed, ')']);
            hax.XLim = [0 1];
        end
        strP(contains(strP, 'n/a')) = "";
        formatAxesAllPop(plotName, charTblVarNames, 'Normalized time', "none", ["any", "any"])
        for kchar = 1 : numChar
            hax = h.a.(plotName)(kchar);
            if abs(hax.YLim(1)) > hax.YLim(2)
                yt = hax.YLim(1);
                ht = text(hax, xt(1), yt, strP(1, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'BackgroundColor','none');
            else                
                yt = hax.YLim(2);
                ht = text(hax, xt(1), yt, strP(1, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            end
            if sigRanH(1, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
            % % % % ht = text(hax, xt(2), yt, strP(2, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            if abs(hax.YLim(1)) > hax.YLim(2)
                yt = hax.YLim(1);
                ht = text(hax, xt(2), yt, strP(2, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'BackgroundColor','none');
            else                
                yt = hax.YLim(2);
                ht = text(hax, xt(2), yt, strP(2, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            end
            if sigRanH(2, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
            % % % % ht = text(hax, xt(3), yt, strP(3, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            if abs(hax.YLim(1)) > hax.YLim(2)
                yt = hax.YLim(1);
                ht = text(hax, xt(3), yt, strP(3, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'BackgroundColor','none');
            else                
                yt = hax.YLim(2);
                ht = text(hax, xt(3), yt, strP(3, kchar), 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'BackgroundColor','none');
            end
            if sigRanH(3, kchar)
                ht.FontWeight = "bold";
            else
                ht.FontWeight = "normal";
            end
        end
        hax.XTick = [];
        hax.XLabel = [];
        text(0 + 0.33/2, hax.YLim(1) - 0.2*diff(hax.YLim), 'Pre-cluster', 'FontSize',  8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
        text(0.34 + 0.33/2, hax.YLim(1) - 0.2*diff(hax.YLim), 'Cluster', 'FontSize',  8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
        text(0.67 + 0.33/2, hax.YLim(1) - 0.2*diff(hax.YLim), 'Post-cluster', 'FontSize',  8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
    end
end
function [cirTbl, ppTbl, rrTbl] = plotSfCharCiFit(ksubj, subjInfo, szCharTbl, siCharTbl)
    global stg
    global h
    numPointsAroundCircle = 97; % Only for smoothness of the circular sections
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    [cirTbl, edSz, ~, rrSzTbl, ~, ppSiTbl, rrSiTbl, ppCirTbl, rrCirTbl, rayCirTbl, omnCirTbl] = fitSaCharCi(szCharTbl, siCharTbl); % Calculate the circular statistics
    % Prepare the output (interpolate to same length)
    ppTbl = table('Size', [0, numChar], 'VariableTypes', repelem("cell", 1, numChar), 'VariableNames', stg.saCharToPlot);
    rrTbl = ppTbl;
    for kchar = 1 : numChar
        if kchar <= numSzChar
            % Plotting circular bar graph
            r = rrSzTbl{1, kchar};
            if kchar == 1 % Here, the size of the individual seizure vectors has no meaning, so we normalize histogram to max of the histogram
                r = r/max(r);
            else
                y1 = szCharTbl{:, stg.szCharToPlot(kchar)};
                r = r/max(y1);
            end
            circleCenter = max(r)/1000; % Pure zero minimum radius makes trouble probably due to rounding errors
            numBetw = (numPointsAroundCircle-1)/(numel(edSz)-1) - 3;
            pap = NaN(1, (numel(edSz)-1)*(numBetw+3) + 1); par = pap;
            pap(1) = edSz(1);
            par(1) = 0;
            for ke = 1 : numel(edSz) - 1
                pap((ke-1)*(numBetw+3) + (2 : 2+numBetw+1)) = linspace(edSz(ke), edSz(ke+1), numBetw+2);
                pap((ke-1)*(numBetw+3) + 2+numBetw+2) = edSz(ke+1);
                par((ke-1)*(numBetw+3) + (2 : 2+numBetw+1)) = r(ke);
                par((ke-1)*(numBetw+3) + 2+numBetw+2) = circleCenter;
            end
            pap = pap*2*pi;
            par(isnan(par)) = circleCenter;
            par(end) = circleCenter;
            ppTbl{1, kchar} = {pap};
            rrTbl{1, kchar} = {par};
        else
            ppTbl{1, kchar} = {ppSiTbl{1, kchar - numSzChar}}; %#ok<CCAT1>
            rrTbl{1, kchar} = {rrSiTbl{1, kchar - numSzChar}}; %#ok<CCAT1>
        end
    end
    ppTbl = varfun(@(x) x{1}, ppTbl);
    rrTbl = varfun(@(x) x{1}, rrTbl);

    % Plot
    [plotTF, plotName] = createFigInd(1.1);
    if plotTF
        numCol = plotSaCharDataCirc(plotName, ksubj, szCharTbl, ppSiTbl, rrSiTbl); % Refer to the created axes using h.f.(thisFcnName)
        % Plot the fit
        for kchar = 1 : numChar
            hax = h.a.(plotName)(ksubj, kchar);
            if kchar <= numSzChar
                pap = ppTbl{1, kchar};
                par = rrTbl{1, kchar};
                [x, y] = pol2cart(-pap - pi/2, par);
                z = -5*ones(size(x));
                h.pa.(plotName)(ksubj, kchar, 2) = patch(hax, x, y, z, 1 - 0.2*(1 - stg.subjColor(ksubj, :)), 'EdgeColor', stg.subjColor(ksubj, :)*stg.subjColorSubjMeanMult);
            end
            % Plotting resultant vector
            resVecWi = 1.5;
            prv = ppCirTbl{1, kchar};
            rrv = rrCirTbl{1, kchar}*hax.XLim(2);
            [x, y] = pol2cart(-prv - pi/2, rrv);
            h.p.(plotName)(ksubj, kchar, 3) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', stg.subjColor(ksubj, :)*0);
            % Arrowhead
            [x(1), y(1)] = pol2cart(-(prv(2)-5/6*pi) - pi/2, rrv(2)/4);
            x(1) = sum(x);
            y(1) = sum(y);
            h.p.(plotName)(ksubj, kchar, 4) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', stg.subjColor(ksubj, :)*0);
            x(1) = 0; y(1) = 0;
            [x(1), y(1)] = pol2cart(-(prv(2)+5/6*pi) - pi/2, rrv(2)/4);
            x(1) = sum(x);
            y(1) = sum(y);
            h.p.(plotName)(ksubj, kchar, 5) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', stg.subjColor(ksubj, :)*0);
            % % % % R scaling
            % % % if kchar == 1
            % % %     factor = ceilToEven1sig(max([r, rrCirTbl{1, kchar}]));
            % % %     for kplot = 1 : 2
            % % %         h.p.(plotName)(ksubj, kchar, kplot).XData = factor*h.p.(plotName)(ksubj, kchar, kplot).XData;
            % % %         h.p.(plotName)(ksubj, kchar, kplot).YData = factor*h.p.(plotName)(ksubj, kchar, kplot).YData;
            % % %     end
            % % %     h.pa.(plotName)(ksubj, kchar, 1).XData = factor*h.pa.(plotName)(ksubj, kchar, 1).XData;
            % % %     h.pa.(plotName)(ksubj, kchar, 1).YData = factor*h.pa.(plotName)(ksubj, kchar, 1).YData;
            % % %     h.a.(plotName)(ksubj, kchar).XLim = factor*h.a.(plotName)(ksubj, kchar).XLim;
            % % %     h.a.(plotName)(ksubj, kchar).YLim = factor*h.a.(plotName)(ksubj, kchar).YLim;
            % % % end
            % Text
            xt = 0;
            yt = hax.YLim(1);
            text(hax, xt, yt, ['r=', num2str(rrv(2), 1)], 'Interpreter', 'tex', 'FontSize', stg.axFontSize, 'FontWeight', 'normal',...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
            % ylbl = getYlbl(plotName, szCharTbl.Properties.VariableNames, kchar);
            ylbl = stg.saCharNameShort(kchar);
            yt = hax.YLim(2);
            text(hax, xt, yt, ylbl, 'Interpreter', 'tex', 'FontSize', stg.axFontSize, 'FontWeight', 'normal',...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        end
        if rem(numCol, 2) == 0
            hax = h.a.(plotName)(ksubj, numCol/2);
            x = hax.XLim(2);
            y = hax.YLim(2) + 0.3*diff(hax.YLim);
        else
            hax = h.a.(plotName)(ksubj, ceil(numCol/2));
            x = hax.XLim(1) + diff(hax.XLim)/2;
            y = hax.YLim(2) + 0.3*diff(hax.YLim);
        end
        text(hax, x, y, subjInfo.subjNm, 'Interpreter', 'none', 'FontSize', stg.axFontSize + 1,...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        for kchar = 1 : numChar + 1
            set(h.a.(plotName)(ksubj, kchar), 'Units', 'normalized');
        end
    end
end
function plotSfCharCiFitAllPop(cirTbl, ppTbl, rrTbl, szCharTbl, siCharTbl)
    global stg
    global h
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    str = strings(1, numChar); resM = NaN(1, numChar); resR = NaN(1, numChar); rayP = NaN(1, numChar); omnP = NaN(1, numChar);
    for kchar = 1 : numChar
        p = cirTbl{:, kchar};
        p = p(~isnan(p));
        resM(kchar) = circ_mean(p);
        if ~stg.andrzejak
            resR(kchar) = circ_r(p);
        else
            resR(kchar) = circ_t(p);
        end
        rayP(kchar) = circ_rtest(p);
        omnP(kchar) = circ_otest(p);
        clear p
        xNm = 'Time of day';
        str(kchar) = ...
            [repelem(' ', 1, 13-numel(xNm)), xNm, ' = ', num2str(mod(resM(kchar)/2/pi*24, 24), '%.2g'), ', L=', num2str(resR(kchar), '%.2g'), 10, ...
             '        Ray p = ', num2str(rayP(kchar), '%.2g'), 10, ...
             '       Omni p = ', num2str(omnP(kchar), '%.2g'), 10];
    end
    [~, rayH] = getSignificanceFromP(rayP, stg.Q, stg.multCompCorr);
    % [~, omnH] = getSignificanceFromP(omnP, stg.Q, stg.multCompCorr);
    positionCm = [10, 25, stg.singleColumnWidth, min(25, stg.saCharHe*numChar*1.1)];
    [plotTF, plotName] = createFigPos(positionCm);
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp('Circadian statistics (animal-wise):')
        disp(['(', plotName, ')'])
        for kchar = 1 : numChar
            ylbl = getYlbl(plotName, [szCharTbl.Properties.VariableNames, siCharTbl.Properties.VariableNames], kchar);
            disp(['    ', ylbl])
            disp(str(kchar))
        end
    end
    if plotTF
        % Calculate axes positions
        % % % % % % % % margGlob = [0 0 0 0]; % Left, bottom, right, top
        % % % % % % % % marg = [0.1 0.5 0.1 0.5];
        [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlobCi, stg.margCi);
        h.f.(plotName).Units = stg.units;
        
        % Plot the data
        numCol = ceil((numChar+1)^0.6);
        numRow = ceil((numChar+1)/numCol);
        axWiHe = 0.7*min(spWi/numCol, spHe/numRow);
        for kchar = 1 : numChar+1
            h.a.(plotName)(1, kchar) = axes('Units', stg.units, 'Position', [...
                spx(mod(1 - 1, numc) + 1) + mod(kchar-1, numCol)*spWi/numCol + max(spWi/numCol - axWiHe)/2,...
                spy(ceil(1/numc)) + (numRow-ceil(kchar/numCol))*spHe/numRow + 0.0*spHe/numRow,...
                axWiHe, ...
                axWiHe],...
                'NextPlot', 'add', 'Visible', 'off', 'XLimMode', 'manual', 'YLimMode', 'manual');
            hax = h.a.(plotName)(1, kchar);
            
            % % % r1 = NaN(stg.numSubj, numel(rrTbl{1, 1}));
            if kchar ~= numChar+1
                r = NaN(stg.numSubj, size(rrTbl{1, kchar}, 2)); p = r;
                for ksubj = 1 : height(cirTbl)
                    % Subject's circular bar graph
                    r(ksubj, :) = rrTbl{ksubj, kchar};
                    p(ksubj, :) = ppTbl{ksubj, kchar};
                    [x, y] = pol2cart(-p(ksubj, :) - pi/2, r(ksubj, :));
                    % h.p.(plotName)(ksubj, kchar, 1) = plot(hax, x, y, 'Color', 1 - 0.5*(1 -stg.subjColor(ksubj, :)));
                    % Subject's mean resultant vector
                    pv = [0, cirTbl{ksubj, kchar}];
                    % r = rrCirTbl{ksubj, kchar};
                    rv = [0 1];
                    [x, y] = pol2cart(-pv - pi/2, rv);
                    h.p.(plotName)(ksubj, kchar, 1) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-',...
                        'LineWidth', stg.lnWiFit, 'Color', stg.subjColor(ksubj, :)); % Not normalized
                end
                % Circular bar graph mean
                rm = mean(r, 1, 'omitmissing');
                pm = p(1, :);
                pm(end+1) = pm(1); %#ok<AGROW>
                rm(end+1) = rm(1); %#ok<AGROW>
                [x, y] = pol2cart(-pm - pi/2, rm);
                h.p.(plotName)(1, kchar, 2) = plot(hax, x, y, 'Color', 'k', 'LineWidth', 1.5);
            end
            % Circle
            circleP = (0 : 1 : 360)/360*2*pi;
            circleR = 1;
            % if ~all(isnan(rrCirTbl{:, kchar}))
            %     circleR = ceilToEven1sig(max(rrCirTbl{:, kchar}, [], 'all'));
            % else
            %     circleR = 1;
            % end
            [x, y] = pol2cart(-circleP - pi/2, circleR);
            h.p.(plotName)(1, kchar, 3) = plot(x, y, 'k');
            hax.XLim = max([x, y])*[-1 1];
            hax.YLim = hax.XLim;
            % Night shading
            circleP = (-90 : 1 : 90)/360*2*pi;
            % colo = [linspace(1, 0, 90), 1, linspace(0, 1, 90)]'*[1 1 1];
            colo = [linspace(1, 0.7, 90), 1, linspace(0.7, 1, 90)]'*[1 1 1];
            colo = permute(colo, [1 3 2]);
            [x, y] = pol2cart(-circleP - pi/2, circleR);
            z = -10*ones(size(x));
            h.pa.(plotName)(1, kchar, 1) = patch(x, y, z, colo, 'EdgeColor', 'none');
            if kchar ~= numChar+1
                % Mean resultant vector
                resVecWi = 1.5;
                prv = [0, resM(kchar)];
                rrv = [0, resR(kchar)];
                [x, y] = pol2cart(-prv - pi/2, rrv);
                h.p.(plotName)(1, kchar, 4) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', 'k');
                % Arrowhead
                [x(1), y(1)] = pol2cart(-(prv(2)-5/6*pi) - pi/2, rrv(2)/4);
                x(1) = sum(x);
                y(1) = sum(y);
                h.p.(plotName)(1, kchar, 5) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', 'k');
                x(1) = 0; y(1) = 0;
                [x(1), y(1)] = pol2cart(-(prv(2)+5/6*pi) - pi/2, rrv(2)/4);
                x(1) = sum(x);
                y(1) = sum(y);
                h.p.(plotName)(1, kchar, 6) = plot(hax, x, y, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', resVecWi, 'Color', 'k');
                x = 0;
                y = hax.YLim(2);
                % text(hax, x, y, num2str(y), 'Interpreter', 'tex', 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
                % ylbl = getYlbl(plotName, szCharTbl.Properties.VariableNames, kchar);
                ylbl = stg.saCharNameShort(kchar);
                text(hax, x, y, ylbl, 'Interpreter', 'tex', 'FontSize', stg.axFontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
                strP = string(['  p=', num2str(rayP(1, kchar), 2)]);
                % ht = text(hax, xrv(2)/2, yrv(2)/2, strP, 'FontSize', stg.statFontSize, 'HorizontalAlignment','center', 'VerticalAlignment', 'middle', 'BackgroundColor','w');
                ht = text(hax, 0, -1, strP, 'FontSize', stg.axFontSize, 'HorizontalAlignment','center', 'VerticalAlignment', 'top', 'BackgroundColor','none');
                if rayH(1, kchar)
                    ht.FontWeight = "bold";
                else
                    ht.FontWeight = "normal";
                end
            end
            % Plot explanatory circle
            if kchar == numChar+1
                plotExplanatoryCircle(hax)
            end
        end
        % Title
        if rem(numCol, 2) == 0
            hax = h.a.(plotName)(1, numCol/2);
            x = hax.XLim(2);
            y = hax.YLim(2) + 0.3*diff(hax.YLim);
        else
            hax = h.a.(plotName)(1, ceil(numCol/2));
            x = hax.XLim(1) + diff(hax.XLim)/2;
            y = hax.YLim(2) + 0.3*diff(hax.YLim);
        end
        text(hax, x, y, 'All mice', 'Interpreter', 'none', 'FontSize', stg.axFontSize + 1,...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        for kchar = 1 : numChar + 1
            set(h.a.(plotName)(1, kchar), 'Units', 'normalized');
        end
    end
end
function plotSsExplainSumOfExp(siCharCurPop)
    global stg
    global h
    positionCm = [10, 25, stg.singleColumnWidth, 6];
    [plotTF, plotName] = createFigPos(positionCm);
    linestyle = {'k--'; 'k-.'; 'k:'};
    if plotTF
        spx = 2;
        spy = 1.3;
        spWi = 6;
        spHe = 4;
        h.f.(plotName).Units = stg.units;
        h.a.(plotName) = axes('Units', stg.units, 'Position', [spx, spy, spWi, spHe]);
        expDurS = 24*3600; % Exponential durations
        tax = (stg.dpBinLenS/10) : (stg.dpBinLenS/10) : expDurS;
        x = tax/3600;
        coeff = coeffvalues(siCharCurPop.fe{1});
        for ke = 1 : numel(coeff)/2
            a = coeff((ke-1)*2 + 1);
            b = coeff(ke*2);
            e(ke, :) = a*exp(b*x); %#ok<AGROW>
            legText{ke} = [num2str(a, 2), '*exp(', num2str(b, 2), '*t)']; %#ok<AGROW>
        end
        % e(2, :) = siCharCurPop.fe{1}.c*exp(siCharCurPop.fe{1}.d*x);
        % e(3, :) = siCharCurPop.fe{1}.e*exp(siCharCurPop.fe{1}.f*x);
        e(ke + 1, :) = sum(e, 1);
        hpesum = plot(x, e(end, :), 'Color', [1 0.8 0], 'LineWidth', 2);
        % % % legText{ke + 1} = strjoin(legText, ' + ');
        legText{ke + 1} = 'sum';
        hold on
        for ke = 1 : numel(coeff)/2
            hpe(ke) = plot(x, e(ke, :), linestyle{ke}, 'LineWidth', 1); %#ok<AGROW>
        end
        xlim([0 24])
        legend([hpe, hpesum], legText, 'Location', 'northeast')
        hax = gca;
        hax.XTick = [0 6 12 18 24];
        if 1/(-coeff(end)) < 1
            xlim([0 2])
            hax.XTick = ([0 0.5 1 1.5 2]);
        end
    end
end
function plotSsExplainConvolution(subjInfo, siCharTbl, siCharCurPop)
    global stg
    global h
    positionCm = [10, 25, 2*stg.singleColumnWidth, 6];
    [plotTF, plotName] = createFigPos(positionCm);
    if plotTF
        if ~strcmpi(subjInfo.subjNm, 'Mouse13')
            return
        end
        spx = 2;
        spy = 1.3;
        spWi = 14;
        spHe = 4;
        h.f.(plotName).Units = stg.units;
        h.a.(plotName) = axes('Units', stg.units, 'Position', [spx, spy, spWi, spHe]);
        % timeStart = 60;
        % timeStop = 61.5;
        timeStart = 83;
        timeStop = 86;
        siCharTbl = siCharTbl((siCharTbl.tax >= subjInfo.dob + timeStart) & (siCharTbl.tax < subjInfo.dob + timeStop), :);
        x = siCharTbl.tax - subjInfo.dob;
        siCharTbl.sz(isnan(siCharTbl.sz)) = 0;
        szSub = find(siCharTbl.sz);
        for ksz = 1 : size(szSub)
            z1 = zeros(size(siCharTbl.tax));
            z1(szSub(ksz)) = 1;
            y = filter(siCharCurPop.expFltB{1, 1}, siCharCurPop.expFltA{1, 1}, z1);
            plot(x, y, 'Marker', 'none', 'LineWidth', 2, 'Color', hsv2rgb(ksz/size(szSub, 1), 1, 1), 'Tag', 'simIedInd');
            hold on
        end
        y = filter(siCharCurPop.expFltB{1, 1}, siCharCurPop.expFltA{1, 1}, siCharTbl.sz*stg.dpBinLenS/3600/24);
        simMovAveLen = stg.simMovAveLen;
        ySm = filter(1/simMovAveLen*ones(1, simMovAveLen), 1, y); % Smoothed
        hpSmSim = plot(x, ySm, 'Marker', 'none', 'LineStyle', ':', 'LineWidth', 2, 'Color', [0.6 0.6 0.6], 'Tag', 'simIed');
        hpSim = plot(x, y, 'Marker', 'none', 'LineWidth', 0.5, 'Color', 'k', 'Tag', 'simIed');
        hsSz = scatter(x(szSub), zeros(size(szSub)), 'Marker', '*', 'MarkerEdgeColor', 'k');
        clear x y
        h.f.(plotName).Units = 'normalized';
        h.a.(plotName).Units = 'normalized';
        legend([hpSim, hpSmSim, hsSz], {'Simulated IED rate', 'Smoothed simulated IED rate', 'Seizures'}, 'Location', 'northeast')
    end
end
function plotSimSim(siCharSimSim)
    global stg
    global h
    numChar = 2*numel(stg.siCharToPlot);
    positionCmViolin = [20, 10, 6, 6];
    [plotTF, plotName] = createFigPos(positionCmViolin);
    if plotTF
        hax(1, 1) = subplot(221);
        y = siCharSimSim.rmsnrmse;
        v = violinMedianConf(hax(1, 1), y, stg.simColor - 0.1, stg.numBootstrapSamples, 0.05);
        hold on
        swarmchart(hax(1, 1), ones(size(y)), y, 3, stg.simColor - 0.1, "filled")
        ylabel('Relative error (-)')
        y1 = y;
        
        hax(1, 2) = subplot(222);
        y = siCharSimSim.rmsnrmsePop;
        v = violinMedianConf(hax(1, 2), y, stg.simPopColor - 0.1, stg.numBootstrapSamples, 0.05);
        hold on
        swarmchart(hax(1, 2), ones(size(y)), y, 3, stg.simPopColor - 0.1, "filled")
        y2 = y;
        y = [y1; y2];
        y(abs(y) > 1e12) = NaN;
        yl = [0 1.25*ceil(max(y, [], 1, 'omitmissing'))];
        set(hax(1, :), 'YLim', yl)

        hax(2, 1) = subplot(223);
        y = siCharSimSim.rho;
        v = violinMedianConf(hax(2, 1), y, stg.simColor - 0.1, stg.numBootstrapSamples, 0.05);
        hold on
        swarmchart(hax(2, 1), ones(size(y)), y, 3, stg.simColor - 0.1, "filled")
        ylabel('Pearson r (-)')
        y1 = y;

        hax(2, 2) = subplot(224);
        y = siCharSimSim.rhoPop;
        v = violinMedianConf(hax(2, 2), y, stg.simPopColor - 0.1, stg.numBootstrapSamples, 0.05);
        hold on
        swarmchart(hax(2, 2), ones(size(y)), y, 3, stg.simPopColor - 0.1, "filled")
        y2 = y;
        y = [y1; y2];
        y(abs(y) > 1e12) = NaN;
        yl = [1.25*floor(2*min(y, [], 1, 'omitmissing'))/2, 1.25*ceil(2*max(y, [], 1, 'omitmissing') + 0.05)/2];
        set(hax(2, :), 'YLim', yl)

        % set(hax, 'YLim', yl)
        set(hax, 'Box', 'on')
        set(hax, 'XTick', [])
    end

                % % for kt = 1 : numcoeff/2
                % %     h.a.([plotName, 'Violin'])(kchar, kt) = subplot(numChar, numcoeff/2, (kchar-1)*numChar + kt);
                % %     hax = h.a.([plotName, 'Violin'])(kchar, kt);
                % %     y = siCharCur.tauH{:, kchar}(:, kt);
                % %     v = violinMedianConf(hax, y, stg.curColor, stg.numBootstrapSamples, 0.05);
                % %     hold on
                % %     swarmchart(hax, ones(size(y)), y, 3, stg.curColor, "filled")
                % %     line(hax, v.DensityWidth*[-1, 1]/2 + 1, tauHPop{kchar}(1, kt)*[1, 1], 'LineWidth', 1, 'Color', 'k')
                % %     ylabel(['tau', num2str(kt), ' (hours)'])
                % % end
                % % yAll = siCharCur.tauH{:, kchar}(:, :);
                % % [yl, ~] = getYLimYTick([2*min(yAll, [], 'all'), 2*max(yAll, [], 'all')]);
                % % yl = yl/2;
                % % set(h.a.([plotName, 'Violin'])(kchar, :), 'YLim', yl)
                % % set(h.a.([plotName, 'Violin'])(kchar, :), 'Box', 'on')
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp(['Signal characteristics simulation (fit length ', num2str(stg.fitSzDurS), ' seconds, smoothing ', num2str(stg.simMovAveLen), ' blocks):'])
        disp(['(', plotName, ')'])

        resultStr = getBasicStats(siCharSimSim.rmsnrmse, 'rmse');
        ss = strsplit(resultStr, char(10)); %#ok<CHARTEN>
        disp(ss{1})
        resultStr = getBasicStats(siCharSimSim.rmsnrmsePop, 'rmsePop');
        ss = strsplit(resultStr, char(10)); %#ok<CHARTEN>
        disp(ss{1})
        resultStr = getBasicStats(siCharSimSim.rho, 'pearson');
        ss = strsplit(resultStr, char(10)); %#ok<CHARTEN>
        disp(ss{1})
        resultStr = getBasicStats(siCharSimSim.rhoPop, 'pearsonPop');
        ss = strsplit(resultStr, char(10)); %#ok<CHARTEN>
        disp(ss{1})
    end
end

% Calculation functions
function [szRate, ed, xOk, yOk] = calcSzRate(subjInfo, szCharTbl, siCharTbl, binlen)
    % szRate ..... number of seizures in the bin divided by the amount of valid signal in the bin, it is in 1/days
    % ed ......... bin edges
    % xOk ........ beginnings of all analysis bins (not the seizure rate analysis bins), used for signal OK plotting
    % yOk ........ NaN where no signal is present, used for signal OK plotting
    % subjInfo ... to get the date of analysis end and date of birth
    % szCharTbl .. to get seizure times
    % siCharTbl .. to get times of valid signal
    % binlen ..... length of bins for seizure rate analysis, specified in the calling function
    requiredProportionValid = 0.25;
    [xOk, yOk, xx] = getSzOkXY(subjInfo, siCharTbl);
    % ed = subjInfo.dob : binlen : subjInfo.anEndN;
    ed = 0 : binlen : subjInfo.anEndN - subjInfo.dob; % Edges for histcounts to get number of seizures in each bin
    sigLen = NaN(1, numel(ed) - 1); % Length of valid signal in that bin (valid for seizure identification, not for IEDs or other signal features)
    for kb = 1 : numel(ed) - 1
        stSub = find(xx(2, :) > ed(kb), 1, 'first'); % Start subscript
        enSub = find(xx(1, :) < ed(kb+1), 1,'last'); % End subscript
        if isempty(stSub) || isempty(enSub) || enSub < stSub
            sigLen(kb) = NaN; % It is, in fact zero, but if we entered zero, it could lead to division by zero below
        else
            dxx = diff(xx, 1, 1);
            sigLen(kb) = sum(dxx(stSub : enSub));
            sigLen(kb) = sigLen(kb) - max((ed(kb) - xx(1, stSub)), 0) - max((xx(2, enSub) - ed(kb+1)), 0);
        end
    end
    hc = histcounts(szCharTbl.szOnsN - subjInfo.dob, ed)./sigLen;
    hc(sigLen/binlen < requiredProportionValid) = NaN; % Data point where less than requiredProportionValid is valid, is considered untrustworthy
    szRate = hc;
end
function [paxTbl, psdTbl, psdciTbl] = psdSiChar(siCharTbl)
    global stg
    charToPlot = ["sz", stg.siCharToPlot];
    numChar = numel(charToPlot);
    Ts = stg.dpBinLenS/3600/24;
    paxpax = cell(1, numChar); psdpsd = cell(1, numChar); psdcipsdci = cell(1, numChar); % Initialization
    for kchar = 1 : numChar
        y = siCharTbl.(charToPlot(kchar));
        maxPer = 12; % In days
        y = fillmissing(y, 'pchip');
        y = y - mean(y);
        [psd, fax, psdci] = pwelch(y, fix(maxPer/Ts), fix(maxPer/2/Ts), [], 1/Ts); % The Fs is in days^-1
        pax = 1./fax(fax >= 1/maxPer & fax <= 2); % Period axis, instead of frequency axis (or in some graphs the more familiar time axis)
        paxpax{kchar} = pax'; % It is the same for all characteristics and all subjects. This is just for consistency.
        psdpsd{kchar} = psd(fax >= 1/maxPer & fax <= 2)';
        psdcipsdci{kchar} = psdci(fax >= 1/maxPer & fax <= 2)';
    end
    paxTbl = cell2table(paxpax, "VariableNames", charToPlot);
    psdTbl = cell2table(psdpsd, "VariableNames", charToPlot);
    psdciTbl = cell2table(psdcipsdci, "VariableNames", charToPlot);
end
function stats = subjectStats(stg, subjInfo, ds, dp, clustStats)
    stats = table;
    stats.Subject = subjInfo.subjNm;
    stats.observPer = subjInfo.anEndDt - subjInfo.anStartDt; % Observation period in days (not accounting for dropouts)
    % Data to stem
    fn = fieldnames(ds);
    for kfn = 1 : numel(fn)
        stats.([fn{kfn}, 'Num']) = height(ds.(fn{kfn}));
        vn = ds.(fn{kfn}).Properties.VariableNames;
        for kvn = 1 : numel(vn)
            if strcmp(vn{kvn}, 'OnsDt') % Mean or median of onset times is irrelevant
                continue
            end
            stats.([fn{kfn}, vn{kvn}]) = stg.withinSubjectStat(ds.(fn{kfn}).(vn{kvn}));
        end
    end
    % Data to plot
    fn = fieldnames(dp);
    for kfn = 1 : numel(fn)
        vn = dp.(fn{kfn}).Properties.VariableNames;
        for kvn = 1 : numel(vn)
            if strcmp(vn{kvn}, 'tax') % Mean or median of onset times is irrelevant
                continue
            end
            stats.([fn{kfn}, vn{kvn}]) = stg.withinSubjectStat(dp.(fn{kfn}).(vn{kvn}), 'omitnan');
        end
    end
% % % % % % % subjInfo
% % % % % % % 
% % % % % % % 
% % % % % % %     % % stats.szNum = length(szCharTbl.szOnsN);
% % % % % % %     % % % stats.szFreq = stats.szNum/stats.observPer; % This is wrong because it does not take into account the dropouts
% % % % % % %     colNames = szCharTbl.Properties.VariableNames;
% % % % % % %     for kc = 1 : numel(colNames)
% % % % % % %         if strcmp(colNames{kc}, 'szOnsN')
% % % % % % %             % % % stats.szIsiH = stg.withinSubjectStat(diff(szCharTbl.szOnsN)*24); % This is wrong because it does not take into account the dropouts
% % % % % % %         else
% % % % % % %             stats.(colNames{kc}) = stg.withinSubjectStat(szCharTbl.(colNames{kc}), 'omitmissing');
% % % % % % %         end
% % % % % % %     end
% % % % % % %     colNames = clustStats.Properties.VariableNames;
% % % % % % %     for kc = 1 : numel(colNames)
% % % % % % %         stats.(colNames{kc}) = clustStats.(colNames{kc});
% % % % % % %     end
% % % % % % %     colNames = siCharTbl.Properties.VariableNames;
% % % % % % %     for kc = 6 : numel(colNames)
% % % % % % %         stats.(colNames{kc}) = stg.withinSubjectStat(siCharTbl.(colNames{kc}), 'omitmissing');
% % % % % % %     end
end
function [fitTbl, ed, xxTbl, yyTbl, xxFitTbl, yyFitTbl] = fitSzCharWh(subjInfo, szCharTbl, siCharTbl)
    % subjInfo ... used to get date of birth (dob) and period of monitoring
    % szCharTbl .. seizure characteristics
    % siCharTbl .. signal characteristics, used to get info on signal dropouts
    % fitTbl ..... table containing linear fits
    % ed ......... edges of bins in which seizure count or mean seizure characteristics are calculated
    % xxTbl ...... centers of the bins (not used in any function so far)
    % yyTbl ...... seizure count or mean seizure characteristics in each bin
    % xxFitTbl ... start and end of the fitted line for each characteristic
    % yyFitTbl ... start and end of the fitted line for each characteristic
    global stg
    numChar = numel(stg.szCharToPlot);
    fitTbl = fitTblInit(1, stg.szCharToPlot);
    
    % Bar graph
    firstSz = szCharTbl.szOnsN(1) - subjInfo.dob;
    lastSz = szCharTbl.szOnsN(end) - subjInfo.dob;
    % ed = floor(siCharTbl.tax(1) - subjInfo.anStartN) : stg.whWinLen : ceil((siCharTbl.tax(end) - subjInfo.anStartN)/7)*7; % From the start to
    % the end of the recording, might be biased since some subjects start seizing later which makes an impression of lower sz rate at the beginning.
    % And, also, having the end of the recording after a long inter-cluster period makes an impression of a decline in sz rate which may induce bias
    % as well (consieder perfectly regular seizures to see the problem). Therefore, I set it from the first to the last seizure.
    % THIS IS NOT ENTIRELY FAIR EITHER SINCE IT OVERESTIMATES THE RATE IN THE FIRST AND LAST BIN. But at least the first-order polynomial should be
    % correct.
    % ed = linspace(szCharTbl.szOnsN(1) - subjInfo.anStartN, szCharTbl.szOnsN(end) - subjInfo.anStartN, 11);
    ed = [firstSz, ceil(firstSz) : stg.whWinLen : floor(lastSz), lastSz];
    xx = cell(1, numChar); yy = cell(1, numChar); xxFit = cell(1, numChar); yyFit = cell(1, numChar);
    for kchar = 1 : numChar
        x = NaN(1, numel(ed)-1); y = NaN(1, numel(ed)-1); binWeights = NaN(1, numel(ed)-1);
        for kb = 1 : numel(ed) - 1
            x(kb) = (ed(kb) + ed(kb+1))/2; % Bin centers are at the mean of the edges
            szSub = find((szCharTbl.szOnsN - subjInfo.dob) >= ed(kb), 1, 'first') :...
                    find((szCharTbl.szOnsN - subjInfo.dob) <=  ed(kb + 1), 1, 'last'); % Subscripts of seizures belonging to the bin
            % Bin weight is computed as the (number of valid points)/(expected number of points in a standard bin)
            [xOK, yOK] = getSzOkXY(subjInfo, siCharTbl);
            xOK = xOK(1 : 3 : end-2);
            yOK = yOK(1 : 3 : end-2);
            numAllDatapoints = (ed(3)-ed(2))/(xOK(2) - xOK(1)); % Number of datapoints in a standard bin (the first and last bin is usually shorter). Might be non-integer.
            numValidDatapoints = sum(~isnan(yOK(xOK > ed(kb) & xOK < ed(kb+1))));
            binWeights(kb) = numValidDatapoints/numAllDatapoints;
            % Compute the value in the given bin
            if kchar == 1
                % y(kb) = numel(szSub)/(ed(kb+1) - ed(kb)); % Does not take into account recording dropouts
                y(kb) = numel(szSub)/(numValidDatapoints*stg.dpBinLenS)*3600*24; % Does not take into account recording dropouts
            else
                if isempty(szSub)
                    y(kb) = NaN;
                else
                    y(kb) = mean(szCharTbl.(stg.szCharToPlot(kchar))(szSub), "omitmissing");
                end
            end
        end
        
        % Line fit. We are fitting the bin values so that the fit is not influenced by the number of seizure occuring in given time period but to the
        % mean value the seizure characteristic in that time period. Imagine having just 2 seizures with high power in the first bin and 10 seizures
        % with low power in the second bin and average number of seizures with average value in the other bins. Fitting the bins gives a negative
        % trend whereas fitting the seizures per se would give a positive trend because the 10 seizures would beat the 2. Or imagine the first and
        % last bin containing different number of seizures each but the same value of the characteristic, e.g. above average power.
        fitTbl = fillInFit(fitTbl, 1, kchar, numChar, x, y, binWeights);
        xx{kchar} = x(:)';
        yy{kchar} = y(:)';
        fitSlope = fitTbl{1, kchar + 0*numChar};
        fitOffset = fitTbl{1, kchar + 3*numChar};
        xxFit{kchar} = [x(1), x(end)];
        yyFit{kchar} = xxFit{kchar}*fitSlope + fitOffset;
    end
    xxTbl = cell2table(xx, "VariableNames", stg.szCharToPlot);
    yyTbl = cell2table(yy, "VariableNames", stg.szCharToPlot);
    xxFitTbl = cell2table(xxFit, "VariableNames", stg.szCharToPlot);
    yyFitTbl = cell2table(yyFit, "VariableNames", stg.szCharToPlot);
end
function [fitTbl, eded, xxTbl, yyTbl, xxFitTbl, yyFitTbl] = fitSzCharCl(subjInfo, siCharTbl, cl)
    % subjInfo ... used to get date of birth (dob) and period of monitoring
    % siCharTbl .. signal characteristics, used to get info on signal dropouts
    % cl ......... to have info on the intracluster seizures
    % fitTbl ..... table containing the slopes, offsets and their confidence intervals
    % eded ....... matrix, each row is bin edges for one cluster
    % xx ......... cell array, each cell contains bin centers for one cluster and one characteristic
    % yy ......... cell array, each cell contains bin means for one cluster and one characteristic
    % xxFitTbl ... each cell contains fitted line at the onset and offset of the cluster for one characteristic
    % yyFitTbl ... each cell contains fitted line at the onset and offset of the cluster for one characteristic
    global stg
    numBins = 4;
    numCl = numel(cl);
    numChar = numel(stg.szCharToPlot);
    fitTbl = fitTblInit(numCl, stg.szCharToPlot);
    eded = NaN(numCl, numBins + 1); xx = cell(numCl, numChar); yy = cell(numCl, numChar); xxFit = cell(numCl, numChar); yyFit = cell(numCl, numChar);
    for kcl = 1 : numCl
        szCharTbl = cl(kcl).szCharTbl;
        firstSzDifficultMethod = szCharTbl.szOnsN(1) - subjInfo.dob;
        lastSzDifficultMethod = szCharTbl.szOnsN(end) - subjInfo.dob;
        firstSz = cl(kcl).szOnsN(1) - subjInfo.dob;
        lastSz = cl(kcl).szOnsN(end) - subjInfo.dob;
        if firstSz ~= firstSzDifficultMethod || lastSz ~= lastSzDifficultMethod
            error('_jk firstSz or lastSz is different when using DifficultMethod')
        end
        % Bar graph
        ed = linspace(firstSz, lastSz, numBins + 1);
        eded(kcl, :) = ed;
        for kchar = 1 : numChar
            x = NaN(1, numel(ed)-1); y = NaN(1, numel(ed)-1); binWeights = NaN(1, numel(ed)-1);
            for kb = 1 : numel(ed) - 1
                x(kb) = (ed(kb) + ed(kb+1))/2; % Bin centers are at the mean of the edges
                szSub = find((szCharTbl.szOnsN - subjInfo.dob) >= ed(kb), 1, 'first') :...
                        find((szCharTbl.szOnsN - subjInfo.dob) <=  ed(kb + 1), 1, 'last'); % Subscripts of seizures belonging to the bin
                % Bin weight is computed as the (number of valid points)/(expected number of points in a standard bin)
                [xOK, yOK] = getSzOkXY(subjInfo, siCharTbl);
                xOK = xOK(1 : 3 : end-2);
                yOK = yOK(1 : 3 : end-2);
                numAllDatapoints = (ed(3)-ed(2))/(xOK(2) - xOK(1)); % Number of datapoints in a standard bin (the first and last bin is usually shorter). Might be non-integer.
                if numAllDatapoints < numBins % If numAllDatapoints in a bin is < 1, then we will always get zero data points in the kb-th bin. The if numAllDatapoints is low, the estimate of binWeights is inaccurate.
                    numValidDatapoints = numAllDatapoints;
                    binWeights(kb) = 1;
                else
                    numValidDatapoints = sum(~isnan(yOK(xOK > ed(kb) & xOK < ed(kb+1))));
                    binWeights(kb) = numValidDatapoints/numAllDatapoints;
                end
                % Compute the value in the given bin
                if kchar == 1
                    y(kb) = numel(szSub)/(ed(kb+1) - ed(kb)); % Does not take into account recording dropouts
                    % % % % % y(kb) = numel(szSub)/(numValidDatapoints*stg.dpBinLenS);
                else
                    if isempty(szSub)
                        y(kb) = NaN;
                    else
                        y(kb) = mean(szCharTbl{szSub, stg.szCharToPlot(kchar)}, "omitmissing");
                    end
                end
            end
            % Line fit. We are fitting the bin values so that the fit is not influenced by the number of seizure occuring in given time period but to the
            % mean value the seizure characteristic in that time period. Imagine having just 2 seizures with high power in the first bin and 10 seizures
            % with low power in the second bin and average number of seizures with average value in the other bins. Fitting the bins gives a negative
            % trend whereas fitting the seizures per se would give a positive trend because the 10 seizures would beat the 2. Or imagine the first and
            % last bin containing different number of seizures each but the same value of the characteristic, e.g. above average power.
            fitTbl = fillInFit(fitTbl, kcl, kchar, numChar, x, y, binWeights);
            xx{kcl, kchar} = x(:)';
            yy{kcl, kchar} = y(:)';
            fitSlope = fitTbl{kcl, kchar + 0*numChar};
            fitOffset = fitTbl{kcl, kchar + 3*numChar};
            xxFit{kcl, kchar} = [x(1), x(end)];
            yyFit{kcl, kchar} = xxFit{kcl, kchar}*fitSlope + fitOffset;
        end
    end
    xxTbl = cell2table(xx, "VariableNames", stg.szCharToPlot);
    yyTbl = cell2table(yy, "VariableNames", stg.szCharToPlot);
    xxFitTbl = cell2table(xxFit, "VariableNames", stg.szCharToPlot);
    yyFitTbl = cell2table(yyFit, "VariableNames", stg.szCharToPlot);
end
function [cirTbl, ed, ppTbl, rrTbl, ppCirTbl, rrCirTbl] = fitSzCharCi(szCharTbl, siCharTbl)
    % szCharTbl .. seizure characteristics
    % siCharTbl .. signal characteristics, used to get info on signal dropouts
    % cirTbl ..... table with circular statistics cirMean (i.e. phase), lower and upper CI (do not know how it works) and resultant length (cirR)
    % ed ......... edges of bins in which seizure count or mean seizure characteristics are calculated in radians
    % ppTbl ...... centers of the bins (not used in any function so far)
    % rrTbl ...... seizure count or mean seizure characteristics in each bin
    % ppCirTbl ... phase of the start and end of the mean resultant vector
    % rrFitTbl ... radius of the start and end of the mean resultant vector
    binlenH = 3; % In hours
    global stg
    numChar = numel(stg.szCharToPlot);
    cirTbl = cirTblInit(1, stg.szCharToPlot);
    
    % Bar graph
    szTod = szCharTbl.szOnsN - floor(szCharTbl.szOnsN); % Seizure time of the day
    siTod = siCharTbl.tax - floor(siCharTbl.tax); % Time of the day of the signal characeteristics time axis
    ed = (0 : binlenH : 24)/24;
    pp = cell(1, numChar); rr = cell(1, numChar); ppCir = cell(1, numChar); rrCir = cell(1, numChar);
    for kchar = 1 : numChar
        p = NaN(1, numel(ed)-1); r = NaN(1, numel(ed)-1); binWeights = NaN(1, numel(ed)-1);
        for kb = 1 : numel(ed) - 1
            p(kb) = (ed(kb) + ed(kb+1))/2*2*pi; % Bin centers are at the mean of the edges
            szInd = szTod >= ed(kb)  &  szTod < ed(kb + 1); % Logical indices of seizures belonging to the bin
            siInd = siTod >= ed(kb)  &  siTod < ed(kb + 1); % Logical indices of signal char time axis belonging to the bin
            % Bin weight is computed as the (number of valid points)/(expected number of points in a standard bin)
            binWeights(kb) = sum(siCharTbl.szValidS(siInd))/3600; % How many hours there was a valid signal in total during this time of day across the whole recording
            % Compute the value in the given bin
            if kchar == 1
                % y(kb) = numel(szSub)/(ed(kb+1) - ed(kb)); % Does not take into account recording dropouts
                r(kb) = sum(szInd)/binWeights(kb); % Expected seizure frequency at this time of the day in the units of sz per hour
            else
                if ~any(szInd)
                    r(kb) = NaN;
                else
                    r(kb) = mean(szCharTbl.(stg.szCharToPlot(kchar))(szInd), "omitmissing");
                end
            end
        end
        % r = r - min(r);
        % r = r/max(r);
        r(isnan(r)) = mean(r, 2, 'omitmissing');
        szOccurrenceTF = true;
        cirTbl = fillInFitResultantVector(cirTbl, 1, kchar, numChar, p, r, szOccurrenceTF);
        pp{kchar} = p(:)';
        rr{kchar} = r(:)';
        ppCir{kchar} = [0, cirTbl{1, kchar + 0*numChar}];
        rrCir{kchar} = [0, cirTbl{1, kchar + 3*numChar}];
    end
    ppTbl = cell2table(pp, "VariableNames", stg.szCharToPlot);
    rrTbl = cell2table(rr, "VariableNames", stg.szCharToPlot);
    ppCirTbl = cell2table(ppCir, "VariableNames", stg.szCharToPlot);
    rrCirTbl = cell2table(rrCir, "VariableNames", stg.szCharToPlot);
end
function [stitched_wt, common_freqs, stitched_coi] = cwtPiecewiseGemini(tax, y, fs)
    % Gemini %
    t_full = tax';
    data = y';
    nan_indices = isnan(data);
    % Find where NaNs start and end (or where non-NaNs start and end)
    % This helps in identifying segments
    diff_nan = diff([1 nan_indices 1]); % Add sentinels at start/end
    start_non_nan = find(diff_nan == -1);
    end_non_nan = find(diff_nan == 1) - 1;
    num_segments = length(start_non_nan);
    % Initialize cells to store results for each segment
    all_cfs = cell(num_segments, 1);
    all_freqs = cell(num_segments, 1);
    all_segment_times = cell(num_segments, 1);
    all_coi = cell(num_segments, 1); % Cone of Influence
    for i = 1:num_segments
        segment_data = data(start_non_nan(i):end_non_nan(i));
        segment_time = t_full(start_non_nan(i):end_non_nan(i));
        if ~isempty(segment_data)
            % Perform CWT on the current segment
            voicesPerOctave = 12;
            minFreq = 1/32;
            maxFreq = 2;
            [wt, freqs, coi] = cwt(segment_data, "amor", fs, VoicesPerOctave=voicesPerOctave, FrequencyLimits=[minFreq maxFreq]); %#ok<ASGLU>
            fb = cwtfilterbank('SignalLength', 1e6, ...
                               'SamplingFrequency', fs, ...
                               'VoicesPerOctave', voicesPerOctave, ...
                               'FrequencyLimits', [minFreq maxFreq]);
            freqs = centerFrequencies(fb);
            num_scales = numel(freqs);
            if size(wt, 1) < num_scales % If fewer rows than expected
                temp_wt = NaN(num_scales, size(wt, 2)); % Create NaN matrix of desired size
                temp_wt(1:size(wt,1), :) = wt; % Copy the available coefficients
                wt = temp_wt; % Use the padded version
                temp_freqs = NaN(num_scales, 1); % Create NaN matrix of desired size
                temp_freqs(1:size(freqs,1), 1) = freqs;
                freqs = temp_freqs;

            end
            % Store the results
            all_wt{i} = wt;
            all_freqs{i} = freqs; % Frequencies should be the same for all segments if fs is constant
            all_segment_times{i} = segment_time;
            all_coi{i} = coi;
        else
            % Handle empty segments (e.g., if there were consecutive NaNs)
            all_wt{i} = [];
            all_freqs{i} = [];
            all_segment_times{i} = [];
            all_coi{i} = [];
        end
    end
    % Note: freqs should be identical for all non-empty segments, so you can pick one
    if any(~cellfun(@isempty, all_freqs))
        common_freqs = all_freqs{find(~cellfun(@isempty, all_freqs), 1)};
    else
        common_freqs = []; % No valid segments
    end
    % Determine the size of the combined CWT matrix
    num_scales = length(common_freqs);
    total_time_points = length(data);
    % Initialize the combined CWT coefficient matrix with NaNs
    stitched_wt = NaN(num_scales, total_time_points);
    stitched_coi = NaN(1, total_time_points); % Cone of Influence
    % Populate the stitched matrix
    for i = 1:num_segments
        if ~isempty(all_wt{i})
            segment_start_idx = start_non_nan(i);
            segment_end_idx = end_non_nan(i);
            stitched_wt(:, segment_start_idx:segment_end_idx) = all_wt{i};
            stitched_coi(segment_start_idx:segment_end_idx) = all_coi{i};
        end
    end
end
function [fitTbl, ed, xxSzTbl, yySzTbl, xxSiTbl, yySiTbl, xxFitTbl, yyFitTbl] = fitSaCharWh(subjInfo, szCharTbl, siCharTbl)
    % subjInfo ... used to get date of birth (dob) and period of monitoring
    % szCharTbl .. seizure characteristics
    % siCharTbl .. signal characteristics, used to get info on signal dropouts
    % fitTbl ..... table containing linear fits
    % ed ......... edges of bins in which seizure count or mean seizure characteristics are calculated
    % xxTbl ...... centers of the bins (not used in any function so far)
    % yyTbl ...... seizure count or mean seizure characteristics in each bin
    % xxFitTbl ... start and end of the fitted line for each characteristic
    % yyFitTbl ... start and end of the fitted line for each characteristic
    global stg
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    fitTbl = fitTblInit(1, [stg.szCharToPlot; stg.siCharToPlot]);
    
    %% Seizure characteristics
    firstSz = szCharTbl.szOnsN(1) - subjInfo.dob;
    lastSz = szCharTbl.szOnsN(end) - subjInfo.dob;
    % ed = floor(siCharTbl.tax(1) - subjInfo.anStartN) : stg.whWinLen : ceil((siCharTbl.tax(end) - subjInfo.anStartN)/7)*7; % From the start to
    % the end of the recording, might be biased since some subjects start seizing later which makes an impression of lower sz rate at the beginning.
    % And, also, having the end of the recording after a long inter-cluster period makes an impression of a decline in sz rate which may induce bias
    % as well (consieder perfectly regular seizures to see the problem). Therefore, I set it from the first to the last seizure.
    % THIS IS NOT ENTIRELY FAIR EITHER SINCE IT OVERESTIMATES THE RATE IN THE FIRST AND LAST BIN. But at least the first-order polynomial should be
    % correct.
    % ed = linspace(szCharTbl.szOnsN(1) - subjInfo.anStartN, szCharTbl.szOnsN(end) - subjInfo.anStartN, 11);
    ed = [firstSz, ceil(firstSz) : stg.whWinLen : floor(lastSz), lastSz];
    xx = cell(1, numSzChar); yy = cell(1, numSzChar); xxFit = cell(1, numChar); yyFit = cell(1, numChar);
    for kchar = 1 : numSzChar 
        x = NaN(1, numel(ed)-1); y = NaN(1, numel(ed)-1); binWeights = NaN(1, numel(ed)-1);
        for kb = 1 : numel(ed) - 1
            x(kb) = (ed(kb) + ed(kb+1))/2; % Bin centers are at the mean of the edges
            szSub = find((szCharTbl.szOnsN - subjInfo.dob) >= ed(kb), 1, 'first') :...
                    find((szCharTbl.szOnsN - subjInfo.dob) <=  ed(kb + 1), 1, 'last'); % Subscripts of seizures belonging to the bin
            % Bin weight is computed as the (number of valid points)/(expected number of points in a standard bin)
            [xOK, yOK] = getSzOkXY(subjInfo, siCharTbl);
            xOK = xOK(1 : 3 : end-2);
            yOK = yOK(1 : 3 : end-2);
            numAllDatapoints = (ed(3)-ed(2))/(xOK(2) - xOK(1)); % Number of datapoints in a standard bin (the first and last bin is usually shorter). Might be non-integer.
            numValidDatapoints = sum(~isnan(yOK(xOK > ed(kb) & xOK < ed(kb+1))));
            binWeights(kb) = numValidDatapoints/numAllDatapoints;
            % Compute the value in the given bin
            if kchar == 1
                % y(kb) = numel(szSub)/(ed(kb+1) - ed(kb)); % Does not take into account recording dropouts
                y(kb) = numel(szSub)/(numValidDatapoints*stg.dpBinLenS)*3600*24;
            else
                if isempty(szSub)
                    y(kb) = NaN;
                else
                    y(kb) = mean(szCharTbl.(stg.szCharToPlot(kchar))(szSub), "omitmissing");
                end
            end
        end
        
        % Line fit. We are fitting the bin values so that the fit is not influenced by the number of seizure occuring in given time period but to the
        % mean value the seizure characteristic in that time period. Imagine having just 2 seizures with high power in the first bin and 10 seizures
        % with low power in the second bin and average number of seizures with average value in the other bins. Fitting the bins gives a negative
        % trend whereas fitting the seizures per se would give a positive trend because the 10 seizures would beat the 2. Or imagine the first and
        % last bin containing different number of seizures each but the same value of the characteristic, e.g. above average power.
        fitTbl = fillInFit(fitTbl, 1, kchar, numChar, x, y, binWeights);
        xx{kchar} = x(:)';
        yy{kchar} = y(:)';
        fitSlope = fitTbl{1, kchar + 0*numChar};
        fitOffset = fitTbl{1, kchar + 3*numChar};
        xxFit{kchar} = [x(1), x(end)];
        yyFit{kchar} = xxFit{kchar}*fitSlope + fitOffset;
    end
    xxSzTbl = cell2table(xx, "VariableNames", stg.szCharToPlot);
    yySzTbl = cell2table(yy, "VariableNames", stg.szCharToPlot);

    %% Signal characteristics
    xx = cell(1, numSiChar); yy = cell(1, numSiChar);
    for kchar = 1 : numSiChar
        x = siCharTbl.tax - subjInfo.dob;
        y = siCharTbl.(stg.siCharToPlot(kchar));
        w = siCharTbl.(stg.siCharBinWeights(kchar));
        fitTbl = fillInFit(fitTbl, 1, kchar+numSzChar, numChar, x, y, w);
        xx{kchar} = x(:)';
        yy{kchar} = y(:)';
        fitSlope = fitTbl{1, kchar+numSzChar + 0*numChar};
        fitOffset = fitTbl{1, kchar+numSzChar + 3*numChar};
        xxFit{kchar+numSzChar} = [x(1), x(end)];
        yyFit{kchar+numSzChar} = xxFit{kchar}*fitSlope + fitOffset;
    end
    xxSiTbl = cell2table(xx, "VariableNames", stg.siCharToPlot);
    yySiTbl = cell2table(yy, "VariableNames", stg.siCharToPlot);
    xxFitTbl = cell2table(xxFit, "VariableNames", [stg.szCharToPlot; stg.siCharToPlot]);
    yyFitTbl = cell2table(yyFit, "VariableNames", [stg.szCharToPlot; stg.siCharToPlot]);
end
function [fitTbl, eded, xxSzTbl, yySzTbl, xxSiTbl, yySiTbl, xxFitTbl, yyFitTbl] = fitSaCharCl(subjInfo, siCharTbl, cl, fitPosition)
    % subjInfo ..... subject info, used to get date of birth (dob)
    % siCharTbl .... table with signal characteristics
    % cl ........... structure with all clusters of the subject
    % fitPosition .. "before" or "during" or "after"
    % fitTbl ....... table containing the slopes, offsets and their confidence intervals
    % xxTbl ........ fitted segments of time axis in seconds from the recording start
    % yyTbl ........ fitted segments of signal characteristics
    % xxFitTbl ..... start and end of the fitted segments of time axis in seconds from the recording start
    % yyFitTbl ..... start and end of the fitted line
    global stg
    numBins = 4;
    numCl = numel(cl);
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    fitTbl = fitTblInit(numCl, [stg.szCharToPlot; stg.siCharToPlot]);
    eded = NaN(numCl, numBins + 1);
    xxSi = cell(numCl, numSiChar); yySi = cell(numCl, numSiChar);
    xxSz = cell(numCl, numSzChar); yySz = cell(numCl, numSzChar);
    xxFit = cell(numCl, numChar); yyFit = cell(numCl, numChar);
    for kcl = 1 : numCl
        clDur = cl(kcl).szOnsN(end) - cl(kcl).szOnsN(1);
        switch fitPosition
            case "before"
                for kchar = 1 : numSzChar
                    xxFit{kcl, kchar} = [NaN, NaN];
                    yyFit{kcl, kchar} = [NaN, NaN];
                    fitTbl = fillInFit(fitTbl, kcl, kchar, numChar, NaN, NaN, NaN);
                end
                st = cl(kcl).szOnsN(1) - stg.clDurMultiple*clDur; % Multiple of the clusters' duration before the cluster onset
                if kcl > 1
                    clOffsets = sort(arrayfun(@(x) x.szOnsN(end), cl));
                    prevClEn = clOffsets(find(clOffsets < cl(kcl).szOnsN(1), 1, 'last'));
                    if st < prevClEn + stg.clDurMultiple*clDur
                        % warning('on', 'all');
                        % warning('_jk Non-nested cluster separated less than 2*stg.clDurMultiple*clDur from the previous one.')
                        % warning('off', 'all');
                        st = NaN;
                    end
                end
                st = max([st; siCharTbl.tax(1)], [], 'includemissing'); % If the the precluster start is before signal start, set it to the signal start
                en = cl(kcl).szOnsN(1);
                clTaxSub = find(siCharTbl.tax >= st, 1, 'first') - 1 : find(siCharTbl.tax > en, 1, 'first') - 1;
                if numel(clTaxSub) == 1 %#ok<ISCL>
                    clTaxSub = [clTaxSub-1, clTaxSub]; %#ok<AGROW> % If there is only one point, we violate the rule for determining the pre-cluster period and add one more point before.
                end
                clTaxSub = clTaxSub(clTaxSub > 0);
                if numel(clTaxSub) == 1 %#ok<ISCL>
                    clTaxSub = [];
                end
            case "during"
                szCharTbl = cl(kcl).szCharTbl;
                st = cl(kcl).szOnsN(1);
                en = cl(kcl).szOnsN(end);
                firstSz = st - subjInfo.dob;
                lastSz = en - subjInfo.dob;
                ed = linspace(firstSz, lastSz, numBins + 1);
                eded(kcl, :) = ed;
                for kchar = 1 : numSzChar
                    x = NaN(1, numel(ed)-1); y = NaN(1, numel(ed)-1); binWeights = NaN(1, numel(ed)-1);
                    for kb = 1 : numel(ed) - 1
                        x(kb) = (ed(kb) + ed(kb+1))/2; % Bin centers are at the mean of the edges
                        szSub = find((szCharTbl.szOnsN - subjInfo.dob) >= ed(kb), 1, 'first') :...
                                find((szCharTbl.szOnsN - subjInfo.dob) <=  ed(kb + 1), 1, 'last'); % Subscripts of seizures belonging to the bin
                        % Bin weight is computed as the (number of valid points)/(expected number of points in a standard bin)
                        [xOK, yOK] = getSzOkXY(subjInfo, siCharTbl);
                        xOK = xOK(1 : 3 : end-2);
                        yOK = yOK(1 : 3 : end-2);
                        numAllDatapoints = (ed(3)-ed(2))/(xOK(2) - xOK(1)); % Number of datapoints in a standard bin (the first and last bin is usually shorter). Might be non-integer.
                        if numAllDatapoints < 3 % If numAllDatapoints in a bin is < 1, then we will always get zero data points in the kb-th bin. The if numAllDatapoints is  low, the estimate of binWeights is inaccurate.
                            numValidDatapoints = numAllDatapoints; %#ok<NASGU>
                            binWeights(kb) = 1;
                        else
                            numValidDatapoints = sum(~isnan(yOK(xOK > ed(kb) & xOK < ed(kb+1))));
                            binWeights(kb) = numValidDatapoints/numAllDatapoints;
                        end
                        % Compute the value in the given bin
                        if kchar == 1
                            y(kb) = numel(szSub)/(ed(kb+1) - ed(kb)); % Does not take into account recording dropouts
                        else
                            if isempty(szSub)
                                y(kb) = NaN;
                            else
                                y(kb) = mean(szCharTbl{szSub, stg.szCharToPlot(kchar)}, "omitmissing");
                            end
                        end
                    end
                    % Line fit. We are fitting the bin values so that the fit is not influenced by the number of seizure occuring in given time period but to the
                    % mean value the seizure characteristic in that time period. Imagine having just 2 seizures with high power in the first bin and 10 seizures
                    % with low power in the second bin and average number of seizures with average value in the other bins. Fitting the bins gives a negative
                    % trend whereas fitting the seizures per se would give a positive trend because the 10 seizures would beat the 2. Or imagine the first and
                    % last bin containing different number of seizures each but the same value of the characteristic, e.g. above average power.
                    fitTbl = fillInFit(fitTbl, kcl, kchar, numChar, x, y, binWeights);
                    xxSz{kcl, kchar} = x(:)';
                    yySz{kcl, kchar} = y(:)';
                    fitSlope = fitTbl{kcl, kchar + 0*numChar};
                    fitOffset = fitTbl{kcl, kchar + 3*numChar};
                    xxFit{kcl, kchar} = [x(1), x(end)];
                    yyFit{kcl, kchar} = xxFit{kcl, kchar}*fitSlope + fitOffset;
                end
                clTaxSub = find(siCharTbl.tax >= st, 1, 'first') - 1 : find(siCharTbl.tax > en, 1, 'first');
            case "after"
                for kchar = 1 : numSzChar
                    xxFit{kcl, kchar} = [NaN, NaN];
                    yyFit{kcl, kchar} = [NaN, NaN];
                    fitTbl = fillInFit(fitTbl, kcl, kchar, numChar, NaN, NaN, NaN);
                end
                st = cl(kcl).szOnsN(end);
                en = cl(kcl).szOnsN(end) + stg.clDurMultiple*clDur; % Multiple of the clusters' duration after the cluster offset
                if kcl < numCl
                    clOnsets = sort(arrayfun(@(x) x.szOnsN(1), cl));
                    nextClSt = clOnsets(find(clOnsets > cl(kcl).szOnsN(end), 1, 'first'));

                    if nextClSt - en < stg.clDurMultiple*clDur
                        % warning('on', 'all');
                        % warning('_jk Non-nested cluster separated less than stg.clDurMultiple*clDur from the next one.')
                        % warning('off', 'all');
                        st = NaN;
                    end
                end
                en = min([en; siCharTbl.tax(end)]); % If the the postcluster end is after signal end, set it to the signal end
                clTaxSub = find(siCharTbl.tax >= st, 1, 'first') : find(siCharTbl.tax > en, 1, 'first');
        end
        if isempty(clTaxSub)
            numCols = width(fitTbl); % Get the number of columns in the table
            nanRow = array2table(NaN(1, numCols)); % Create a row of NaNs
            nanRow.Properties.VariableNames = fitTbl.Properties.VariableNames;
            fitTbl = [fitTbl; nanRow]; %#ok<AGROW> % Append the row to the table
            for kchar = 1 : numSiChar
                xxFit{kcl, kchar+numSzChar} = [NaN, NaN];
                yyFit{kcl, kchar+numSzChar} = [NaN, NaN];
            end
        else
            x = siCharTbl.tax(clTaxSub) - subjInfo.dob;
            for kchar = 1 : numSiChar
                y = siCharTbl.(stg.siCharToPlot(kchar))(clTaxSub);
                w = siCharTbl.(stg.siCharBinWeights(kchar))(clTaxSub);
                fitTbl = fillInFit(fitTbl, kcl, kchar+numSzChar, numChar, x, y, w);
                xxSi{kcl, kchar} = x(:)';
                yySi{kcl, kchar} = y(:)';
                fitSlope = fitTbl{kcl, kchar+numSzChar + 0*numChar};
                fitOffset = fitTbl{kcl, kchar+numSzChar + 3*numChar};
                xxFit{kcl, kchar+numSzChar} = [x(1), x(end)];
                yyFit{kcl, kchar+numSzChar} = xxFit{kcl, kchar+numSzChar}*fitSlope + fitOffset;
            end
        end
    end
    % % % % % % xxSz = cellfun(@(x) {x}, xxSz, 'UniformOutput', false)
    xxSzTbl = cell2table(xxSz, "VariableNames", stg.szCharToPlot); % Data are in cells because each row has different length
    % % % % % % yySz = cellfun(@(x) {x}, yySz, 'UniformOutput', false);
    yySzTbl = cell2table(yySz, "VariableNames", stg.szCharToPlot);
    xxSi = cellfun(@(x) {x}, xxSi, 'UniformOutput', false);
    xxSiTbl = cell2table(xxSi, "VariableNames", stg.siCharToPlot); % Data are in cells because each row has different length
    yySi = cellfun(@(x) {x}, yySi, 'UniformOutput', false);
    yySiTbl = cell2table(yySi, "VariableNames", stg.siCharToPlot);
    xxFitTbl = cell2table(xxFit, "VariableNames", [stg.szCharToPlot; stg.siCharToPlot]);
    yyFitTbl = cell2table(yyFit, "VariableNames", [stg.szCharToPlot; stg.siCharToPlot]);
end
function [cirTbl, edSz, ppSzTbl, rrSzTbl, edSi, ppSiTbl, rrSiTbl, ppCirTbl, rrCirTbl, rayCirTbl, omnCirTbl] = fitSaCharCi(szCharTbl, siCharTbl)
    % szCharTbl .. seizure characteristics
    % siCharTbl .. signal characteristics
    % cirTbl ..... table with circular statistics cirMean (i.e. phase), lower and upper CI (do not know how it works) and resultant length (cirR)
    % ed ......... edges of bins in which seizure count or mean seizure characteristics are calculated in radians
    % ppTbl ...... centers of the bins (not used in any function so far)
    % rrTbl ...... seizure rate or mean seizure characteristics in each bin
    % ppCirTbl ... phase of the start and end of the mean resultant vector
    % rrFitTbl ... radius of the start and end of the mean resultant vector
    global stg
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    cirTbl = cirTblInit(1, stg.saCharToPlot);
    
    %% Circular - seizures
    binlenH = 3; % In hours
    szTod = szCharTbl.szOnsN - floor(szCharTbl.szOnsN); % Seizure time of the day
    siTod = siCharTbl.tax - floor(siCharTbl.tax); % Time of the day of the signal characeteristics time axis
    edSz = (0 : binlenH : 24)/24;
    ppSz = cell(1, numSzChar); rrSz = cell(1, numSzChar); ppCir = cell(1, numChar); rrCir = cell(1, numChar); rayCir = cell(1, numChar); omnCir = cell(1, numChar);
    for kchar = 1 : numSzChar
        p = NaN(1, numel(edSz)-1); r = NaN(1, numel(edSz)-1); binWeights = NaN(1, numel(edSz)-1);
        for kb = 1 : numel(edSz) - 1
            p(kb) = (edSz(kb) + edSz(kb+1))/2*2*pi; % Bin centers are at the mean of the edges
            szInd = szTod >= edSz(kb)  &  szTod < edSz(kb + 1); % Logical indices of seizures belonging to the bin
            siInd = siTod >= edSz(kb)  &  siTod < edSz(kb + 1); % Logical indices of signal char time axis belonging to the bin
            % Bin weight is computed as the (number of valid points)/(expected number of points in a standard bin)
            binWeights(kb) = sum(siCharTbl.szValidS(siInd), 'omitmissing')/3600; % How many hours there was a valid signal in total during this time of day across the whole recording
            nominalBinWeight = ceil(siCharTbl.tax(end) - siCharTbl.tax(1))*binlenH; % How many hours there would be in each time-of-day bin if there were no dropouts.
            binWeightsNorm = binWeights/nominalBinWeight;
            % Compute the value in the given bin
            if kchar == 1
                % y(kb) = numel(szSub)/(ed(kb+1) - ed(kb)); % Does not take into account recording dropouts
                % r(kb) = sum(szInd)/binWeights(kb); % Expected seizure frequency at this time of the day in the units of sz per hour
                r(kb) = sum(szInd)/binWeightsNorm(kb);
                szOccurrenceTF = true;
            else
                if ~any(szInd)
                    r(kb) = NaN;
                else
                    r(kb) = mean(szCharTbl.(stg.szCharToPlot(kchar))(szInd), "omitmissing");
                end
                szOccurrenceTF = false;
            end
        end
        % r = r - min(r);
        % r = r/max(r);
        r(isnan(r)) = mean(r, 2, 'omitmissing'); % Replace bins containing NaN (no data were available) by mean of other bins
        cirTbl = fillInFitResultantVector(cirTbl, 1, kchar, numChar, p, r, szOccurrenceTF);
        ppSz{kchar} = p(:)';
        rrSz{kchar} = r(:)';
        ppCir{kchar} = [0, cirTbl{1, kchar + 0*numChar}];
        rrCir{kchar} = [0, cirTbl{1, kchar + 3*numChar}];
        if kchar == 1
            rayCir{kchar} = circ_rtest(2*pi*szTod);
            omnCir{kchar} = circ_otest(2*pi*szTod);
        else
            % % % rayCir{kchar} = circ_rtest(p', r'); % Probably incorrect use
            rayCir{kchar} = NaN;
            omnCir{kchar} = NaN; % Omnibus test not suitable for binned data
        end

    end
    
    %% Circular - signals
    binlenH = stg.dpBinLenS/3600; % In hours
    siTod = siCharTbl.tax - floor(siCharTbl.tax); % Time of the day of the signal characeteristics time axis
    edSi = (0 : binlenH : 24)/24; % From 0 to 1
    ppSi = cell(1, numSiChar); rrSi = cell(1, numSiChar);
    for kchar = 1 : numSiChar
        p = NaN(1, numel(edSi)-1); r = NaN(1, numel(edSi)-1);
        for kb = 1 : numel(edSi) - 1 % Over bins
            p(kb) = (edSi(kb) + edSi(kb+1))/2*2*pi; % Bin centers are at the mean of the edges
            siInd = siTod >= edSi(kb)  &  siTod < edSi(kb + 1); % Logical indices of signal char time axis belonging to the bin
            if ~any(siInd)
                r(kb) = NaN;
            else
                r(kb) = mean(siCharTbl.(stg.siCharToPlot(kchar))(siInd), "omitmissing");
            end
        end
        % r = r - min(r); % Increases the length of the resultant vector but does not change the direction
        rNorm = r/max(r); % Normalization so that maximum is 1
        szOccurrenceTF = false;
        cirTbl = fillInFitResultantVector(cirTbl, 1, kchar+numSzChar, numChar, p, r, szOccurrenceTF);
        ppSi{kchar} = p(:)';
        rrSi{kchar} = rNorm(:)';
        ppCir{kchar+numSzChar} = [0, cirTbl{1, kchar+numSzChar + 0*numChar}];
        rrCir{kchar+numSzChar} = [0, cirTbl{1, kchar+numSzChar + 3*numChar}];
        % % % rayCir{kchar+numSzChar} = circ_rtest(p', r'); % Probably incorrect use
        if stg.saCharToPlot(kchar+numSzChar) == "ied"
            r = r*binlenH/(stg.dpBinLenS/3600);
            rayCir{kchar+numSzChar} = circ_rtest(p', r'); % Should be OK since IED occurrence is in fact binned data.
        else
            rayCir{kchar+numSzChar} = NaN;
        end
        omnCir{kchar+numSzChar} = NaN; % Omnibus test not suitable for binned data
    end
    ppSzTbl = cell2table(ppSz, "VariableNames", stg.szCharToPlot);
    rrSzTbl = cell2table(rrSz, "VariableNames", stg.szCharToPlot);
    ppSiTbl = cell2table(ppSi, "VariableNames", stg.siCharToPlot);
    rrSiTbl = cell2table(rrSi, "VariableNames", stg.siCharToPlot);
    ppCirTbl = cell2table(ppCir, "VariableNames", stg.saCharToPlot);
    rrCirTbl = cell2table(rrCir, "VariableNames", stg.saCharToPlot);
    rayCirTbl = cell2table(rayCir, "VariableNames", stg.saCharToPlot);
    omnCirTbl = cell2table(omnCir, "VariableNames", stg.saCharToPlot);
end
function [fitTbl, xxTbl, yyTbl, xxFitTbl, yyFitTbl] = fitSiCharWh(subjInfo, siCharTbl)
    % subjInfo ..... subject info, used to get date of birth (dob)
    % siCharTbl .... table with signal characteristics
    % fitTbl ....... table containing the slopes, offsets and their confidence intervals
    % xxTbl ........ fitted segments of time axis
    % yyTbl ........ fitted segments of signal characteristics
    % xxFitTbl ..... start and end of the fitted line
    % yyFitTbl ..... start and end of the fitted line
    global stg
    numChar = numel(stg.siCharToPlot);
    fitTbl = fitTblInit(1, stg.siCharToPlot);
    xx = cell(1, numChar); yy = cell(1, numChar); xxFit = cell(1, numChar); yyFit = cell(1, numChar);
    for kchar = 1 : numChar
        x = siCharTbl.tax - subjInfo.dob;
        y = siCharTbl.(stg.siCharToPlot(kchar));
        w = siCharTbl.(stg.siCharBinWeights(kchar));
        fitTbl = fillInFit(fitTbl, 1, kchar, numChar, x, y, w);
        xx{kchar} = x(:)';
        yy{kchar} = y(:)';
        fitSlope = fitTbl{1, kchar + 0*numChar};
        fitOffset = fitTbl{1, kchar + 3*numChar};
        xxFit{kchar} = [x(1), x(end)];
        yyFit{kchar} = xxFit{kchar}*fitSlope + fitOffset;
    end
    xxTbl = cell2table(xx, "VariableNames", stg.siCharToPlot);
    yyTbl = cell2table(yy, "VariableNames", stg.siCharToPlot);
    xxFitTbl = cell2table(xxFit, "VariableNames", stg.siCharToPlot);
    yyFitTbl = cell2table(yyFit, "VariableNames", stg.siCharToPlot);
end
function [fitTbl, xxTbl, yyTbl, xxFitTbl, yyFitTbl] = fitSiCharCl(subjInfo, siCharTbl, cl, fitPosition)
    % subjInfo ..... subject info, used to get date of birth (dob)
    % siCharTbl .... table with signal characteristics
    % cl ........... structure with all clusters of the subject
    % fitPosition .. "before" or "during" or "after"
    % fitTbl ....... table containing the slopes, offsets and their confidence intervals
    % xxTbl ........ fitted segments of time axis in seconds from the recording start
    % yyTbl ........ fitted segments of signal characteristics
    % xxFitTbl ..... start and end of the fitted segments of time axis in seconds from the recording start
    % yyFitTbl ..... start and end of the fitted line
    global stg
    numCl = numel(cl);
    numChar = numel(stg.siCharToPlot);
    fitTbl = fitTblInit(numCl, stg.siCharToPlot);
    xx = cell(numCl, numChar); yy = cell(numCl, numChar); xxFit = cell(numCl, numChar); yyFit = cell(numCl, numChar);
    for kcl = 1 : numCl
        clDur = cl(kcl).szOnsN(end) - cl(kcl).szOnsN(1);
        switch fitPosition
            case "before"
                st = cl(kcl).szOnsN(1) - stg.clDurMultiple*clDur; % Multiple of the clusters' duration before the cluster onset
                if kcl > 1
                    clOffsets = sort(arrayfun(@(x) x.szOnsN(end), cl));
                    prevClEn = clOffsets(find(clOffsets < cl(kcl).szOnsN(1), 1, 'last'));
                    if st < prevClEn + stg.clDurMultiple*clDur
                        % warning('on', 'all');
                        % warning('_jk Non-nested cluster separated less than 2*stg.clDurMultiple*clDur from the previous one.')
                        % warning('off', 'all');
                        st = NaN;
                    end
                end
                st = max([st; siCharTbl.tax(1)], [], 'includemissing'); % If the the precluster start is before signal start, set it to the signal start
                en = cl(kcl).szOnsN(1);
                clTaxSub = find(siCharTbl.tax >= st, 1, 'first') - 1 : find(siCharTbl.tax > en, 1, 'first') - 1;
                if isscalar(clTaxSub)
                    clTaxSub = [clTaxSub-1, clTaxSub]; %#ok<AGROW> % If there is only one point, we violate the rule for determining the pre-cluster period and add one more point before.
                end
                clTaxSub = clTaxSub(clTaxSub > 0);
                if isscalar(clTaxSub)
                    clTaxSub = [];
                end
            case "during"
                st = cl(kcl).szOnsN(1);
                en = cl(kcl).szOnsN(end);
                clTaxSub = find(siCharTbl.tax >= st, 1, 'first') - 1 : find(siCharTbl.tax > en, 1, 'first');
            case "after"
                st = cl(kcl).szOnsN(end);
                en = cl(kcl).szOnsN(end) + stg.clDurMultiple*clDur; % Multiple of the clusters' duration after the cluster offset
                if kcl < numCl
                    clOnsets = sort(arrayfun(@(x) x.szOnsN(1), cl));
                    nextClSt = clOnsets(find(clOnsets > cl(kcl).szOnsN(end), 1, 'first'));
                    if nextClSt - en < stg.clDurMultiple*clDur
                        % warning('on', 'all');
                        % warning('_jk Non-nested cluster separated less than stg.clDurMultiple*clDur from the next one.')
                        % warning('off', 'all');
                        st = NaN;
                    end
                end
                en = min([en; siCharTbl.tax(end)]); % If the the postcluster end is after signal end, set it to the signal end
                clTaxSub = find(siCharTbl.tax >= st, 1, 'first') : find(siCharTbl.tax > en, 1, 'first');
        end
        if isempty(clTaxSub)
            numCols = width(fitTbl); % Get the number of columns in the table
            nanRow = array2table(NaN(1, numCols)); % Create a row of NaNs
            nanRow.Properties.VariableNames = fitTbl.Properties.VariableNames;
            fitTbl = [fitTbl; nanRow]; %#ok<AGROW> % Append the row to the table
            for kchar = 1 : numChar
                xxFit{kcl, kchar} = [NaN, NaN];
                yyFit{kcl, kchar} = [NaN, NaN];
            end
        else
            x = siCharTbl.tax(clTaxSub) - subjInfo.dob;
            for kchar = 1 : numChar
                y = siCharTbl.(stg.siCharToPlot(kchar))(clTaxSub);
                w = siCharTbl.(stg.siCharBinWeights(kchar))(clTaxSub);
                fitTbl = fillInFit(fitTbl, kcl, kchar, numChar, x, y, w);
                xx{kcl, kchar} = x(:)';
                yy{kcl, kchar} = y(:)';
                fitSlope = fitTbl{kcl, kchar + 0*numChar};
                fitOffset = fitTbl{kcl, kchar + 3*numChar};
                xxFit{kcl, kchar} = [x(1), x(end)];
                yyFit{kcl, kchar} = xxFit{kcl, kchar}*fitSlope + fitOffset;
            end
        end
    end
    xx = cellfun(@(x) {x}, xx, 'UniformOutput', false);
    xxTbl = cell2table(xx, "VariableNames", stg.siCharToPlot); % Data are in cells because each row has different length
    yy = cellfun(@(x) {x}, yy, 'UniformOutput', false);
    yyTbl = cell2table(yy, "VariableNames", stg.siCharToPlot);
    xxFitTbl = cell2table(xxFit, "VariableNames", stg.siCharToPlot);
    yyFitTbl = cell2table(yyFit, "VariableNames", stg.siCharToPlot);
end
function [cirTbl, ed, ppTbl, rrTbl, ppCirTbl, rrCirTbl] = fitSiCharCi(siCharTbl)
    % szCharTbl .. seizure characteristics
    % siCharTbl .. signal characteristics, used to get info on signal dropouts
    % cirTbl ..... table with circular statistics cirMean (i.e. phase), lower and upper CI (do not know how it works) and resultant length (cirR)
    % ed ......... edges of bins in which seizure count or mean seizure characteristics are calculated in radians
    % ppTbl ...... centers of the bins (not used in any function so far)
    % rrTbl ...... seizure count or mean seizure characteristics in each bin
    % ppCirTbl ... phase of the start and end of the mean resultant vector
    % rrFitTbl ... radius of the start and end of the mean resultant vector
    global stg
    binlenH = 1; % In hours
    numChar = numel(stg.siCharToPlot);
    cirTbl = cirTblInit(1, stg.siCharToPlot);
    siTod = siCharTbl.tax - floor(siCharTbl.tax); % Time of the day of the signal characeteristics time axis
    ed = (0 : binlenH : 24)/24;
    pp = cell(1, numChar); rr = cell(1, numChar); ppCir = cell(1, numChar); rrCir = cell(1, numChar);
    for kchar = 1 : numChar
        p = NaN(1, numel(ed)-1); r = NaN(1, numel(ed)-1);
        for kb = 1 : numel(ed) - 1 % Over bins
            p(kb) = (ed(kb) + ed(kb+1))/2*2*pi; % Bin centers are at the mean of the edges
            siInd = siTod >= ed(kb)  &  siTod < ed(kb + 1); % Logical indices of signal char time axis belonging to the bin
            if ~any(siInd)
                r(kb) = NaN;
            else
                r(kb) = mean(siCharTbl.(stg.siCharToPlot(kchar))(siInd), "omitmissing");
            end
        end
        r = r - min(r); % Increases the length of the resultant vector but does not change the direction
        r = r/max(r); % Normalization so that maximum is 1
        szOccurrenceTF = false;
        cirTbl = fillInFitResultantVector(cirTbl, 1, kchar, numChar, p, r, szOccurrenceTF);
        pp{kchar} = p(:)';
        rr{kchar} = r(:)';
        ppCir{kchar} = [0, cirTbl{1, kchar + 0*numChar}];
        rrCir{kchar} = [0, cirTbl{1, kchar + 3*numChar}];
    end
    ppTbl = cell2table(pp, "VariableNames", stg.siCharToPlot);
    rrTbl = cell2table(rr, "VariableNames", stg.siCharToPlot);
    ppCirTbl = cell2table(ppCir, "VariableNames", stg.siCharToPlot);
    rrCirTbl = cell2table(rrCir, "VariableNames", stg.siCharToPlot);
end
function [fitTbl, xxTbl, yyTbl, xxFitTbl, yyFitTbl, xxOtherTbl, yyOtherTbl] = aroundSzVsOther(subjInfo, szCharTbl, siCharTbl, fitPosition)
    % subjInfo ..... subject info, used to get date of birth (dob)
    % szCharTbl .... table with seizure characteristics
    % siCharTbl .... table with signal characteristics
    % fitPosition .. "before" or "after"
    % fitTbl ....... table containing the slopes, offsets and their confidence intervals
    % xxTbl ........ fitted segments of time axis in seconds from the recording start
    % yyTbl ........ fitted segments of signal characteristics
    % xxFitTbl ..... start and end of the fitted segments of time axis in seconds from the recording start
    % yyFitTbl ..... start and end of the fitted line
    global stg
    % Determine which seizures will be used
    ons = szCharTbl.szOnsN;
    offs = szCharTbl.szOnsN + szCharTbl.szDurS/3600/24;
    switch fitPosition
        case "before"
            toKeepInd = [true; ons(2:end) - offs(1:end-1) >= stg.leadSzTimeS/3600/24]; % We keep the first sz and seizures that are at least leadTimeS after the previous seizure offset
            szCharTbl = szCharTbl(toKeepInd, :); % Keeping szCharTbl might be useful for more filtering, e.g. only seizures of certain duration or severity
            toKeepInd = ~(szCharTbl.szOnsN - stg.fitSzDurS/3600/24 < siCharTbl.tax(1) + stg.leadSzTimeS/3600/24);
            szCharTbl = szCharTbl(toKeepInd, :); % Keeping szCharTbl might be useful for more filtering, e.g. only seizures of certain duration or severity
        case "after"
            toKeepInd = [ons(2:end) - offs(1:end-1) >= stg.leadSzTimeS/3600/24; true]; % We keep the seizures that are followed by an ISI at least leadTimeS and the last seizure
            szCharTbl = szCharTbl(toKeepInd, :); % Keeping szCharTbl might be useful for more filtering, e.g. only seizures of certain duration or severity
            toKeepInd = ~(szCharTbl.szOnsN + szCharTbl.szDurS/3600/24 + stg.fitSzDurS/3600/24 > siCharTbl.tax(end)); % If the fit goes beyond the recording
            szCharTbl = szCharTbl(toKeepInd, :); % Keeping szCharTbl might be useful for more filtering, e.g. only seizures of certain duration or severity
    end
    % Initialize variables
    numSz = size(szCharTbl, 1);
    numChar = numel(stg.siCharToPlot);
    fitTbl = fitTblInit(numSz, stg.siCharToPlot);
    taxSub = cell(numSz, numChar);
    xx = cell(numSz, numChar); yy = cell(numSz, numChar);
    xxFit = cell(numSz, numChar); yyFit = cell(numSz, numChar);
    xxOther = cell(1, numChar); yyOther = cell(1, numChar);
    if numSz == 0
        for kchar = 1 : numChar
            len = stg.fitSzDurS/stg.dpBinLenS + 1;
            x = NaN(len, 1);
            y = NaN(len, 1);
            fitTbl = fillInFit(fitTbl, 1, kchar, numChar, x, y);
            xx{1, kchar} = x(:)';
            yy{1, kchar} = y(:)';
            xxOther{1, kchar} = siCharTbl.tax - subjInfo.dob;
            yyOther{1, kchar} = siCharTbl.(stg.siCharToPlot(kchar));
            fitSlope = fitTbl{1, kchar + 0*numChar};
            fitOffset = fitTbl{1, kchar + 3*numChar};
            xxFit{1, kchar} = [x(1), x(end)];
            yyFit{1, kchar} = xxFit{1, kchar}*fitSlope + fitOffset;
        end
    end
    for ksz = 1 : numSz
        for kchar = 1 : numChar
            switch fitPosition
                case "before"
                    st = szCharTbl.szOnsN(ksz) - stg.fitSzDurS/3600/24;
                    en = szCharTbl.szOnsN(ksz);
                    taxSub{ksz, kchar} = find(siCharTbl.tax >= st, 1, 'first') - 1 : find(siCharTbl.tax > en, 1, 'first') - 1;
                case "after"
                    st = szCharTbl.szOnsN(ksz) + szCharTbl.szDurS(ksz)/3600/24;
                    en = szCharTbl.szOnsN(ksz) + szCharTbl.szDurS(ksz)/3600/24 + stg.fitSzDurS/3600/24;
                    taxSub{ksz, kchar} = find(siCharTbl.tax > st, 1, 'first') : find(siCharTbl.tax >= en, 1, 'first');
            end
        end
    end
    sigLen = min(cellfun(@(x)numel(x), taxSub), [], 'all'); % Get the length of the shortest signal snippet
    for k = 1 : numel(taxSub)
        tSub = taxSub{k}(:)';
        taxSub{k} = tSub(1 + (size(tSub, 2)-sigLen) : end); % Remove a few samples at the beginning of the snippet if needed so that all snippets will be the same length
    end
    for ksz = 1 : numSz
        for kchar = 1 : numChar
            x = siCharTbl.tax(taxSub{ksz, kchar}') - subjInfo.dob;
            y = siCharTbl.(stg.siCharToPlot(kchar))(taxSub{ksz, kchar}');
            % w = siCharTbl.(stg.siCharBinWeights(kchar))(taxSub{ksz, kchar}');
            if sum(~isnan(y))/numel(y) < 0.5 % If less than half of the data points are non-NaN, make them NaN all so that they are not further considered
                y = NaN(size(y));
            end
            fitTbl = fillInFit(fitTbl, ksz, kchar, numChar, x, y);
            xx{ksz, kchar} = x(:)';
            yy{ksz, kchar} = y(:)';
            fitSlope = fitTbl{ksz, kchar + 0*numChar};
            fitOffset = fitTbl{ksz, kchar + 3*numChar};
            xxFit{ksz, kchar} = [x(1), x(end)];
            yyFit{ksz, kchar} = xxFit{ksz, kchar}*fitSlope + fitOffset;
        end
    end
    for kchar = 1 : numChar
        taxSubOneColumn = cell2mat(taxSub(:, kchar)'); % Select only the desired char, put all cells in one row of cells, convert to matrix (one row)
        taxSubOther = setdiff(1 : numel(siCharTbl.tax), taxSubOneColumn)'; % Get the subscripts into tax which are not contained in the taxSubOneColumn and transpose to get it in a column
        xxOther{1, kchar} = siCharTbl.tax(taxSubOther) - subjInfo.dob; % Probably will not be used but it is here for symmetry and cleanliness
        yyOther{1, kchar} = siCharTbl.(stg.siCharToPlot(kchar))(taxSubOther);
    end
    xxTbl = cell2table(xx, "VariableNames", stg.siCharToPlot);
    yyTbl = cell2table(yy, "VariableNames", stg.siCharToPlot);
    xxOtherTbl = cell2table(xxOther, "VariableNames", stg.siCharToPlot);
    yyOtherTbl = cell2table(yyOther, "VariableNames", stg.siCharToPlot);
    xxFitTbl = cell2table(xxFit, "VariableNames", stg.siCharToPlot);
    yyFitTbl = cell2table(yyFit, "VariableNames", stg.siCharToPlot);
end
function [fitTbl, xxTbl, yyTbl, xxFitTbl, yyFitTbl] = fitSiCharSz(subjInfo, szCharTbl, siCharTbl, fitPosition)
    % subjInfo ..... subject info, used to get date of birth (dob)
    % szCharTbl .... table with seizure characteristics
    % siCharTbl .... table with signal characteristics
    % fitPosition .. "before" or "after"
    % fitTbl ....... table containing the slopes, offsets and their confidence intervals
    % xxTbl ........ fitted segments of time axis in seconds from the recording start
    % yyTbl ........ fitted segments of signal characteristics
    % xxFitTbl ..... start and end of the fitted segments of time axis in seconds from the recording start
    % yyFitTbl ..... start and end of the fitted line
    global stg
    ons = szCharTbl.szOnsN;
    offs = szCharTbl.szOnsN + szCharTbl.szDurS/3600/24;
    switch fitPosition
        case "before"
            toKeepInd = [true; ons(2:end) - offs(1:end-1) >= stg.leadSzTimeS/3600/24]; % We keep the first sz and seizures that are at least leadTimeS after the previous seizure offset
            szCharTbl = szCharTbl(toKeepInd, :); % Keeping szCharTbl might be useful for more filtering, e.g. only seizures of certain duration or severity
            toKeepInd = ~(szCharTbl.szOnsN - stg.fitSzDurS/3600/24 < siCharTbl.tax(1) + stg.leadSzTimeS/3600/24);
            szCharTbl = szCharTbl(toKeepInd, :); % Keeping szCharTbl might be useful for more filtering, e.g. only seizures of certain duration or severity
        case "after"
            toKeepInd = [ons(2:end) - offs(1:end-1) >= stg.leadSzTimeS/3600/24; true]; % We keep the seizures that are followed by an ISI at least leadTimeS and the last seizure
            szCharTbl = szCharTbl(toKeepInd, :); % Keeping szCharTbl might be useful for more filtering, e.g. only seizures of certain duration or severity
            toKeepInd = ~(szCharTbl.szOnsN + szCharTbl.szDurS/3600/24 + stg.fitSzDurS/3600/24 > siCharTbl.tax(end)); % If the fit goes beyond the recording
            szCharTbl = szCharTbl(toKeepInd, :); % Keeping szCharTbl might be useful for more filtering, e.g. only seizures of certain duration or severity
    end
    numSz = size(szCharTbl, 1);
    numChar = numel(stg.siCharToPlot);
    fitTbl = fitTblInit(numSz, stg.siCharToPlot);
    taxSub = cell(numSz, numChar); xx = cell(numSz, numChar); yy = cell(numSz, numChar); xxFit = cell(numSz, numChar); yyFit = cell(numSz, numChar);
    if numSz == 0
        for kchar = 1 : numChar
            len = stg.fitSzDurS/stg.dpBinLenS + 1;
            x = NaN(len, 1);
            y = NaN(len, 1);
            fitTbl = fillInFit(fitTbl, 1, kchar, numChar, x, y);
            xx{1, kchar} = x(:)';
            yy{1, kchar} = y(:)';
            fitSlope = fitTbl{1, kchar + 0*numChar};
            fitOffset = fitTbl{1, kchar + 3*numChar};
            xxFit{1, kchar} = [x(1), x(end)];
            yyFit{1, kchar} = xxFit{1, kchar}*fitSlope + fitOffset;
        end
    end
    for ksz = 1 : numSz
        for kchar = 1 : numChar
            switch fitPosition
                case "before"
                    st = szCharTbl.szOnsN(ksz) - stg.fitSzDurS/3600/24;
                    en = szCharTbl.szOnsN(ksz);
                    taxSub{ksz, kchar} = find(siCharTbl.tax >= st, 1, 'first') - 1 : find(siCharTbl.tax > en, 1, 'first') - 1;
                case "after"
                    st = szCharTbl.szOnsN(ksz) + szCharTbl.szDurS(ksz)/3600/24;
                    en = szCharTbl.szOnsN(ksz) + szCharTbl.szDurS(ksz)/3600/24 + stg.fitSzDurS/3600/24;
                    taxSub{ksz, kchar} = find(siCharTbl.tax > st, 1, 'first') : find(siCharTbl.tax >= en, 1, 'first');
            end
        end
    end
    sigLen = min(cellfun(@(x)numel(x), taxSub), [], 'all'); % Get the length of the shortest signal snippet
    for k = 1 : numel(taxSub)
        tSub = taxSub{k}(:)';
        taxSub{k} = tSub(1 + (size(tSub, 2)-sigLen) : end);
    end
    for ksz = 1 : numSz
        for kchar = 1 : numChar
            x = siCharTbl.tax(taxSub{ksz, kchar}') - subjInfo.dob;
            y = siCharTbl.(stg.siCharToPlot(kchar))(taxSub{ksz, kchar}');
            % w = siCharTbl.(stg.siCharBinWeights(kchar))(taxSub{ksz, kchar}');
            if sum(~isnan(y))/numel(y) < 0.5 % If less than half of the data points are non-NaN, make them NaN all so that they are not further considered
                y = NaN(size(y));
            end
            fitTbl = fillInFit(fitTbl, ksz, kchar, numChar, x, y);
            xx{ksz, kchar} = x(:)';
            yy{ksz, kchar} = y(:)';
            fitSlope = fitTbl{ksz, kchar + 0*numChar};
            fitOffset = fitTbl{ksz, kchar + 3*numChar};
            xxFit{ksz, kchar} = [x(1), x(end)];
            yyFit{ksz, kchar} = xxFit{ksz, kchar}*fitSlope + fitOffset;
        end
    end
    xxTbl = cell2table(xx, "VariableNames", stg.siCharToPlot);
    yyTbl = cell2table(yy, "VariableNames", stg.siCharToPlot);
    xxFitTbl = cell2table(xxFit, "VariableNames", stg.siCharToPlot);
    yyFitTbl = cell2table(yyFit, "VariableNames", stg.siCharToPlot);
end
function [resultStr, m, sem, sdv, med, iql, sigRanP] = getBasicStats(x, xNm)
    % Based on Reviewer #2 I added the option of reporting standard deviation instead of SEM.
    global stg
    m = mean(x, 1, 'omitmissing');
    sem = std(x, 1, 'omitmissing')/sqrt(sum(~isnan(x)));
    sdv = std(x, 1, 'omitmissing');
    med = median(x, 1, 'omitmissing');
    iql = iqr(x, 1);
    if any(~isnan(x))
        sigRanP = signrank(x);
    else
        sigRanP = NaN;
    end
    if stg.sdTF
        resultStr = ...
            [repelem(' ', 1, 13-numel(xNm)), xNm, ' = ', num2str(m, '%.2g'), '±', num2str(sdv, '%.2g'), ' (', num2str(med, '%.2g'), '±', num2str(iql, '%.2g'), ')', 10, ...
             '            p = ', num2str(sigRanP, '%.2g')];
    else
        resultStr = ...
            [repelem(' ', 1, 13-numel(xNm)), xNm, ' = ', num2str(m, '%.2g'), '±', num2str(sem, '%.2g'), ' (', num2str(med, '%.2g'), '±', num2str(iql, '%.2g'), ')', 10, ...
             '            p = ', num2str(sigRanP, '%.2g')];
    end
end
function [fecoeff, amplitudes, decayFactors, tauH, fe, gofe] = fitExponentials(x, y, baseline, varargin)
    global stg
    if isempty(varargin)
        stSub = max(2, (60/stg.dpBinLenS + 1)); % Start subscript
    else
        stSub = varargin{1};
    end
    e = ''; % The expression to be fitted
    stPt = NaN(1, 2*stg.simSiCharFltOrder);
    if stg.fitSzDurS < 6*3600
        decayRateMult = 10;
    else
        decayRateMult = 1;
    end
    AstPt = max(y);
    for k = 1 : stg.simSiCharFltOrder
        if rem(stg.simSiCharFltOrder - k, 2) == 0
            % % % sgn(k) = 1;
            % % % s = ['+', char(97+(k-1)), '*'];
            stPt(2*(k-1) + 1) = AstPt;
        else
            % % % sgn(k) = -1;
            % % % s = ['+', char(97+(k-1)), '*'];
            stPt(2*(k-1) + 1) = -AstPt;
        end
        s = ['+', char(97+2*(k-1)), '*'];
        e = [e, s, 'exp(', 98+2*(k-1), '*x)']; %#ok<AGROW>
    
        stPt(2*k) = -decayRateMult/(5^(k-1)); % Starting point. First exponential is the fastest.
    end
    e = [e, ''];
    % Here we substitute first few (supposedly invalid) samples by either baseline or the first valid value (at position stSub)
    if rem(stg.simSiCharFltOrder, 2) == 0 % If even number of exponential, start at zero
        xf = x;
        yf = [zeros(stSub-1, 1); y(stSub:end)' - baseline];
    else % If odd number of exponentials, start at the top
        xf = x;
        yf = [y(stSub)*ones(stSub-1, 1)  - baseline; y(stSub:end)' - baseline];
    end
    % Fitting proper
    [fe, gofe] = fit(xf', yf, e,...
        'StartPoint', stPt, 'Lower', stPt./20.^sign(stPt), 'Upper', stPt.*20.^sign(stPt)); % Fit Exponential, Goodness Of Fit Exponential
    fecoeff = [coeffvalues(fe), gofe.adjrsquare, gofe.rmse];
    amplitudes = fecoeff(1 : 2 : end - 2); % Last two are R^2 and RMSE and not coefficients of the exponentials
    decayFactors = exp(stg.dpBinLenS/3600*fecoeff(2 : 2 : end - 2)); % Convert rates of decay to decay factors (essentially poles of the filter)
    tauH = stg.dpBinLenS/3600./(-log(decayFactors)); % The same as below
    % % % tauH = 1./-fecoeff(2 : 2 : end - 2); % The same as above. The fecoeff(2 : 2 : end - 2) are in the same units as x (probably hours)
end
function [b, a] = designFiltFromExponentials(decayFactors, amplitudes)
    pol = zeros(numel(amplitudes), numel(amplitudes)); % Polynomials for the numerator
    for kp = 1 : numel(amplitudes)
        allExceptThis = setdiff(1 : numel(decayFactors), kp); % Used for partial fractions
        decayFactorsToUse = decayFactors(allExceptThis);
        pol(kp, :) = amplitudes(kp)*poly(decayFactorsToUse);
    end
    b = sum(pol, 1);
    a = poly(decayFactors); % Get denominator from the poles 
end
function [fpcoeff, fp, gofp] = fitPowerLaw(x, y, baseline)
    global stg
    stSub = max(2, round((60/stg.dpBinLenS + 1))); % Start subscript
    [fp, gofp] = fit(x(stSub:end)', y(stSub:end)' - baseline, 'power1');
    fpcoeff = [coeffvalues(fp), gofp.adjrsquare, gofp.rmse];
end
function siCharTbl = addSimulatedData(siCharTbl, siCharCur, siCharCurPop, ksubj)
    global stg
    sz = siCharTbl.sz/24;
    sz(isnan(sz)) = 0;
    % Subject specific IED prediction filter
    for kchar = 1 : numel(stg.siCharToPlot)
        if any(isnan(siCharCur.expFltB{ksubj, kchar})) % If filter coefficients contain NaN, insert NaN in the output
            sim = NaN(size(siCharTbl, 1), 1);
            simPop = NaN(size(siCharTbl, 1), 1);
        else
            % Subject-specific simulation
            ori = siCharTbl.(stg.siCharToPlot(kchar));
            oriSm = filter(1/stg.simMovAveLen*ones(1, stg.simMovAveLen), 1, ori); % Smoothed
            maOriSm = max(oriSm, [], 1, 'omitmissing');
            sim = filter(siCharCur.expFltB{ksubj, kchar}, siCharCur.expFltA{ksubj, kchar}, sz);
            simSm = filter(1/stg.simMovAveLen*ones(1, stg.simMovAveLen), 1, sim); % Smoothed
            maSimSm = max(abs(simSm), [], 1, 'omitmissing'); % The abs helps when the whole sim signal is negative. This rarely happens in case of erroneous fitting. I admit it is a sloppy solution. But the fit is inaccurate anyway.
            bsl = siCharCur.baseline(ksubj, kchar);
            sim = sim/maSimSm*(maOriSm - bsl) + bsl; % We normalize using the smoothed signals but save the non-smoothed version in siCharTbl. Maybe weird, I admit.
            % Population-based simulation
            simPop = filter(siCharCurPop.expFltB{1, kchar}, siCharCurPop.expFltA{1, kchar}, sz);
            simPopSm = filter(1/stg.simMovAveLen*ones(1, stg.simMovAveLen), 1, simPop); % Smoothed
            maSimPopSm = max(simPopSm, [], 1, 'omitmissing');
            bslPop = siCharCurPop.baseline(1, kchar); % Here, the baseline is in the units of the filter coefficients which may be normalized to gain 1
            simPop = simPop + bslPop;
            simPop = simPop/maSimPopSm*maOriSm;
        end
        t = table(sim);
        t.Properties.VariableNames = stg.siCharToPlot(kchar) + "Sim";
        tPop = table(simPop);
        tPop.Properties.VariableNames = stg.siCharToPlot(kchar) + "SimPop";
        siCharTbl = [siCharTbl, t, tPop]; %#ok<AGROW>
    end
end
function [rmsnrmse, rho, pval, rmsnrmsePop, rhoPop, pvalPop] = simulatedDataSimilarity(siCharTbl)
% Try limitting sz to e.g. 2
    global stg
    rmsnrmse = NaN(1, numel(stg.siCharToPlot));
    rmsnrmsePop = NaN(1, numel(stg.siCharToPlot));
    rho = NaN(1, numel(stg.siCharToPlot));
    rhoPop = NaN(1, numel(stg.siCharToPlot));
    pval = NaN(1, numel(stg.siCharToPlot));
    pvalPop = NaN(1, numel(stg.siCharToPlot));
    for kchar = 1 : numel(stg.siCharToPlot)
        % Get data
        ori = fillmissing(siCharTbl.(stg.siCharToPlot(kchar)), 'linear', 1, 'MaxGap', ceil(4*3600/stg.dpBinLenS));
        sim = fillmissing(siCharTbl.(stg.siCharToPlot(kchar) + "Sim"), 'linear', 1, 'MaxGap', ceil(4*3600/stg.dpBinLenS));
        simPop = fillmissing(siCharTbl.(stg.siCharToPlot(kchar) + "SimPop"), 'linear', 1, 'MaxGap', ceil(4*3600/stg.dpBinLenS));
        % Smooth
        ori = filter(1/stg.simMovAveLen*ones(1, stg.simMovAveLen), 1, ori);
        sim = filter(1/stg.simMovAveLen*ones(1, stg.simMovAveLen), 1, sim);
        simPop = filter(1/stg.simMovAveLen*ones(1, stg.simMovAveLen), 1, simPop);
        % Remove NaN
        toKeep = ~isnan(ori);
        ori = ori(toKeep);
        sim = sim(toKeep);
        simPop = simPop(toKeep);
        % Compute the statistics
        
% citatel = sqrt(mean((sim - ori).*(sim - ori)))
% jmenovatel = sqrt(mean(ori.*ori))
        rmsnrmse(kchar) = sqrt(mean((sim - ori).*(sim - ori)))/sqrt(mean(ori.*ori)); % RMS-normalized root mean squared error
% pause
        rmsnrmsePop(kchar) = sqrt(mean((simPop - ori).*(simPop - ori)))/sqrt(mean(ori.*ori));
        [rho(kchar), pval(kchar)] = corr(sim, ori);
        [rhoPop(kchar), pvalPop(kchar)] = corr(simPop, ori);
    end
end

% Plotting functions
% % % % % % % % % function [plotTF, plotName] = createFigInd(multHe)
function h = createFigInd(stg, h, figName, multHe)
    % Create figure for plots of individual subjects
    % If stg.plotXXX related to the calling function is true, creates figure and stores it in h.f structure
    % Based on the number of subjects, it decides on the optimal size of the figure
    % multHe .... multiplier of height, provided by the calling function
    % ret ....... true if stg.plotXXX is false so the calling function will know it should not procede
    % plotName .. char array derived from the name of the calling function
    % % % % % % % % % % global stg
    % % % % % % % % global h
    st = dbstack;
    callerName = st(2).name;
    plotName = [lower(callerName(5)), callerName(6 : end)];

    % Decide on whether we do the calculations and the plot
    if stg.(['plot', upper(plotName(1)), plotName(2 : end)])
        plotTF = true;
        if ~isfield(h.f, plotName)
            if stg.sbNCol == 1
                wiCm = stg.singleColumnWidth;
            elseif stg.sbNCol == 2
                wiCm = 12.75;
            else
                wiCm = 17;
            end
            type = plotName(1 : 2);
            charHe = stg.([type, 'CharHe']);
            if numel(plotName) >= 6 
                if strcmp(plotName(3 : 6), 'Char')
                    numChar = numel(stg.([type, 'CharToPlot']));
                else
                    numChar = 1;
                end
            else
                numChar = 1;
            end
            heCm = min(25, charHe*numChar*multHe*stg.sbNRow);
            positionCm = [20, 1, wiCm, heCm];
            h.f.(plotName) = figure('Units', stg.units, 'Position', positionCm, 'Tag', ['fig', upper(plotName(1)), plotName(2 : end)]);
            h.f.(plotName).Units = 'pixels';
            h.f.(plotName).Name = plotName;
            h.f.(plotName).Color = [1 1 1];
        else
            figure(h.f.(plotName))
        end
    else
        plotTF = false;
    end
end
function [plotTF, plotName] = createFigPos(positionCm, varargin)
    % Create figure, the calling function decides on the position (including size)
    % If stg.plotXXX related to the calling function is true, creates figure and stores it in h.f structure
    % Based on the number of subjects, it decides on the optimal size of the figure
    % positionCm .. position in centimeters
    % ret ......... true if stg.plotXXX is false so the calling function will know it should not procede
    % plotName .... char array derived from the name of the calling function
    % varargin .... optionally specify figure name. Can be used if a figure with the standard name was already created within given caller function and you need another figure.
    global h
    global stg
    if isempty(varargin)
        st = dbstack;
        callerName = st(2).name;
        plotName = [lower(callerName(5)), callerName(6 : end)];
    else
        plotName = varargin{1};
        plotTF = true;
    end
    
    % Decide whether we do the calculations and the plot
    if ~exist('plotTF', 'var')
        plotTF = logical(stg.(['plot', upper(plotName(1)), plotName(2 : end)]));
    end
    if plotTF
        if ~isfield(h.f, plotName)
            % h.f.(plotName) = figure('Position', [300 300 800 400], 'Tag', ['fig', upper(plotName(1)), plotName(2 : end)]);
            h.f.(plotName) = figure('Units', stg.units, 'Position', positionCm, 'Tag', ['fig', upper(plotName(1)), plotName(2 : end)]);
            h.f.(plotName).Units = 'pixels';
            h.f.(plotName).Name = plotName;
            h.f.(plotName).Color = [1 1 1];
        else
            figure(h.f.(plotName))
        end
    end
end
function createAxesAllPopSlopeBox(plotName)
    % Prepares axes for all characteristics for the plots where all subjects are in a single plot
    global stg
    global h
    % Calculate axes positions
    type = plotName(1 : 2);
    numChar = numel(stg.([type, 'CharToPlot']));
    if contains(plotName, 'CharClFitAllPop')
        numChar = numChar + 2*numel(stg.siCharToPlot);
    end
    if contains(plotName, 'CharSzFitAllPop')
        numChar = numChar + numel(stg.siCharToPlot);
    end
    numCol = 2;
    numRow = ceil(numChar/2);
    [spx, spy, spWi, spHe, ~, ~] = getSubplotXYWH(plotName, stg.margGlobSlopeBox, stg.margSlopeBox);
    rowHeShift = spHe/numRow;
    rowHeOff = 0.8; % Centimeters
    rowHe = spHe/numRow - rowHeOff - 0.3;
    colWiShift = spWi/numCol;
    colWiOff1 = 1.2; % Centimeters
    colWi1 = (spWi/numCol - colWiOff1)*0.55;
    colWiOff2 = 0.2; % Centimeters
    colWi2 = (spWi/numCol - colWiOff1)*0.3;
    h.f.(plotName).Units = stg.units;
    for kchar = 1 : numChar
        nCol = rem(kchar-1, numCol) + 1;
        nRow = ceil(kchar/numCol);
        h.a.(plotName)(1, kchar) = axes('Units', stg.units, 'Position', ...
            [spx + (nCol-1)*colWiShift + colWiOff1, spy + (numRow - nRow)*rowHeShift + rowHeOff, colWi1, rowHe], 'NextPlot', 'add');
        h.a.(plotName)(2, kchar) = axes('Units', stg.units, 'Position', ...
            [spx + (nCol-1)*colWiShift + colWiOff1 + colWi1 + colWiOff2, spy + (numRow - nRow)*rowHeShift + rowHeOff, colWi2, rowHe], 'NextPlot', 'add');
    end
end
function createAxesAllPopStack(plotName)
    % Prepares axes for all characteristics for the plots where all subjects are in a single plot
    global stg
    global h

    % Calculate axes positions
    type = plotName(1 : 2);
    numChar = numel(stg.([type, 'CharToPlot']));
    % % % % % % % % margGlob = [2.2 0.5 0 0.5]; % Left, bottom, right, top
    % % % % % % % % marg = [0.7 0.7 0.5 0.5]; % Left, bottom, right, top
    [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlob, stg.marg);
    h.f.(plotName).Units = stg.units;
    for kchar = 1 : numChar
        h.a.(plotName)(kchar) = axes('Units', stg.units, 'Position', ...
            [spx, spy(ceil(1/numc)) - (kchar - numChar)*spHe/numChar, spWi, spHe/numChar], 'NextPlot', 'add');
    end
end
function plotSzCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl)
    global stg
    global h
    
    % Calculate axes positions
    numChar = numel(stg.szCharToPlot);
    % % % % % % % % margGlob = [2.2 0.5 0 0.5]; % Left, bottom, right, top
    % % % % % % % % marg = [0.7 0.7 0.5 0.5]; % Left, bottom, right, top
    [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlob, stg.marg);
    h.f.(plotName).Units = stg.units;
    
    % Plot the data
    for kchar = 1 : numChar
        h.a.(plotName)(ksubj, kchar) = axes('Units', stg.units, 'Position', [...
            spx(mod(ksubj - 1, numc) + 1), spy(ceil(ksubj/numc)) - (kchar - numChar)*spHe/numChar, spWi, spHe/numChar]);
        
        % Signal OK marker
        [x, y] = getSzOkXY(subjInfo, siCharTbl);
        h.p.(plotName)(ksubj, kchar, 1) = plot(x, y, 'Marker', 'none', 'LineWidth', 3, 'Color', stg.subjColor(ksubj, :));
        clear x y
        hold on
        
        % Plot seizures
        % x = (szCharTbl.szOnsN - subjInfo.anStartN); % X data represent days since recording start
        x = (szCharTbl.szOnsN - subjInfo.dob); % X data represent age in days
        % x = datetime(szCharTbl.szOnsN, 'ConvertFrom', 'datenum'); % X data in absolute time
        x = repelem(x, 3);
        y1 = szCharTbl{:, stg.szCharToPlot(kchar)};
        
        y(1 : 3 : 3*size(szCharTbl, 1)) = 0;
        if kchar == 1
            y(2 : 3 : 3*size(szCharTbl, 1)) = 1;
        else
            y(2 : 3 : 3*size(szCharTbl, 1)) = y1;
        end
        y(3 : 3 : 3*size(szCharTbl, 1)) = NaN;
        h.p.(plotName)(ksubj, kchar, 2) = plot(x, y, 'Marker', 'none', 'LineWidth', 0.5, 'Color', stg.subjColor(ksubj, :), 'Tag', 'sz');
        clear x y
    end
end
function plotSzCharDataCirc(plotName, ksubj, szCharTbl)
    global stg
    global h
    % Calculate axes positions
    numChar = numel(stg.szCharToPlot);
    % % % % % % % % margGlob = [0 0 0 0];
    % % % % % % % % marg = [0.1 0.5 0.1 0.5];
    [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlobCi, stg.margCi);
    h.f.(plotName).Units = stg.units;
    
    % Plot the data
    numCol = ceil((numChar+1)^0.8);
    numRow = ceil((numChar+1)/numCol);
    axWiHe = 0.8*min(spWi/numCol, spHe/numRow);
    for kchar = 1 : numChar+1
        h.a.(plotName)(ksubj, kchar) = axes('Units', stg.units, 'Position', [...
            spx(mod(ksubj - 1, numc) + 1) + mod(kchar-1, numCol)*spWi/numCol + max(spWi/numCol - axWiHe, 0)/2,...
            spy(ceil(ksubj/numc)) + (numRow-ceil(kchar/numCol))*spHe/numRow + 0.0*spHe/numRow,...
            axWiHe, ...
            axWiHe],...
            'NextPlot', 'add', 'Visible', 'off', 'XLimMode', 'manual', 'YLimMode', 'manual');
        hax = h.a.(plotName)(ksubj, kchar);
        if kchar ~= numChar+1
            % Plot seizures
            p = (szCharTbl.szOnsN - floor(szCharTbl.szOnsN))*2*pi; % Circadian phase
            p = repelem(p, 3)';
            y1 = szCharTbl{:, stg.szCharToPlot(kchar)};
            if kchar == 1
                r(1 : 3 : 3*size(szCharTbl, 1)) = 0;
                r(2 : 3 : 3*size(szCharTbl, 1)) = 1;
            else
                r(1 : 3 : 3*size(szCharTbl, 1)) = 0;
                % r(2 : 3 : 3*size(szCharTbl, 1)) = (y1-min(y1))/range(y1); % Normalize by range, this effectively moves the origin of the circle
                r(2 : 3 : 3*size(szCharTbl, 1)) = y1/max(y1); % Normalize by max, origin corresponds to 0
            end
            r(3 : 3 : 3*size(szCharTbl, 1)) = NaN;
            [x, y] = pol2cart(-p - pi/2, r);
            h.p.(plotName)(ksubj, kchar, 1) = plot(x, y, 'Marker', 'none', 'LineWidth', 0.5, 'Color', stg.subjColor(ksubj, :), 'Tag', 'sz');
        end
        % Circle
        circleP = (0 : 1 : 360)/360*2*pi;
        circleR = 1;
        % if ~all(isnan(y1))
        %     circleR = ceilToEven1sig(max(r));
        % else
        %     circleR = 1;
        % end
        [x, y] = pol2cart(-circleP - pi/2, circleR);
        h.p.(plotName)(ksubj, kchar, 2) = plot(x, y, 'k');
        hax.XLim = max([x, y])*[-1 1];
        hax.YLim = hax.XLim;
        % Night shading
        circleP = (-90 : 1 : 90)/360*2*pi;
        % colo = [linspace(1, 0, 90), 1, linspace(0, 1, 90)]'*[1 1 1];
        colo = [linspace(1, 0.7, 90), 1, linspace(0.7, 1, 90)]'*[1 1 1];
        colo = permute(colo, [1 3 2]);
        [x, y] = pol2cart(-circleP - pi/2, circleR);
        z = -10*ones(size(x));
        h.pa.(plotName)(ksubj, kchar, 1) = patch(x, y, z, colo, 'EdgeColor', 'none');
        if kchar == numChar+1 % Explanatory circle
            text(hax, 0, -0.65, 'midnight', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', stg.statFontSize)
            % text(hax, -0.65, 0, 'morning', 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', stg.statFontSize)
            text(hax, 0, 0.65, 'noon', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', stg.statFontSize)
            % text(hax, 0.65, 0, 'afternoon', 'Rotation', 270, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', stg.statFontSize)
            % Left arrow
            p = (15 : 35)/100;
            p = p*2*pi;
            r = 0.7*ones(size(p));
            [x, y] = pol2cart(-p - pi/2, r);
            plot(hax, x, y, 'k', 'LineWidth', 1)
            p1 = p(end) - 3/4*pi;
            r1 = 0.2;
            [x1, y1] = pol2cart(-p1 - pi/2, r1);
            x1 = x(end) + [0, x1];
            y1 = y(end) + [0, y1];
            plot(hax, x1, y1, 'k', 'LineWidth', 1)
            p2 = p(end) - 1/3*pi;
            r2 = 0.2;
            [x2, y2] = pol2cart(-p2 - pi/2, r2);
            x2 = x(end) + [0, x2];
            y2 = y(end) + [0, y2];
            plot(hax, x2, y2, 'k', 'LineWidth', 1)
            % Right arrow
            p = (65 : 85)/100;
            p = p*2*pi;
            r = 0.7*ones(size(p));
            [x, y] = pol2cart(-p - pi/2, r);
            plot(hax, x, y, 'k', 'LineWidth', 1)
            p1 = p(end) - 3/4*pi;
            r1 = 0.2;
            [x1, y1] = pol2cart(-p1 - pi/2, r1);
            x1 = x(end) + [0, x1];
            y1 = y(end) + [0, y1];
            plot(hax, x1, y1, 'k', 'LineWidth', 1)
            p2 = p(end) - 1/3*pi;
            r2 = 0.2;
            [x2, y2] = pol2cart(-p2 - pi/2, r2);
            x2 = x(end) + [0, x2];
            y2 = y(end) + [0, y2];
            plot(hax, x2, y2, 'k', 'LineWidth', 1)
        end
    end
end
function plotSaCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl)
    global stg
    global h
    % Calculate axes positions
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    % % % % % % % % margGlob = [2.2 0.5 0 0.5]; % Left, bottom, right, top
    % % % % % % % % marg = [0.7 0.7 0.5 0.5]; % Left, bottom, right, top
    [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlob, stg.marg); % Provides the position of the space for one subject (not one char of the subject)
    h.f.(plotName).Units = stg.units;
    % Plot the data
    for kchar = 1 : numChar
        h.a.(plotName)(ksubj, kchar) = axes('Units', stg.units, 'Position', [...
            spx(mod(ksubj - 1, numc) + 1), spy(ceil(ksubj/numc)) - (kchar - numChar)*spHe/numChar, spWi, spHe/numChar]);

        % Signal OK marker
        [x, y] = getSzOkXY(subjInfo, siCharTbl);
        h.p.(plotName)(ksubj, kchar, 1) = plot(x, y, 'Marker', 'none', 'LineWidth', 3, 'Color', stg.subjColor(ksubj, :));
        clear x y
        hold on
        
        % Plot seizures
        % x = (szCharTbl.szOnsN - subjInfo.anStartN); % X data represent days since recording start
        x = (szCharTbl.szOnsN - subjInfo.dob); % X data represent age in days
        % x = datetime(szCharTbl.szOnsN, 'ConvertFrom', 'datenum'); % X data in absolute time
        x = repelem(x, 3);
        if kchar <= numSzChar
            y1 = szCharTbl{:, stg.szCharToPlot(kchar)};
        else
            y1 = szCharTbl{:, stg.szCharToPlot(1)};
        end
        y(1 : 3 : 3*size(szCharTbl, 1)) = 0;
        if kchar == 1 || kchar > numSzChar
            y(2 : 3 : 3*size(szCharTbl, 1)) = 1;
        else
            y(2 : 3 : 3*size(szCharTbl, 1)) = y1;
        end
        y(3 : 3 : 3*size(szCharTbl, 1)) = NaN;
        h.p.(plotName)(ksubj, kchar, 2) = plot(x, y, 'Marker', 'none', 'LineWidth', 0.5, 'Color', stg.subjColor(ksubj, :), 'Tag', 'sz');
        clear x y
        
        if kchar > numSzChar
            % Plot signal characteristic
            y = siCharTbl{:, stg.siCharToPlot(kchar - numSzChar)};
            x = siCharTbl.tax - subjInfo.dob;
            h.p.(plotName)(ksubj, kchar, 3) = plot(x, y, 'Marker', 'none', 'LineWidth', 0.5, 'Color', 'k');
            yl = getYLimYTick(y, stg.siCharYLim(kchar - numSzChar, :));
    
            % Format the plot
            h.p.(plotName)(ksubj, kchar, 1).YData(~isnan(h.p.(plotName)(ksubj, kchar, 1).YData)) = yl(1) + 1e-12;
            h.p.(plotName)(ksubj, kchar, 2).YData(1 : 3 : end-2) = yl(1) + 1e-12; % Probably for the vertical lines to not go over the axes. Try removing it.
            h.p.(plotName)(ksubj, kchar, 2).YData(2 : 3 : end-1) = yl(2) - 1e-12;
        end
    end
end
function plotSaCharDataCWT(plotName, ksubj, subjInfo, szCharTbl, siCharTbl) %#ok<INUSD>
    global stg
    global h
    % Calculate axes positions
    numSzChar = 1; % We will only plot seizure rate, not other characteristics
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar + numSiChar; % Seizure rate + number of signal chars + number of coherences between signal chars and seizures
    [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlob, stg.marg); % Provides the position of the space for one subject (not one char of the subject)
    h.f.(plotName).Units = stg.units;
    tax = (siCharTbl.tax - subjInfo.dob); % Time axis
    voicesPerOctave = 12;
    minFreq = 1/16;
    maxFreq = 2;
    % Plot the data
    for kchar = 1 : numChar
        h.a.(plotName)(ksubj, kchar) = axes('Units', stg.units, 'Position', [...
            spx(mod(ksubj - 1, numc) + 1), spy(ceil(ksubj/numc)) - (kchar - numChar)*spHe/numChar, spWi, spHe/numChar]);
        hax = h.a.(plotName)(ksubj, kchar);
        fs = 1/stg.dpBinLenS*3600*24; % Frequency is in cycles per day
        if kchar <= numSzChar
            % Plot seizures
            y = siCharTbl.sz;
            % y = fillmissing(y, 'linear', 'MaxGap', 6);
            y = fillmissing(y, 'linear');
            [wt, fax, coi] = cwtPiecewiseGemini(tax, y, fs); % Frequency is in cycles per day
            % [cfs, fax] = cwt(y, 1/stg.dpBinLenS*3600*24); % Frequency is in cycles per day
            % % % h.p.(plotName)(ksubj, kchar, 1) = imagesc("XData",tax,"YData",fax,"CData",abs(wt),"CDataMapping","scaled");
            % h.p.(plotName)(ksubj, kchar, 1) = imagesc(tax, fax, abs(wt), 'Parent', hax);
            h.p.(plotName)(ksubj, kchar, 1) = pcolor(tax, fax, abs(wt));
            shading flat;
            % ylabel(['Sz rate', 10, 'Cyc/day'], 'FontSize', stg.axFontSize)
            ylabel('Sz rate', 'FontSize', stg.axFontSize)
            hold on;
            plot(tax, coi, 'w--', 'LineWidth', 1.5);
        end
        if kchar > numSzChar && kchar <= numSzChar + numSiChar
            % Plot signal characteristic
            y = siCharTbl{:, stg.siCharToPlot(kchar - numSzChar)};
            % y = fillmissing(y, 'linear', 'MaxGap', 6);
            y = fillmissing(y, 'linear');
            % % % y = 0.5*sin(2*pi*0.25*tax) + 0.5*sin(2*pi*1*tax) + 0.5*sin(2*pi*2*tax);
            [wt, fax] = cwtPiecewiseGemini(tax, y, fs); % Frequency is in cycles per day
            % [cfs, fax] = cwt(y, 1/stg.dpBinLenS*3600*24); % Frequency is in cycles per day
            % % % h.p.(plotName)(ksubj, kchar, 1) = image("XData",tax,"YData",fax,"CData",abs(wt),"CDataMapping","scaled");
            % h.p.(plotName)(ksubj, kchar, 1) = imagesc(tax, fax, abs(wt), 'Parent', hax);
            h.p.(plotName)(ksubj, kchar, 1) = pcolor(tax, fax, abs(wt));
            shading flat;
            % ylabel(['IED rate', 10, 'Cyc/day'], 'FontSize', stg.axFontSize);
            ylabel('IED rate', 'FontSize', stg.axFontSize);
            hold on;
            plot(tax, coi, 'w--', 'LineWidth', 1.5);
        end
        if kchar > numSzChar + numSiChar
            y1 = siCharTbl.sz;
            y1 = fillmissing(y1, 'linear') + 0.1*randn(size(y1));
            y2 = siCharTbl{:, stg.siCharToPlot(kchar - numSzChar - numSiChar)};
            y2 = fillmissing(y2, 'linear') + 0.1*randn(size(y2));
            [wcoh, wcs, freqs, coi] = wcoherence(y1, y2, fs, VoicesPerOctave=voicesPerOctave, FrequencyLimits=[minFreq maxFreq]); %#ok<ASGLU>
            fb = cwtfilterbank('SignalLength', 1e6, ...
                   'SamplingFrequency', fs, ...
                   'VoicesPerOctave', voicesPerOctave, ...
                   'FrequencyLimits', [minFreq maxFreq]);
            freqs = centerFrequencies(fb);
            num_scales = numel(freqs);
            if size(wcoh, 1) > num_scales
                pause
            end
            if size(wcoh, 1) < num_scales % If fewer rows than expected
                temp_wcoh = NaN(num_scales, size(wcoh, 2)); % Create NaN matrix of desired size
                temp_wcoh(1:size(wcoh,1), :) = wcoh; % Copy the available coefficients
                wcoh = temp_wcoh; % Use the padded version
                temp_wcs = NaN(num_scales, size(wcs, 2)); % Create NaN matrix of desired size
                temp_wcs(1:size(wcs,1), :) = wcs;
                wcs = temp_wcs;
                temp_freqs = NaN(num_scales, 1); % Create NaN matrix of desired size
                temp_freqs(1:size(freqs,1), 1) = freqs;
                freqs = temp_freqs;
            end
            wt = wcoh;
            fax = freqs;
            h.p.(plotName)(ksubj, kchar, 1) = pcolor(tax, fax, abs(wt));
            shading flat;
            clim([0 1])
            hold on;
            plotWaveletCoherenceArrows(wcoh, wcs, freqs, coi, hax)
            % ylabel(['Sz vs. IED', 10, 'Cyc/day'], 'FontSize', stg.axFontSize);
            ylabel('Sz vs. IED', 'FontSize', stg.axFontSize);
        end
        hold on
        yline(1, ':w')
        yline(0.25, ':w')
        hax.Layer = 'top';
        hax.XLim = [tax(1), tax(end)];
        hax.YLim = [minFreq maxFreq];
        hax.YTick = [0 0.25 0.5 1 2];
        hax.YDir = 'normal';
        hax.YScale = 'log';
        hax.FontSize = stg.axFontSize;
        hax.Box = stg.box;
        hax.Units = 'normalized';
        plot(tax, coi, 'w--', 'LineWidth', 1.5);
        patch('XData', [tax; flipud(tax)], 'YData', [0.01*ones(size(tax)); coi'], 'ZData', 1000*ones(size([tax; flipud(tax)])), 'FaceColor', 0.0*[1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3)
        if kchar == 1
            hax.Title.String = subjInfo.subjNm;
            hax.Title.Interpreter = 'none';
            % hax.Title.Color = stg.subjColor(ksubj, :);
            hax.Title.Color = 'k';
            if ~subjInfo.sex
                hax.Title.FontAngle = 'italic';
            end
        end
        if kchar < numChar
            hax.XTick = [];
        end
        if kchar == numChar
            % Signal OK marker
            [x, y] = getSzOkXY(subjInfo, siCharTbl);
            h.p.(plotName)(ksubj, kchar, 4) = plot(x, y+1/16, 'Marker', 'none', 'LineWidth', 3, 'Color', stg.subjColor(ksubj, :));
            clear x y
            hold on
            if ksubj > stg.numSubj - 3
                xlabel('Age (days)');
            end
        end
    end
end
function plotWaveletCoherenceArrows(wcoh, wcs, freqs, coi, hax)
% plotWaveletCoherenceArrows: Plots phase arrows on a wavelet coherence scalogram
%   using manually drawn lines and controlling for logarithmic Y-axis distortion.
%
%   This function overlays phase arrows on an existing time-frequency plot
%   (e.g., created with pcolor or imagesc) to indicate the phase relationship
%   between two signals. Arrows are only plotted in regions of high coherence
%   and outside the Cone of Influence (COI). It manually draws arrow shafts
%   and heads to ensure a consistent visual length and shape across a
%   logarithmic Y-axis.
%
% Inputs:
%   wcoh  - Magnitude-squared wavelet coherence matrix (real, 0 to 1).
%           (Output from wcoherence, first argument)
%   wcs   - Smoothed wavelet cross-spectrum matrix (complex).
%           (Output from wcoherence, second argument)
%   freqs - Vector of frequencies (Hz) corresponding to rows of wcoh/wcs.
%           (Output from wcoherence, third argument 'f')
%   coi   - Vector of Cone of Influence periods, same length as time dimension.
%           (Output from wcoherence, fourth argument 'coi')
%   hax   - Axes object handle where the scalogram is already plotted.
%           Crucially, the Y-axis of these axes should be set to 'log' scale.

    % --- Input Validation ---
    if ~ishandle(hax) || ~strcmp(get(hax, 'Type'), 'axes')
        error('plotWaveletCoherenceArrows:InvalidAxesHandle', ...
              'Invalid axes handle provided. Please ensure hax is a valid axes object.');
    end

    % --- Parameters for arrow plotting (adjust these for desired visual density and length) ---
    % Default vertical density (number of rows to skip for arrows)
    freq_subsample_factor       = 10;    
    coherence_threshold         = 0.7;   % Min wcoh value to display an arrow (0 to 1)
    
    % These factors control the *visual length* of the arrow in X and Y (log) directions.
    % You'll likely need to fine-tune 'arrow_visual_length_factor' to get the overall size
    % and 'arrowhead_angle_degrees' and 'arrowhead_length_ratio' for the head's appearance.
    arrow_visual_length_factor  = 0.15;  % Overall length relative to a typical plot dimension.
                                         % Larger value means longer arrows.
    arrowhead_angle_degrees     = 25;    % Angle of arrowhead barbs relative to main shaft (in degrees)
    arrowhead_length_ratio      = 0.6;   % Length of arrowhead barbs relative to total arrow length (0 to 1)
    
    arrow_color                 = 'k';   % Color of the arrows (e.g., 'k' for black, 'w' for white)
    arrow_line_width            = 0.3;     % Line width of the arrow shafts and heads

    % --- 1. Get the time axis from the existing plot ---
    plot_obj = findobj(hax, 'Type', 'Surface', '-or', 'Type', 'Image');
    if isempty(plot_obj)
        warning('plotWaveletCoherenceArrows:NoPlotFound', ...
                'No pcolor/imagesc plot found in the provided axes. Using axes X-limits as fallback for time.');
        times = linspace(hax.XLim(1), hax.XLim(2), size(wcoh, 2));
    else
        times = get(plot_obj(1), 'XData'); 
    end

    % --- Dynamic adjustment of time_subsample_factor ---
    % Ensure a maximum of 20 arrows in a row for long datasets
    max_arrows_in_row = 20;
    num_total_time_points = length(times); % Total number of time points in the data

    % Calculate the required subsample factor to meet the max_arrows_in_row constraint
    calculated_time_subsample_factor = ceil(num_total_time_points / max_arrows_in_row);

    % Use a default factor for shorter data, or the calculated factor for longer data
    default_min_time_subsample_factor = 25; % Ensures at least some reasonable sparsity for short data
    time_subsample_factor = max(default_min_time_subsample_factor, calculated_time_subsample_factor);
    
    % Ensure time_subsample_factor is at least 1 (to avoid errors with extremely short data)
    time_subsample_factor = max(1, time_subsample_factor);

    % --- 2. Calculate phase from the cross-spectrum ---
    phase = angle(wcs); % Phase in radians (-pi to pi)

    % --- 3. Create subsampled grids for arrow positions ---
    X_arrows = times(1:time_subsample_factor:end);
    Y_arrows = freqs(1:freq_subsample_factor:end);
    [X_grid, Y_grid] = meshgrid(X_arrows, Y_arrows);

    % --- 4. Subsample coherence and phase data for the arrow grid ---
    phase_subsampled = phase(1:freq_subsample_factor:end, 1:time_subsample_factor:end);
    coherence_subsampled = wcoh(1:freq_subsample_factor:end, 1:time_subsample_factor:end);

    % --- 5. Determine which arrows to plot (coherence threshold & COI) ---
    plot_mask = (coherence_subsampled > coherence_threshold); % Coherence thresholding

    coi_periods_interp = interp1(times, coi, X_arrows, 'linear', 'extrap'); % Interpolate COI
    coi_periods_interp = coi_periods_interp(:)'; % Ensure it's a row vector for repmat
    
    periods_Y_grid = 1 ./ Y_grid; % Convert Y_grid frequencies to periods for COI comparison
    coi_mask_matrix = periods_Y_grid < repmat(coi_periods_interp, size(periods_Y_grid, 1), 1);
    
    plot_mask = plot_mask & ~coi_mask_matrix; % Combine masks: high coherence AND outside COI

    % --- 6. Calculate arrow components in (X_data, log_Y_data) space ---
    % First, define the 'unit' length for arrows in this transformed space.
    % We'll use a fraction of the average spacing in X and the log-Y range.
    x_range = max(hax.XLim) - min(hax.XLim);
    log_y_range = max(log(hax.YLim)) - min(log(hax.YLim)); % Calculate range in log Y space

    % This 'reference_length' determines the overall visual size of the arrows.
    % It's a balance between X and logY axis dimensions.
    reference_length = min(x_range, log_y_range) * arrow_visual_length_factor; 

    % U component (horizontal): In linear X-data units
    U_comp = reference_length * cos(phase_subsampled) .* plot_mask;
    
    % V component (vertical): In log(Y_data) units
    V_log_comp = reference_length * sin(phase_subsampled) .* plot_mask;

    % --- 7. Calculate actual start and end points for drawing ---
    X_start = X_grid;
    log_Y_start = log(Y_grid); % Convert Y_grid to its log value for calculation

    X_end = X_start + U_comp;
    log_Y_end = log_Y_start + V_log_comp;
    Y_end = exp(log_Y_end); % Convert log_Y_end back to linear Y-data units for plotting

    % --- 8. Loop through and draw each arrow manually ---
    axes(hax); % Make sure the target axes is current
    hold on;   % Add to the existing plot

    arrowhead_angle_rad = deg2rad(arrowhead_angle_degrees);

    for i = 1:numel(X_grid)
        if plot_mask(i) % Only draw if the mask is true for this arrow
            x1 = X_start(i);
            y1 = Y_grid(i); % CORRECTED: Use Y_grid for the starting Y-coordinate
            x2 = X_end(i);
            y2 = Y_end(i);

            % Draw the arrow shaft
            line([x1, x2], [y1, y2], 'Color', arrow_color, 'LineWidth', arrow_line_width);

            % --- Draw the arrowhead ---
            % Angle of the main arrow shaft in the transformed (X, logY) space
            current_arrow_angle = atan2(log_Y_end(i) - log_Y_start(i), X_end(i) - X_start(i));
            
            % Head length in the transformed (X, logY) space
            head_len = reference_length * arrowhead_length_ratio;

            % Calculate arrowhead barb points in (X_data, log_Y_data) space
            % Barb 1
            x_barb1_transformed = X_end(i) - head_len * cos(current_arrow_angle - arrowhead_angle_rad);
            log_y_barb1_transformed = log_Y_end(i) - head_len * sin(current_arrow_angle - arrowhead_angle_rad);
            y_barb1 = exp(log_y_barb1_transformed); % Convert back to linear Y data

            % Barb 2
            x_barb2_transformed = X_end(i) - head_len * cos(current_arrow_angle + arrowhead_angle_rad);
            log_y_barb2_transformed = log_Y_end(i) - head_len * sin(current_arrow_angle + arrowhead_angle_rad);
            y_barb2 = exp(log_y_barb2_transformed); % Convert back to linear Y data
            
            % Draw the arrowhead barbs
            line([x2, x_barb1_transformed], [y2, y_barb1], 'Color', arrow_color, 'LineWidth', arrow_line_width);
            line([x2, x_barb2_transformed], [y2, y_barb2], 'Color', arrow_color, 'LineWidth', arrow_line_width);
        end
    end
    
    hold off; % Release the hold on the axes
end
function plotSsCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl)
    global stg
    global h
    % Calculate axes positions
    numSiChar = numel(stg.siCharToPlot);
    numChar = 2*numSiChar;
    [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlob, stg.marg); % Provides the position of the space for one subject (not one char of the subject)
    h.f.(plotName).Units = stg.units;
    % Plot the data
    for kchar = 1 : numChar
        h.a.(plotName)(ksubj, kchar) = axes('Units', stg.units, 'Position', [...
            spx(mod(ksubj - 1, numc) + 1), spy(ceil(ksubj/numc)) - (kchar - numChar)*spHe/numChar, spWi, spHe/numChar]);

        % Signal OK marker
        [x, y] = getSzOkXY(subjInfo, siCharTbl);
        h.p.(plotName)(ksubj, kchar, 1) = plot(x, y, 'Marker', 'none', 'LineWidth', 3, 'Color', stg.subjColor(ksubj, :));
        clear x y
        hold on

        % Plot seizures
        x = (szCharTbl.szOnsN - subjInfo.dob); % X data represent age in days
        x = repelem(x, 3);
        y = NaN(1, 3*size(szCharTbl, 1));
        y(1 : 3 : 3*size(szCharTbl, 1)) = 0;
        y(2 : 3 : 3*size(szCharTbl, 1)) = 1;
        h.p.(plotName)(ksubj, kchar, 2) = plot(x, y, 'Marker', 'none', 'LineWidth', 0.5, 'Color', stg.subjColor(ksubj, :), 'Tag', 'sz');
        clear x y
        if kchar <= numSiChar
            % Plot simulated signal characteristic
            y = siCharTbl{:, stg.siCharToPlot(kchar) + "Sim"};
            y = filter(1/stg.simMovAveLen*ones(1, stg.simMovAveLen), 1, y);
            x = siCharTbl.tax - subjInfo.dob;
            h.p.(plotName)(kchar, ksubj, 3) = plot(x, y, 'Marker', 'none', 'LineWidth', 2, 'Color', stg.simColor);
            % Plot signal characteristic
            y = siCharTbl{:, stg.siCharToPlot(kchar)};
            y = fillmissing(y, 'linear', 1, 'MaxGap', ceil(4*3600/stg.dpBinLenS));
            y = filter(1/stg.simMovAveLen*ones(1, stg.simMovAveLen), 1, y);
            x = siCharTbl.tax - subjInfo.dob;
            h.p.(plotName)(kchar, ksubj, 3) = plot(x, y, 'Marker', 'none', 'LineWidth', 0.5, 'Color', 'k');
            % Format the plot
            yl = getYLimYTick(y, stg.siCharYLim(kchar, :));
            h.p.(plotName)(ksubj, kchar, 1).YData(~isnan(h.p.(plotName)(ksubj, kchar, 1).YData)) = yl(1) + 1e-12;
            h.p.(plotName)(ksubj, kchar, 2).YData(1 : 3 : end-2) = yl(1) + 1e-12; % Probably for the vertical lines to not go over the axes. Try removing it.
            h.p.(plotName)(ksubj, kchar, 2).YData(2 : 3 : end-1) = yl(2) - 1e-12;
        end

        if kchar > numSiChar
            % Plot simulated signal characteristic
            y = siCharTbl{:, stg.siCharToPlot(kchar - numSiChar) + "SimPop"};
            y = filter(1/stg.simMovAveLen*ones(1, stg.simMovAveLen), 1, y);
            x = siCharTbl.tax - subjInfo.dob;
            h.p.(plotName)(kchar, ksubj, 3) = plot(x, y, 'Marker', 'none', 'LineWidth', 2, 'Color', stg.simPopColor);
            % Plot artificial signal characteristic
            y = siCharTbl{:, stg.siCharToPlot(kchar  - numSiChar)};
            y = fillmissing(y, 'linear', 1, 'MaxGap', ceil(4*3600/stg.dpBinLenS));
            y = filter(1/stg.simMovAveLen*ones(1, stg.simMovAveLen), 1, y);
            x = siCharTbl.tax - subjInfo.dob;
            h.p.(plotName)(kchar, ksubj, 3) = plot(x, y, 'Marker', 'none', 'LineWidth', 0.5, 'Color', 'k');
            % Format the plot
            yl = getYLimYTick(y, stg.siCharYLim(kchar - numSiChar, :));
            h.p.(plotName)(ksubj, kchar, 1).YData(~isnan(h.p.(plotName)(ksubj, kchar, 1).YData)) = yl(1) + 1e-12;
            h.p.(plotName)(ksubj, kchar, 2).YData(1 : 3 : end-2) = yl(1) + 1e-12; % Probably for the vertical lines to not go over the axes. Try removing it.
            h.p.(plotName)(ksubj, kchar, 2).YData(2 : 3 : end-1) = yl(2) - 1e-12;
        end
    end
end
function numCol = plotSaCharDataCirc(plotName, ksubj, szCharTbl, ppSiTbl, rrSiTbl)
    global stg
    global h
    % Calculate axes positions
    numSzChar = numel(stg.szCharToPlot);
    numSiChar = numel(stg.siCharToPlot);
    numChar = numSzChar + numSiChar;
    % % % % % % % % margGlob = [0 0 0 0];
    % % % % % % % % marg = [0.1 0.5 0.1 0.5];
    [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlobCi, stg.margCi);
    h.f.(plotName).Units = stg.units;
    
    % Plot the data
    axWiHeSizeCoef = 0.7;
    numCol = ceil((numChar+1)^0.6); % Within subject's space
    numRow = ceil((numChar+1)/numCol);
    axWiHe = axWiHeSizeCoef*min(spWi/numCol, spHe/numRow);
    for kchar = 1 : numChar+1
        h.a.(plotName)(ksubj, kchar) = axes('Units', stg.units, 'Position', [...
            spx(mod(ksubj - 1, numc) + 1) + mod(kchar-1, numCol)*spWi/numCol + max(spWi/numCol - axWiHe, 0)/2,...
            spy(ceil(ksubj/numc)) + (numRow-ceil(kchar/numCol))*spHe/numRow + 0.0*spHe/numRow,...
            axWiHe, ...
            axWiHe],...
            'NextPlot', 'add', 'Visible', 'off', 'XLimMode', 'manual', 'YLimMode', 'manual');
        hax = h.a.(plotName)(ksubj, kchar);

        % Plot seizure characteristics
        if kchar ~= numChar+1
            % Plot seizures
            p = (szCharTbl.szOnsN - floor(szCharTbl.szOnsN))*2*pi; % Circadian phase
            p = repelem(p, 3)';
            if kchar == 1 || kchar > numSzChar
                r(1 : 3 : 3*size(szCharTbl, 1)) = 0;
                r(2 : 3 : 3*size(szCharTbl, 1)) = 1;
            else
                y1 = szCharTbl{:, stg.szCharToPlot(kchar)};
                r(1 : 3 : 3*size(szCharTbl, 1)) = 0;
                % r(2 : 3 : 3*size(szCharTbl, 1)) = (y1-min(y1))/range(y1); % Normalize by range, this effectively moves the origin of the circle
                r(2 : 3 : 3*size(szCharTbl, 1)) = y1/max(y1); % Normalize by max, origin corresponds to 0
            end
            r(3 : 3 : 3*size(szCharTbl, 1)) = NaN;
            [x, y] = pol2cart(-p - pi/2, r);
            h.p.(plotName)(ksubj, kchar, 1) = plot(x, y, 'Marker', 'none', 'LineWidth', 0.5, 'Color', stg.subjColor(ksubj, :), 'Tag', 'sz');
        end

        % Plot signal characteristics
        if kchar > numSzChar && kchar < numChar + 1
            p = [ppSiTbl{1, kchar - numSzChar}, ppSiTbl{1, kchar - numSzChar}(1)];
            r = [rrSiTbl{1, kchar - numSzChar}, rrSiTbl{1, kchar - numSzChar}(1)]/max(rrSiTbl{1, kchar - numSzChar});
            [x, y] = pol2cart(-p - pi/2, r);
            h.p.(plotName)(ksubj, kchar, 2) = plot(x, y, 'k', 'LineWidth', 1);
        end

        % Plot circle
        circleP = (0 : 1 : 360)/360*2*pi;
        circleR = 1;
        % if ~all(isnan(y1))
        %     circleR = ceilToEven1sig(max(r));
        % else
        %     circleR = 1;
        % end
        [x, y] = pol2cart(-circleP - pi/2, circleR);
        h.p.(plotName)(ksubj, kchar, 2) = plot(x, y, 'k');
        hax.XLim = max([x, y])*[-1 1];
        hax.YLim = hax.XLim;
        % Night shading
        circleP = (-90 : 1 : 90)/360*2*pi;
        % colo = [linspace(1, 0, 90), 1, linspace(0, 1, 90)]'*[1 1 1];
        colo = [linspace(1, 0.7, 90), 1, linspace(0.7, 1, 90)]'*[1 1 1];
        colo = permute(colo, [1 3 2]);
        [x, y] = pol2cart(-circleP - pi/2, circleR);
        z = -10*ones(size(x));
        h.pa.(plotName)(ksubj, kchar, 1) = patch(x, y, z, colo, 'EdgeColor', 'none');

        % Plot explanatory circle
        if kchar == numChar+1
            plotExplanatoryCircle(hax)
        end
    end
end
function plotSiCharData(plotName, ksubj, subjInfo, szCharTbl, siCharTbl)
    % Plots seizures and signal characteristics on top of them
    global stg
    global h
    % Calculate axes positions
    numChar = numel(stg.siCharToPlot);
    % % % % % % % % margGlob = [1.8 0.5 0 0.5]; % Left, bottom, right, top
    % % % % % % % % marg = [0.7 0.7 0.5 0.5]; % Left, bottom, right, top
    [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlob, stg.marg);
    h.f.(plotName).Units = stg.units;
    % Plot the data
    for kchar = 1 : numChar
        h.a.(plotName)(ksubj, kchar) = axes('Units', stg.units, 'Position', ...
            [spx(mod(ksubj - 1, numc) + 1), spy(ceil(ksubj/numc)) - (kchar - numChar)*spHe/numChar, spWi, spHe/numChar]);
        h.a.(plotName)(ksubj, kchar).Units = 'normalized';

        % Signal OK marker
        [x, y] = getSzOkXY(subjInfo, siCharTbl);
        h.p.(plotName)(ksubj, kchar, 1) = plot(x, y, 'Marker', 'none', 'LineWidth', 3, 'Color', stg.subjColor(ksubj, :));
        clear x y
        hold on
        
        % Plot seizures
        % x = (szCharTbl.szOnsN - subjInfo.anStartN); % X data common for polynomial fitting and plotting
        x = (szCharTbl.szOnsN - subjInfo.dob); % X data common for polynomial fitting and plotting
        % x = datetime(szCharTbl.szOnsN, 'ConvertFrom', 'datenum');
        x = repelem(x, 3);
        y(1 : 3 : 3*size(szCharTbl, 1)) = 0;
        y(2 : 3 : 3*size(szCharTbl, 1)) = 1;
        y(3 : 3 : 3*size(szCharTbl, 1)) = NaN;
        h.p.(plotName)(ksubj, kchar, 2) = plot(x, y, 'Marker', 'none', 'LineWidth', 0.5, 'Color', stg.subjColor(ksubj, :), 'Tag', 'sz');
        hold on
        
        % Plot signal characteristic
        y = siCharTbl{:, stg.siCharToPlot(kchar)};
        x = siCharTbl.tax - subjInfo.dob;
        h.p.(plotName)(kchar, ksubj, 3) = plot(x, y, 'Marker', 'none', 'LineWidth', 0.5, 'Color', 'k');
        yl = getYLimYTick(y, stg.siCharYLim(kchar, :));

        % Format the plot
        h.p.(plotName)(ksubj, kchar, 1).YData(~isnan(h.p.(plotName)(ksubj, kchar, 1).YData)) = yl(1) + 1e-12;
        h.p.(plotName)(ksubj, kchar, 2).YData(1 : 3 : end-2) = yl(1) + 1e-12;
        h.p.(plotName)(ksubj, kchar, 2).YData(2 : 3 : end-1) = yl(2) - 1e-12;
        
        clear x y
    end
end
function plotSiCharDataCirc(plotName, ksubj, szCharTbl, ppTbl, rrTbl)
    global stg
    global h
    % Calculate axes positions
    numChar = numel(stg.siCharToPlot);
    [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlobCi, stg.margCi);
    h.f.(plotName).Units = stg.units;
    
    % Plot the data
    numCol = ceil(numChar^0.7);
    numRow = ceil(numChar/numCol);
    axWiHe = 0.8*min(spWi/numCol, spHe/numRow);
    for kchar = 1 : numChar
        h.a.(plotName)(ksubj, kchar) = axes('Units', stg.units, 'Position', [...
            spx(mod(ksubj - 1, numc) + 1) + mod(kchar-1, numCol)*spWi/numCol + max(spWi/numCol - axWiHe, 0)/2,...
            spy(ceil(ksubj/numc)) + (numRow-ceil(kchar/numCol))*spHe/numRow + 0.0*spHe/numRow,...
            axWiHe, ...
            axWiHe],...
            'NextPlot', 'add', 'Visible', 'off', 'XLimMode', 'manual', 'YLimMode', 'manual');
        hax = h.a.(plotName)(ksubj, kchar);
        % circleR = ceilToEven1sig(max(rrTbl{1, kchar}));
        circleR = 1;
        % Plot seizures
        p = (szCharTbl.szOnsN - floor(szCharTbl.szOnsN))*2*pi; % Circadian phase
        p = repelem(p, 3)';
        r(1 : 3 : 3*size(szCharTbl, 1)) = 0; % Resultant length (aka modulus)
        r(2 : 3 : 3*size(szCharTbl, 1)) = circleR;
        r(3 : 3 : 3*size(szCharTbl, 1)) = NaN;
        [x, y] = pol2cart(-p - pi/2, r);
        h.p.(plotName)(ksubj, kchar, 2) = plot(x, y, 'Marker', 'none', 'LineWidth', 0.5, 'Color', stg.subjColor(ksubj, :), 'Tag', 'sz');
        % % Plot signals
        % p = mod(siCharTbl{:, "tax"}, 1)*2*pi;
        % r = siCharTbl{:, stg.siCharToPlot(kchar)};
        % [x, y] = pol2cart(-p - pi/2, r);
        % h.p.(plotName)(ksubj, kchar, 2) = plot(x, y, 'Color', 0.7*[1 1 1], 'LineWidth', 0.5);
        % Plot signals' means
        p = [ppTbl{1, kchar}, ppTbl{1, kchar}(1)];
        r = [rrTbl{1, kchar}, rrTbl{1, kchar}(1)]/max(rrTbl{1, kchar});
        [x, y] = pol2cart(-p - pi/2, r);
        h.p.(plotName)(ksubj, kchar, 2) = plot(x, y, 'k', 'LineWidth', 1);
        % Circle
        circleP = (0 : 1 : 360)/360*2*pi;
        [x, y] = pol2cart(-circleP - pi/2, circleR);
        h.p.(plotName)(ksubj, kchar, 3) = plot(x, y, 'k');
        hax.XLim = max([x, y])*[-1 1];
        hax.YLim = hax.XLim;
        % Night shading
        circleP = (-90 : 1 : 90)/360*2*pi;
        % colo = [linspace(1, 0, 90), 1, linspace(0, 1, 90)]'*[1 1 1];
        colo = [linspace(1, 0.7, 90), 1, linspace(0.7, 1, 90)]'*[1 1 1];
        colo = permute(colo, [1 3 2]);
        [x, y] = pol2cart(-circleP - pi/2, circleR);
        z = -10*ones(size(x));
        h.pa.(plotName)(ksubj, kchar, 1) = patch(x, y, z, colo, 'EdgeColor', 'none');
    end
end
function numRealOut = plotPeriEventData(plotName, ksubj, x, yy, myy, plotColor)
    % plotName .. plot name derived from the name of the calling function
    % ksubj ..... index of subject
    % x ......... vector of x data
    % yy ........ table with y data, each cell of the table contains a vector (not a cell array of size 1x1)
    % myy ....... table similar to yy but with only one row containing the mean (or another statistic)
    global stg
    global h
    % Calculate axes positions
    [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlob, stg.marg);
    h.f.(plotName).Units = stg.units;
    % Get some basic variables
    numChar = numel(stg.([plotName([1, 2]), 'CharToPlot'])); % Usually either siCharToPlot or saCharToPlot
    numReal = size(yy, 1); % Number of realizations (e.g. number of seizures or number of subjects)
    % Create axes if needed
    if ~isfield(h.a, plotName)
        for kchar = 1 : numChar
            h.a.(plotName)(ksubj, kchar) = axes('Units', stg.units, 'NextPlot', 'add', 'Position', ...
                [spx(mod(ksubj - 1, numc) + 1), spy(ceil(ksubj/numc)) - (kchar - numChar)*spHe/numChar, spWi, spHe/numChar]);
            h.a.(plotName)(ksubj, kchar).Units = 'normalized';
        end
    else
        if size(h.a.(plotName), 1) < ksubj
            for kchar = 1 : numChar
                h.a.(plotName)(ksubj, kchar) = axes('Units', stg.units, 'NextPlot', 'add', 'Position', ...
                    [spx(mod(ksubj - 1, numc) + 1), spy(ceil(ksubj/numc)) - (kchar - numChar)*spHe/numChar, spWi, spHe/numChar]);
                h.a.(plotName)(ksubj, kchar).Units = 'normalized';
            end
        end
    end
    % Plot proper
    numRealOut = NaN(1, numChar);
    for kchar = 1 : numChar
        hax = h.a.(plotName)(ksubj, kchar);
        if isempty(yy)
            numRealOut(kchar) = 0;
        else
            numRealOut(kchar) = sum(~all(isnan(yy{:, kchar}), 2));
        end
        for kreal = 1 : numReal
            plot(hax, x, yy{kreal, kchar}, 'LineWidth', 0.5, 'Color', plotColor(kreal, :))
        end
        plot(hax, x, myy{1, kchar}, 'LineWidth', 1, 'Color', 'k')
        % % % % Plot mean of other
        % % % plot(hax, [x(1), x(end)], [yOther{1, kchar}, yOther{1, kchar}], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
        % % % hax.YScale = 'log';
    end
end
function numRealOut = plotPeriEventDataWithOther(plotName, ksubj, x, yy, myy, yOther, plotColor)
    % plotName .. plot name derived from the name of the calling function
    % ksubj ..... index of subject
    % x ......... vector of x data
    % yy ........ table with y data, each cell of the table contains a vector (not a cell array of size 1x1)
    % myy ....... table similar to yy but with only one row containing the mean (or another statistic)
    global stg
    global h
    % Calculate axes positions
    [spx, spy, spWi, spHe, ~, numc] = getSubplotXYWH(plotName, stg.margGlob, stg.marg);
    h.f.(plotName).Units = stg.units;
    % Get some basic variables
    % numChar = numel(stg.siCharToPlot);
    numChar = numel(stg.siCharToPlot);
    numReal = size(yy, 1); % Number of realizations (e.g. number of seizures or number of subjects)
    % Create axes if needed
    if ~isfield(h.a, plotName)
        for kchar = 1 : numChar
            h.a.(plotName)(ksubj, kchar) = axes('Units', stg.units, 'NextPlot', 'add', 'Position', ...
                [spx(mod(ksubj - 1, numc) + 1), spy(ceil(ksubj/numc)) - (kchar - numChar)*spHe/numChar, spWi, spHe/numChar]);
            h.a.(plotName)(ksubj, kchar).Units = 'normalized';
        end
    else
        if size(h.a.(plotName), 1) < ksubj
            for kchar = 1 : numChar
                h.a.(plotName)(ksubj, kchar) = axes('Units', stg.units, 'NextPlot', 'add', 'Position', ...
                    [spx(mod(ksubj - 1, numc) + 1), spy(ceil(ksubj/numc)) - (kchar - numChar)*spHe/numChar, spWi, spHe/numChar]);
                h.a.(plotName)(ksubj, kchar).Units = 'normalized';
            end
        end
    end
    % Plot proper
    numRealOut = NaN(1, numChar);
    for kchar = 1 : numChar
        hax = h.a.(plotName)(ksubj, kchar);
        if isempty(yy)
            numRealOut(kchar) = 0;
        else
            numRealOut(kchar) = sum(~all(isnan(yy{:, kchar}), 2));
        end
        for kreal = 1 : numReal
            plot(hax, x, yy{kreal, kchar}, 'LineWidth', 0.5, 'Color', plotColor(kreal, :))
        end
        plot(hax, x, myy{1, kchar}, 'LineWidth', 1, 'Color', 'k')
        % Plot mean of other
        plot(hax, [x(1), x(end)], [yOther{1, kchar}, yOther{1, kchar}], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
        hax.YScale = 'log';
    end
end
function formatAxesInd(plotName, ksubj, subjInfo, charTblVarNames, xlbl, xlimMethod , varargin)
    % Format axes in plots of individual subjects
    % plotName .... so that we know which axes to format
    % ksubj .......... index of subject
    % subjInfo ....... tells the analysis start and end, used for x-axis limits
    % charTbl ..... seizure characteristics or signal characteristics table, used to get variable names
    % xlbl ........ x-axis label
    % xlimMethod .. "none", "analysisPeriod", "XData"
    % ylimMethod .. optional, matrix of strings, e.g. ["nonpos", "nonneg"], the height can be 1 if it is the same for all characteristics
    % varargin .... ylimMethod and extension of the y-axis label (e.g. PSD (dB))
    global stg
    global h
    drawnow
    extension = '';
    if ~isempty(varargin)
        extension = varargin{1};
    end
    type = plotName(1 : 2);
    if strcmp(plotName(3 : 6), 'Char')
        numChar = numel(stg.([type, 'CharToPlot']));
        if strcmp(plotName, 'siCharPsd')
            numChar = 2;
        end
    else
        disp(plotName)
        stackInfo = dbstack;
        % Display the line number
        fprintf('This code is at line %d\n', stackInfo(1).line);
        pause
        numChar = 1; % Not used so far. Probably could be used for the ISI histograms, Karoly plot, etc.
    end
    for kchar = 1 : numChar
        hax = h.a.(plotName)(ksubj, kchar);
        % Title
        if kchar == 1
            hax.Title.String = subjInfo.subjNm;
            hax.Title.Interpreter = 'none';
            % hax.Title.Color = stg.subjColor(ksubj, :);
            hax.Title.Color = 'k';
            if ~subjInfo.sex
                hax.Title.FontAngle = 'italic';
            end
        end
        % x-axis
        switch xlimMethod
            case "none"
            case "analysisPeriod"
                xmi = subjInfo.anStartN - subjInfo.dob;
                % xmi = round(xmi - 10^floor(log10(xmi))/2, 1, "significant", TieBreaker = "tozero");
                xmi = floor(xmi);
                xma = subjInfo.anEndN - subjInfo.dob;
                % xma = round(xma + 10^floor(log10(xma))/2, 1, "significant", TieBreaker = "tozero");
                xma = ceil(xma);
                hax.XLim = [xmi, xma];
            case "xdata"
                chldrn = hax.Children;
                chldrn = chldrn(arrayfun(@(obj) isa(obj, 'matlab.graphics.chart.primitive.Line'), chldrn));
                xmi = min([chldrn.XData]);
                xma = max([chldrn.XData]);
                if xmi ~= 0
                    xmi = round(xmi, 2, "significant");
                end
                if xma ~= 0
                    xma = round(xma, 2, "significant");
                end
                hax.XLim = [xmi, xma];
                % hax.XTick = [xmi, xmi/2, 0, xma/2, xma];
                hax.XTick = xmi : (xma-xmi)/4 : xma;
            otherwise
                disp(['xlimMethod is', char(xlimMethod)])
                error('_jk formatAxesAllPop Unknown xlimMethod')
        end
        if kchar ~= numChar
            hax.XTickLabel = [];
        end
        % y-axis
        numChildren = numel(hax.Children);
        ymi = NaN(numChildren, 1); yma = NaN(numChildren, 1);
        for kchild = 1 : numChildren
            if any(contains(properties(hax.Children(kchild)), 'YData'))
                ymi(kchild) = min(hax.Children(kchild).YData(:));
                yma(kchild) = max(hax.Children(kchild).YData(:));
            else
                ymi(kchild) = NaN;
                yma(kchild) = NaN;
            end
        end
        ymi = min(ymi, [], 'omitmissing');
        yma = max(yma, [], 'omitmissing');
        if ymi == -Inf
            ymi = 0;
        end
        if ~contains(plotName, 'Psd')
            [yl, yt] = getYLimYTick([ymi yma], stg.([type, 'CharYLim'])(kchar, :)); % Input could be the whole signal but getYLimYTick takes the min and max of it anyway
        else
            [yl, yt] = getYLimYTick([ymi yma], ["any", "any"]); % Input could be the whole signal but getYLimYTick takes the min and max of it anyway
        end
        hax.YLim = yl;
        hax.YTick = yt;
        hax.TickLength = [0.01 0.01];
        if mod(ksubj - 1, stg.sbNCol) == 0
            ylbl = getYlbl(plotName, charTblVarNames, kchar);
            hax.YLabel.String = [ylbl, extension];
            if numChar > 1
                hax.YLabel.Rotation = -35;
                hax.YLabel.Position(1) = hax.YLabel.Position(1) - (hax.XLim(1) - hax.YLabel.Position(1))/3;
                hax.YLabel.HorizontalAlignment = 'right';
                hax.YLabel.VerticalAlignment = 'middle';
            end
        end
        hax.FontSize = stg.axFontSize;
        hax.Layer = 'top';
        hax.Box = stg.box;
        hax.Units = 'normalized';
    end
    if ksubj > stg.numSubj - stg.sbNCol
        if strcmpi(xlbl, 'Time (hours)') && hax.XLim(2) <= 1
            hax.XTickLabel = cellfun(@(x) num2str(str2double(x)*60), hax.XTickLabel, 'UniformOutput', false);
            xlabel('Time (minutes)');
        else
            xlabel(xlbl)
        end
    end
end
function formatAxesAllPop(plotName, charTblVarNames, xlbl, xlimMethod, varargin)
    % Format axes in plots of all subjects in one axes and population mean (or other statistic)
    % plotName .... so that we know which axes to format
    % charTbl ..... seizure characteristics or signal characteristics table, used to get variable names
    % xlbl ........ x-axis label
    % xlimMethod .. "none", "analysisPeriod", "XData"
    % ylimMethod .. optional, matrix of strings, e.g. ["nonpos", "nonneg"], the height can be 1 if it is the same for all characteristics
    % varargin .... ylimMethod and extension of the y-axis label (e.g. PSD (dB))
    global stg
    global h
    drawnow
    type = plotName(1 : 2); % si, sz or sa, signal char, seizure char or all
    switch numel(varargin)
        case 0
            ylimMethod = [];
            extension = [];
        case 1
            ylimMethod = varargin{1};
            extension = [];
        case 2
            ylimMethod = varargin{1};
            extension = varargin{2};
    end
    if isempty(ylimMethod)
        ylimMethod = stg.([type, 'CharYLim']);
    end
    if size(ylimMethod, 1) == 1
        ylimMethod = repelem(ylimMethod, size(stg.([type, 'CharYLim']), 1), 1);
    end
    if isempty(extension)
        extension = '';
    end

    if strcmp(plotName(3 : 6), 'Char')
        numChar = numel(stg.([type, 'CharToPlot'])); % It is a plot of characteristics
    else
        numChar = 1; % Not used so far. Probably could be used e.g. for the ISI histogram
    end
    for kchar = 1 : numChar
        hax = h.a.(plotName)(1, kchar);
        if kchar == 1
            hax.Title.String = 'All mice';
            hax.Title.Interpreter = 'none';
            % hax.Title.Color = stg.subjColor(ksubj, :);
            hax.Title.Color = 'k';
        end

        % x-axis
        hax.XLimMode = "manual";
        switch xlimMethod
            case "none"
            case "analysisPeriod"
                xmi = subjInfo.anStartN - subjInfo.dob;
                % xmi = round(xmi - 10^floor(log10(xmi))/2, 1, "significant", TieBreaker = "tozero");
                xmi = floor(xmi);
                xma = subjInfo.anEndN - subjInfo.dob;
                % xma = round(xma + 10^floor(log10(xma))/2, 1, "significant", TieBreaker = "tozero");
                xma = ceil(xma);
                hax.XLim = [xmi, xma];
            case "xdata"
                xmi = min([hax.Children.XData]);
                xma = max([hax.Children.XData]);
                if xmi ~= 0
                    xmi = round(xmi, 2, "significant");
                end
                if xma ~= 0
                    xma = round(xma, 2, "significant");
                end
                hax.XLim = [xmi, xma];
                hax.XTick = xmi : (xma-xmi)/4 : xma;
            otherwise
                disp(['xlimMethod is', char(xlimMethod)])
                error('_jk formatAxesAllPop Unknown xlimMethod')
        end
        if kchar ~= numChar
            hax.XTickLabel = [];
        end
        % y-axis
        hax.YLimMode = "manual";
        numChildren = numel(hax.Children);
        ymi = NaN(numChildren, 1); yma = NaN(numChildren, 1);
        for kchild = 1 : numChildren
            if any(contains(properties(hax.Children(kchild)), 'YData'))
                ymi(kchild) = min(hax.Children(kchild).YData(:));
                yma(kchild) = max(hax.Children(kchild).YData(:));
            else
                ymi(kchild) = NaN;
                yma(kchild) = NaN;
            end
        end
        ymi = min(ymi, [], 'omitmissing');
        yma = max(yma, [], 'omitmissing');
        if all(~isnan([ymi, yma]))
            [yl, yt] = getYLimYTick([ymi yma], ylimMethod(kchar, :)); % Input could be the whole signal but getYLimYTick takes the min and max of it anyway
        else
            yl = [0 1];
            yt = [0 0.5];
        end
        hax.YLim = yl;
        hax.YTick = yt;
        hax.TickLength = [0.01 0.01];
        ylbl = getYlbl(plotName, charTblVarNames, kchar);
        hax.YLabel.String = [ylbl, extension];
        if numChar > 1
            hax.YLabel.Rotation = -35;
            hax.YLabel.Position(1) = hax.YLabel.Position(1) - (hax.XLim(1) - hax.YLabel.Position(1))/3;
            hax.YLabel.HorizontalAlignment = 'right';
            hax.YLabel.VerticalAlignment = 'middle';
        end
        hax.FontSize = stg.axFontSize;
        hax.Layer = 'top';
        hax.Box = stg.box;
        hax.Units = 'normalized';
    end
    if strcmpi(xlbl, 'Time (hours)') && hax.XLim(2) <= 1
        hax.XTickLabel = cellfun(@(x) num2str(str2double(x)*60), hax.XTickLabel, 'UniformOutput', false);
        xlabel(hax, 'Time (minutes)');
    else
        xlabel(hax, xlbl)
    end
    hax.Units = 'normalized';
end
function v = violinMedianConf(hax, y, facecolor, numBootstrapSamples, alpha)
    global stg
    v = violinplot(hax, y);
    v.FaceColor = facecolor;
    line(hax, v.DensityWidth*[-1, 1]/2 + 1, median(y, 'omitmissing')*[1, 1], 'LineWidth', 1, 'Color', facecolor)
    bootstrapMedians = NaN(stg.numBootstrapSamples, 1);
    % Generate bootstrap samples and calculate medians
    for kb = 1 : numBootstrapSamples
        bootstrapSample = y(randi(length(y), length(y), 1));
        bootstrapMedians(kb) = median(bootstrapSample);
    end
    ciLow = prctile(bootstrapMedians, 100*alpha / 2);
    ciHigh = prctile(bootstrapMedians, 100*(1 - alpha / 2));
    line(hax, v.DensityWidth*[-1, 1]/2 + 1, ciLow*[1, 1], 'LineWidth', 0.5, 'Color', facecolor, 'LineStyle', '--')
    line(hax, v.DensityWidth*[-1, 1]/2 + 1, ciHigh*[1, 1], 'LineWidth', 0.5, 'Color', facecolor, 'LineStyle', '--')
end
function printFigures
    global stg
    global h
    if ~stg.printFigures
        return
    end
    figNm = fieldnames(h.f);
    for k = 1 : numel(figNm)
        % Get figure number
        fn = fieldnames(stg);
        fn = fn(startsWith(fn, 'plot')); % Select only the fields of stg determining if a figure is plotted
        nFig = find(endsWith(fn, figNm{k}, 'IgnoreCase', true));
        figure(h.f.(figNm{k}))
        % filen = [char(datetime('now', 'Format', 'yyMMdd_HHmmss')), ' ', num2str(nFig, '%02d'), ...
        %    figNm{k}, ' N=', num2str(stg.numSubj, '%02d'), ' T=', num2str(stg.dpBinLenS, '%06d'),];
        filen = [figNm{k}, ' N=', num2str(stg.numSubj, '%02d'), ' T=', num2str(stg.dpBinLenS, '%06d'),];
        savefig(gcf, ['./_fig/', filen, '.fig'])
        print(['./_eps/', filen, '.eps'], '-depsc', '-vector')
        print(['./_jpg/', filen, '.jpg'], '-djpeg', '-r1800')
        print(['./_pdf/', filen, '.pdf'], '-dpdf', '-vector')
    end
end
function setFormat
    set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
    set(0, 'DefaultTextFontName', 'Arial');
    set(0, 'DefaultAxesFontName', 'Arial');
    set(0, 'DefaultUicontrolFontName', 'Arial');
    format compact
    datetime.setDefaultFormats('default', 'yyyy-MM-dd HH:mm:ss')
end

% Helper functions
function [yl, yt] = getYLimYTick(y, varargin)
    % y ......... signal to accommodate within the axes
    % varargin .. two-element string array indicating requirements for the the y-axis limits
    cemaab = ceilToEven1sig(max(abs(y(~isnan(y))))); % Ceiling of maximum of abs of the y
    if cemaab == 0 || isnan(cemaab)
        yl(1) = 0;
        yl(2) = 1;
    else
        if mod(log10(cemaab), 1) == 0
            numDec = -1*(floor(log10(cemaab)) - 1);
            additive = 0.99*10^(-numDec)/2; % The 0.9 ensures that the graph can go partially outside axes. If it is undesired, set it to 1.
        else
            numDec = -1*(floor(log10(cemaab)));
            additive = 0.99*10^(-numDec)/2;
        end
        mi = min(y);
        ma = max(y);
        yl(1) = ceilToEven1sig(ma);
        yl(1) = round(mi/2 - additive, numDec, "decimals", "TieBreaker", "tozero")*2;
        yl(2) = round(ma/2 + additive, numDec, "decimals", "TieBreaker", "tozero")*2;
        exponent = floor(log10(abs(yl(1))));
        mantissa = yl(1)/(10^exponent);
        if mantissa == -2 || mi >= yl(1)/2
            yl(1) = yl(1)/2;
        end
    end
    % Force limits if they are provided by the varargin
    if ~isempty(varargin)
        for kl = 1 : 2 % Lower and upper limit
            if ~isnan(str2double(varargin{1}(kl)))
                yl(kl) = str2double(varargin{1}(kl));
            else
                switch varargin{1}(kl)
                    case "nonpos"
                        yl(kl) = min(yl(kl), 0);
                    case "zero"
                        yl(kl) = 0;
                    case "nonneg"
                        yl(kl) = max(yl(kl), 0);
                    case "any"
                    otherwise
                        error('_jk getYLimYTick: Unknown axis limit requirement in varargin.')
                end
            end
        end
    end
    if yl(1) == yl(2)
        yl(2) = yl(1) + 1;
    end
    yt(1) = yl(1);
    yt(2) = mean(yl);
end
function y = ceilToEven1sig(x)
    % Ceiling to even number on 1 significant digit
    if numel(x) ~= 1
        y = NaN;
    else
        y = round((x/2 + 10^floor(log10(x/2))/2), 1, "significant", TieBreaker = "tozero")*2;
    end
end
function y = floorToEven1sig(x) %#ok<DEFNU> % DOES NOT WORK FOR NEGATIVE NUMBERS
    % Floor to even number on 1 significant digit
    if numel(x) ~= 1
        y = NaN;
    else
        y = round((x/2 - 10^floor(log10(x/2))/2), 1, "significant", TieBreaker = "tozero")*2;
    end
end
function fitTbl = fitTblInit(numRows, charToPlot)
    numChar = numel(charToPlot);
    % Declare and initialize the stats table
    fitTbl = table('Size', [numRows, 6*numChar], 'VariableTypes', repelem("double", 6*numChar));
    fitTbl{1, :} = deal(NaN);
    for kchar = 1 : numChar
        fitTbl.Properties.VariableNames{1, kchar + 0*numChar} = char(charToPlot(kchar) + "fitSlope");
        fitTbl.Properties.VariableNames{1, kchar + 1*numChar} = char(charToPlot(kchar) + "fitSlopeCIL");
        fitTbl.Properties.VariableNames{1, kchar + 2*numChar} = char(charToPlot(kchar) + "fitSlopeCIH");
        fitTbl.Properties.VariableNames{1, kchar + 3*numChar} = char(charToPlot(kchar) + "fitOffset");
        fitTbl.Properties.VariableNames{1, kchar + 4*numChar} = char(charToPlot(kchar) + "fitOffsetCIL");
        fitTbl.Properties.VariableNames{1, kchar + 5*numChar} = char(charToPlot(kchar) + "fitOffsetCIH");
    end        
end
function cirTbl = cirTblInit(numRows, charToPlot)
    numChar = numel(charToPlot);
    % Declare and initialize the stats table
    cirTbl = table('Size', [numRows, 6*numChar], 'VariableTypes', repelem("double", 6*numChar));
    cirTbl{1, :} = deal(NaN);
    for kchar = 1 : numChar
        cirTbl.Properties.VariableNames{1, kchar + 0*numChar} = char(charToPlot(kchar) + "cirMean"); % Direction
        cirTbl.Properties.VariableNames{1, kchar + 1*numChar} = char(charToPlot(kchar) + "cirMeanCIL");
        cirTbl.Properties.VariableNames{1, kchar + 2*numChar} = char(charToPlot(kchar) + "cirMeanCIH");
        cirTbl.Properties.VariableNames{1, kchar + 3*numChar} = char(charToPlot(kchar) + "cirR"); % Length
        cirTbl.Properties.VariableNames{1, kchar + 4*numChar} = char(charToPlot(kchar) + "cirRCIL");
        cirTbl.Properties.VariableNames{1, kchar + 5*numChar} = char(charToPlot(kchar) + "cirRCIH");
    end        
end
function stats = fillInFit(stats, k, kchar, numChar, x, y, varargin)
    % stats ..... statistics table, usually contains the fit and other data
    % k ......... row in the table, usually index of the cluster or seizure
    % kchar ..... index of the characteristic fitted
    % numChar ... total number of characteristics analyzed
    % x ......... x data to fit
    % y ......... y data to fit
    % varargin .. weigths of each data point (if the data point comes from the data with a lot of dropouts, we assign it lower weight)
    if ~isempty(varargin)
        binWeights = varargin{1};
        binWeights = binWeights(~isnan(y));
        binWeights = binWeights(:);
    end
    xf = x(~isnan(y)); % x to fit
    yf = y(~isnan(y)); % y to fit
    xf = xf(:);
    yf = yf(:);
    if numel(yf) >= 3
        if ~isempty(varargin)
            xf(isnan(binWeights)) = [];
            yf(isnan(binWeights)) = [];
            binWeights(isnan(binWeights)) = [];
            ff = fit(xf, yf, "poly1", "Weights", binWeights);
        else
            ff = fit(xf, yf, "poly1");
        end
        ci = confint(ff);
        fitSlope = ff.p1;
        fitOffset = ff.p2;
    elseif numel(yf) == 2
        if ~isempty(varargin)
            ff = fit(xf, yf, "poly1", "Weights", binWeights);
        else
            ff = fit(xf, yf, "poly1");
        end
        ci = NaN(2, 2);
        fitSlope = ff.p1;
        fitOffset = ff.p2;
    else
        ci = NaN(2, 2);
        fitSlope = NaN;
        fitOffset = NaN;
    end
    stats{k, kchar + 0*numChar} = fitSlope;
    stats{k, kchar + 1*numChar} = ci(1, 1);
    stats{k, kchar + 2*numChar} = ci(2, 1);
    stats{k, kchar + 3*numChar} = fitOffset;
    stats{k, kchar + 4*numChar} = ci(1, 2);
    stats{k, kchar + 5*numChar} = ci(2, 2);
end
function stats = fillInFitResultantVector(stats, k, kchar, numChar, p, r, szOccurrenceTF)
    global stg
    % stats ..... statistics table, usually contains the fit and other data
    % k ......... row in the table, usually index of the cluster or seizure or just 1 it the whole recording is analyzed
    % kchar ..... index of the characteristic fitted
    % numChar ... total number of characteristics analyzed
    % p ......... 
    % f ......... 
    pf = p(~isnan(r)); % x to fit
    rf = r(~isnan(r)); % y to fit
    pf = pf(:);
    rf = rf(:);
    if all(isnan(rf))
        pf = NaN;
        rf = NaN;
    end
    % My method
    % % % resultantVector = rf'*exp(1i*pf); % It sums the vector due to matrix-style multiplication
    % % % p0 = angle(resultantVector)
    % % % r0 = abs(resultantVector)
    % % % r00 = r0/sum(rf)
    
    % Circular statistics toolbox method (should be equal to my method but provides also confidence limits)
    [p, plu, pll] = circ_mean(pf, rf);
        % % % % % % R = circ_r(pf, rf)
        % % % % % % T = circ_t(pf, rf)
        % % % % % % pause

    if stg.andrzejak && szOccurrenceTF
        r = circ_t(pf, rf);
    else
        r = circ_r(pf, rf);
    end
    stats{k, kchar + 0*numChar} = p;
    stats{k, kchar + 1*numChar} = pll;
    stats{k, kchar + 2*numChar} = plu;
    stats{k, kchar + 3*numChar} = r;
    stats{k, kchar + 4*numChar} = NaN;
    stats{k, kchar + 5*numChar} = NaN;
end
function ylbl = getYlbl(plotName, charTblVarNames, kchar)
    global stg
    type = plotName(1 : 2); % si, sz, sa or ss
    if strcmpi(charTblVarNames(1), 'sz') && kchar == 1
        ylbl = 'Sz rate';
        return
    elseif strcmpi(charTblVarNames(1), 'sz') && kchar == 2
        ylbl = 'IED rate';
        return
    end
    ylbl = stg.([type, 'CharYLbl']){strcmp(charTblVarNames, stg.([type, 'CharToPlot']){kchar})};
    if (contains(plotName, 'szChar') && ~strcmp(plotName, 'szChar'))   ||   (contains(plotName, 'saChar') && ~strcmp(plotName, 'saChar'))
        if strcmp(ylbl, 'Sz time')
            ylbl = 'Sz rate (day^{-1})';
        end
    end
end
function plotExplanatoryCircle(hax)
    global stg
    text(hax, 0, -0.65, 'midnight', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', stg.statFontSize)
    % text(hax, -0.65, 0, 'morning', 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', stg.statFontSize)
    text(hax, 0, 0.65, 'noon', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', stg.statFontSize)
    % text(hax, 0.65, 0, 'afternoon', 'Rotation', 270, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', stg.statFontSize)
    % Left arrow
    p = (15 : 35)/100;
    p = p*2*pi;
    r = 0.7*ones(size(p));
    [x, y] = pol2cart(-p - pi/2, r);
    plot(hax, x, y, 'k', 'LineWidth', 1)
    p1 = p(end) - 3/4*pi;
    r1 = 0.2;
    [x1, y1] = pol2cart(-p1 - pi/2, r1);
    x1 = x(end) + [0, x1];
    y1 = y(end) + [0, y1];
    plot(hax, x1, y1, 'k', 'LineWidth', 1)
    p2 = p(end) - 1/3*pi;
    r2 = 0.2;
    [x2, y2] = pol2cart(-p2 - pi/2, r2);
    x2 = x(end) + [0, x2];
    y2 = y(end) + [0, y2];
    plot(hax, x2, y2, 'k', 'LineWidth', 1)
    % Right arrow
    p = (65 : 85)/100;
    p = p*2*pi;
    r = 0.7*ones(size(p));
    [x, y] = pol2cart(-p - pi/2, r);
    plot(hax, x, y, 'k', 'LineWidth', 1)
    p1 = p(end) - 3/4*pi;
    r1 = 0.2;
    [x1, y1] = pol2cart(-p1 - pi/2, r1);
    x1 = x(end) + [0, x1];
    y1 = y(end) + [0, y1];
    plot(hax, x1, y1, 'k', 'LineWidth', 1)
    p2 = p(end) - 1/3*pi;
    r2 = 0.2;
    [x2, y2] = pol2cart(-p2 - pi/2, r2);
    x2 = x(end) + [0, x2];
    y2 = y(end) + [0, y2];
    plot(hax, x2, y2, 'k', 'LineWidth', 1)
end
function IMGD = dilateN(IMG, MASK)
    % IMG ... n-dim binary image (logical)
    % MASK ... dilate n-dim mask (logical) or size of dilatation mask (pixels)
    %
    % Example:
    % IMGD=dilateN(IMG,10) extends boundaries of 5 pixels 

    if ~islogical(IMG)
       error('only for logical image') 
    end    

    if numel(MASK)==1
        MASK=single(ones(MASK*ones(1,length(size(IMG)))));
    end

    IMG=single(IMG);

    IMGD=convn(IMG,MASK,'same')>0;
end
function IMGE = erodeNjk(IMG, MASK)
    % IMG ... n-dim binary image (logical)
    % MASK ... erode n-dim mask (logical) or size of dilatation mask (pixels)
    %
    % Example:
    % IMGE=dilateN(IMG,10) crops boundaries of 5 pixels 

    if ~islogical(IMG)
       error('only for logical image') 
    end    

    if numel(MASK)==1
        MASK=true(MASK*ones(1,length(size(IMG))));
    end

    % IMG=single(IMG);

    % IMGE=convn(IMG,MASK,'same')==sum(MASK(:));
    IMGE=~rot90(logical(convn(~rot90(IMG, 2), MASK, 'same')), 2);
end
function baseline = getBaseline(y)
    baseline = mean(y(end-min(numel(y)-1, 5) : end));
end
function [pValues, isSignificant] = getSignificanceFromP(p, alpha, varargin)
    % pValues ........ corrected p-values
    % isSignificant .. decision whether given p-values means a significant effect
    % p .............. column vector or matrix of p-values
    % alpha .......... significance level or false discovery rate (usually called Q)
    % varargin ....... can be empty, "none", "Bonferroni", "BenjaminiHochberg"
    % SETTINGS
    actOnEachRowSeparately = true;
    %
    numVarargin = numel(varargin);
    switch numVarargin
        case 0
            multCompCorrection = "none";
            actOnEachRowSeparately = true;
        case 1
            multCompCorrection = varargin{1};
            actOnEachRowSeparately = true;
        case 2
            multCompCorrection = varargin{1};
            actOnEachRowSeparately = varargin{2};
    end
    switch multCompCorrection
        case "none"
            pValues = p;
            isSignificant = pValues < alpha;
        case "Bonferroni"
            if actOnEachRowSeparately
                pValues = p*size(p, 2);
                isSignificant = pValues < alpha;
            else
                pValues = p*numel(p);
                isSignificant = pValues < alpha;
            end
        case "BenjaminiHochberg"
            % Transpose, computes everything in columns instead of rows, transpose back
            pValues = p';
            % % % sigRanP_builtin_mafdr = mafdr(sigRanP', 'BHFDR', true)
            if actOnEachRowSeparately
                for kc = 1 : size(pValues, 2)
                    p = pValues(:, kc);
                    nump = numel(p);
                    pfdr = sortrows([(1 : nump)', p], 2); % Matrix with ascending column and column of p-values. The first column will be used to sort out p-values after the whole process.
                    pfdr = [pfdr, (1 : nump)'/nump*alpha, false(nump, 1)]; %#ok<AGROW> % Add two more columns - the Benjamini-Hochberg limit and a column of false
                    lastSignificant = find(pfdr(:, 2) < pfdr(:, 3), 1, 'last');
                    pfdr(1 : lastSignificant, 4) = true;
                    pfdr = sortrows(pfdr, 1);
                    isSignificant(:, kc) = pfdr(:, 4)'; %#ok<AGROW>
                end
            else
                origShape = size(p); % Original shape of the matrix of p-values
                p = p(:);
                nump = numel(p);
                pfdr = sortrows([(1 : nump)', p], 2); % Matrix with ascending column and column of p-values. The first column will be used to sort out p-values after the whole process.
                pfdr = [pfdr, (1 : nump)'/nump*alpha, false(nump, 1)]; % Add two more columns - the Benjamini-Hochberg limit and a column of false
                lastSignificant = find(pfdr(:, 2) < pfdr(:, 3), 1, 'last');
                pfdr(1 : lastSignificant, 4) = true;
                pfdr = sortrows(pfdr, 1);
                isSignificant = reshape(pfdr(:, 4)', origShape);
            end
            pValues = pValues';
            isSignificant = logical(isSignificant)';
    end
end
function signChar = getSignChar(x)
    if isnan(sign(x))
        signChar = 'n/a';
    else
        switch sign(x)
            case -1
                signChar = '-';
            case 0
                signChar = '0';
            case 1
                signChar = '+';
        end
    end
end

% Text display functions
function descriptiveStats(tbl)
    global stg
    if stg.showStat
        disp([newline, '-------------------------------------------------------------'])
        disp('Basic descriptive statistics:')
        disp('(descriptiveStats)')
        if stg.sdTF
            disp('mean±SD (median±IQR)')
        else
            disp('mean±SEM (median±IQR)')
        end
        disp(['Total number of subjects: ', num2str(size(tbl, 1))])
        for k = 2 : size(tbl, 2)
            varNm = tbl.Properties.VariableNames{k};
            m = mean(tbl{:, k}, 1, 'omitmissing');
            sem = std(tbl{:, k}, 1, 'omitmissing')/sqrt(sum(~isnan(tbl{:, k})));
            sdv = std(tbl{:, k}, 1, 'omitmissing');
            med = median(tbl{:, k}, 1, 'omitmissing');
            iql = iqr(tbl{:, k}, 1);
            if stg.sdTF
                resultStr = ...
                    [repelem(' ', 1, 20-numel(varNm)), varNm, ' = ', num2str(m, '%.2g'), '±', num2str(sdv, '%.2g'), ' (', num2str(med, '%.2g'), '±', num2str(iql, '%.2g'), ')'];
            else
                resultStr = ...
                    [repelem(' ', 1, 20-numel(varNm)), varNm, ' = ', num2str(m, '%.2g'), '±', num2str(sem, '%.2g'), ' (', num2str(med, '%.2g'), '±', num2str(iql, '%.2g'), ')'];
            end
            disp(resultStr)
        end
    end
end

% Functions from the Internet
function [alpha, xmin, L] = plfit(x, varargin)
    % PLFIT fits a power-law distributional model to data.
    %    Source: http://www.santafe.edu/~aaronc/powerlaws/
    % 
    %    PLFIT(x) estimates x_min and alpha according to the goodness-of-fit
    %    based method described in Clauset, Shalizi, Newman (2007). x is a 
    %    vector of observations of some quantity to which we wish to fit the 
    %    power-law distribution p(x) ~ x^-alpha for x >= xmin.
    %    PLFIT automatically detects whether x is composed of real or integer
    %    values, and applies the appropriate method. For discrete data, if
    %    min(x) > 1000, PLFIT uses the continuous approximation, which is 
    %    a reliable in this regime.
    %   
    %    The fitting procedure works as follows:
    %    1) For each possible choice of x_min, we estimate alpha via the 
    %       method of maximum likelihood, and calculate the Kolmogorov-Smirnov
    %       goodness-of-fit statistic D.
    %    2) We then select as our estimate of x_min, the value that gives the
    %       minimum value D over all values of x_min.
    %
    %    Note that this procedure gives no estimate of the uncertainty of the 
    %    fitted parameters, nor of the validity of the fit.
    %
    %    Example:
    %       x = (1-rand(10000,1)).^(-1/(2.5-1));
    %       [alpha, xmin, L] = plfit(x);
    %
    %    The output 'alpha' is the maximum likelihood estimate of the scaling
    %    exponent, 'xmin' is the estimate of the lower bound of the power-law
    %    behavior, and L is the log-likelihood of the data x>=xmin under the
    %    fitted power law.
    %    
    %    For more information, try 'type plfit'
    %
    %    See also PLVAR, PLPVA
    
    % Version 1.0   (2007 May)
    % Version 1.0.2 (2007 September)
    % Version 1.0.3 (2007 September)
    % Version 1.0.4 (2008 January)
    % Version 1.0.5 (2008 March)
    % Version 1.0.6 (2008 July)
    % Version 1.0.7 (2008 October)
    % Version 1.0.8 (2009 February)
    % Version 1.0.9 (2009 October)
    % Copyright (C) 2008-2009 Aaron Clauset (Santa Fe Institute)
    % Distributed under GPL 2.0
    % http://www.gnu.org/copyleft/gpl.html
    % PLFIT comes with ABSOLUTELY NO WARRANTY
    % 
    % Notes:
    % 
    % 1. In order to implement the integer-based methods in Matlab, the numeric
    %    maximization of the log-likelihood function was used. This requires
    %    that we specify the range of scaling parameters considered. We set
    %    this range to be [1.50 : 0.01 : 3.50] by default. This vector can be
    %    set by the user like so,
    %    
    %       a = plfit(x,'range',[1.001:0.001:5.001]);
    %    
    % 2. PLFIT can be told to limit the range of values considered as estimates
    %    for xmin in three ways. First, it can be instructed to sample these
    %    possible values like so,
    %    
    %       a = plfit(x,'sample',100);
    %    
    %    which uses 100 uniformly distributed values on the sorted list of
    %    unique values in the data set. Second, it can simply omit all
    %    candidates above a hard limit, like so
    %    
    %       a = plfit(x,'limit',3.4);
    %    
    %    Finally, it can be forced to use a fixed value, like so
    %    
    %       a = plfit(x,'xmin',3.4);
    %    
    %    In the case of discrete data, it rounds the limit to the nearest
    %    integer.
    % 
    % 3. When the input sample size is small (e.g., < 100), the continuous 
    %    estimator is slightly biased (toward larger values of alpha). To
    %    explicitly use an experimental finite-size correction, call PLFIT like
    %    so
    %    
    %       a = plfit(x,'finite');
    %    
    %    which does a small-size correction to alpha.
    %
    % 4. For continuous data, PLFIT can return erroneously large estimates of 
    %    alpha when xmin is so large that the number of obs x >= xmin is very 
    %    small. To prevent this, we can truncate the search over xmin values 
    %    before the finite-size bias becomes significant by calling PLFIT as
    %    
    %       a = plfit(x,'nosmall');
    %    
    %    which skips values xmin with finite size bias > 0.1.
    
    vec     = [];
    sample  = [];
    xminx   = [];
    limit   = [];
    finite  = false;
    nosmall = false;
    nowarn  = false;
    
    % parse command-line parameters; trap for bad input
    i=1; 
    while i<=length(varargin)
      argok = 1; 
      if ischar(varargin{i}) 
        switch varargin{i}
            case 'range',        vec     = varargin{i+1}; i = i + 1;
            case 'sample',       sample  = varargin{i+1}; i = i + 1;
            case 'limit',        limit   = varargin{i+1}; i = i + 1;
            case 'xmin',         xminx   = varargin{i+1}; i = i + 1;
            case 'finite',       finite  = true;    i = i + 1;
            case 'nowarn',       nowarn  = true;    i = i + 1;
            case 'nosmall',      nosmall = true;    i = i + 1;
            otherwise, argok=0; 
        end
      end
      if ~argok
        disp(['(PLFIT) Ignoring invalid argument #' num2str(i+1)]); 
      end
      i = i+1; 
    end
    if ~isempty(vec) && (~isvector(vec) || min(vec)<=1)
	    fprintf('(PLFIT) Error: ''range'' argument must contain a vector; using default.\n');
        vec = [];
    end
    if ~isempty(sample) && (~isscalar(sample) || sample<2)
	    fprintf('(PLFIT) Error: ''sample'' argument must be a positive integer > 1; using default.\n');
        sample = [];
    end
    if ~isempty(limit) && (~isscalar(limit) || limit<min(x))
	    fprintf('(PLFIT) Error: ''limit'' argument must be a positive value >= 1; using default.\n');
        limit = [];
    end
    if ~isempty(xminx) && (~isscalar(xminx) || xminx>=max(x))
	    fprintf('(PLFIT) Error: ''xmin'' argument must be a positive value < max(x); using default behavior.\n');
        xminx = [];
    end
    
    % reshape input vector
    x = reshape(x,numel(x),1);
    
    % select method (discrete or continuous) for fitting
    if     isempty(setdiff(x,floor(x))), f_dattype = 'INTS';
    elseif isreal(x),    f_dattype = 'REAL';
    else                 f_dattype = 'UNKN'; %#ok<SEPEX>
    end
    if strcmp(f_dattype,'INTS') && min(x) > 1000 && length(x)>100
        f_dattype = 'REAL';
    end
    
    % estimate xmin and alpha, accordingly
    switch f_dattype
        case 'REAL'
            xmins = unique(x);
            xmins = xmins(1:end-1);
            if ~isempty(xminx)
                xmins = xmins(find(xmins>=xminx,1,'first'));
            end
            if ~isempty(limit)
                xmins(xmins>limit) = [];
            end
            if ~isempty(sample)
                xmins = xmins(unique(round(linspace(1,length(xmins),sample))));
            end
            dat   = zeros(size(xmins));
            z     = sort(x);
            for xm=1:length(xmins)
                xmin = xmins(xm);
                z    = z(z>=xmin); 
                n    = length(z);
                % estimate alpha using direct MLE
                a    = n ./ sum( log(z./xmin) );
                if nosmall
                    if (a-1)/sqrt(n) > 0.1
                        dat(xm:end) = [];
                        xm = length(xmins)+1; %#ok<FXSET,NASGU>
                        break;
                    end
                end
                % compute KS statistic
                cx   = (0:n-1)'./n;
                cf   = 1-(xmin./z).^a;
                dat(xm) = max( abs(cf-cx) );
            end
            D     = min(dat);
            xmin  = xmins(find(dat<=D,1,'first'));
            z     = x(x>=xmin);
            n     = length(z); 
            alpha = 1 + n ./ sum( log(z./xmin) );
            if finite, alpha = alpha*(n-1)/n+1/n; end % finite-size correction
            if n < 50 && ~finite && ~nowarn
                fprintf('(PLFIT) Warning: finite-size bias may be present.\n');
            end
            L = n*log((alpha-1)/xmin) - alpha.*sum(log(z./xmin));
    
        case 'INTS'
            if isempty(vec)
                vec  = (1.50:0.01:3.50);    % covers range of most practical 
            end                             % scaling parameters
            zvec = zeta(vec);
    
            xmins = unique(x);
            xmins = xmins(1:end-1);
            if ~isempty(xminx)
                xmins = xmins(find(xmins>=xminx,1,'first'));
            end
            if ~isempty(limit)
                limit = round(limit);
                xmins(xmins>limit) = [];
            end
            if ~isempty(sample)
                xmins = xmins(unique(round(linspace(1,length(xmins),sample))));
            end
            if isempty(xmins)
                fprintf('(PLFIT) Error: x must contain at least two unique values.\n');
                alpha = NaN; xmin = x(1); D = NaN; %#ok<NASGU>
                return;
            end
            xmax   = max(x);
            dat    = zeros(length(xmins),2);
            z      = x;
            fcatch = 0;
    
            for xm=1:length(xmins)
                xmin = xmins(xm);
                z    = z(z>=xmin);
                n    = length(z);
                % estimate alpha via direct maximization of likelihood function
                if fcatch==0
                    try
                        % vectorized version of numerical calculation
                        zdiff = sum( repmat((1:xmin-1)',1,length(vec)).^-repmat(vec,xmin-1,1) ,1);
                        L = -vec.*sum(log(z)) - n.*log(zvec - zdiff);
                    catch
                        % catch: force loop to default to iterative version for
                        % remainder of the search
                        fcatch = 1;
                    end
                end
                if fcatch==1
                    % force iterative calculation (more memory efficient, but 
                    % can be slower)
                    L       = -Inf*ones(size(vec));
                    slogz   = sum(log(z));
                    xminvec = (1:xmin-1);
                    for k=1:length(vec)
                        L(k) = -vec(k)*slogz - n*log(zvec(k) - sum(xminvec.^-vec(k)));
                    end
                end
                [Y,I] = max(L);
                % compute KS statistic
                fit = cumsum((((xmin:xmax).^-vec(I)))./ (zvec(I) - sum((1:xmin-1).^-vec(I))));
                cdi = cumsum(hist(z,xmin:xmax)./n); %#ok<HIST>
                dat(xm,:) = [max(abs( fit - cdi )) vec(I)];
            end
            % select the index for the minimum value of D
            [D,I] = min(dat(:,1)); %#ok<ASGLU>
            xmin  = xmins(I);
            n     = sum(x>=xmin);
            alpha = dat(I,2);
            if finite, alpha = alpha*(n-1)/n+1/n; end % finite-size correction
            if n < 50 && ~finite && ~nowarn
                fprintf('(PLFIT) Warning: finite-size bias may be present.\n');
            end
            L = Y;
    
        otherwise
            fprintf('(PLFIT) Error: x must contain only reals or only integers.\n');
            alpha = [];
            xmin  = [];
            L     = [];
            return;
    end
end
function [p, gof] = plpva(x, xmin, varargin)
% PLPVA calculates the p-value for the given power-law fit to some data.
%    Source: http://www.santafe.edu/~aaronc/powerlaws/
% 
%    PLPVA(x, xmin) takes data x and given lower cutoff for the power-law
%    behavior xmin and computes the corresponding p-value for the
%    Kolmogorov-Smirnov test, according to the method described in 
%    Clauset, Shalizi, Newman (2007).
%    PLPVA automatically detects whether x is composed of real or integer
%    values, and applies the appropriate method. For discrete data, if
%    min(x) > 1000, PLPVA uses the continuous approximation, which is 
%    a reliable in this regime.
%   
%    The fitting procedure works as follows:
%    1) For each possible choice of x_min, we estimate alpha via the 
%       method of maximum likelihood, and calculate the Kolmogorov-Smirnov
%       goodness-of-fit statistic D.
%    2) We then select as our estimate of x_min, the value that gives the
%       minimum value D over all values of x_min.
%
%    Note that this procedure gives no estimate of the uncertainty of the 
%    fitted parameters.
%
%    Example:
%       x = (1-rand(10000,1)).^(-1/(2.5-1));
%       [p, gof] = plpva(x, 1);
%
%    For more information, try 'type plpva'
%
%    See also PLFIT, PLVAR

% Version 1.0   (2007 May)
% Version 1.0.2 (2007 September)
% Version 1.0.3 (2007 September)
% Version 1.0.4 (2008 January)
% Version 1.0.5 (2008 March)
% Version 1.0.6 (2008 April)
% Version 1.0.7 (2009 October)
% Version 1.0.8 (2012 January)
% Copyright (C) 2008-2012 Aaron Clauset (Santa Fe Institute)
% Distributed under GPL 2.0
% http://www.gnu.org/copyleft/gpl.html
% PLPVA comes with ABSOLUTELY NO WARRANTY
% 
% Notes:
% 
% 1. In order to implement the integer-based methods in Matlab, the numeric
%    maximization of the log-likelihood function was used. This requires
%    that we specify the range of scaling parameters considered. We set
%    this range to be [1.50 : 0.01 : 3.50] by default. This vector can be
%    set by the user like so,
%    
%       p = plpva(x, 1,'range',[1.001:0.001:5.001]);
%    
% 2. PLPVA can be told to limit the range of values considered as estimates
%    for xmin in two ways. First, it can be instructed to sample these
%    possible values like so,
%    
%       a = plpva(x,1,'sample',100);
%    
%    which uses 100 uniformly distributed values on the sorted list of
%    unique values in the data set. Second, it can simply omit all
%    candidates above a hard limit, like so
%    
%       a = plpva(x,1,'limit',3.4);
%    
%    Finally, it can be forced to use a fixed value, like so
%    
%       a = plpva(x,1,'xmin',1);
%    
%    In the case of discrete data, it rounds the limit to the nearest
%    integer.
% 
% 3. The default number of semiparametric repetitions of the fitting
% procedure is 1000. This number can be changed like so
%    
%       p = plpva(x, 1,'reps',10000);
% 
% 4. To silence the textual output to the screen, do this
%    
%       p = plpva(x, 1,'reps',10000,'silent');
% 

vec    = [];
sample = [];
limit  = [];
xminx  = [];
Bt     = [];
quiet  = false;
persistent rand_state;

% parse command-line parameters; trap for bad input
i=1; 
while i<=length(varargin) 
  argok = 1; 
  if ischar(varargin{i}) 
    switch varargin{i}
        case 'range',        vec    = varargin{i+1}; i = i + 1;
        case 'sample',       sample = varargin{i+1}; i = i + 1;
        case 'limit',        limit  = varargin{i+1}; i = i + 1;
        case 'xmin',         xminx  = varargin{i+1}; i = i + 1;
        case 'reps',         Bt     = varargin{i+1}; i = i + 1;
        case 'silent',       quiet  = true;
        otherwise, argok=0; 
    end
  end
  if ~argok
    disp(['(PLPVA) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end
if ~isempty(vec) && (~isvector(vec) || min(vec)<=1)
	fprintf('(PLPVA) Error: ''range'' argument must contain a vector; using default.\n');
    vec = [];
end
if ~isempty(sample) && (~isscalar(sample) || sample<2)
	fprintf('(PLPVA) Error: ''sample'' argument must be a positive integer > 1; using default.\n');
    sample = [];
end
if ~isempty(limit) && (~isscalar(limit) || limit<1)
	fprintf('(PLPVA) Error: ''limit'' argument must be a positive value >= 1; using default.\n');
    limit = [];
end
if ~isempty(Bt) && (~isscalar(Bt) || Bt<2)
	fprintf('(PLPVA) Error: ''reps'' argument must be a positive value > 1; using default.\n');
    Bt = [];
end
if ~isempty(xminx) && (~isscalar(xminx) || xminx>=max(x))
	fprintf('(PLPVA) Error: ''xmin'' argument must be a positive value < max(x); using default behavior.\n');
    xminx = [];
end

% reshape input vector
x = reshape(x,numel(x),1);

% select method (discrete or continuous) for fitting
if     isempty(setdiff(x,floor(x))), f_dattype = 'INTS';
elseif isreal(x),    f_dattype = 'REAL';
else                 f_dattype = 'UNKN'; %#ok<SEPEX>
end
if strcmp(f_dattype,'INTS') && min(x) > 1000 && length(x)>100
    f_dattype = 'REAL';
end
N = length(x);
x = reshape(x,N,1); % guarantee x is a column vector
if isempty(rand_state)
    rand_state = cputime;
    rand('twister',sum(100*clock)); %#ok<CLOCK,RAND>
end
if isempty(Bt)
    Bt = 1000;
end
nof = zeros(Bt,1);

if ~quiet
    fprintf('Power-law Distribution, p-value calculation\n');
    fprintf('   Copyright 2007-2010 Aaron Clauset\n');
    fprintf('   Warning: This can be a slow calculation; please be patient.\n');
    fprintf('   n    = %i\n   xmin = %6.4f\n   reps = %i\n',length(x),xmin,length(nof));
end
tic;
% estimate xmin and alpha, accordingly
switch f_dattype
    
    case 'REAL'
        
        % compute D for the empirical distribution
        z     = x(x>=xmin);	nz   = length(z);
        y     = x(x<xmin); 	ny   = length(y);
        alpha = 1 + nz ./ sum( log(z./xmin) );
        cz    = (0:nz-1)'./nz;
        cf    = 1-(xmin./sort(z)).^(alpha-1);
        gof   = max( abs(cz - cf) );
        pz    = nz/N;

        % compute distribution of gofs from semi-parametric bootstrap
        % of entire data set with fit
        for B=1:length(nof)
            % semi-parametric bootstrap of data
            n1 = sum(rand(N,1)>pz);
            q1 = y(ceil(ny.*rand(n1,1)));
            n2 = N-n1;
            q2 = xmin*(1-rand(n2,1)).^(-1/(alpha-1));
            q  = sort([q1; q2]);

            % estimate xmin and alpha via GoF-method
            qmins = unique(q);
            qmins = qmins(1:end-1);
            if ~isempty(xminx)
                qmins = qmins(find(qmins>=xminx,1,'first'));
            end
            if ~isempty(limit)
                qmins(qmins>limit) = [];
                if isempty(qmins), qmins = min(q); end
            end
            if ~isempty(sample)
                qmins = qmins(unique(round(linspace(1,length(qmins),sample))));
            end
            dat   = zeros(size(qmins));
            for qm=1:length(qmins)
                  qmin = qmins(qm);
                  zq   = q(q>=qmin);
                  nq   = length(zq);
                  a    = nq ./ sum( log(zq./qmin) );
                  cq   = (0:nq-1)'./nq;
                  cf   = 1-(qmin./zq).^a;
                  dat(qm) = max( abs(cq - cf) );
            end
            if ~quiet
                fprintf('[%i]\tp = %6.4f\t[%4.2fm]\n',B,sum(nof(1:B)>=gof)./B,toc/60);
            end
            % store distribution of estimated gof values
            nof(B) = min(dat);
        end
        p = sum(nof>=gof)./length(nof);

    case 'INTS'

        if isempty(vec)
            vec  = (1.50:0.01:3.50);    % covers range of most practical 
        end                            % scaling parameters
        zvec = zeta(vec);

        % compute D for the empirical distribution
        z     = x(x>=xmin);	nz   = length(z);	xmax = max(z);
        y     = x(x<xmin); 	ny   = length(y);

        L  = -Inf*ones(size(vec));
        for k=1:length(vec)
            L(k) = -vec(k)*sum(log(z)) - nz*log(zvec(k) - sum((1:xmin-1).^-vec(k)));
        end
        [Y,I] = max(L); %#ok<ASGLU>
        alpha = vec(I);

        fit = cumsum((((xmin:xmax).^-alpha))./ (zvec(I) - sum((1:xmin-1).^-alpha)));
        cdi = cumsum(hist(z,(xmin:xmax))./nz); %#ok<HIST>
        gof = max(abs( fit - cdi ));
        pz  = nz/N;

        mmax = 20*xmax;
        pdf = [zeros(xmin-1,1); (((xmin:mmax).^-alpha))'./ (zvec(I) - sum((1:xmin-1).^-alpha))];
        cdf = [(1:mmax+1)' [cumsum(pdf); 1]];

        % compute distribution of gofs from semi-parametric bootstrap
        % of entire data set with fit
        for B=1:length(nof)
            % semi-parametric bootstrap of data
            n1 = sum(rand(N,1)>pz);
            q1 = y(ceil(ny.*rand(n1,1)));
            n2 = N-n1;

            % simple discrete zeta generator
            r2 = sort(rand(n2,1));  c = 1;
            q2 = zeros(n2,1);	    k = 1;
            for i=xmin:mmax+1
                while c<=length(r2) && r2(c)<=cdf(i,2), c=c+1; end
                q2(k:c-1) = i;
                k = c;
                if k>n2, break; end
            end
            q = [q1; q2];

            % estimate xmin and alpha via GoF-method
            qmins = unique(q);
            qmins = qmins(1:end-1);
            if ~isempty(xminx)
                qmins = qmins(find(qmins>=xminx,1,'first'));
            end
            if ~isempty(limit)
                qmins(qmins>limit) = [];
                if isempty(qmins), qmins = min(q); end
            end
            if ~isempty(sample)
                qmins = qmins(unique(round(linspace(1,length(qmins),sample))));
            end
            dat   = zeros(size(qmins));
            qmax  = max(q); zq = q;
            for qm=1:length(qmins)
                qmin = qmins(qm);
                zq   = zq(zq>=qmin);
                nq   = length(zq);
                if nq>1
                    try
                        % vectorized version of numerical calculation
                        zdiff = sum( repmat((1:qmin-1)',1,length(vec)).^-repmat(vec,qmin-1,1) ,1);
                        L = -vec.*sum(log(zq)) - nq.*log(zvec - zdiff);
                    catch
                       % iterative version (more memory efficient, but slower)
                       L       = -Inf*ones(size(vec));
                       slogzq  = sum(log(zq));
                       qminvec = (1:qmin-1);
                       for k=1:length(vec)
                           L(k) = -vec(k)*slogzq - nq*log(zvec(k) - sum(qminvec.^-vec(k)));
                       end
                    end
                    [Y,I] = max(L); %#ok<ASGLU>

                    fit = cumsum((((qmin:qmax).^-vec(I)))./ (zvec(I) - sum((1:qmin-1).^-vec(I))));
                    cdi = cumsum(hist(zq,(qmin:qmax))./nq); %#ok<HIST>
                    dat(qm) = max(abs( fit - cdi ));
                else
                    dat(qm) = -Inf;
                end

            end
            if ~quiet
                fprintf('[%i]\tp = %6.4f\t[%4.2fm]\n',B,sum(nof(1:B)>=gof)./B,toc/60);
            end
            % -- store distribution of estimated gof values
            nof(B) = min(dat);
        end
        p = sum(nof>=gof)./length(nof);

    otherwise
        fprintf('(PLPVA) Error: x must contain only reals or only integers.\n');
        p   = [];
        gof = [];
        return;
end
end
function [mu, ul, ll] = circ_mean(alpha, w)
    %
    % mu = circ_mean(alpha, w)
    %   Computes the mean direction for circular data.
    %
    %   Input:
    %     alpha	sample of angles in radians
    %     [w		weightings in case of binned angle data]
    %
    %   Output:
    %     mu		mean direction
    %     ul    upper 95% confidence limit
    %     ll    lower 95% confidence limit 
    %
    % PHB 7/6/2008
    %
    % References:
    %   Statistical analysis of circular data, N. I. Fisher
    %   Topics in circular statistics, S. R. Jammalamadaka et al. 
    %   Biostatistical Analysis, J. H. Zar
    %
    % Circular Statistics Toolbox for Matlab
    
    % By Philipp Berens, 2009
    % berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
    

    % check vector size
    if size(alpha,2) > size(alpha,1)
	    alpha = alpha';
    end
    
    if nargin<2
      % if no specific weighting has been specified
      % assume no binning has taken place
      w = ones(size(alpha));
    else
      if size(w,2) > size(w,1)
        w = w';
      end
    end
    
    % % % w(isnan(w)) = mean(w, 1, 'omitmissing');
    
    if isempty(alpha) || isempty(w)
        alpha = NaN;
        w = NaN;
    end

    % compute weighted sum of cos and sin of angles
    r = w'*exp(1i*alpha);
    
    % obtain mean by
    mu = angle(r);
    
    % confidence limits if desired
    if nargout > 1
      t = circ_confmean(alpha,0.05,w);
      ul = mu + t;
      ll = mu - t;
    end
end
function t = circ_t(alpha, w, d)
    % r = circ_t(alpha, w, d)
    %   Computes mean resultant vector length for circular data.
    %
    %   Input:
    %     alpha	sample of angles in radians
    %     [w		number of incidences in case of binned angle data]
    %     [d    spacing of bin centers for binned data, if supplied 
    %           correction factor is used to correct for bias in 
    %           estimation of r, in radians (!)]
    %
    %   Output:
    %     r		mean resultant length
    %
    % PHB 7/6/2008
    %
    % References:
    %   Statistical analysis of circular data, N.I. Fisher
    %   Topics in circular statistics, S.R. Jammalamadaka et al. 
    %   Biostatistical Analysis, J. H. Zar
    %
    % Circular Statistics Toolbox for Matlab
    
    % By Philipp Berens, 2009
    % berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
    
    % Jan Kudlacek modification 2025-06-09
    % Addition of the normalized version according to
    % Andrzejak et al. 2023: High expectations on phase locking: Better quantifying the concentration of circular data. Chaos 33,  091106
    % https://doi.org/10.1063/5.0166468
    % Full text available here
    % https://www.upf.edu/web/ntsa/publications-featured
    % Requires global variable stg to be defined with a field "andrzejak". If true, returns the T instead of R, see equation (5).
    
    
    % check vector size
    if size(alpha,2) > size(alpha,1)
	    alpha = alpha';
    end
    
    if nargin<2
      % if no specific weighting has been specified
      % assume no binning has taken place
      w = ones(size(alpha));
      n = numel(alpha); % Added by Jan. Used for Andrzejak's re-normalization later
    else
      if size(w,2) > size(w,1)
        w = w';
      end
      n = sum(w); % Added by Jan. Used for Andrzejak's re-normalization later
    end
    if nargin<3
      % per default do not apply correct for binned data
      d = 0;
    end
    % % % isValid = ~isnan(alpha) & ~isnan(w);
    % % % alpha = alpha(isValid);
    % % % w = w(isValid);

    % % % w(isnan(w)) = mean(w, 1, 'omitmissing');

    if isempty(alpha) || isempty(w)
        alpha = NaN;
        w = NaN;
    end
    
    % compute weighted sum of cos and sin of angles
    r = w'*exp(1i*alpha);
    
    % obtain length 
    r = abs(r)/sum(w);
    
    % for data with known spacing, apply correction factor to correct for bias
    % in the estimation of r (see Zar, p. 601, equ. 26.16)
    if d ~= 0
      c = d/2/sin(d/2);
      r = c*r;
    end

    % Andrzejak's re-normalization
    gamma = 1/2*sqrt(pi/n); % Equation (4)
    t = (r - gamma)/(1 - gamma);
end
function r = circ_r(alpha, w, d)
    % r = circ_r(alpha, w, d)
    %   Computes mean resultant vector length for circular data.
    %
    %   Input:
    %     alpha	sample of angles in radians
    %     [w		number of incidences in case of binned angle data]
    %     [d    spacing of bin centers for binned data, if supplied 
    %           correction factor is used to correct for bias in 
    %           estimation of r, in radians (!)]
    %
    %   Output:
    %     r		mean resultant length
    %
    % PHB 7/6/2008
    %
    % References:
    %   Statistical analysis of circular data, N.I. Fisher
    %   Topics in circular statistics, S.R. Jammalamadaka et al. 
    %   Biostatistical Analysis, J. H. Zar
    %
    % Circular Statistics Toolbox for Matlab
    
    % By Philipp Berens, 2009
    % berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
    
    % check vector size
    if size(alpha,2) > size(alpha,1)
	    alpha = alpha';
    end
    
    if nargin<2
      % if no specific weighting has been specified
      % assume no binning has taken place
      w = ones(size(alpha));
    else
      if size(w,2) > size(w,1)
        w = w';
      end
    end
    if nargin<3
      % per default do not apply correct for binned data
      d = 0;
    end
    % % % isValid = ~isnan(alpha) & ~isnan(w);
    % % % alpha = alpha(isValid);
    % % % w = w(isValid);

    % % % w(isnan(w)) = mean(w, 1, 'omitmissing');

    if isempty(alpha) || isempty(w)
        alpha = NaN;
        w = NaN;
    end
    
    % compute weighted sum of cos and sin of angles
    r = w'*exp(1i*alpha);
    
    % obtain length 
    r = abs(r)/sum(w);
    
    % for data with known spacing, apply correction factor to correct for bias
    % in the estimation of r (see Zar, p. 601, equ. 26.16)
    if d ~= 0
      c = d/2/sin(d/2);
      r = c*r;
    end
end
function t = circ_confmean(alpha, xi, w, d)
    %
    % t = circ_mean(alpha, xi, w, d)
    %   Computes the confidence limits on the mean for circular data.
    %
    %   Input:
    %     alpha	sample of angles in radians
    %     [xi   (1-xi)-confidence limits are computed, default 0.05]
    %     [w		number of incidences in case of binned angle data]
    %     [d    spacing of bin centers for binned data, if supplied 
    %           correction factor is used to correct for bias in 
    %           estimation of r, in radians (!)]
    
    %
    %   Output:
    %     t     mean +- d yields upper/lower (1-xi)% confidence limit
    %
    % PHB 7/6/2008
    %
    % References:
    %   Statistical analysis of circular data, N. I. Fisher
    %   Topics in circular statistics, S. R. Jammalamadaka et al. 
    %   Biostatistical Analysis, J. H. Zar
    %
    % Circular Statistics Toolbox for Matlab
    
    % By Philipp Berens, 2009
    % berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
    
    
    % check vector size
    if size(alpha,2) > size(alpha,1)
	    alpha = alpha';
    end
    
    % set confidence limit size to default
    if nargin<2 || isempty(xi)
      xi = 0.05;
    end
    
    if nargin<2
      % if no specific weighting has been specified
      % assume no binning has taken place
      w = ones(size(alpha));
    else
      if size(w,2) > size(w,1)
        w = w';
      end
      if length(alpha)~=length(w)
        error('Input dimensions do not match.')
      end
    end
    % % % isValid = ~isnan(alpha) & ~isnan(w);
    % % % alpha = alpha(isValid);
    % % % w = w(isValid);
    
    % % % w(isnan(w)) = mean(w, 1, 'omitmissing');
    
    if isempty(alpha) || isempty(w)
        t = NaN;
        return
    end

    if nargin<4
      % per default do not apply correct for binned data
      d = 0;
    end
    
    % compute ingredients for conf. lim.
    r = circ_r(alpha,w,d);
    n = sum(w);
    R = n*r;
    c2 = chi2inv((1-xi),1);
    
    % check for resultant vector length and select appropriate formula
    if r < .9 && r > sqrt(c2/2/n)
      t = sqrt((2*n*(2*R^2-n*c2))/(4*n-c2));  % equ. 26.24
    elseif r >= .9
      t = sqrt(n^2-(n^2-R^2)*exp(c2/n));      % equ. 26.25
    else 
      t = NaN;
      % warning('Resultant vector does not allow to specify confidence limits on mean. \nResults may be wrong or inaccurate.'); 
    end
    
    % apply final transform
    t = acos(t/R);
end
function [pval, z] = circ_rtest(alpha, w, d)
    %
    % [pval, z] = circ_rtest(alpha,w)
    %   Computes Rayleigh test for non-uniformity of circular data.
    %   H0: the population is uniformly distributed around the circle
    %   HA: the populatoin is not distributed uniformly around the circle
    %   Assumption: the distribution has maximally one mode and the data is 
    %   sampled from a von Mises distribution!
    %
    %   Input:
    %     alpha	sample of angles in radians
    %     [w		number of incidences in case of binned angle data]
    %     [d    spacing of bin centers for binned data, if supplied 
    %           correction factor is used to correct for bias in 
    %           estimation of r, in radians (!)]
    %
    %   Output:
    %     pval  p-value of Rayleigh's test
    %     z     value of the z-statistic
    %
    % PHB 7/6/2008
    %
    % References:
    %   Statistical analysis of circular data, N. I. Fisher
    %   Topics in circular statistics, S. R. Jammalamadaka et al. 
    %   Biostatistical Analysis, J. H. Zar
    %
    % Circular Statistics Toolbox for Matlab
    
    % By Philipp Berens, 2009
    % berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
    
    if size(alpha,2) > size(alpha,1)
	    alpha = alpha';
    end

    % % % isValid = ~isnan(alpha);
    % % % if nargin > 1
    % % %   if size(w,2) > size(w,1)
	% % %     w = w';
    % % %   end
    % % %   isValid = isVAlid & ~isnan(w);
    % % %   w = w(isValid);
    % % % end
    % % % alpha = alpha(isValid);
    
    % % % w(isnan(w)) = mean(w, 1, 'omitmissing');
    
    if nargin < 2
	  r =  circ_r(alpha);
      n = length(alpha);
    else
      if length(alpha)~=length(w)
        error('Input dimensions do not match.')
      end
      if nargin < 3
        d = 0;
      end
      r =  circ_r(alpha,w(:),d);
      n = sum(w);
    end
    
    % compute Rayleigh's R (equ. 27.1)
    R = n*r;
    
    % compute Rayleigh's z (equ. 27.2)
    z = R^2 / n;
    
    % compute p value using approxation in Zar, p. 617
    pval = exp(sqrt(1+4*n+4*(n^2-R^2))-(1+2*n));
    
    % outdated version:
    % compute the p value using an approximation from Fisher, p. 70
    % pval = exp(-z);
    % if n < 50
    %   pval = pval * (1 + (2*z - z^2) / (4*n) - ...
    %    (24*z - 132*z^2 + 76*z^3 - 9*z^4) / (288*n^2));
    % end
end
function [pval, m] = circ_otest(alpha, sz, w)
    %
    % [pval, m] = circ_otest(alpha,sz,w)
    %   Computes Omnibus or Hodges-Ajne test for non-uniformity of circular data.
    %   H0: the population is uniformly distributed around the circle
    %   HA: the population is not distributed uniformly around the circle
    %
    %   Alternative to the Rayleigh and Rao's test. Works well for unimodal,
    %   bimodal or multimodal data. If requirements of the Rayleigh test are 
    %   met, the latter is more powerful.
    %
    %   Input:
    %     alpha	sample of angles in radians
    %     [sz   step size for evaluating distribution, default 1 degree
    %     [w		number of incidences in case of binned angle data]
    
    %   Output:
    %     pval  p-value 
    %     m     minimum number of samples falling in one half of the circle
    %
    % PHB 3/16/2009
    %
    % References:
    %   Biostatistical Analysis, J. H. Zar
    %   A bivariate sign test, J. L. Hodges et al., 1955
    %   A simple test for uniformity of a circular distribution, B. Ajne, 1968
    %
    % Circular Statistics Toolbox for Matlab
    
    % By Philipp Berens, 2009
    % berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
    
    if size(alpha,2) > size(alpha,1)
	  alpha = alpha';
    end

    % % % isValid = ~isnan(alpha);
    % % % if nargin > 1
    % % %   if size(w,2) > size(w,1)
	% % %     w = w';
    % % %   end
    % % %   isValid = isVAlid & ~isnan(w);
    % % %   w = w(isValid);
    % % % end
    % % % alpha = alpha(isValid);
    
    % % % w(isnan(w)) = mean(w, 1, 'omitmissing');
    
    if nargin < 2 || isempty(sz)
      sz = circ_ang2rad(1);
    end
    
    if nargin < 3
      w = ones(size(alpha));
    else
      if length(alpha)~=length(w)
        error('Input length does not match.')
      end
      w =w(:);  
    end
    
    alpha = mod(alpha,2*pi);
    n = sum(w);
    dg = 0:sz:pi;
    
    m = zeros(size(dg));
    for i=1:length(dg)
      m(i) = sum((alpha > dg(i) & alpha < pi + dg(i)).*w);    
    end
    m = min(m);
    
    if n > 50
      % approximation by Ajne (1968)
      A = pi*sqrt(n) / 2 / (n-2*m);
      pval = sqrt(2*pi) / A * exp(-pi^2/8/A^2);
    else
      % exact formula by Hodges (1955)
      pval = 2^(1-n) * (n-2*m) * nchoosek(n,m);  
    end
end
function alpha = circ_ang2rad(alpha)
        % alpha = circ_ang2rad(alpha)
    %   converts values in degree to radians
    %
    % Circular Statistics Toolbox for Matlab
    
    % By Philipp Berens, 2009
    % berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
    
    alpha = alpha * pi /180;
end
function ht = textbp(string,varargin)
    % TEXTBP  implements 'best' location for text, a la legend
    %    TEXTBP uses a modified LSCAN algorithm from the old MATLAB
    %    LEGEND command to place text such that it minimizes the
    %    obscuration of data points.
    %
    %    TEXTBP(STRING) is the simplest use of this function.  Any text
    %    properties can be passed in by the same methods implemented in
    %    the MATLAB TEXT builtin function. ie, following the STRING
    %    with (PropertyName,PropertyValue) pairs. 
    %
    %    HT = TEXTBP(STRING) returns the handle to the text object
    TOL = 50; % Max # of data points we are allowed to obscure
    % first get the size of the text in plot-normalized units
    h_temp = text(0,0,string,'units','normalized',varargin{:});
    extent = get(h_temp,'Extent');
    width = extent(3);
    height = extent(4);
    delete(h_temp);
    % do the hard work
    pos = tscan(gca,width,height,TOL);
    % if everything went fine, then put the text onto the plot
    if (pos ~= -1)
      ht_local = text(pos(1),pos(2),string,'units','normalized',...
		      'Vert','bottom',varargin{:});
    end
    % export the text object handle, if requested.
    if nargout > 0
      ht = ht_local;
    end
end
function [Pos]=tscan(ha,wdt,hgt,tol,stickytol,hl) %#ok<INUSD>
    %TSCAN  Scan for good text location.
    %   TSCAN is used by TEXTBP to determine a "good" place for
    %   the text to appear. TSCAN returns either the
    %   position (in figure normalized units) the text should
    %   appear, or a -1 if no "good" place is found.
    %
    %   TSCAN searches for the best place on the graph according
    %   to the following rules.
    %       1. Text must obscure as few data points as possible.
    %          Number of data points the text may cover before plot
    %          is "squeezed" can be set with TOL. The default is a 
    %          large number so to enable squeezing, TOL must 
    %          be set. A negative TOL will force squeezing.
    %       2. Regions with neighboring empty space are better.
    %       3. Bottom and Left are better than Top and Right.
    %      x 4. If a legend already exists and has been manually placed,
    %      x    then try to put new legend "close" to old one. 
    %
    %   TSCAN(HA,WDT,HGT,TOL,STICKYTOL,HL) returns a 2 element
    %   position vector. WDT and HGT are the Width and Height of
    %   the legend object in figure normalized units. TOL
    %   and STICKYTOL are tolerances for covering up data points.
    %   HL is the handle of the current legend or -1 if none exist. 
    %
    %   TSCAN(HA,WDT,HGT,TOL) allows up to TOL data
    %   points to be covered when selecting the best
    %   text location.
    %
    %changes from LSCAN
    %1. existing text bracketing fixed
    %2. sticky references removed for clarity
    %3. returns position in plot normalized units, not figure
    %normalized (0,0 is LL axis, not LL of figure)
    % data point extraction for-loop modified to handle histograms
    % properly (Peter Mao, 6/21/11).
    % modified from LSCAN by Peter Mao 6/16/06.
    %   Drea Thomas     5/7/93
    %   Copyright 1984-2005 The MathWorks, Inc.
    %   $Revision: 1.1 $  $Date: 2006/06/19 21:08:48 $
    %   $Revision: 1.2 $  $Date: 2011/06/21 09:41 $
    % Defaults
    debug=0;
    if debug>1
      holdstatOFF = ~ishold;%%
      hold on;%%
    end
    % Calculate tile size
    % save old units
    %% this part makes text walk across screen with repeated use!!
    %axoldunits = get(ha,'units');
    %set(ha,'units','normalized')
    %cap=get(ha,'Position'); %[fig]
    %set(ha,'units',axoldunits);
    cap = [0 0 1 1]; %'position' in [norm] units
    xlim=get(ha,'Xlim'); %[data]
    ylim=get(ha,'Ylim'); %[data]
    H=ylim(2)-ylim(1);
    W=xlim(2)-xlim(1);
    dh=.03*H;
    dw=.03*W;   % Scale so legend is away from edge of plot
    H=.94*H;
    W=.94*W;
    xlim=xlim+[dw -dw];
    ylim=ylim+[dh -dh];
    Hgt=hgt/cap(4)*H; %[data]
    Wdt=wdt/cap(3)*W;
    Thgt=H/round(-.5+H/Hgt);
    Twdt=W/round(-.5+W/Wdt);
    % Get data, points and text
    % legend, not included here, is a child of gcf with 'tag' 'legend'
    Kids=get(ha,'children');
    Xdata=[];Ydata=[];
    for i=1:size(Kids,1)
      Xtemp = [];
      Ytemp = [];
      if strcmp(get(Kids(i),'type'),'line')
        Xtemp = get(Kids(i),'Xdata');
        Ytemp = get(Kids(i),'Ydata');
      elseif strcmp(get(Kids(i),'type'),'patch') % for histograms
        % X/Ydata from patch are LL,UL,UR,LR for each bar.  
        % loop below fills in the edges of the bar with fake data
        Xtemp0 = get(Kids(i),'Xdata');
        Ytemp0 = get(Kids(i),'Ydata');
        Xstart = Xtemp0(1,:);
        Yend   = Ytemp0(2,:);
        for jj=1:length(Xstart)
          thisY = Yend(jj):-Thgt:0;
          thisX = repmat(Xstart(jj),1,length(thisY));
          Xtemp = [Xtemp, thisX]; %#ok<AGROW>
          Ytemp = [Ytemp, thisY]; %#ok<AGROW>
        end
      elseif strcmp(get(Kids(i),'type'),'text')
        tmpunits = get(Kids(i),'units');
        set(Kids(i),'units','data')
        %        tmp=get(Kids(i),'Position');
        ext=get(Kids(i),'Extent');
        set(Kids(i),'units',tmpunits);
        %        Xdata=[Xdata,[tmp(1) tmp(1)+ext(3)]];
        %        Ydata=[Ydata,[tmp(2) tmp(2)+ext(4)]];
        Xtemp = [ext(1) ext(1) ext(1)+ext(3) ext(1)+ext(3)*.5 ext(1)+ext(3)];
        Ytemp = [ext(2) ext(2)+ext(4) ext(2) ext(2)+ext(4)*.5 ext(2)+ext(4)];
      end
      Xdata=[Xdata; Xtemp(:)];
      Ydata=[Ydata; Ytemp(:)];
    end
    if debug>1, plot(Xdata,Ydata,'r.'); end
    %   Determine # of data points under each "tile"
    i=1;j=1;
    for yp=ylim(1):Thgt/2:(ylim(2)-Thgt)
        i=1;
        for xp=xlim(1):Twdt/2:(xlim(2)-Twdt)
           pop(j,i) = ...
               sum(sum((Xdata >= xp).*(Xdata<=xp+Twdt).*(Ydata>=yp).*(Ydata<=yp+Thgt)));    
    %       line([xp xp],[ylim(1) ylim(2)]);
           i=i+1;   
        end
    %    line([xlim(1) xlim(2)],[yp yp]);
        j=j+1;
    end
    % Cover up fewest points.
    minpop = min(min(pop));
    if debug, disp(sprintf('minimally covering tile convers %d points',minpop)); end %#ok<UNRCH,DSPSP>
    if minpop > tol
        Pos=-1;
        warning('Raise TOL in calling function to %d',minpop);
        return
    end
    %%%%%%%%%%%%%%%%%%%%%%
    %                    %
    %sticky stuff removed%
    %                    %
    %%%%%%%%%%%%%%%%%%%%%%
    popmin = pop == min(min(pop));
    if sum(sum(popmin))>1     % Multiple minima in # of points
      [a,b]=size(pop);
      if min(a,b)>1 %check over all tiles
        
        % Look at adjacent tiles and see if they are empty
        % adds in h/v nearest neighbors, double add if on an edge
        pop=[pop(2,:)',pop(1:(a-1),:)']'+[pop(2:(a),:)',pop((a-1),:)']'+...
	    [pop(:,2),pop(:,1:(b-1))]+[pop(:,2:b),pop(:,(b-1))] + pop;
        % LSCAN had two calls to the line above w/o the trailing "+ pop"
        popx=popmin.*(pop==min(pop(popmin)));
        if sum(sum(popx))>1 % prefer bottom left to top right
          flag=1;i=1;j=1;
          while flag
	    if flag == 2
	      if popx(i,j) == 1
	        popx=popx*0;popx(i,j)=1;
	        flag = 0;
	      else
	        popx=popx*0;popx(i,j+1)=1;
	        flag = 0;
	      end
	    else
	      if popx(i,j)==1
	        flag = 2;
	        popx=popx*0;popx(i,j)=1;
	      else
	        j=j+1;
	        if j==b+1
	          j=1;i=i+1;
	        end
	        if i==a+1 % my add'n
	          i=1;
	          flag=2;
	        end
	      end
	    end
          end
        end
      else % only one tile
        popx=popmin*0;popx(1,1)=1;
      end
    else   % Only 1 minima in # covered points
      popx=popmin;
    end
   %recover i,j location that we want to use
   i=find(max(popx));i=i(1); 
   j=find(max(popx'));j=j(1); %#ok<UDIM>
    Pos=[((i-1)/(W/Twdt*2/.94)+.03)*cap(3)+cap(1),((j-1)/(H/Thgt*2/.94)+.03)*cap(4)+cap(2)];
    if debug, disp(sprintf('(i,j) = (%d,%d)',i,j)); end %#ok<UNRCH,DSPSP>
    if debug>1
      if holdstatOFF
        hold off
      end
    end
end

% end


% % % % % % % % % % % % % % % % %% OLD SETTINGS
% % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % %% Select plots
% % % % % % % % % % % % % % % % % Seizure occurrence
% % % % % % % % % % % % % % % % stg.plotSzRaster            =  0; % Raster plot of seizures
% % % % % % % % % % % % % % % % stg.plotSzKaroly            =  0; % Plot according to Karoly et al., Brain 2016
% % % % % % % % % % % % % % % % stg.plotSzRate              =  0; % Seizure rate which is used for the PSD computation
% % % % % % % % % % % % % % % % stg.plotSzPsd               =  0; % Dropouts accounted for in szRate but it is impossible to compensate for them in the PSD
% % % % % % % % % % % % % % % % stg.plotSzPsdAllPop         =  0; % Dropouts accounted for in szRate but it is impossible to compensate for them in the PSD
% % % % % % % % % % % % % % % % stg.plotSzIsiHist           =  0; % Dropouts not accounted for
% % % % % % % % % % % % % % % % stg.plotSzIsiHistAll        =  0; % Dropouts not accounted for
% % % % % % % % % % % % % % % % stg.plotSzIsiHistPop        =  0; % Dropouts not accounted for
% % % % % % % % % % % % % % % % % Seizure characteristics
% % % % % % % % % % % % % % % % stg.plotSzChar              =  0; % Just plot the data
% % % % % % % % % % % % % % % % stg.plotSzCharWhFit         =  0; % Fit whole recording
% % % % % % % % % % % % % % % % stg.plotSzCharWhFitAllPop   =  0; % Fit whole recording
% % % % % % % % % % % % % % % % stg.plotSzCharCl            =  0; % Fit during cluster
% % % % % % % % % % % % % % % % stg.plotSzCharClFit         =  0; % Fit during cluster
% % % % % % % % % % % % % % % % stg.plotSzCharClFitAllPop   =  0; % Fit during cluster
% % % % % % % % % % % % % % % % stg.plotSzCharCiFit         =  0; % Circadian profile
% % % % % % % % % % % % % % % % stg.plotSzCharCiFitAllPop   =  0; % Circadian profile
% % % % % % % % % % % % % % % % % Seizure and signal characteristics in one figure
% % % % % % % % % % % % % % % % stg.plotSaChar              =  0; % Just plot the data
% % % % % % % % % % % % % % % % stg.plotSaCharCWT           =  0; % Continuous wavelet transform and wavelet coherence
% % % % % % % % % % % % % % % % stg.plotSaCharWhFit         =  0; % Fit whole recording
% % % % % % % % % % % % % % % % stg.plotSaCharWhFitAllPop   =  0; % Fit whole recording
% % % % % % % % % % % % % % % % stg.plotSaCharCl            =  0; % Data during cluster
% % % % % % % % % % % % % % % % stg.plotSaCharClFit         =  0; % Fit during cluster
% % % % % % % % % % % % % % % % stg.plotSaCharClFitAllPop   =  0; % Fit during cluster
% % % % % % % % % % % % % % % % % % % % % stg.clusterExampleMouseJc20190509_2 = 0;
% % % % % % % % % % % % % % % % stg.plotSaCharCiFit         =  0; % Circadian profile
% % % % % % % % % % % % % % % % stg.plotSaCharCiFitAllPop   =  0; % Circadian profile
% % % % % % % % % % % % % % % % % Signal characteristics
% % % % % % % % % % % % % % % % stg.plotSiChar              =  0; % Just plot the data
% % % % % % % % % % % % % % % % stg.plotSiCharPsd           =  0; % Power spectral density
% % % % % % % % % % % % % % % % stg.plotSiCharPsdAllPop     =  0; % Power spectral density
% % % % % % % % % % % % % % % % stg.plotSiCharWhFit         =  0; % Fit whole recording
% % % % % % % % % % % % % % % % stg.plotSiCharWhFitAllPop   =  0; % Fit whole recording
% % % % % % % % % % % % % % % % stg.plotSiCharCl            =  0; % Raw data before, during and after the cluster
% % % % % % % % % % % % % % % % stg.plotSiCharClAllPop      =  0; % Raw data before, during and after the cluster
% % % % % % % % % % % % % % % % stg.plotSiCharClFit         =  0; % Fit before, during and after the cluster
% % % % % % % % % % % % % % % % stg.plotSiCharClFitAllPop   =  0; % Fit before, during and after the cluster
% % % % % % % % % % % % % % % % stg.plotSiCharCiFit         =  0; % Circadian profile
% % % % % % % % % % % % % % % % stg.plotSiCharCiFitAllPop   =  0; % Circadian profile
% % % % % % % % % % % % % % % % stg.plotSiCharSzBeAfVsOther =  0; % Compare the IED rate around seizure (before or after) vs. at other times (added in rev01)
% % % % % % % % % % % % % % % % stg.plotSiCharSz            =  0; % Raw data before and after the seizure
% % % % % % % % % % % % % % % % stg.plotSiCharSzAllPop      =  0; % Raw data before and after the seizure
% % % % % % % % % % % % % % % % stg.plotSiCharSzFit         =  0; % Line fit before and after the seizure
% % % % % % % % % % % % % % % % stg.plotSiCharSzFitAllPop   =  0; % Line fit before and after the seizure
% % % % % % % % % % % % % % % % stg.plotSiCharSzCur         =  0; % Curve fit after the seizure
% % % % % % % % % % % % % % % % stg.plotSiCharSzCurAllPop   =  0; % Curve fit after the seizure
% % % % % % % % % % % % % % % % % Seizure and signal characteristics and filter-derived IED rate in one figure
% % % % % % % % % % % % % % % % stg.plotSsChar              =  0; % Just plot the data
% % % % % % % % % % % % % % % % stg.plotSfCharWhFit         =  0; % Fit whole recording
% % % % % % % % % % % % % % % % stg.plotSfCharWhFitAllPop   =  0; % Fit whole recording
% % % % % % % % % % % % % % % % stg.plotSfCharClFit         =  0; % Fit during cluster
% % % % % % % % % % % % % % % % stg.plotSfCharClFitAllPop   =  0; % Fit during cluster
% % % % % % % % % % % % % % % % stg.plotSfCharCiFit         =  0; % Circadian profile
% % % % % % % % % % % % % % % % stg.plotSfCharCiFitAllPop   =  0; % Circadian profile
% % % % % % % % % % % % % % % % stg.plotSimSim              =  0;
% % % % % % % % % % % % % % % % stg.plotSsExplainConvolution=  0; % Explanation of convolution with individual responses
% % % % % % % % % % % % % % % % stg.plotSsExplainSumOfExp   =  0;
% % % % % % % % % % % % % % % % % General
% % % % % % % % % % % % % % % % stg.showStat                =  0;
% % % % % % % % % % % % % % % % stg.printFigures            =  0;


% FINISHED Add tauH to printed stats output
% FINISHED Problem in szChar circular. Sharp spike at midnight.
% FINISHED Sort out percentages
% FINISHED Name handles automatically by function names in all plotting functions
% FINISHED Function for stats table initialization
% FINISHED Deal with the warning of the polyfit (compute it only if it makes sense)
% FINISHED Deal with dropouts in szChar plots
% FINISHED "If the pre-cluster start is before previous cluster's end, set it to the previous cluster's end" might be not the best solution since this way the post-cluster phenomena confound the results
% FINISHED Fix the seizure rate in plotSzRate so that the dropouts are taken into account
% FINISHED, added an option to do either. szRate should maybe show the same bins as is used for the PSD and by plot and not patch
% FINISHED Add number of clusters evaluated in each subject to the stats tables
% FINISHED Order of plots in siCharTrend
% FINISHED Sort out sbStats
% FINISHED Re-run IED and EMG detection on jc20190327. Time stamp discrepancy in the label files was found.
% FINISHED Sort out YLim in szChar plots
% FINISHED Automatic figure size according to the number of subjects: If one subject, make it A5 or A7.
% FINISHED Add figure names
% FINISHED Test whether selecting just some chars works
% FINISHED Is PSD in dB or linear? Answer: pwelch returns linear, so I convert it to dB in the plotting function
% FINISHED Subjects should have their color irrespective of which subjects are plotted. Commenting out certain subjects might be not the best way.
% FINISHED Plot signal ok mark in all figures
% FINISHED All: cluster-wise x ANIMAL-WISE
% FINISHED All: normalized x ABSOLUTE
% FINISHED All: color coding: subject x TIME WINDOW
% FINISHED What are the units of seizure rate?
% FINISHED Box: ON x off
% FINISHED szRaster fix margins
% FINISHED forceZeroYLim1TF change to char system (e.g. 'nonneg' or 'zero')
% FINISHED All: in absolute time x-axis limits should be same in all chars
% FINISHED Change y-labels so that they fit better (possibly everything in one row)
% FINISHED szCharClFit fix labels
% FINISHED Subject names black
% FINISHED createAxesAllPop should use plotName to get ylbl
% FINISHED Smaller markers in PSDs
% FINISHED All: orig sig char: overlay x NOT OVERLAY
% FINISHED Should clust be a table instead of structure? NO, IT WOULD BE NOT PRACTICAL
% FINISHED markClusterInRaster
% FINISHED Depict clusters
% FINISHED Put figure creation to the plotting function (remove createFigures in the future)
% FINISHED Solve the histograms with logarithmic vs. linear bins.
% FINISHED Use formatAxesInd to all plots
% FINISHED Auto XLim in plotSzRaster
% FINISHED ISI non-exponential - provide stats
% FINISHED Call also szChar by name and not column number
% FINISHED siCharClFitAll use it or delete it
% FINISHED Use the term subject instead of animal or mouse
% FINISHED YLabel should only be tilted if numChar is > something
% FINISHED Auto formatting of the seizure occurrence plots
% FINISHED Cluster-wise or SUBJECT-WISE? Unify figures and stats
% FINISHED Each figure should have a separate function
% FINISHED Plot the histogram of all ISI. Is it exponential? Is it power-law?
% FINISHED Mark clusters in raster plots, szChar plots and siChar plots
% FINISHED Take dropouts into consideration when finding clusters (cluster does not begin or end by a dropout)
% FINISHED Does seizure dur and rac decrease after inter-cluster period? Compute dur(startOfCluster(N+1)) - dur(endOfCluster(N)). (I would expect it to be negative most of the time)
% FINISHED Intercluster periods durations
% FINISHED siCharPsd should also return something and we should make a population figure from it
% FINISHED Find a better way of determining the method of getting xlim in formatAxes
% FINISHED Rename statsTableInit and statsTable to fitTableInit and fitTable
% FINISHED Get y-label from the getData (and thus from the saved files)
% FINISHED Sz time is probably a wrong label in some cases, should be Sz rate
% FINISHED Test different number of characteristics
% FINISHED Change ret to plotTF (which will equal ~ret). Might be easier to read.
% FINISHED stg.szCharToPlot - call them by name not by number, similarly to stg.siCharToPlot
% FINISHED Design a better subplot function and possibly override the Matlab subplot
% FINISHED Substitute stg.beColor etc. by stg.fitColor which is a matrix
% FINISHED TO DO: showStat
% FINISHED If it plots all subjects in one figure, it is All. If it adds mean of all subjects, it is AllPop. If there is only the mean, it is Pop.
% FINISHED Sort out sbStats
% FINISHED Stats table should be named more descriptively such as fitTable and psdTable. And the functions as well.
% FINISHED Why is stg.siCharBinWeights not used? How do I compute bin weights? PROBABLY NOT NEEDED, IT IS TAKEN INTO ACCOUNT ALREADY IN getData FUNCTION
% FINISHED Use median and mean consistently across functions. E.g. mean within subject, median over subjects.
% FINISHED plotName should be always the first argument
% FINISHED In plotSzCharClFitAllPop, second and third argument are not used
% FINISHED Change {1, kchar} to {kchar}
% FINISHED stg should not appear in the input argument of functions
% FINISHED If trend is significant, type the p-value in bold
% FINISHED Add cluster duration and cluster cycle period to basic stats.
% FINISHED stats.szFreq is wrong since it does not take into account dropouts
% FINISHED Add mean IED rate to basic stats? And mean signal characteristics?
% FINISHED Add Bonferroni correction to Pop plots and stats print
% FINISHED Longer pre- and post-cluster trend should be analyzed. Take rather whole cluster length and not only half-cluster length.
% FINISHED Plot raw pre and post data to see the possibilities of curve fitting (power-law or exponential?)
% FINIShED Normalize circadian data so that maximum is 1. Then, plot mean resultant vector with the length depicting PLV.
% FINISHED Circadian distribution of signal characteristics
% FINISHED Circular add an explanatory circle.
% FINISHED Normalize pre-sz and post-sz trends in AllPop before averaging
% FINISHED Put szChar and siChar in one plot
% FINISHED Add mean circular histogram of szChar to AllPop
% FINISHED Add box plot, violin plot or swarm chart
% FINISHED Show longer peri-seizure data (maybe 3 days)
% FINISHED Fix the bugs in plotting long peri-seizure data
% FINISHED Sort out units of slopes
% FINISHED % Indicate subject sex in the raster plot or in the PSD and comment on catamenial epilepsy


% MAYBE LATER Plot signal characteristics under the trends
% MAYBE LATER Convert char arrays to strings??? Not sure what would be the advantage.
% MAYBE LATER Rename char to feat??? It would be consistent with Isa's paper, otherwise it seems ok to me.
% MAYBE LATER Change the y-labels in the data files - Variance should have a unit
% MAYBE LATER If the fit is made mostly on NaNs it should be also NaN
% MAYBE LATER Make use of Richard labels
% MAYBE LATER Solve the PSD plot labels (now outside the figure)
% MAYBE LATER Signal OK line goes slightly below the x-axis which looks ugly
% MAYBE LATER Circular statistic (subject-wise and sample-wise) should be selected similarly to the other statistics (@mean or @median)

% PROBABLY NOT Do similar to the sz trend with raw chars also for clusters
% PROBABLY NOT Peri-seizure and peri-lead-seizure. Lead seizures are the ones at the cluster onset and outside of clusters

