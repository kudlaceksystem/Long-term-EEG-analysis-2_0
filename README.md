# ALPACA Analyzer of long-term profiles and circadian arrangement

Use fcdHfoLong06.m. We may rename the function to ALPACA once it can plot at least one figure :-)

I am also splitting the whole program into multiple files, some functions are in +gd already.

The getData function is now finished for the work with label files. If some characteristic needs loading the signal files we have to implement it. It
is prepared but not finished.

I am also changing from datenum to datetime since Matlab recommends it.

## General naming conventions
I use the following naming conventions. I suggest we all follow them in this program or change them (I am open to discussion).
- In agreement with Matlab convention, variable names start with lower case. In case of a multi-word variable name, I use camel case.
- If something has an abbreviation written in all caps (e.g. SEM, PSD, IED), in variable names, I spell it with only first capital (i.e. Sem, Psd, Ied)
- Field names in structures and column names in tables may start with capitals.
- In the following cases, I do not use camel case:
    - Iteration variable starts with k. E.g. kch ... iteration over channels, kf ... iteration over files, k ... just iteration over anything.
    - num... stands for number of ... . Typically numbin ... number of bins, numch ... number of channels.
    - Similarly ...len stands for ... length, such as binlen ... bin length.
    - ...p stands for ... path, e.g. lblp ... label path (e.g. "\\neurodata\Lab Neurophysiology root\Mouse01\"). Caution: As of now, we sometimes use the "\\" at the end and sometimes not. We should standardize this one or the other way. Which one do you prefer?
    - ...n stands for ... name, e.g. lbln ... label name (e.g. "Mouse01-241224_200000-lbl3.mat")
    - ...pn stands for ... path and name (e.g. "\\neurodata\Lab Neurophysiology root\Mouse01\Mouse01-241224_200000-lbl3.mat")
- Signal is usually abbreviated snl except for some old compatibility requirements (sig could cause confusion with the signum function)
- Suffixes may indicate data class
    - N .... datenum
    - Dt ... datetime
    - Du ... duration
    - TF ... true of false, i.e. logical
    - S .... seconds
    - Ind .. index (logical)
    - Sub .. subscript (sequential number)

## +fcn folder
Functions called directly by alpaca. E.g. getData, extractClusters, etc.

## getData
The function getData creates three variables: subjectInfo, ds, dp.

The philosophy is to first get all the data that might be relevant in the ds and dp sturctures using dsDesc and dpDesc.
The data might be saved so that we do not have to load all the individual data files again and again when just improving analyses and figures.
Then, there will be figDesc structure, which will describe the figures. Figures may combine data from both ds and dp.

### Input
#### dsDesc
Data to stem description. Structure.
The first field is "Name" and contains a column vector or strings indicating names of all other fields,
i.e. names of the analyzed phenomena (e.g. seizures, epochs of drinking).
The other fields contain structure arrays.
The length of the structure array is given by the number of parameters of the phenomenon that we may wish to analyze, e.g. onset time, duration, severity.
The structure will have several fields with information required for the computation of given characteristic, see the commented example below.
Different characteristics may require different fields.
If a characteristic does not require certain field which is needed for other characteristics, just do not initialize it.
Matlab will automatically initialize it with some default value, which will be, however, ignored by the rest of the program.
##### Example
```matlab
dsDesc.Name = ["Seizure; Drink"];
% It is often useful to declare some of the variables in advance when they are the same in many elements of the structure array.
mainLbl = ["Seizure", "seizure", "SEIZURE", "S", "s"]; % Names of the main labels to analyze (here any name for a seizure label class)
exLblAn = ["Noise", "Artifact"]; % Labels to exclude in all channels if present in any.
minSepSzS = 60; % Minimum separation of seizures in seconds
dsDesc.Seizure(1).VarName    = "OnsDt"; % Variable name, i.e. name of the characteristic. Here onset in datetime.
dsDesc.Seizure(1).VarType    = "datetime"; % Variable type (i.e. Matlab class). We may consider using the word class instead of type but it could be confused the the label class of OSEL labels.
dsDesc.Seizure(1).CalcFcn    = "gd.dsfGetOnsDt"; % Calculating function. Name of the function which will be called to calculate the data. The function may be in a +folder.
dsDesc.Seizure(1).SrcData    = "Lbl"; % Is the data computed from label files, signal files or both? Permitted values are "Lbl", "Snl" or "LblSnl".
dsDesc.Seizure(1).MainLbl    = mainLbl; % See the declaration above
dsDesc.Seizure(1).ExLblAn    = exLblAn; % See the declaration above
dsDesc.Seizure(1).MinSepS    = minSepSzS; % See the declaration above
dsDesc.Seizure(1).PlotTitle  = "Seizure occurrence"; % Title used in the final plots
dsDesc.Seizure(1).YAxisLabel = ""; % y-axis label used in the final plots (should include units, btw units should be in parentheses and not brackets)
dsDesc.Seizure(2).VarName    = "DurDu";
dsDesc.Seizure(2).VarType    = "duration";
dsDesc.Seizure(2).CalcFcn    = "gd.dsfGetDurDu";
dsDesc.Seizure(2).SrcData    = "Lbl";
dsDesc.Seizure(2).MainLbl    = mainLbl;
dsDesc.Seizure(2).ExLblAn    = exLblAn;
dsDesc.Seizure(2).MinSepS    = minSepSzS;
dsDesc.Seizure(2).PlotTitle  = "Seizure duration";
dsDesc.Seizure(2).YAxisLabel = "Sz dur (s)";
dsDesc.Seizure(3).VarName    = "Pow";
dsDesc.Seizure(3).VarType    = "double";
dsDesc.Seizure(3).CalcFcn    = "gd.dsfGetPow";
dsDesc.Seizure(3).SrcData    = "Lbl";
dsDesc.Seizure(3).MainLbl    = mainLbl;
dsDesc.Seizure(3).ExLblAn    = exLblAn;
dsDesc.Seizure(3).MinSepS    = minSepSzS;
dsDesc.Seizure(3).PlotTitle  = "Seizure signal power";
dsDesc.Seizure(3).YAxisLabel = "Sz power (a.u.)";
... and similarly for the Drink field.
```

#### dpDesc
Data to plot description. Similar to dsDesc.
There may be additional fields, most importantly CalcLvl, which indicates calculation level (on the level of files or bins).
It indicates in which loop the calculation takes place (over bins or over files within the bin).
The reason is that sometimes, we first calculate some characteristics for each file, then aggregate them (e.g. sum them or take the mean) and then,
once we have the whole bin ready, we calculate some other characteristics. A typical example is the calculation of IED rate.
First, we get the IED count per file and the amount of valid signal in each file (i.e. not contaminated by artifacts preventing IED counting).
This represents calculation at the level of files.
Then, we sum the IED count over all files of the block and we also sum the amount of the valid signal.
Then, we calculate at the level of bins the IED rate by dividing IED count by the amount of valid signal.

### Output
#### subjectInfo
Structure containing subject name (code), date of birth, and start and end of the analyzed period. Subject name is a string, others are datetime objects.
#### ds
Data for stem. These data will be shown by the stem type of graph. Structure containing dynamically named fields according to what is analyzed.
In each table, each row belongs to one event, each column to one characteristic of the event.
The first column should always be the onset time of the event in datetime.
If we analyze each channel separately, then each column of the data tables contains a row vector.
##### Example
- ds.Seizure ... table where the first column is "onsDt" (onset, Dt indicates datetime format), then there could be columns like "durDu" (duration in seconds), "pow" (signal power), etc.
- ds.Drink ... columns could be "onsN", "durS", "volM" (volume drank in milliliters)
#### dp
Data for plot. These data will be plotted by line graphs. The data were analyzed in time bins (windows, blocks), e.g. 1-hour bins. The structure has 3 fields, each containing a table.
##### dp.tax
Time axis. Times of ends of the bins in datenum. We will use bins so that the analysis is causal.
E.g. if the result would be seizure risk estimated from the data in given bin,
we have to time-stamp the calculated seizure risk to the end of the time bin when the whole bin is finally recorded
so that we can analyze it and provide the seizure risk estimate for the future (e.g. next bin).
##### Other fields
Other fileds have a similar structure to ds just each row corresponds to one time bin.
If we analyze each channel separately, then each column of the data tables contains a row vector.

### Computation
There are two main parts of the function. The first fills in the ds.
It loads label files sequentially and extracts data.
The other part fills in dp.
This one splits time into bins and loads the files that it needs.
The code is vastly commented so no description is needed here.
It uses a functions from the +gd folder.

#### +gd folder
+gd (stand for getData) folder contains functions called from within getData function. The function names should follow the following convention:
- Functions for computations related to stem data should have the prefix ds
- Functions for computations related to plot data should have the prefix dp
- Functions for computations on file basis should have the prefix f
- Functions for computations on bin basis should have the prefix b
- Functions for computations related to both should have the prefix db (b stands for both)
- Functions directly called by getData function should be named ..Get...
- Functions for manipulation of labels not directly called by getData should be named ..Lbl...
- Functions for manipulation of signal files not directly called by getData should be named ..Snl...

## figMain
### Current state
Defines figure size internally.

Calls a function to create a figure at a specified position.


###
New philosophy:

figDesc.Name ... names of the figures

figDesc.(figDesc.Name(k)) ... structure with fields:
MainFcn ... function to be called
then some description which the MainFcn will use