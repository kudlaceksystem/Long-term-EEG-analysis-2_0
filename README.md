# ALPACA Analyzer of long-term profiles and circadian arrangement

I started implementing it in fcdHfoLong04.m. We may rename the function afterwards. We may also split it into multiple files.

## getData
The function getData creates three variables: subjectInfo, ds, dp.

### Variables
#### subjectInfo
Structure containing subject name (code), date of birth, and start and end of the analyzed period.
#### ds
Data for stem. These data will be shown by the stem type of graph. Structure containing dynamically named fields according to what is analyzed. Each field contains a table. In the table, each row belongs to one event, each column to one characteristic of the event. The first column will always be the onset time of the event in datenum.
##### Example
- ds.Seizure ... table where the first column is "onsN" (onset, N indicates datenum format), then there could be columns like "durS" (duration in seconds), "pow" (signal power), etc.
- ds.Drink ... columns could be "onsN", "durS", "volM" (volume drank in milliliters)
#### dp
Data for plot. These data will be plotted by line graphs. The data were analyzed in time bins (windows, blocks), e.g. 1-hour bins. The structure has 3 fields, each containing a table.
##### dp.tax
Time axis. Times of ends of the bins in datenum. We will use bins so that the analysis is causal. E.g. if the result would be seizure risk estimated from the data in given bin, we have to time-stamp the calculated seizure risk to the end of the time bin when the whole bin is finally recorded so that we can analyze it and provide the seizure risk estimate for the future (e.g. next bin).
##### dp.timeS
Columns contain total time in seconds spent when given label was present during that bin. E.g. during a bin there were three 5-minute epochs of slow-wave sleep, which could be depicted in a column "sws" by number 900 (3\*5=15 minutes is 900 seconds). There could also be e.g. a column "swsPow" showing mean slow-wave sleep EEG signal power during the bin (only the slow-wave sleep epochs would be considered and other signal would be ingored for this calculation). There can be separate columns also for various logical operations between different label classes (e.g. "emgOrArt" for EMG or artifact).
##### dp.count
Counts of various patterns in the given bin (e.g. seizures, interictal discharges, action potentials).
##### dp.rate
Rates of the patterns computed as count/timeOfUsableSignal. timeOfUsableSignal definition may be different for different patterns but it always has to possible to calculated from the dp.timeS and bin length.
##### dp.char
The user may also define some characteristics of the EEG patterns (e.g. spike amplitude or polarity). The computation of the characteristic should be implemented for each pattern (e.g. mean, median or maximum across all realization in given bin).

### Computation
The main part of the function will manage the division of the data into time bins.
Function for specific calculations will be implemented for given application.

### Input
Input to getData function includes: dsName, dsDescr.

#### dsName
A string column vector containing the names of the event types (often OSEL label classes).

#### dsDescr
A structure with fields named by the strings in dsName.
Each field contains a 2D string array.

Columns: Characteristics of the event (must include onsN, i.e. onset in datenum).

Rows:
1. Characteristic name
2. Characteristic data type (double, logical, ...)
3. Name of the function which will calculate the values.
4. Name of the label class which will be used to get the values.




