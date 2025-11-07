# ALPACA Analyzer of long-term profiles and circadian arrangement
# LAMA Long-term analysis of medical annotations
## getData
The function getData creates three variables:
### subjectInfo
Structure containing subject name (code), date of birth, and start and end of the analyzed period.
### ds
Data for stem. These data will be shown by the stem type of graph. Structure containing dynamically named fields according to what is analyzed. Each field contains a table. In the table, each row belongs to one event, each column to one characteristic of the event. The first column will always be the onset time of the event in datenum.
#### Example
- ds.Seiz ... table where the first column is "onsN" (onset, N indicates datenum format), then there could be columns like "durS" (duration in seconds), "pow" (signal power), etc.
- ds.Drink ... columns could be "onsN", "durS", "volM" (volume drank in milliliters)

### dp
Data for plot. These data will be plotted by line graphs. Again, a structure with tables. Each table contains first column "tax"