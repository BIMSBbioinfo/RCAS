# RCAS 1.1.1 

## New features 
  - New function introduced: plotFeatureBoundaryCoverage that creates interative plotly graphs for coverage profiles at transcript feature boundaries
  - Changes in the HTML output of runReport function:
    - The interactive plots now can be exported to not only PNG but also SVG for higher quality figures.
    - The interactive tables for GO and GSEA analyses now have additional column for Benjamini-Hochberg multiple testing corrected p-values. 
    - The design of the coverage profiles have changed:
      - The coverage profiles are represented not only by a line but also a confidence envelope surrounding the line. The solid-line in the middle represents the mean coverage score (output of ScoreMatrix function from genomation package), and the thickness of the surrounding confidence envelope is as large as the standard error of the mean multiplied by 1.96 (95% confidence interval). 
    - The font sizes in all plots have been increased to 14. 
  - calculateCoverageProfile function returns a Score Matrix object rather than a data.frame. 
  - calculateCoverageProfileList object returns a data.frame rather than a list of data.frame objects. 

## Bug fixes
  - In inst/report.Rmd script (called by the runReport function), the motif analysis chunk has been updated to account for the number of input intervals before down-sampling to 10000. 
  - Coverage profile calculations for coding exon-intron junctions has changed:
    - Now, firstly all internal exons are found such that all exons have at least two neighboring introns on both 5' and 3' directions. Then, the coverage profiles on the 5' and 3' directions are calculated. Previously, this filtering for the criteria that the exon **must** be in the middle of two introns was not applied and only a 50 bp region was plotted. Now, the profiled region is the area that spans 1000 bps both upstream and downstream of the internal exons.  

## Changes in Dependencies 
  - Due to a bug fix in older version of genomation package, RCAS now depends on genomation releases later than 1.5.5. 
  - plotrix library is added to dependencies. The std.error function from this library is use to calculate standard error of the mean of a given vector of real numbers. 


# RCAS 0.99.7
- Added selfContained argument to runReport function so that the generated html report can be arranged to be standalone (slow to load but self-contained) or not (fast to load but external dependencies have to be shipped with the html file). 
- Updated motif analysis chunk in report.Rmd script. Query regions shorter than 15 bp are resized to 15bp to enable motif search in datasets with short interval queries. 

# RCAS 0.99.6
fixed repository issue - version bump to 0.99.6

# RCAS 0.99.5
New features:
- DT::datatables in the report now contain buttons to export data tables to CSV, Excel, PDF, or
print or copy the data. 
- with a new option in runReport function 'printProcessedTables', when set to TRUE, the raw data generated to make all tables and figures can be printed to text files in the working directory. 
- A new table is added at the beginning of the runReport output, that shows the input settings
used to run the document
- A sessionInfo() is added to the end of the runReport output
- Added 'quiet' option to runReport function to optionally suppress progress bars and messages during report generation.

Bug fixes
- With the new plotly version release, the download button for plotly generated figures in the 
report html works. 


# RCAS 0.99.4
New features: 
- Additional coverage profile plots for Transcription Start/End Sites
- Better ranking of GO and MSIGDB terms using both FDR and fold change values
- Minor bug fix for motif analysis in the reporting script

# RCAS 0.99.3
Made one minor revision 
- Cleaned up the directory inst/

# RCAS 0.99.2 
Made two minor revisions asked by the reviewer:
- manual loading of some packages manually in the vignette is cancelled
- a NOTE flagged by R CMD Check about a global variable has been fixed 

# RCAS 0.99.1
RCAS has been built on Linux, Mac, and Windows.
The package built without errors or warnings on all platforms.

http://bioconductor.org/spb_reports/RCAS_0.99.1_buildreport_20160419074009.html

# RCAS 0.99.0
Pre-release version of RCAS (RNA Centric Annotation System)



