# StAMPP

# *News*

# StAMPP 1.6.2 _(2020-04-22)_
## Bug fixes
* removed use of temp environment and associated unintentional binding to base environment 

# StAMPP 1.6.1 _(2020-03-18)_
## Bug fixes
* fixed error and example data documentation

# StAMPP 1.6.0 _(2020-03-12)_
## Other changes
* updated to reflect the new R 4.0.0 default of stringsAsFactors = FALSE in data.frame() and read.table()

# StAMPP 1.5.1 _(2017-10-10)_
## Bug fixes
* bug fixes during bootstrapping Fst calculation if Fst values are NaN

# StAMPP 1.5.0 _(2016-10-10)_
## Bug fixes
* bug fixes during importation of genlight objects

# StAMPP 1.4.0 _(2015-06-30)_
## Bug fixes
* fixed format compatibility issues with stamppPhylip and the updated software program Darwin v6
* bug fixes during importation of polyploid genotype data

# StAMPP 1.3.0 _(2014-08-25)_
## Other changes
* updated for compatibility with pegas 0.6

# StAMPP 1.2.0 _(2014-02-26)_
## Enhancements
* speed improvements to stamppConvert, stamppNeisD and stamppGmatrix for large datasets

## Bug fixes
* fixed an error in calculating Nei's genetic distance between individuals using genotype data imported in 'freq' format; -- Prior to this update StAMPP had a bug that when missing data was ecountered it would print an error to the screen and would not calculate the genetic distances. This error was only associated with data imported in 'freq' format when calculating genetic distances between individuals. Any results from data analysis in StAMPP were and are still accurate as this error prevented the calculation of genetic distance and therefore no incorrect values would have been generated.
* corrected an error in calculating Nei's genetic distance where StAMPP did not split the calculation up when using large datasets and consequently produced a memory allocation error when memory was limiting.

# StAMPP 1.1.0 _(2013-06-06)_
## Other changes
* updated for compatibility with adegenet 1.3-9

# StAMPP 1.0.0 _(2013-05-29)_
## First public release.


