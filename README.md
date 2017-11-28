# EEMindices
Automates the calculation of indices derived from excitation-emission matrix (EEM) spectra of dissolved organic matter (DOM).

## Overview

This function was developed and tested for peatland DOM samples run on a Jobin Yvon SPEX Fluoromax-4 spectrometer (Florida State University, Dept. of Chemistry and Biochemistry, Tallahassee, FL), with excitation wavelengths scanned from 240-500nm in 5nm increments, emission wavelengths from 290-600nm in 2nm increments, and other instrument settings used by Hodgkins et al. (2016). It will most likely work with other environmental DOM samples, but has not been tested on non-environmental samples, on data run with other instruments, or with other wavelength ranges.

This function works with the dilution and scatter-corrected ASCII matrix data files (.txt) exported by FLToolbox 2.10 (Zepp et al., 2004). It assumes that these files have not been edited since they were created.

All of your datafiles must have a common string in the filename, e.g. "porewater_[sample name].txt". This string (in the previous example, "porewater_") is used as the `filepattern` argument. Strings interpretable as regular expressions also work, with the argument fixed=FALSE; however, **note that anything matching `filepattern` will be dropped from the sample names in the final output file.** In addition, filenames with special characters (except underscores, dashes, and periods) might produce an error. 

**Required packages:** This function uses the `colorRamps` package to make EEM plots with the familiar color scheme used by FLToolbox.

An example dataset (27 EEM datafiles) is included in the folder "example_data". Results for these spectra are published in the following manuscript:
Hodgkins S. B., Tfaily M. M., Podgorski D. C., McCalley C. K., Saleska S. R., Crill P. M., Rich V. I., Chanton J. P. and Cooper W. T. (2016) Elemental composition and optical properties reveal changes in dissolved organic matter along a permafrost thaw chronosequence in a subarctic peatland. Geochim. Cosmochim. Acta 187, 123–140.

## Output data

The function outputs a .csv file with the following indices:

  * Biological index (BIX): Ratio of emission intensities at (380 nm)/(430 nm) at excitation = 310 nm. This increases with the proportion of microbial exudates in the DOM pool (Huguet et al., 2009).

  * Fluorescence index (FI): Ratio of emission intensities at (470 nm)/(520 nm) at excitation = 370 nm. This increases with the proportion of microbially-derived DOM relative to plant-derived DOM (McKnight et al., 2001) 

  * Humification index (HIX): Defined as H/L, where H and L are the integrated emission intensities from 434–480 nm and 300–344 nm (respectively) measured at excitation = 255 nm (Zsolnay et al., 1999). This index essentially measures the ratio of peaks A/T (as defined by Coble, 1996), and is thus positively correlated with the proportion of aromatic, humic-like DOM.

The output file also includes the following data on peak C (defined by Coble, 1996):

  * C.intensity: Fluorescence intensity of peak C. This intensity normalized to the absorbance at 340 nm (from separate measurements) is inversely related to molecular weight (Stewart and Wetzel, 1981; Baker et al., 2008).

  * C.em.max: Emission wavelength of peak C maximum. This wavelength increases with DOM unsaturation (Senesi, 1990).

  * C.ex.max: Excitation wavelength of peak C maximum.

## Instructions 

1. Put all your EEM data files into the same folder with a copy of EEMindices.R.

2. If not already, make sure your data filenames are in the correct format, with a common string in all the filenames (see above).

3. Open and source EEMindices.R.

4. Call the function using the following arguments:

   `filepattern`: String common to all filenames.

   `outputfile` (optional): Filename for output file, including the extension (.csv) (default='EEMS_results.csv').
        
   `show.plots` (optional): Whether to show plots with the C peak location marked (default=TRUE).
        
   `fixed` (optional): Whether to interpret `filepattern` as a literal string instead of a regular expression (default=TRUE).

5. Unless you specified show.plots=FALSE, check each plot for correct C peak positioning. If it appears incorrect, type any notes to appear in the final output file for that sample.

6. The function outputs a .csv file with the results (BIX, FI, HIX, C.intensity, C.em.max), along with other intermediate variables used to locate peak C (used for troubleshooting) and any notes you typed during step 5.

## Notes on accuracy

As they are based on fixed wavelengths, BIX, FI, and HIX should all be calculated correctly with no problems.

Since the location of peak C varies slightly between samples, this function uses a recursive peak-finding algorithm to locate peak C. Therefore, the peak location (C.em.max and C.ex.max) and intensity (C.intensity) have a slight chance of being inaccurate, especially for samples with a small (peak C)/(peak A) ratio or a small signal/noise ratio. I tested it on a broad range of samples, but it still sometimes locates the peak incorrectly. For example, a common problem is the misidentification of noise on peak A as peak C. This problem and others are easily visible as an incorrect peak C location on the plots, which you can see with show.plots=TRUE.


The peak-finding algorithm uses the following steps:

1. Determines the emission wavelength of maximum fluorescence at an excitation wavelength of 325nm. This emission wavelength is ini.C.em.max.
  
2. Using an emission wavelength of ini.C.em.max, locates the excitation wavelength of the trough between peaks A and C. This excitation wavelength is min.C.ex.
  
3. Finds the maximum fluorescence within the window:
  
   excitation wavelengths:  min.C.ex to 400nm,   and
      
   emission wavelengths:    ini.C.em.max - 26nm to ini.C.em.max + 26nm.
      
4. Makes sure the maximum fluorescence from step 3 does not have a neighboring fluorescence value that is greater. This can (but does not always) happen if its excitation wavelength = min.C.ex.
  
5. If there are no greater neighboring fluorescence values, the maximum fluorescence from step 3 is assigned as C.intensity, and its emission and excitation wavelengths as C.em.max and C.ex.max, respectively.
  
6. If there is a greater neighboring fluorescence value, the location of the fluorescence value from step 3 is eliminated from further consideration. The process then re-starts at step 3, not including previously-eliminated values, with no more than 100 tests performed.


To aid in checking the results, this function does the following:

* Shows plots of the EEMs with the final peak C location indicated with a point (if `show.plots=TRUE`), and prompts for notes on the accuracy of each C peak assignment (see Instructions).

* Outputs intermediate variables used to locate peak C, including ini.C.em.max, min.C.ex, and the number of tests (testcount). If the outputted C.ex.max is still the same as min.C.ex, you should check it manually.

## Licensing

Copyright © 2017 Suzanne Hodgkins and Florida State University.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

## References

Baker A., Tipping E., Thacker S. A. and Gondar D. (2008) Relating dissolved organic matter fluorescence and functional properties. Chemosphere 73, 1765–1772.

Coble P. G. (1996) Characterization of marine and terrestrial DOM in seawater using excitation-emission matrix spectroscopy. Mar. Chem. 51, 325–346.

Hodgkins S. B., Tfaily M. M., Podgorski D. C., McCalley C. K., Saleska S. R., Crill P. M., Rich V. I., Chanton J. P. and Cooper W. T. (2016) Elemental composition and optical properties reveal changes in dissolved organic matter along a permafrost thaw chronosequence in a subarctic peatland. Geochim. Cosmochim. Acta 187, 123–140.

Huguet A., Vacher L., Relexans S., Saubusse S., Froidefond J. M. and Parlanti E. (2009) Properties of fluorescent dissolved organic matter in the Gironde Estuary. Org. Geochem. 40, 706–719.

McKnight D. M., Boyer E. W., Westerhoff P. K., Doran P. T., Kulbe T. and Andersen D. T. (2001) Spectrofluorometric characterization of dissolved organic matter for indication of precursor organic material and aromaticity. Limnol. Oceanogr. 46, 38–48.

Senesi N. (1990) Molecular and quantitative aspects of the chemistry of fulvic acid and its interactions with metal ions and organic chemicals: Part II. The fluorescence spectroscopy approach. Anal. Chim. Acta 232, 77–106.

Stewart A. J. and Wetzel R. G. (1981) Asymmetrical relationships between absorbance, fluorescence, and dissolved organic carbon. Limnol. Oceanogr. 26, 590–597.

Zepp R. G., Sheldon W. M. and Moran M. A. (2004) Dissolved organic fluorophores in southeastern US coastal waters: correction method for eliminating Rayleigh and Raman scattering peaks in excitation-emission matrices. Mar. Chem. 89, 15–36.

Zsolnay A., Baigar E., Jimenez M., Steinweg B. and Saccomandi F. (1999) Differentiating with fluorescence spectroscopy the sources of dissolved organic matter in soils subjected to drying. Chemosphere 38, 45–50.
