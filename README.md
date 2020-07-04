# Vedelaar-et-al-2020
Matlab script to determine the fraction of persisters as described in: "A robust method for generating, quantifying and testing large amounts of Escherichia coli persisters"

The file "fit.m" is the actual Matlab script. In the subfolder, there are three versions of this file (called 'goodfit.m', 'badfit1.m', 'badfit2.m'). By running those files, the input data (files 1.csv until 15.csv) will be processed and fitted with different variable bounds. As described in the paper,
this will either generate good or bad fits. 

## Prerequisites

* Matlab R2019A or higher
* CSV files with flowcytometry files as described above
* Script

## Getting Started

The following steps describe the data acquisition and data format requirements. 

* Generate tolerant cells as described in section 1: “Generation of large fractions of tolerant cells” and to stain them as described in the paragraphs. 
* To obtain data that can be efficiently and easily used with the Matlab script, the measurements must be done at specific times, see paper.
* The data needs to be gathered over the growth of the culture.
* The input data needs to be provided as a csv file containing the flow cytometry measurements for each cell at each time point, without a header (File 1). 
* In the Matlab script, the following input needs to be provided: 
	The number identifying the column of the CSV file containing the relevant fluorescence intensity (FI) data (in case of our data, the number is 3 – variable SC). 
	The scaling factor. Different flow cytometers have different sensitivity and numerical output values. (variable maxval_new_FC). 
	The times when the samples were taken (variable tt) 
	The absolute cell concentration in the culture at the respective time points (variable cc) 
	The number of cells (or rows) in each CSV data file (variable g) 
	A specification of which of the data files should be used for the fitting (variables indx and cc_indx). 

```
%% TO SPECIFY 1:
% select column containing FI data and scale the data to fit the histogram
% provide file names to load data files, log10 transformation of data
% and fitting into histogram bins.

%% TO SPECIFY 2:
% enter input data for cell concentration, number of data
% points in each file and relative time at which the sample was meaured

%% TO SPECIFY 3
% specify which samples should be used to perform the fit

%% TO SPECIFY 4
% specify the way data is pre-treated and some initial parameters

%% TO SPECIFY 5
% initial guesses for parameters

```

## The script will provide you with the following:
* Calculated alpha
* Display output, alpha, mu growing population, mu non-growing population, sgma_grow, sgma_nong, I_0, EI, dye_degr, BG cutoff, BG magnit, BG width, CC weight, Bigaus pt, CC pts
* Plot with all bigaussian fits in one figure
* Plot with bigaussian fits, one per subplot
* Plot with all time points in one figure, without bigaussian fit
* Plot with cellcount curve
