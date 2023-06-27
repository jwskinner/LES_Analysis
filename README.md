# MATLAB LES Analysis

This Repo contains MATLAB scripts that read Large-Eddy Simulation (LES) data, perform analysis and make plots.

## Usage

1. Clone or download the MATLAB folder.
2. Open MATLAB and navigate to the folder.

## Script Description

The main script in this folder is `main.m`. It performs the following tasks:

1. Reads in data from specified locations and time averages the 1D profiles.
2. Plots out either: plots the fields (e.g., wind speed, temp, and moisture); plots moiture variance and liquid water portential temperature budget terms; computes and plot Liquid Water Path; computes and plots 1D horizontal mean profiles in z; plots the vertical structure of variables in height; or, computes and plots kinetic energy spectra in height.
3. main.m can also loop over all files in the LES simulation directory to produce movies of the aforementioned plots. 

In addition, more specific plotting scripts, which are used to produce plots for papers, are included in 'plotting', these are:  
* plot_CKE_budgets.m -- plots the CKE budgets from the Kardunov & Randall Paper.
* plot_tavg_1D_profs.m -- Reads in LES data and time averages the fields, plots 1D time- and horizontally- averaged profiles of field parameters and fluxes.
* plot_tavg_budgets.m -- Reads in LES data in, computes the liquid water potential temperature and moisture budgets, time averages the budgets.   
* plot_LWP_panels.m -- Reads in LES data at specified times, computes LWP and plots panels at those times (used for producing paper figures).
* plot_compstd_spectra.m -- Reads in LES data, computes and plots the compensated 2D kinetic energy spectrum from Matheou & Teixiera 2019 (figure 14)
* plot_KE_spectra.m -- Reads in LES data, computes and plots a simple 2D kinetic energy spectrum.
* plot_tavg_hist.m -- Reads in LES data, computes and plots time averaged histograms for  U, V, W, Theta and TKE.
* write_to_paraview.m -- Reads in LES data, converts and outputs to paraview format for 3D visualisation.

   
## Prerequisites

To run the script, you need the following:

- MATLAB software installed on your computer.

## Instructions

1. Modify the script according to your specific data locations and parameters.

    - Set the `data_loc` variable to the location of your data files.
    - Modify `folder` variables to specify the folders containing data for different cases.
    - Adjust the `i_start` and `i_end` variables to set the range of files for averaging or specify the time frames for plotting in 'files'.
    - Customize the plot settings as needed.

2. Run the script in MATLAB.

## File Structure

The MATLAB folder has the following structure:

- `LES_Data_Plot.m`: The main script for reading and plotting LES data.
- Other auxiliary functions and scripts.

## Additional Notes

- Make sure the required physical constants are loaded using the `get_constants` function.
- The script uses various auxiliary functions to process the data and generate the plots.
- The resulting plots will be displayed in separate figures.

Feel free to explore the code and customize it to suit your needs!

For more information or assistance, please contact J.W. Skinner (author of the script).

