# BIRD
Blink-Identified Robust Detrending for EEG / ERP Data in MATLAB

## Requirements
### MATLAB
The method has been developed for use in the MATLAB computing environment. Matlab can be downloaded [here](https://www.mathworks.com/products/matlab.html).
### EEGLAB / ERPLAB
The method has been developed for analyzing EEG data with EEGLAB and ERPLAB. EEGLAB can be downloaded [here](https://sccn.ucsd.edu/eeglab/download.php). ERPLAB can be downloaded [here](https://github.com/lucklab/erplab/releases).

## Installation

Simply download a copy of the BIRD.m and nt_detrend.m files, and add them to your MATLAB path. If you already have a folder added to the MATLAB path, you can place the files there and MATLAB will read them.

## Getting Started

An example script and dataset is provided to outline the steps for running BIRD. The current version of BIRD is designed for ERP analyses. It thus assumes that a file has been loaded into EEGLAB, that an EVENTLIST has been created, and that bins have been assigned with a BDF file. It uses a VEOG channel to detect blinks, and so a VEOG channel must be specified.
1. Do stuff
