# TRANSECT

Repository associated with *"Relief dynamics of the Vosges through morphometric and cosmogenic nuclides analyses: a case of transient landscapes in a low deformation context"* (PhD thesis, Mathieux, 2025).  

## Overview  

The **TRANSECT** toolbox provides a framework to build, analyze, and visualize cross-transects perpendicular to a user-defined baseline in a Digital Elevation Model (DEM).  

TRANSECT constructs a baseline (e.g., drainage divide or river reach) and builds connected paths perpendicular to it using either geometric buffers or flow-routing along the topography. Each transect is stored in a structured object, with detailed fields for interpolated paths (`int`) and connected nodes (`conn`). The object supports flexible resampling, statistics computation, and visualization.  

The toolbox leverages **TopoToolbox** for DEM handling, flow routing, and coordinate transformations.  

## Structure  

- **TRANSECT functions**  
  - `TRANSECT.m` — class constructor to build transects along a baseline.  
  - `plot.m` — visualize baseline and transects.  
  - `resample.m` — resample TRANSECT object with another DEM resolution.  
  - `extractT.m` — accessor to extract arrays (x, y, z, d, indices) by transect index.
  - `pairing.m` — pair transect nodes with either coordinates or a logical mask (e.g., floodplain or basin outline).  

- **Helpers**  
  - `SelectDivide.m` — Interactive function to select drainage divide segments from a `DIVIDEobj` (TopoToolbox). 
  - `shortpath.m` — Simplifies and reconnects unordered `(x,y)` coordinates into a valid shortest path. 
  - `ProgressBar.m` — optional utility for tracking iterations.
  - `maskWidth.m` — Computes optimal transect half-width (`w`) from a logical mask (e.g., floodplain or basin outline) for use with `TRANSECT` and `pairing`.

## Installation

1. Clone the repository: `git clone <repository-url>` (replace `<repository-url>` with your repository URL) or download the ZIP file from GitHub and extract it.
2. In MATLAB, navigate to the `TRANSECTs-main` folder and run:
   ```matlab
   addpath(genpath('TRANSECT-main'));
3. Ensure TopoToolbox v3 (https://github.com/TopoToolbox/topotoolbox3) is installed and added to your MATLAB path. MATLAB R2023b or later is recommended (tested on R2025b). Statistics and Machine Learning Toolbox is required. The Mapping Toolbox is recommended, and the Parallel Computing Toolbox is optional for parallel processing.

## Status
This toolbox is still under development and will continue to evolve.

## Contact

For questions, contact bastien.mathieux@gmail.com or open an issue on GitHub. 

Bastien Mathieux - University of Strasbourg
