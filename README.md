# EenSMITEscripts

MATLAB scripts and functions for microscopy image analysis using SMITE (Single Molecule Imaging Toolbox Extended).

## Overview

This repository provides streamlined workflows for:
- **SPT (Single Particle Tracking)**: Track single molecules across two color channels and analyze diffusion
- **SMLM (Single Molecule Localization Microscopy)**: High-density single molecule imaging and super-resolution reconstruction
- **SiMPull (Single Molecule Pull-down)**: Multi-channel colocalization analysis

Primarily designed for Olympus IX83 microscope data with EMCCD cameras.

---

## Repository Structure

```
EenSMITEscripts/
├── conversion/              # Olympus .vsi to .mat file conversion
├── tracking/                # Single particle tracking algorithms
├── analysis/                # Analysis workflows
│   ├── diffusion/          # Diffusion coefficient analysis (MSD, CDF)
│   ├── smlm/               # SMLM processing and super-resolution
│   ├── advanced/           # HMM, channel registration, clustering
│   └── simpull/            # SiMPull colocalization workflow
├── utilities/               # Helper functions
│   ├── export/             # Data export utilities (text, CSV, HDF5)
│   ├── plotting/           # Visualization and movie generation
│   ├── image_processing/   # Image conversion and processing
│   └── helpers/            # Miscellaneous utilities
├── workflows/               # High-level workflow entry points (future)
└── archive/                 # Outdated scripts for reference
```

---

## Quick Start

### Prerequisites

- MATLAB R2020b or later
- SMITE toolbox (Single Molecule Imaging Toolbox Extended)
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox

### Setup

1. Clone this repository:
   ```bash
   git clone https://github.com/IMAdams/EenSMITEscripts.git
   cd EenSMITEscripts
   ```

2. Add SMITE to your MATLAB path

3. Add this repository to your MATLAB path:
   ```matlab
   addpath(genpath('/path/to/EenSMITEscripts'))
   ```

---

## Workflows

### WORKFLOW 1: SPT (Single Particle Tracking) - 2 Channel

Track single molecules across two color channels and analyze diffusion behavior.

#### Step 1: File Conversion (Olympus .vsi → .mat)

**Script**: `conversion/Olympus_file_separation_EB2.m`

Open the script and use **Section A** for 2-channel sequence data with alternating frame acquisition.

**Parameters to set**:
- `fiducial = 0` or `1` (fiducial markers present?)
- `gain = 0` or `1` (gain calibration file?)
- `background = 0` or `1` (background file?)
- File path to your .vsi data

**Output**:
- `Channel_1.mat` - First color channel
- `Channel_2.mat` - Second color channel
- Preview images (.tif and .fig)

#### Step 2: Channel Registration (if needed)

**Script**: `analysis/advanced/IMA_channelreg_script.m`

Use fiducial images to create registration transforms for aligning channels.

**Output**:
- `RegistrationTransform_*.mat`

#### Step 3: Particle Tracking

**Script**: `tracking/IMA_trackingloop_ix83_11052025.m` ⭐ **CURRENT PRODUCTION VERSION**

**Key features**:
- EMCCD camera calibration (avoids sCMOS SPT bug)
- Dual-channel processing
- AnalysisID: '1106'

**Parameters to configure**:
- `rawDataDir` - Directory with Channel_1.mat and Channel_2.mat
- ROI settings
- Tracking parameters

**Output**:
- `SMD_*.mat` - Single molecule detections
- `SMF_*.mat` - Single molecule fits with tracking

**For batch multi-day processing**: Use `tracking/IMA_trackingloop072225.m`

#### Step 4: Diffusion Analysis

**Script**: `analysis/diffusion/IMA_diffusion_script.m` ⭐ **PRIMARY**

**Analysis includes**:
- MSD (Mean Square Displacement) ensemble fitting
- CDF (Cumulative Distribution Function) of jumps with 3-component Brownian fitting
- Per-cell and per-trajectory analysis
- Sliding window analysis

**Output**:
- Diffusion coefficients
- Summary statistics
- Publication-quality plots

**Optional**: Use `analysis/diffusion/IMA_diffusion_script_photons.m` for photon-filtered analysis

#### Step 5: Advanced Analysis (Optional)

- **Dimer Analysis**: `analysis/advanced/IMA_HMM_script.m` - Hidden Markov Model for protein dimer states
- **Clustering**: `analysis/smlm/IMA_BaGoL.m` - Bayesian Grouping of Localizations

---

### WORKFLOW 2: SMLM (Single Molecule Localization Microscopy)

High-density single molecule imaging for super-resolution reconstruction.

#### Step 1: File Conversion

**Script**: `conversion/Olympus_file_separation_EB2.m`

Use **Section B** for 2-channel snapshot data (SiMPull 488/647 format).

**Output**:
- Channel-specific .mat files
- Merged RGB images

#### Step 2: Batch SMLM Processing

**Script**: `analysis/smlm/IMA_pub_SMLM.m` ⭐ **PRODUCTION**

**Features**:
- Iterates over multiple SMF parameter sets
- Processes coverslip directories automatically
- Calls `smi.Publish.performFullAnalysis()`
- Converts cell images to PNG

**Output**:
- SMD/SMF files
- Analysis results
- PNG images of cells

#### Step 3: Advanced SMLM Analysis (Optional)

- **Super-Resolution Fitting**: `analysis/smlm/IMA_SR_analysis.m`
- **Bayesian Clustering**: `analysis/smlm/IMA_BaGoL.m`
  - Hierarchical RJMCMC clustering
  - Super-resolution image reconstruction
  - ROI-based processing

---

### WORKFLOW 3: SiMPull (Single Molecule Pull-down)

3-channel colocalization analysis (488/561/647 nm).

#### Steps

**Script**: `analysis/simpull/SiMPullMain.m`

This comprehensive 367-line workflow handles:
1. File conversion (use Section B of Olympus_file_separation_EB2.m)
2. QuadView data format processing
3. Channel registration
4. Gain/offset calibration
5. Overlap percentage calculations
6. Background correction

**Output**:
- Colocalization statistics
- Multi-channel visualizations

---

## Utilities Reference

### Export Tools (`utilities/export/`)

- `IMA_exportStructToText.m` - Export SMF parameters to text files
- `IMA_exportD_est_Nvals.m` - Export diffusion summary statistics (3 formats)
- `IMA_extractDataToCSV.m` - General data extraction to CSV
- `IMA_extractCDFJumpsToh5.m` - Export CDF jump data to HDF5

### Plotting Tools (`utilities/plotting/`)

- `IMA_plotTRAJ2.m` - Advanced trajectory visualization
- `IMA_TwoChannelMovieGenerator.m` - Generate 2-channel movies
- `IMA_plotROIDriver.m` / `MJW_plotROIDriver.m` - ROI-based plotting
- `IMA_makeFrame.m` / `IMA_makeFrameTwoChannels.m` - Frame rendering
- `MJW_cluster_overlay.m` - Cluster overlay visualization

### Image Processing (`utilities/image_processing/`)

- `IMA_cell_images_to_png.m` - Convert cell images to normalized PNG
- `writeTiff16.m` - Write 16-bit TIFF files
- `compute_lim.m` - Compute image limits

### Helper Functions (`utilities/helpers/`)

- `IMA_Cell_id.m` - Extract cell ID from filename
- `IMA_loadResultsFile.m` - Load analysis results
- `IMA_prepareAxes.m` - Prepare axes for plotting

---

## Key Scripts Summary

### Current Production Versions

| Script | Purpose | Location |
|--------|---------|----------|
| `Olympus_file_separation_EB2.m` | File conversion (all modes) | `conversion/` |
| `IMA_trackingloop_ix83_11052025.m` | SPT for IX83 (Nov 5, 2025) | `tracking/` |
| `IMA_diffusion_script.m` | Primary diffusion analysis | `analysis/diffusion/` |
| `IMA_pub_SMLM.m` | SMLM batch processing | `analysis/smlm/` |
| `IMA_HMM_script.m` | Dimer analysis | `analysis/advanced/` |
| `IMA_BaGoL.m` | Bayesian clustering | `analysis/smlm/` |
| `SiMPullMain.m` | SiMPull workflow | `analysis/simpull/` |

### Archived Scripts

Outdated scripts are in `archive/` for reference:
- `tracking_script_ix83_121222_EB.m` - 2022 version (use newer scripts instead)
- `Olympus_file_separation.m` - Basic converter (superseded by EB2 version)
- `Olympus_RegistrationFileSeparationRejoinQuadview.m` - Specialized converter

---

## Common Issues & Solutions

### Issue: sCMOS camera causing tracking errors
**Solution**: Use `IMA_trackingloop_ix83_11052025.m` which specifies EMCCD camera type

### Issue: Two-channel data not aligned
**Solution**: Run `analysis/advanced/IMA_channelreg_script.m` with fiducial images before tracking

### Issue: File paths not found
**Solution**: Ensure all subdirectories are in MATLAB path:
```matlab
addpath(genpath('/path/to/EenSMITEscripts'))
```

---

## File Naming Conventions

- `IMA_*` - Ian Adams scripts
- `MJW_*` - Michael Wester scripts
- `Olympus_*` - File conversion utilities
- Dated scripts (e.g., `_20240912_`) indicate specific versions

---

## Contributing

This is a research codebase under active development. Scripts are functional but may require parameter adjustments for your specific experimental setup.

---

## Notes

- **SPT IX83 workflow**: Uses EMCCD gain calibration to avoid sCMOS bug
- **Two-channel SPT**: Requires alternating frame acquisition mode in Olympus .vsi files
- **SMLM**: Best results with proper calibration files (gain, background)
- **Batch processing**: Multi-day tracking available via `IMA_trackingloop072225.m`

---

## Contact

For questions or issues, please contact the repository maintainer or open an issue on GitHub.

---

## License

Research use only. See repository license for details.
