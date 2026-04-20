# **Python Implementation of an Automated Algorithm for Mapping Regions of Cold-Air Pooling in Complex Terrain**

Cramblitt, J. & Lundquist, J.

Cramblitt, J. & Lundquist, J.

Version: 2.0.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17541537.svg)](https://doi.org/10.5281/zenodo.17541537)

Version 2.0.0 corrects a typo in the curvature formula and improves handling of coordinate systems, including transitioning from the geographic WGS84 coordinate sytem (EPSG:4326) to the projected NAD83/Conus Albers coordinate system for accurate distance-dependent calculations of slope, rank and curvature. These corrections resulted in minor changes to terrain classifications, but are more faithful to the original algorithm of Lundquist et al., [2008](#references). Additionally, performance is now drastically improved through targeted improvements to the **get_rank_curvature()** function and calculations are now saved iteratively to a zarr directory. Some unecessary functions have been removed, and one helper function for rank calculations was added. 

Originally developed by Lundquist et al., [2008](#references), and implemented in Python in this repository by John Cramblitt.

## CONTACT

John Cramblitt \
jcramblitt@berkeley.edu \
University of Washington, Seattle \
University of California, Berkeley

## REQUIREMENTS

Python (3.11 or later) \
numpy, matplotlib, cartopy, scipy, xarray, rioxarray, numba, dask, zarr

## INSTALLATION

Clone the repository and install dependencies. 

With Anaconda (recommended):
```
conda env create -f environment.yml

conda activate cap-algorithm
```
With pip:
```
pip install numpy matplotlib cartopy scipy xarray rioxarray numba dask zarr
```
## FILE STRUCTURE

- `data/`
    - `DEM_30m.tif` (accessed and resampled from [USGS 3D Elevation Program](#references) using the [py3dep](#references) module)
- `figures/` (empty, for plots)
- `output/` (empty, for zarr directory)
- `src/`
    - `capfuncs.py`
- `environment.yml`
- `notebook.ipynb`
- `README.md`

## QUICKSTART

Step by step instructions are provided in the example notebook, but the core functionality is contained in two functions:

```python
from src.capfuncs import read_dem, cap_analysis 

# 1. load DEM (.tif format). Replace with your filepath if DEM is not in data/.
terrain_data = read_dem()

# 2. run CAP analysis. 
#    - xv_dist = list of cross-valley distances in meters (e.g., [3000]).
ds = cap_analysis(terrain_data=terrain_data, xv_dist=[3000])
```

## ALGORITHM FUNCTION

Lundquist et al., [2008](#references) developed an automated algorithm to map regions of cold-air pooling (CAP) in mountain terrain from a digital elevation model (DEM). Terrain analysis identifies areas of concave, low-lying topography where cold air is likely to accumulate. Terrain-based classifications provide a strong first-order guess for regions of CAP, though other meteorological controls (particularly valley winds and cloudiness) may influence local CAP dynamics. Algorithm output should be considered alongside a variety of factors, and ideally evaluated against local observations. 

The algorithm was designed with narrow, mountain valley terrain features in mind. It is not explicitly intended to describe large-scale cold-air trapping, such as what occurs in the Columbia Basin in Eastern Washington or Salt Lake City in Utah, so proceed with caution in large-scale applications.

This implementation is not fully optimized for performance, and system memory limitations may arise when analyzing large regions or high-resolution DEMs. Please contact jcramblitt@berkeley.edu if you are facing significant performance issues. 

The Python implementation in this repository functions as outlined: 

### **1. Upload Digital Elevation Model (DEM)**

The **read_dem()** function automatically reads and formats .tif files for processing.  

The domain of the file must be larger than the domain of interest to accomodate calculations involving nearby pixels. Specifically, the domain of interest must be padded on all sides by a length equal to half the user-chosen cross-valley distance (see discussion on cross-valley distance and critical radii below). This can be expressed by number of pixels as follows:

$\text{Number of Margin Pixels} \geq \frac{r}{\text{DEM resolution}}$

### **2. Calculate Terrain Stats**

The **cap_analysis()** function integrates a series of calculations, and outputs CAP classifications across the DEM. 

For each pixel, the algorithm calculates slope, rank, curvature, and finally the CAP signal. The 'xv_dist' argument is critical to the calculation of rank and curvature. The internal functions used to calculate each of these variables are summarized at the bottom of this file, and exact code can be viewed in the capfuncs.py source code file. 

*Slope* is calculated with numpy's gradient function. A critical radius *r*, is user-defined as half the approximate cross-valley distance. For each pixel, *rank* is the proportion of pixels within the surrounding square of radius *r* which are lower in elevation than that pixel. *Curvature* is calculated at each pixel according to:

$cv = \frac{1}{4} (\frac{1}{2r} (- \frac{z_w+z_e}{2}+z-\frac{z_n+z_s}{2})) + \frac{1}{2\sqrt{2}r}(z-\frac{z_{sw}+z_{ne}}{2}+z-\frac{z_{nw}+z_{se}}{2})$

where *z* is the elevation at that pixel, *r* is the critical radius, $z_{n/s/e/w}$ is the elevation of the pixel a distance *r* to the north/south/east/west (i.e., the elevation at the edge centers of the local square used to calculate *rank*), and $z_{ne/se/sw/nw}$ is the elevation of the pixel a distance *r* in both specified directions (i.e., the elevation at the corners of the local square used to calculate *rank*). This formula was first introduced in Liston & Elder, 2006. 

CAP classifications are determined by the following conditions:
1) No CAP anywhere where slope is greater than 30°
2) No CAP anywhere where curvature is positive
3) For the remaining pixels, the choice between CAP and no signal is determined by the cutoff: $rank =-\frac{1}{60}slope+\frac{1}{2}$. \
    3a. CAP where $rank < -\frac{1}{60}slope+\frac{1}{2}$ \
    3b. no signal where $rank > -\frac{1}{60}slope+\frac{1}{2}$ 

Lundquist et al., 2008 provide further detail on algorithm development. Corrections for topographic amplification factor are not yet included in this repository. 

DEM size, resolution, and number of iterations through distance values (i.e., the length of the array passed to xv_dist) all influence algorithm speed. The algorithm is quick to run over small scales at a reasonable resolution, but is not fully optimized for performance and may fail under system memory constraints for large or high resolution DEMs. Future versions may improved speed or memory handling if this poses a barrier for users.   

### **3. Plotting**

The notebook provides a few crude plotting examples, as well as one custom plotting function to visualize CAP distributions across terrain. This function is provided for illustrative purposes only. If you use or adapt figures generated by this code for publication, please acknowledge this repository.

## ACRONYMS

DEM - digital elevation model \
CAP - cold-air pooling

## CITATION

If you use this code, please cite this repository via its DOI, and Lundquist et al., [2008](#references) for algorithm design. 

## LICENSE

This work is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## REFERENCES

Liston, G. E., and K. Elder, 2006: A Meteorological Distribution System for High-Resolution Terrestrial Modeling (MicroMet). J. Hydrometeor., 7, 217–234, https://doi.org/10.1175/JHM486.1.

Lundquist, J. D., N. Pepin, and C. Rochford (2008), Automated algorithm for mapping regions of CAP in complex terrain, J. Geophys. Res., 113, D22107, doi:10.1029/2008JD009879.

U.S. Geological Survey. (2024). 1 Arc-second Digital Elevation Models (DEMs) - USGS National Map 3DEP Downloadable Data Collection. Retrieved from https://www.usgs.gov/3d-elevation-program.

Chegini, T., Li, H.-Y., & Leung, L. R. (2021). HyRiver: Hydroclimate Data Retriever. Journal of Open Source Software, 6(66), 1–3. https://doi.org/10.21105/joss.03175.

## FUNCTIONS

**read_dem** (filepath='data/DEM_30m.tif', bounding_box=None):

    Opens digital elevation model (DEM) file. 

    Args:
        filepath: str. REQUIRED. The path to the digital elevation model file, relative to current directory.
        bounding_box: list. Default = None. A lat/lon bounding box for the domain of analysis in the form [xmin, ymin, xmax, ymax]. Use WGS84. 

    Returns:
        xarray.Dataset containing elevation variable and spatial y, x coordinates. 

**terrainplot**(ds, bounding_box=None, cmap=None, contour=False, linewidths=0.05, grid=True, output='show', dpi=300):

    Create terrain map of the DEM Dataset with hillshade and optional contour lines. 

    Args:
        ds: xarray.Dataset. REQUIRED. xarray.Dataset containing elevation variable and spatial coordinates y, x.
        bounding_box: list. Default = None. Lat/lon bounding box in the form [xmin, ymin, xmax, ymax]. Use WGS84. 
        cmap: string. Default = None. Optional custom colormap, plots with tan/blue shading by default.
        contour: bool. Default = False. Option to overlay contour lines.
        linewidths: float. Default = 0.05. Linewidth for contour lines.
        grid: bool. Default = True. Draw lat/lon gridlines and lat/lon labels.
        output: str. 'show' or 'save'. Default = 'show.' Save figure or display it.
        dpi: int or float. Default = 300. Resolution for saved figure.

    Returns: 
        None. Displays plot or saves .png to figures/.

**cap_analysis** (terrain_data, resolution, xv_dist):

    Integrates all necessary calculations to determine CAP signal at each pixel across the DEM. 

    Args:
        terrain_data: xarray.Dataset. REQUIRED. xarray.Dataset with variable 'elevation' with coordinates 'y' and 'x'. 
        xv_dist: list. REQUIRED. list or array of cross-valley distances (meters).

    Returns:
        None. Saves elevation and calculated slope, rank, curvature, and CAP iteratively to a zarr directory in output/. 

**plotcap**(terrain_stats, xv_dist, bounding_box=None, depth_plot=False, shade=False, contour=False, linewidths=0.05, grid=True, output='show', dpi=300):

    Create contour map of the digital elevation model and shade topography prone to CAP. 

    Args:
        terrain_stats: xarray.Dataset. REQUIRED. xarray.Dataset containing elevation and CAP variables with spatial coordinates y, x. Must be in the NAD83 / Conus Albers coordinate system (EPSG:5070).
        xv_dist: coordinate value. Default = None. REQUIRED if depth_plot=False. The cross-valley distance (meters) for which to plot the CAP signal.
        bounding_box: list. Default = None. Lat/lon bounding box in the form [xmin, xmax, ymin, ymax]. Use WGS84.
        depth_plot: bool. Default = False. Option to shade regions based on the number of iterations for which each pixel was classified as CAP.
        shade: bool. Default = False. Option to shade regions with no signal.
        contour: bool. Default = False. Option to overlay elevation contour lines.
        linewidths: float. Default = 0.05. Linewidth for contour lines. Also scales CAP region contours. 
        output: str. Default = 'show'. Choose from 'save' or 'show'. Option to output the figure or save it to PNG.
        dpi: int or float. If output='save', dpi is the resolution of the saved figure.

    Returns:
        None. Displays plot or saves .png to figures/.

## INTERNAL FUNCTIONS

**get_slope** (subset_elevation, resolution):

    Calculate slope (degrees) at each pixel within the digital elevation model.  

    Args:
        elevation: array. REQUIRED. 2D array of elevations (meters).
        resolution_x: int or float. REQUIRED. x-grid spacing in meters (DEM resolution).
        resolution_y: int or float. REQUIRED. y-grid spacing in meters (DEM resolution).

    Returns:
        slope: 2D numpy array of slope (degrees) for each pixel. 

**rolling_rank** (elevation, hop_x, hop_y):

    Calculates rank at each pixel within the digital elevation model.

    Args:
        elevation: array. REQUIRED. 2D array of elevations (meters).
        hop_x: int. REQUIRED. Hop distance in x-direction (pixels).
        hop_y: int. REQUIRED. Hop distance in y-direction (pixels).

    Returns:
        rank: 2D numpy array of rank values (fraction) for each pixel.

**get_rank_curvature** (terrain_data, subset_elevation, margin_pixels, resolution, r):

    Calculates rank and curvature at each pixel within the digital elevation model following Lundquist et al., 2008. 

    Args:
        elevation: array. REQUIRED. A 2D numpy array of elevation including margin pixels.
        resolution_x: float or int. REQUIRED. Grid spacing of the DEM in the x-direction (meters).
        resolution_y: float or int. REQUIRED. Grid spacing of the DEM in the y-direction (meters).
        r: float or int. REQUIRED. Critical radius used in rank and curvature calculations (meters).
        max_margin_x: int. REQUIRED. Maximum margin in the x-direction (pixels) corresponding to the largest r-value for the iterative set.
        max_margin_y: int. REQUIRED. Maximum margin in the y-direction (pixels) corresponding to the largest r-value for the iterative set.

    Returns:
        rank, curvature: Tuple of 2D numpy arrays (rank, curvature).

**classify_CAP** (subset_elevation, slope, curvature, rank):

    Classifies cold-air pooling signal at pixel in the digital elevation model. 

    Args:
        slope: REQUIRED. 2D numpy array of slope (degrees).
        curvature: REQUIRED. 2D numpy array of curvature.
        rank: REQUIRED. 2D numpy array of rank. 

    Returns:
        2D numpy array of cold-air pooling signal for each pixel in the digital elevation model.
            -1 = no-CAP
            0 = no-signal
            1 = CAP

**add_colorbar** (ax, mappable, label,ticks=None, position='right', orientation='vertical'):
    Adds an external colorbar to a plot.

    Args:
        ax: matplotlib axis object. REQUIRED. The axis to which the colorbar will be added
        mappable: matplotlib mappable object. REQUIRED. The mappable object to which the colorbar will be linked.
        label: str. REQUIRED. The label for the colorbar.
        ticks: list. Default = None. A list of tick values for the colorbar.
        position: str. Default = 'right'. The position of the colorbar. Choose from 'left', 'right', 'top', or 'bottom'.
        orientation: str. Default = 'vertical'. The orientation of the colorbar, either 'vertical' or 'horizontal'.

    Returns:
        None. Adds colorbar to the specified axis.
