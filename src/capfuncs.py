import os
import sys
import warnings
#
import numpy as np
import math
#
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, LightSource, ListedColormap
from matplotlib.patches import Patch
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable  
import matplotlib as mpl
mpl.rcParams['hatch.linewidth'] = 0.03
plt.rcParams['figure.constrained_layout.use'] = True
#
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
#
import xarray as xr
import rioxarray as rx
from rasterio.enums import Resampling
#
from scipy.ndimage import generic_filter
#
import numba as nb
import gc
#
########### ---------------------begin functions--------------------- ###########

def read_dem(filepath, bounding_box=None):
    """
    Opens digital elevation model (DEM) file. 
    Args:
        filepath: str. REQUIRED. The path to the digital elevation model file, relative to current directory.
        bounding_box: list. Default = None. A lat/lon bounding box for the domain of analysis in the form [xmin, ymin, xmax, ymax]. Use WGS84. 
    Returns:
        xarray.Dataset containing elevation variable and spatial y, x coordinates. 
    """
    print("reading file...")
    dem = rx.open_rasterio(filepath, chunks="auto") # open DEM with rioxarray

    # Convert to NAD83/Conus Albers if not already in that projection
    if dem.rio.crs is None:
        raise AttributeError("DEM has no CRS.")
    elif dem.rio.crs.to_epsg() != 5070:
        print("\033[1m" + "File Note: " + "\033[0m" + f"Reprojected from {dem.rio.crs} to EPSG:5070 (NAD83/Conus Albers)")
        dem = dem.rio.reproject("EPSG:5070")
    else:
        print("\033[1m" + "File Note: " + "\033[0m" + "DEM is already in NAD83/Conus Albers (EPSG:5070)")

    # clip to bounding box if provided
    if bounding_box:
        dem = dem.rio.clip_box(*bounding_box, crs="EPSG:4326")

    # return DEM as dataset
    if 'band' in dem.dims:
        if len(dem.band) != 1:
            raise ValueError("DEM contains multiple bands. Please ensure the DEM is a single-band raster .tif file.")
        ds = dem.squeeze(dim='band', drop=True).to_dataset(name='elevation')
    else:
        raise KeyError("DEM missing 'band' coordinate. Please ensure the DEM is a valid raster .tif file.")

    ds = ds.rio.set_spatial_dims(x_dim='x', y_dim='y')

    print('DEM successfully loaded, size',ds.elevation.shape)

    return ds

########### ---------------------next function------------------------ ###########

def add_colorbar(ax,mappable,label,ticks=None,position='right',orientation='vertical'):
    """
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
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(position, size="5%", pad=0.1, axes_class=plt.Axes) 
    cbar = plt.colorbar(mappable, cax=cax,orientation=orientation)
    if ticks is not None:
        cbar.set_ticks(np.arange(1,len(ticks)+1))
    cbar.set_label(label)
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.tick_params(labelsize=10)

########### ---------------------next function------------------------ ###########

def terrainplot(ds, bounding_box=None, cmap=None, contour=False, linewidths=0.05, grid=True, output='show', dpi=300):
    """
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
    """
    # define projection
    albers5070 = ccrs.AlbersEqualArea(central_longitude=-96,central_latitude=23,standard_parallels=(29.5, 45.5))

    # load values to plot
    elev = ds.elevation.values
    x = ds.x.values
    y = ds.y.values
    x2d, y2d = np.meshgrid(x, y)

    fig, ax = plt.subplots(subplot_kw={'projection': albers5070},figsize=(6.5,5))

    bounds = ds.rio.bounds()
    extent = [bounds[0], bounds[2], bounds[1], bounds[3]]

    if bounding_box:
        ax.set_extent(bounding_box, crs=ccrs.PlateCarree())
    else:
        ax.set_extent(extent, crs=albers5070)

    # establish colormap
    if cmap is None:
        hex_color1 = "#e8e7d3"
        hex_color2 = "#374157"
        cmap1 = LinearSegmentedColormap.from_list("custom_cmap", [hex_color1, hex_color2])
    else: 
        cmap1 = plt.get_cmap(cmap)

    # create and plot hillshade
    ls = LightSource(azdeg=335, altdeg=30)
    hillshade = ls.hillshade(elev, vert_exag=.04, dx=ds.rio.resolution()[0], dy=ds.rio.resolution()[0])
    ax.imshow(hillshade, cmap='gray', extent=extent,transform=albers5070) # plot hillshade

    # Plot elevation data
    cfill = ax.pcolormesh(x,y,elev, cmap=cmap, transform=albers5070, alpha=0.65)
    add_colorbar(ax, cfill, 'Elevation (m)') # elevation data colorbar

    # plot contour lines if contour=True
    if contour:
        contour_levels = np.arange(math.floor(elev.min() / 100) * 100, math.ceil(elev.max() / 100) * 100 + 100, 100) # define contour levels for every 100m
        ax.contour(x, y, elev, levels=contour_levels, colors='black', transform=albers5070, linewidths=linewidths) # contour every 100m
        contour_levels = np.arange(math.floor(elev.min() / 1000) * 1000, math.ceil(elev.max() / 1000) * 1000 + 1000, 1000) # define contour levels for every 1000m
        contour = ax.contour(x, y, elev, levels=contour_levels, colors='black', transform=albers5070, linewidths=linewidths*1.5) # contour every 1000m
        ax.clabel(contour, levels=contour_levels, inline=True, fontsize=4) # label 1000m contours

    # plot gridlines if grid=True
    if grid:
        gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.4, color="white", alpha=0.6)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {"size": 8}
        gl.ylabel_style = {"size": 8}
        gl.xformatter = LongitudeFormatter()
        gl.yformatter = LatitudeFormatter()

    # save or show figure
    if output == 'show':
        plt.show()
    elif output == 'save':
        plt.savefig('figures/terrain.png', dpi=dpi, bbox_inches='tight')
        print('saved figure to figures/terrain.png')

########### ---------------------next function------------------------ ###########
def get_slope(elevation, resolution_x,resolution_y):
    """
    Calculate slope (degrees) at each pixel within the digital elevation model.  
    Args:
        elevation: array. REQUIRED. 2D array of elevations (meters).
        resolution_x: int or float. REQUIRED. x-grid spacing in meters (DEM resolution).
        resolution_y: int or float. REQUIRED. y-grid spacing in meters (DEM resolution).
    Returns:
        slope: 2D numpy array of slope (degrees) for each pixel. 
    """ 
    dz_dy, dz_dx = np.gradient(elevation,np.abs(resolution_y),np.abs(resolution_x)) # calculate gradient in z (elevation) with respect to x and y
    slope_radians = np.sqrt(dz_dx**2 + dz_dy**2) # convert to radians
    slope = np.degrees(np.arctan(slope_radians)) # convert to degrees

    return slope

########### ---------------------next function------------------------ ###########

@nb.njit(parallel=True, fastmath=True) # use numba for improved speed
def rolling_rank(elevation, hop_x, hop_y):
    """
    Calculates rank at each pixel within the digital elevation model.
    Args:
        elevation: array. REQUIRED. 2D array of elevations (meters).
        hop_x: int. REQUIRED. Hop distance in x-direction (pixels).
        hop_y: int. REQUIRED. Hop distance in y-direction (pixels).
    Returns:
        rank: 2D numpy array of rank values (fraction) for each pixel.
    """
    ny, nx = elevation.shape 

    rank = np.full((ny - 2*hop_y, nx - 2*hop_x), np.nan, dtype=np.float32) # initialize rank array
    
    # iterate through each pixel in elevation array, grab window of neighbors, and calculate rank
    for i in nb.prange(hop_y, ny - hop_y):
        for j in range(hop_x, nx - hop_x):
            window = elevation[i-hop_y:i+hop_y+1, j-hop_x:j+hop_x+1] # create window of neighbors around center pixel
            center_val = elevation[i, j]
            total_count = window.size - 1
            less_count = np.count_nonzero(window < center_val) # count # of pixels of lower elevation than center pixel
            rank[i-hop_y, j-hop_x] = less_count / total_count # rank as fraction
    
    return rank

########### ---------------------next function------------------------ ###########

def get_rank_curvature(elevation, resolution_x, resolution_y, r, max_margin_x, max_margin_y):
    """
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
    """
    # calculate hop distances (distance to travel in each direction from center pixel for rank/curvature calculations)
    hop_x = int(np.ceil(r/np.abs(resolution_x))) # hop distance in pixels
    hop_y = int(np.ceil(r/np.abs(resolution_y))) # hop distance in pixels

    # calculate curvature
    # this loop is for handling edge case where hop distance equals the max margin
    if (hop_x == max_margin_x) & (hop_y == max_margin_y):
        center = elevation[max_margin_y:-max_margin_y,max_margin_x:-max_margin_x]

        north = elevation[:-(max_margin_y+hop_y),max_margin_x : -max_margin_x]  # north neighbor
        south = elevation[(max_margin_y+hop_y):,max_margin_x : -max_margin_x] # south neighbor
        west  = elevation[max_margin_y : -max_margin_y,:-(max_margin_x + hop_x)] # west neighbor
        east  = elevation[max_margin_y : -max_margin_y,(max_margin_x + hop_x):] # east neighbor
        
        curvature_term_a = (2 * center) - (0.5 * (north + south + west + east))

        nw = elevation[:-(max_margin_y + hop_y),:-(max_margin_x + hop_x)] # northwest neighbor
        ne = elevation[:-(max_margin_y + hop_y),max_margin_x + hop_x :] # northeast neighbor
        sw = elevation[max_margin_y + hop_y:,:-(max_margin_x + hop_x)] # southwest neighbor
        se = elevation[max_margin_y + hop_y:,max_margin_x + hop_x :] # southeast neighbor

        curvature_term_b = (2*center) - (0.5*(nw + ne + sw + se))

        curvature = 0.25 * (
            curvature_term_a / (2*r) +
            curvature_term_b / (2*r*np.sqrt(2))
        )

    # this loop is for all other cases where hop distance is less than the max margin
    else:
        center = elevation[max_margin_y:-max_margin_y,max_margin_x:-max_margin_x]

        north = elevation[(max_margin_y-hop_y):-(max_margin_y+hop_y),max_margin_x : -max_margin_x]  
        south = elevation[(max_margin_y+hop_y):-(max_margin_y-hop_y),max_margin_x : -max_margin_x]
        west  = elevation[max_margin_y:-max_margin_y,(max_margin_x-hop_x):-(max_margin_x + hop_x)]
        east  = elevation[max_margin_y:-max_margin_y,(max_margin_x + hop_x):-(max_margin_x - hop_x)]

        curvature_term_a = (2 * center) - (0.5 * (north + south + west + east))

        nw = elevation[(max_margin_y-hop_y):-(max_margin_y + hop_y),(max_margin_x-hop_x):-(max_margin_x + hop_x)]
        ne = elevation[(max_margin_y-hop_y):-(max_margin_y + hop_y),max_margin_x + hop_x:-(max_margin_x - hop_x)]
        sw = elevation[max_margin_y+hop_y:-(max_margin_y-hop_y),(max_margin_x-hop_x):-(max_margin_x+hop_x)]
        se = elevation[max_margin_y+hop_y:-(max_margin_y-hop_y),(max_margin_x+hop_x):-(max_margin_x-hop_x)]

        curvature_term_b = (2*center) - (0.5*(nw + ne + sw + se))

        curvature = 0.25 * (
            curvature_term_a / (2*r) +
            curvature_term_b / (2*r*np.sqrt(2))
        )

    # calculate rank
    rank_full = rolling_rank(elevation, hop_x, hop_y)

    # trim rank to match curvature
    rank = rank_full[max_margin_y-hop_y:rank_full.shape[0]-(max_margin_y-hop_y),
                     max_margin_x-hop_x:rank_full.shape[1]-(max_margin_x-hop_x)]

    return rank, curvature

########### ---------------------next function------------------------ ###########

def classify_CAP(slope, curvature, rank):
    """
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
    """
    CAP = np.zeros(slope.shape) # initialize array for CAP values
    percentile = -(1/60)*slope+0.5 # calculate perventile threshold, given by Lundquist et al., 2008
    CAP[rank<percentile] = 1 # CAP (CAP = 1) when rank is less than percentile
    CAP[slope>30] = -1 # No CAP (CAP = -1) when slope is greater than 30 degrees
    CAP[curvature>0] = -1 # or No CAP (CAP = -1) when curvature is positive

    return CAP

########### ---------------------next function------------------------ ###########

def cap_analysis(terrain_data, xv_dist):
    """
    Integrates all necessary calculations to determine CAP signal at each pixel across the DEM. 
    Args:
        terrain_data: xarray.Dataset. REQUIRED. xarray.Dataset with variable 'elevation' with coordinates 'y' and 'x'. 
        xv_dist: list. REQUIRED. list or array of cross-valley distances (meters).
    Returns:
        None. Saves elevation and calculated slope, rank, curvature, and CAP iteratively to a zarr directory in output/. 
    """
    if xv_dist is None:
        raise ValueError("the xv_dist arg is required. Please input a list or array of values with at least one element.")
    
    if len(xv_dist) == 0:
        raise ValueError("xv_dist cannot be empty. Must be a list or array of values with at least one element.")

    rvals = np.divide(xv_dist,2) # calculate critical radius from cross-valley distance
    rvals = sorted(rvals, reverse=True) # sort from largest to smallest

    crs = terrain_data.rio.crs # grab crs from original DEM
    attrs = terrain_data.attrs # grab attributes from original DEM
    resolution_x, resolution_y = terrain_data.rio.resolution() # grab resolution in x and y from original DEM
    
    x_1d = terrain_data.x.to_numpy()            
    y_1d = terrain_data.y.to_numpy()           

    # SUBSET THE DEM
    max_margin_x = int(math.ceil(rvals[0]/np.abs(resolution_x)))
    max_margin_y = int(math.ceil(rvals[0]/np.abs(resolution_y)))

    elevation = np.ascontiguousarray(terrain_data['elevation'].values) 
    subset_elevation = elevation[max_margin_y:-max_margin_y,max_margin_x:-max_margin_x] # subset elevation to the smallest necessary size corresponding to the largest rval
    subset_x = x_1d[max_margin_x:-max_margin_x] # subset x to the smallest necessary size corresponding to the largest rval
    subset_y = y_1d[max_margin_y:-max_margin_y] # subset y to the smallest necessary size corresponding to the largest rval

    # CALCULATE SLOPE
    slope = get_slope(elevation=subset_elevation,resolution_x=resolution_x,resolution_y=resolution_y)

    # SAVE ELEVATION AND SLOPE TO ZARR (independent of xv_dist)
    store = "output/CAP_layers.zarr"

    ds_layer = xr.Dataset(
        data_vars=dict(
            elevation=(["y", "x"], subset_elevation),
            slope=(["y", "x"], slope),
        ),
        coords=dict(
            x=("x", subset_x),  
            y=("y", subset_y), 
            xv_dist=[] # empty dimension for appending later
        ),
        attrs = attrs)

    ds_layer.rio.set_spatial_dims(x_dim='x', y_dim='y')
    ds_layer.rio.write_crs(crs,inplace=True)

    ds_layer.to_zarr(store,mode="w")    

    # prepare empty lists to store CAP/rank/curvature arrays
    CAP_list = []
    rank_list = []
    curvature_list = []

    # ITERATE THROUGH CROSS-VALLEY DISTANCES
    print('Iterating through critical values...')
    sys.stdout.write("a progress indicator will appear shortly")
    sys.stdout.flush()
    for index,r in enumerate(rvals):
        xv = r*2
        print("\n processing xv_dist:",xv,"...")

        # calculate rank, curvature
        rank, curvature = get_rank_curvature(elevation=elevation, 
                                                resolution_x=resolution_x,resolution_y=resolution_y, 
                                                max_margin_x=max_margin_x, max_margin_y=max_margin_y,
                                                r=r,
                                                )
        # classify CAP
        CAP = classify_CAP(slope=slope, curvature=curvature, rank=rank)

        rank = np.expand_dims(rank,axis=2) # expand from (mxn) to (mxnx1)
        curvature = np.expand_dims(curvature,axis=2) # expand from (mxn) to (mxnx1)
        CAP = np.expand_dims(CAP,axis=2) # expand from (mxn) to (mxnx1)

        ds_layer = xr.Dataset(
            data_vars=dict(
                CAP=(["y", "x","xv_dist"], CAP),
                rank=(["y", "x","xv_dist"],rank),
                curvature=(["y", "x","xv_dist"],curvature)
            ),
            coords=dict(
                x=("x", subset_x),  
                y=("y", subset_y), 
                xv_dist=("xv_dist",[xv])
            ),
            attrs = attrs) 
        ds_layer = ds_layer.assign_attrs({'CAP Signal': '-1=no-CAP, 0=no-Signal, 1=CAP'}) # assign attribute for future reference

        ds_layer.rio.set_spatial_dims(x_dim='x', y_dim='y')
        ds_layer.rio.write_crs(crs,inplace=True)
        ds_layer.to_zarr(store, mode="a", append_dim="xv_dist") # append to zarr store

        percent_complete = ((index + 1) / len(rvals)) * 100  
        sys.stdout.write(f"\r{percent_complete:.1f}% complete")  
        sys.stdout.flush()  

        del rank, curvature, CAP # remove variables from memory before next loop
        gc.collect() # collect trash      

    percent_complete = 100
    sys.stdout.write(f"\r{percent_complete:.1f}% complete \n")  
    sys.stdout.flush()  
    
    print('CAP analysis complete, saved to zarr.')

########### ---------------------next function------------------------ ###########

def plotcap(terrain_stats,xv_dist=None,bounding_box=None,depth_plot=False,shade=False,contour=False,linewidths=0.05,grid=True,output='show',dpi=300):
    """
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
    """
    elev = terrain_stats['elevation'].values # grab elevation values as array
    x = terrain_stats['x']
    y = terrain_stats['y']

    albers_proj = ccrs.epsg(5070) # define projection
                        
    fig, ax = plt.subplots(subplot_kw={'projection': albers_proj})

    hex_color1 = "#E9E8D4"  
    hex_color2 = "#3B465C"   
    cmap1 = LinearSegmentedColormap.from_list("custom_cmap", [hex_color1, hex_color2])

    extent = [x.min(), x.max(), y.min(), y.max()]

    # set extent based on bounding box, if provided (otherwise use full extent)
    if bounding_box:
        ax.set_extent(bounding_box, crs=ccrs.PlateCarree())
    else:
        ax.set_extent(extent, crs=albers_proj)
        bounding_box = extent

    # create and plot hillshade
    ls = LightSource(azdeg=335, altdeg=30)
    hillshade = ls.hillshade(elev, vert_exag=.8, dx=terrain_stats.rio.resolution()[0], dy=terrain_stats.rio.resolution()[0])
    ax.imshow(hillshade, cmap='gray', extent=[x.min(), x.max(), y.min(), y.max()],transform=albers_proj)

    # plot elevation data
    cfill = plt.pcolormesh(x, y, elev, cmap=cmap1, transform=albers_proj, alpha=0.5)
    add_colorbar(ax, cfill, 'Elevation (m)')  # Adds external colorbar

    if contour:
        contour_levels = np.arange(math.floor(elev.min() / 100) * 100, math.ceil(elev.max() / 100) * 100 + 100, 100) # define contour levels for every 100m
        ax.contour(x, y, elev, levels=contour_levels, colors='black', transform=albers_proj, linewidths=linewidths) # contour every 100m
        contour_levels = np.arange(math.floor(elev.min() / 1000) * 1000, math.ceil(elev.max() / 1000) * 1000 + 1000, 1000) # define contour levels for every 1000m
        contour = ax.contour(x, y, elev, levels=contour_levels, colors='black', transform=albers_proj, linewidths=linewidths*1.5) # contour every 1000m, thicker lines
        ax.clabel(contour, levels=contour_levels, inline=True, fontsize=4) # label 1000m contours
    
    # if depth_plot=True, shade based on number of iterations classified as CAP.
    # consecutively, only shade where previous (larger) xv_dist also classified as CAP.
    if depth_plot:
        cmap2 = plt.get_cmap('Blues')

        num_layers = len(terrain_stats['xv_dist'])
        mask = xr.ones_like(terrain_stats['CAP'].isel(xv_dist=0),dtype='bool')  # Start with a mask of ones
        CAP_masked = xr.zeros_like(terrain_stats['CAP'].isel(xv_dist=0))  

        for i in range(num_layers):
            CAP = terrain_stats['CAP'].isel(xv_dist=i)
            mask &= (CAP == 1)  # retain only where current and previous layers at least 1
            CAP_masked = CAP_masked + mask.where(mask, other=0) 
        CAP_contours = CAP_masked.copy()
        CAP_masked = CAP_masked.where(CAP_masked != 0, np.nan)

        cfill = ax.pcolormesh(x, y, CAP_masked, cmap=cmap2, transform=albers_proj, alpha=0.75)
        ax.contour(x, y, CAP_contours, levels=[CAP_contours.min()],colors='black', transform=albers_proj,linewidths=linewidths*4)
        add_colorbar(ax, cfill, '# of iterations',orientation='horizontal',position='top',ticks=np.arange(num_layers))  # Adds external colorbar

    else: 
        cap_mask = np.where(terrain_stats['CAP'].sel(xv_dist=xv_dist) == 1, 1, np.nan)
        no_signal_mask = np.where(terrain_stats['CAP'].sel(xv_dist=xv_dist) == 0, 1, np.nan)

        ax.contour(x,y,terrain_stats['CAP'].sel(xv_dist=xv_dist), levels=[0,1,2],colors='black', transform=albers_proj, linewidths=0.07)
        legend_patch = Patch(facecolor='skyblue', label='Cold-air pool', alpha=.6)
        ax.legend(handles=[legend_patch], loc='upper right')

        cfill = ax.pcolormesh(
            x, y, cap_mask,
            cmap=ListedColormap(['skyblue']),
            transform=albers_proj,
            alpha=.6,
            shading='auto',
            label='Cold-Air Pool'
        )    

        if shade:
            ax.contourf(terrain_stats.x,terrain_stats.y,no_signal_mask,hatches=[6*'.'],colors='none',transform=albers_proj)

    # plot gridlines and tick labels if grid=True
    if grid:
        gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.4, color="white", alpha=0.6)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {"size": 8}
        gl.ylabel_style = {"size": 8}
        gl.xformatter = LongitudeFormatter()
        gl.yformatter = LatitudeFormatter()

        gl.xlabel_style = {"size": 8, "rotation": 30,"ha":'right'}   # rotate 30°
        gl.ylabel_style = {"size": 8, "rotation": -7,"ha":'right'}
    
    if output == 'show':
        plt.show()
    elif output == 'save':
        plt.savefig('figures/CAP_map.png', dpi=dpi, bbox_inches='tight')
        print('saved figure to CAP_map_test.png')