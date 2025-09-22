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

########### ---------------------begin functions--------------------- ###########

def read_dem(filepath='data/DEM_30m.tif', bounding_box = None):
    """
    Opens digital elevation model (DEM) file. 
    Args:
        filepath: str. Default='data/DEM_30m.tif'. The path to the digital elevation model file, relative to current directory.
        bounding_box: list. Default = None. A lat/lon bounding box for the domain of analysis in the form [xmin, xmax, ymin, ymax]. Use WGS84. 
    Returns:
        xarray.Dataset containing elevation variable and lat, lon coordinates. 
    """
    print("reading file...")
    dem = rx.open_rasterio(filepath, chunks="auto") # open DEM file with rioxarray

    # Convert to WGS84 if not already in that projection
    if dem.rio.crs is None:
        warnings.warn("DEM has no CRS info. Assumed file is in WGS84 (EPSG:4326).", UserWarning)
        dem = dem.rio.write_crs("EPSG:4326")
    elif dem.rio.crs.to_epsg() != 4326:
        print("\033[1m" + "File Note: " + "\033[0m" + f"Reprojected from {dem.rio.crs} to EPSG:4326 (WGS84)")
        dem = dem.rio.reproject("EPSG:4326")
    else:
        print("\033[1m" + "File Note: " + "\033[0m" + "DEM File is already in WGS84 (EPSG:4326)")

    if bounding_box:
        xmin, xmax, ymin, ymax = bounding_box
        dem = dem.rio.clip_box(minx=xmin, miny=ymin, maxx=xmax, maxy=ymax)

    # return DEM as dataset
    ds = dem.to_dataset(dim="band").rename({1: "elevation"})
    ds = ds.rename({"x": "lon", "y": "lat"})
    ds = ds.rio.set_spatial_dims(x_dim='lon', y_dim='lat')

    print('DEM succesfully loaded, size',ds.elevation.shape)

    return ds

########### ---------------------next function------------------------ ###########

def terrainplot(ds,bounding_box=None,cmap=None,contour=False,linewidths=0.05,grid=True,output='show',dpi=300):
    """
    Create terrain map of the DEM Dataset with hillshade and optional contour lines.  
    Args:
        ds: xarray.Dataset. REQUIRED. xarray.Dataset containing elevation variable and coordinates lat, lon.
        bounding_box: list. Default = None. Lat/lon bounding box in the form [xmin, xmax, ymin, ymax]. Use WGS84. 
        cmap: string. Default = None. Optional custom colormap, plots with tan/blue shading by default.
        contour: bool. Default = False. Option to overlay contour lines.
        linewidths: float. Default = 0.05. Linewidth for contour lines.
        grid: bool. Default = True. Draw lat/lon gridlines.
        output: str. 'show' or 'save'. Default = 'show.' Save figure or display it.
        dpi: int or float. Default = 300. Resolution for saved figure.
    Returns: 
        None. Displays plot or saves .png to examples/figures/.
    """
    # set bounding box and load subset data into memory.  
    if bounding_box:
        ds_subset = ds.sel(
            lon=slice(bounding_box[0], bounding_box[1]),
            lat=slice(bounding_box[3], bounding_box[2])  
        )
    else:
        bounds = ds.rio.bounds()
        bounding_box = [bounds[0], bounds[2], bounds[1], bounds[3]] 
        ds_subset = ds
    ds_subset = ds_subset.load()

    # load values to plot
    elev = ds_subset.elevation.values
    lon = ds_subset.lon.values
    lat = ds_subset.lat.values
    lon2d, lat2d = np.meshgrid(lon, lat)

    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()},figsize=(6.5,5))

    # establish colormap
    if cmap is None:
        hex_color1 = "#e8e7d3"
        hex_color2 = "#374157"
        cmap1 = LinearSegmentedColormap.from_list("custom_cmap", [hex_color1, hex_color2])
    else: 
        cmap1 = plt.get_cmap(cmap)

    # set extent based on bounding box
    extent = [lon.min(), lon.max(), lat.min(), lat.max()]
    ax.set_extent(bounding_box, crs=ccrs.PlateCarree())

    # create and plot hillshade
    ls = LightSource(azdeg=335, altdeg=30)
    hillshade = ls.hillshade(elev, vert_exag=.04, dx=1, dy=1)
    ax.imshow(hillshade, cmap='gray', extent=[lon.min(), lon.max(), lat.min(), lat.max()])

    # plot elevation shading
    cfill = plt.pcolormesh(lon, lat, elev, cmap=cmap1, transform=ccrs.PlateCarree(), alpha=0.55)
    add_colorbar(ax, cfill, 'Elevation (m)')  # Adds external colorbar

    # plot contour lines if contour=True
    if contour:
        contour_levels = np.arange(math.floor(elev.min() / 100) * 100, math.ceil(elev.max() / 100) * 100 + 100, 100)
        ax.contour(lon, lat, elev, levels=contour_levels, colors='black', transform=ccrs.PlateCarree(), linewidths=linewidths)
        contour_levels = np.arange(math.floor(elev.min() / 1000) * 1000, math.ceil(elev.max() / 1000) * 1000 + 1000, 1000)
        contour = ax.contour(lon, lat, elev, levels=contour_levels, colors='black', transform=ccrs.PlateCarree(), linewidths=linewidths*2)
        ax.clabel(contour, levels=contour_levels, inline=True, fontsize=4)

    # create gridlines if grid=True
    if grid:
        longitudes = np.arange(math.floor(bounding_box[0] * 10) / 10, math.ceil(bounding_box[1] * 10) / 10 + 0.1, 0.1)
        latitudes = np.arange(math.floor(bounding_box[2] * 10) / 10, math.ceil(bounding_box[3] * 10) / 10 + 0.1, 0.1)

        gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.4, color="gray", alpha=0.6)
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
def get_slope(subset_elevation, resolution):
    """
    Calculate slope (degrees) at each pixel within the digital elevation model.  
    Args:
        subset_elevation: array. REQUIRED. 2D numpy array of elevations (meters).
        resolution: int or float. REQUIRED. Grid spacing in meters (DEM resolution).
    Returns:
        Numpy array of slope (degrees) for each pixel. 
    """ 
    dz_dx, dz_dy = np.gradient(subset_elevation,resolution) # calculate slope at each point
    slope_radians = np.sqrt(dz_dx**2 + dz_dy**2) # convert to radians
    slope = np.degrees(np.arctan(slope_radians)) # convert to degrees

    return slope

########### ---------------------next function------------------------ ###########

def get_rank_curvature(terrain_data, subset_elevation, margin_pixels, resolution, r):
    """
    Calculates rank and curvature at each pixel within the digital elevation model following Lundquist et al., 2008. 
    Args:
        terrain_data: dict. REQUIRED. A dictionary with key 'elevation' containing 2D numpy array of DEM elevations. 
        subset_elevation: array. REQUIRED. A 2D numpy array of elevation with margin pixels removed.
        margin_pixels: int. REQUIRED. The number of pixels to ignore at the edge of the DEM.
        resolution: float or int. REQUIRED. Grid spacing of the DEM (meters). 
        r: float or int. REQUIRED. Critical radius used in rank and curvature calculations (meters). 
    Returns:
        Tuple of 2D numpy arrays (rank, curvature).
    """
    elevation = terrain_data['elevation']

    rank = np.full_like(subset_elevation, fill_value=np.nan, dtype=float) # initialize array of nans
    curvature = np.full_like(subset_elevation, fill_value=np.nan, dtype=float) # initialize array of nans

    start_row, start_col = margin_pixels, margin_pixels # start iterations at margin pixels
    rows, cols = subset_elevation.shape # get number of rows/columns in DEM

    # iterate through each pixel, calculating rank and curvature
    for subset_row in range(rows):
        for subset_col in range(cols):

            row = start_row + subset_row
            col = start_col + subset_col

            hop = round(r/resolution) # hop distance in pixels

            square_matrix = elevation[row-hop:row+hop,col-hop:col+hop] # matrix of surrounding pixels within hop distance

            total_count = square_matrix.size-1
            less_count = np.count_nonzero(square_matrix<elevation[row,col])
            rank[subset_row,subset_col] = less_count/total_count 

            # curvature term a
            curvature_term_a = (2*elevation[row,col]-
                                (elevation[row,col-hop]+elevation[row,col+hop]+
                                elevation[row+hop,col]+elevation[row-hop,col])/2)
                
            # curvature terb b
            curvature_term_b = (2*elevation[row,col]-
                                (elevation[row+hop,col-hop]+elevation[row+hop,col+hop]+
                                elevation[row-hop,col-hop]+elevation[row+hop,col+hop])/2)
            curvature[subset_row,subset_col] = 0.25*((2*r)**(-1)*curvature_term_a+((2*r*np.sqrt(r))**(-1))*curvature_term_b)

    return rank, curvature

########### ---------------------next function------------------------ ###########

def classify_CAP(subset_elevation, slope, curvature, rank):
    """
    Classifies cold-air pooling signal at pixel in the digital elevation model. 
    Args:
        subset_elevation: REQUIRED. 2D numpy array of elevation with margin pixels removed.
        slope: array. REQUIRED. 2D numpy array of slope.
        curvature: array. REQUIRED. 2D numpy array of curvature.
        rank: array. REQUIRED. 2D numpy array of rank. 
    Returns:
        2-dimensional array. Array of cold-air pooling signal for each pixel in the digital elevation model.
            1 = CAP
            0 = no signal
           -1 = no CAP
    """
    CAP = np.zeros(subset_elevation.shape) # initialize array for CAP values
    percentile = -(1/60)*slope+0.5 # calculate perventile threshold, given by Lundquist et al., 2008
    CAP[rank<percentile] = 1 # CAP (CAP = 1) when rank is less than percentile
    CAP[slope>30] = -1 # No CAP (CAP = -1) when slope is greater than 30 degrees
    CAP[curvature>0] = -1 # or No CAP (CAP = -1) when curvature is positive

    return CAP
########### ---------------------next function------------------------ ###########
def to_ds(terrain_stats,rvals):
    """
    Convert numpy array to xarray.Dataset.
    Args:
        terrain_stats: dict. A dictionary with keys 'elevation', 'lon', 'lat', 'slope', 'rank', 'curvature', and 'CAP'.
        rvals: array or list. REQUIRED. An array or list of critical radius values (meters).
    Returns:
        xarray.Dataset containing elevation, slope, rank, curvature, and CAP variables with coordinates lat, lon. 
    """
    ds = xr.Dataset(
        data_vars=dict(
            elevation=(["lat", "lon"], terrain_stats['elevation']),
            slope=(["lat", "lon"], terrain_stats['slope']),
            CAP=(["lat", "lon", "xv_dist"], terrain_stats['CAP']),
            rank=(["lat", "lon", "xv_dist"], terrain_stats['rank']),
            curvature=(["lat", "lon", "xv_dist"], terrain_stats['curvature'])
        ),
        coords=dict(
            lon=("lon", terrain_stats['lon'][0]),  
            lat=("lat", terrain_stats['lat'][:,0]), 
            xv_dist=("xv_dist",np.multiply(rvals,2))
        ),
        attrs=dict(description="Terrain stats."),
    )
    ds.rio.set_spatial_dims(x_dim='lon', y_dim='lat')

    return ds

########### ---------------------next function------------------------ ###########

def cap_analysis(terrain_data, resolution, xv_dist):
    """
    Integrates all necessary calculations to determine CAP signal at each pixel across the DEM. 
    Args:
        terrain_data: xarray.Dataset. REQUIRED. xarray.Dataset with variable 'elevation' with coordinates 'lat' and 'lon'. 
        resolution: int or float. REQUIRED. Grid spacing of DEM (meters).
        xv_dist: list. REQUIRED. list or array of cross-valley distances (meters).
    Returns:
        xarray.Dataset with variables elevation, slope, rank, curvature, and CAP across lat/lon for each xv_dist value. Also writes dataset to .nc file. 
    """
    if xv_dist is None:
        raise ValueError("the xv_dist arg is required. Please input a list or array of values with at least one element.")
    
    if len(xv_dist) == 0:
        raise ValueError("xv_dist cannot be empty. Must be a list or array of values with at least one element.")

    rvals = np.divide(xv_dist,2) # calculate critical radius from cross-valley distance
    rvals = sorted(rvals, reverse=True) # sort from largest to smallest
    
    lon_1d = terrain_data.lon.to_numpy()            
    lat_1d = terrain_data.lat.to_numpy()           
    lon_2d, lat_2d = np.meshgrid(lon_1d, lat_1d)
    crs = terrain_data.rio.crs
    terrain_data = {'elevation':terrain_data.elevation.to_numpy(),'lon':lon_2d,'lat':lat_2d}

    CAP_list = []
    rank_list = []
    curvature_list = []

    # SUBSET THE DEM
    margin_pixels = int(math.ceil(rvals[0]/resolution))
    subset_elevation = terrain_data['elevation'][margin_pixels:-margin_pixels, margin_pixels:-margin_pixels]
    subset_lon =  terrain_data['lon'][margin_pixels:-margin_pixels, margin_pixels:-margin_pixels]
    subset_lat =  terrain_data['lat'][margin_pixels:-margin_pixels, margin_pixels:-margin_pixels]

    # CALCULATE SLOPE
    slope = get_slope(subset_elevation=subset_elevation,resolution=resolution)

    if len(rvals) > 1: 
        print('Iterating through critical values...')
        sys.stdout.write("a progress indicator will appear shortly")  
        sys.stdout.flush()
        for index,r in enumerate(rvals):
            rank, curvature = get_rank_curvature(terrain_data=terrain_data, subset_elevation=subset_elevation, margin_pixels=margin_pixels, resolution=resolution, r=r)
            CAP = classify_CAP(subset_elevation=subset_elevation, slope=slope, curvature=curvature, rank=rank)

            rank_list.append(rank)
            curvature_list.append(curvature)
            CAP_list.append(CAP)

            percent_complete = (index / len(rvals)) * 100  
            sys.stdout.write(f"\r{percent_complete:.1f}% complete")  
            sys.stdout.flush()  

    elif len(rvals) == 1:
        print('Calculating for a single critical value...')
        rank, curvature = get_rank_curvature(terrain_data=terrain_data, subset_elevation=subset_elevation, margin_pixels=margin_pixels, resolution=resolution, r=rvals[0])
        CAP = classify_CAP(subset_elevation=subset_elevation, slope=slope, curvature=curvature, rank=rank)

        rank_list.append(rank)
        curvature_list.append(curvature)
        CAP_list.append(CAP)        

    rank = np.stack(rank_list,axis=-1)
    curvature = np.stack(curvature_list,axis=-1)
    CAP = np.stack(CAP_list,axis=-1)

    terrain_stats = {'elevation':subset_elevation,'lon':subset_lon,'lat':subset_lat,'slope':slope,'rank':rank,'curvature':curvature,'CAP':CAP}
    
    percent_complete = 100
    sys.stdout.write(f"\r{percent_complete:.1f}% complete \n")  
    sys.stdout.flush()  
    
    ds = to_ds(terrain_stats=terrain_stats,rvals=rvals)
    ds.rio.write_crs(crs,inplace=True)

    if len(rvals)==1:
        filepath = f"output/terrain_stats_{int(rvals[0]*2)}.nc"
    else:
        filepath = f"output/terrain_stats_{int(rvals[0]*2)}_to_{int(rvals[-1]*2)}.nc"
    ds.to_netcdf(filepath)
    print(filepath)

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
        cbar.ax.tick_params(labelsize=7)
    cbar.set_label(label)
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')

########### ---------------------next function------------------------ ###########

def plotcap(terrain_stats,xv_dist=None,bounding_box=None,depth_plot=False,shade=False,contour=False,linewidths=0.05,grid=True,output='show',dpi=300):
    """
    Create contour map of the digital elevation model and shade topography prone to CAP. 
    Args:
        terrain_stats: xarray.Dataset. REQUIRED. xarray.Dataset containing elevation and CAP variables with coordinates lat, lon.
        xv_dist: coordinate value. Default = None. REQUIRED if depth_plot=False. The cross-valley distance (meters) for which to plot CAP pooling signal.
        bounding_box: list. Default = None. Lat/lon bounding box in the form [xmin, xmax, ymin, ymax]. Use WGS84.
        depth_plot: bool. Default = False. Option to shade regions based on the number of iterations for which each pixel was classified as CAP.
        shade: bool. Default = False. Option to shade regions with no signal.
        contour: bool. Default = False. Option to overlay elevation contour lines.
        linewidths: float. Default = 0.05. Linewidth for contour lines. Also scales CAP region contours. 
        output: str. Default = 'show'. Choose from 'save' or 'show'. Option to output the figure or save it to PNG.
        dpi: int or float. If output='save', dpi is the resolution of the saved figure.
    Returns:
        None. Displays plot or saves .png to examples/figures/.
    """
    elev = terrain_stats['elevation'].values
    lon = terrain_stats['lon']
    lat = terrain_stats['lat']
    
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()},figsize=(6.5,5))

    # establish colormap
    hex_color1 = "#e8e7d3"
    hex_color2 = "#374157"
    cmap1 = LinearSegmentedColormap.from_list("custom_cmap", [hex_color1, hex_color2])

    extent = [lon.min(), lon.max(), lat.min(), lat.max()]

    # set extent based on bounding box
    if bounding_box:
        ax.set_extent(bounding_box, crs=ccrs.PlateCarree())
    else:
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        bounding_box = extent

    # generate hillshade
    ls = LightSource(azdeg=335, altdeg=30)
    hillshade = ls.hillshade(elev, vert_exag=.04, dx=1, dy=1)
    ax.imshow(hillshade, cmap='gray', extent=[lon.min(), lon.max(), lat.min(), lat.max()])

    # plot elevation shading
    cfill = plt.pcolormesh(lon, lat, elev, cmap=cmap1, transform=ccrs.PlateCarree(), alpha=0.5)
    add_colorbar(ax, cfill, 'Elevation (m)')  # Adds external colorbar

    # plot contour lines if contour=True
    if contour:
        contour_levels = np.arange(math.floor(elev.min() / 100) * 100, math.ceil(elev.max() / 100) * 100 + 100, 100)
        ax.contour(lon, lat, elev, levels=contour_levels, colors='black', transform=ccrs.PlateCarree(), linewidths=linewidths)
        contour_levels = np.arange(math.floor(elev.min() / 1000) * 1000, math.ceil(elev.max() / 1000) * 1000 + 1000, 1000)
        contour = ax.contour(lon, lat, elev, levels=contour_levels, colors='black', transform=ccrs.PlateCarree(), linewidths=linewidths*2)
        ax.clabel(contour, levels=contour_levels, inline=True, fontsize=4)
    
    # if depth_plot=True, shade regions based on the number of iterations for which each pixel was classified as CAP
    if depth_plot:
        cmap2 = plt.get_cmap('Blues')

        num_layers = len(terrain_stats['xv_dist'])
        mask = np.ones_like(terrain_stats['CAP'].isel(xv_dist=0),dtype='bool')  # Start with a mask of ones
        CAP_masked = np.zeros_like(terrain_stats['CAP'].isel(xv_dist=0))  

        for i in range(num_layers):
            CAP = terrain_stats['CAP'].isel(xv_dist=i)
            mask &= (CAP == 1)  # Retain only where current and previous layers at least 1
            CAP_masked = CAP_masked + np.where(mask, 1, 0)
        CAP_contours = CAP_masked.copy()
        CAP_masked[CAP_masked==0] = np.nan

        cfill = ax.pcolormesh(lon, lat, CAP_masked, cmap=cmap2, transform=ccrs.PlateCarree(), alpha=0.75)
        ax.contour(lon, lat, CAP_contours, levels=[CAP_contours.min()],colors='black', transform=ccrs.PlateCarree(),linewidths=linewidths*4)
        add_colorbar(ax, cfill, '# of iterations',orientation='horizontal',position='top',ticks=np.arange(num_layers))  # Adds external colorbar

    else: 
        if xv_dist is None:
            raise ValueError("the xv_dist arg is required if depth_plot=False.")
        cap_mask = np.where(terrain_stats['CAP'].sel(xv_dist=xv_dist) == 1, 1, np.nan)
        no_signal_mask = np.where(terrain_stats['CAP'].sel(xv_dist=xv_dist) == 0, 1, np.nan)

        ax.contour(lon,lat,terrain_stats['CAP'].sel(xv_dist=xv_dist), levels=[0,1,2],colors='black', transform=ccrs.PlateCarree(), linewidths=0.07)
        legend_patch = Patch(facecolor='skyblue', label='Cold-air pool', alpha=.6)
        ax.legend(handles=[legend_patch], loc='upper right')

        cfill = ax.pcolormesh(
            lon, lat, cap_mask,
            cmap=ListedColormap(['skyblue']),
            transform=ccrs.PlateCarree(),
            alpha=.6,
            shading='auto',
            label='Cold-Air Pool'
        )    

        # if shade=True, shade regions with no signal
        if shade:
            ax.contourf(terrain_stats.lon,terrain_stats.lat,no_signal_mask,hatches=[6*'.'],colors='none',transform=ccrs.PlateCarree(),)

    # create gridlines if grid=True
    if grid:
        longitudes = np.arange(math.floor(bounding_box[0] * 10) / 10, math.ceil(bounding_box[1] * 10) / 10 + 0.1, 0.1)
        latitudes = np.arange(math.floor(bounding_box[2] * 10) / 10, math.ceil(bounding_box[3] * 10) / 10 + 0.1, 0.1)

        gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.4, color="gray", alpha=0.6)
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
        plt.savefig('figures/cap_map.png', dpi=dpi, bbox_inches='tight')
        print('saved figure to figures/cap_map.png')