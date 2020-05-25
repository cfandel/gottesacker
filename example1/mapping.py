'''Mapping module for GemPy.
Chloé Fandel 2019.
Functions to visualize and export GemPy model outputs.
2D maps, custom diagonal cross-sections, export to GSLIB and VTK, import from GSLIB.
Based on original code from Elisa Heim.

See accompanying notebook for examples of how to use.'''

#Imports:
import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import gdal
import pyevtk
from copy import copy

#sys.path.append("../../..")   #optional: if gempy has been downloaded from GitHub rather than installed normally, look for it in the folders above the current folder
import gempy as gp


#############################################################################
def importDEM(filename, show=True):
    '''Import DEM from a tif file using gdal package.
    Return a dem object, and xyz extent and resolution.
    (this can be used to set the model extent)
    NOTE: vertical (z) resolution can't be extracted from the raster!
    
    filename: string indicating the filename (must be a rectangular tif)
    show:     option to show a plot of the DEM or not.
    
    Returns:     grid_info = [xmin,xmax,xres,dx,ymin,ymax,yres,dy,zmin,zmax,zres,dz]  

    dem:      gdal dem objecct
    dema:     array of elevation values of dim: xres,yres
    xmin:     minimum x value (same for ymin, zmin)
    xmax:     maximum x value (same for ymax, zmax)
    xres:     x resolution, aka number of columns, aka number of cells along x axis (NOT pixel width)
    dx:       pixel width in x direction 
    etc.
    '''
    
    dem = gdal.Open(filename)    #DEM must be rectangular tif 
    dema = dem.ReadAsArray()     #copy of DEM as a numpy array (defaults to integers)
    dema = dema.astype(float)    #convert integer array to float array
    dema[dema==0] = np.nan       #replace zeros with NaNs (have to convert array to float first)

    ulx, pixelwidthx, xskew, uly, yskew, pixelheighty = dem.GetGeoTransform() #get resolution and coordinate info (for some reason the order of skew and pixel size is flipped for y axis?!)
    ncol = dem.RasterXSize            #number of columns (aka number of cells along x axis)
    nrow = dem.RasterYSize            #number of rows (aka number of cells along y axis)
    lrx = ulx + (ncol * pixelwidthx)  #lower right x coord = upper left x coord + (width of raster cells in x direction * number of raster cells in x direction)
    lry = uly + (nrow * pixelheighty)

    #Get min and max elevations (z):
    #note: gdal's built-in GetRasterBand and GetStatistics return an incorrect zmin (WHY?!)
    zmin = np.nanmin(dema)
    zmax = np.nanmax(dema)
    
    #Assign useful names:
    xmin = ulx
    xmax = lrx
    xres = ncol
    dx =   abs(pixelwidthx)
    ymin = lry
    ymax = uly
    dy =   abs(pixelheighty)
    yres = nrow
    zres = 'na'     #can't be extracted from raster

    #Print results & display raster:
    if show==True:
        print('Raster dimensions: \nxmin: {:<12} xmax: {:<12} xres: {} \nymin: {:<12} ymax: {:<12} yres: {} \nzmin: {:<12} zmax: {:<12} zres: {}'.format(
            xmin,xmax,xres,ymin,ymax,yres,zmin,zmax,zres))
        plt.imshow(dema, extent=(xmin,xmax,ymin,ymax), vmin=zmin, vmax=zmax) #plot raster as image
        #print(gdal.Info(dem))  #for more detailed file info, uncomment this line
        
    return dem,dema,xmin,xmax,xres,dx,ymin,ymax,yres,dy,zmin,zmax,zres




#############################################################################
def crop2topo(sol, dem, a=None):
    '''Crop gempy lith_block to land surface using an imported DEM raster file.
    
    Inputs:
    sol:    gempy solutions object containing lith_block & grid, OR, list [xmin,xmax,xres, ymin,ymax,yres, zmin,zmax,zres]
    dem:    path to DEM raster file, OR array of elevation values, of same xy dimensions as model, to be used for cropping.
    a:      optional: array to be cropped if not using gempy (e.g. as provided by mapping.crop2raster())
            must already be of dimensions (xres,yres,zres)
              
    Output:
    gcrop:    cropped lith_block array'''

    #Get coordinate info from grid & create VTK cells info in coord system instead of in cell indices:
    try:
        xmin = sol.grid.regular_grid.extent[0]        #min (left) coordinate
        xmax = sol.grid.regular_grid.extent[1]        #max (right) coordinate
        xres = sol.grid.regular_grid.resolution[0]    #number of pixels
        ymin = sol.grid.regular_grid.extent[2]
        ymax = sol.grid.regular_grid.extent[3]
        yres = sol.grid.regular_grid.resolution[1]
        zmin = sol.grid.regular_grid.extent[4]      #important: need zmin of model, NOT just DEM, since model may extend below lowest point at land surface
        zmax = sol.grid.regular_grid.extent[5]
        zres = sol.grid.regular_grid.resolution[2]
        
    except:
        xmin = sol[0]   #min (left) coordinate
        xmax = sol[1]   #max (right) coordinate
        xres = sol[2]   #number of pixels
        ymin = sol[3]
        ymax = sol[4]
        yres = sol[5]
        zmin = sol[6]   #important: need zmin of model, NOT just DEM: model may extend below lowest point at land surface
        zmax = sol[7]
        zres = sol[8]
        
    dx   = (xmax-xmin)/xres                       #pixel width
    xvals = np.arange(xmin,xmax+dx,dx)            #coordinates of column edges (borders between pixels)
    dy   = (ymax-ymin)/yres
    yvals = np.arange(ymin,ymax+dy,dy)
    dz   = (zmax-zmin)/zres
    zvals = np.arange(zmin,zmax+dz,dz)

    #Get DEM & convert to topo array:
    try:
        dem = gdal.Open(dem_path)    #DEM must be rectangular tif of same xy resolution as model grid
        dema = dem.ReadAsArray()     #copy of DEM as a numpy array (defaults to integers)
        dema = dema.astype(float)    #convert integer array to float array
        dema[dema==0] = np.nan       #replace zeros with NaNs (must first convert array to float)        
        t = dema.copy()              #make a copy of elevation values directly from dem (do not use topo grid - different indexing!)
        t = np.rot90(t)              #rotate & flip to be same orientation as lith block
        t = np.fliplr(t)                 
        t = np.flipud(t)
    except:
        dema = dem
        dema = dema.astype(float)    #convert integer array to float array
        dema[dema==0] = np.nan       #replace zeros with NaNs (must first convert array to float)        
        t = dema.copy()              #make a copy of elevation values directly from dem (do not use topo grid - different indexing!)
        #t = np.rot90(t)              #rotate & flip to be same orientation as lith block
        #t = np.fliplr(t)                 
        #t = np.flipud(t)

    #Get lith array:
    try:
        g = a.copy()                         #if array already provided, use array (make copy) (must already be 3D)
    except:                                  #otherwise use lith_block
        g = sol.lith_block.copy()            #make a copy of lith block to avoid messing up original
        g = np.reshape(g, (xres,yres,zres))  #reshape lith block to 3D
    g = g.round()                            #round to integers

    #Check dimensions for debugging:
    #print('topo shape:\t', t.shape)
    #print('lith shape:\t', g.shape)
    #print('topo min:\t',   t.min(), np.unravel_index(t.argmin(), t.shape)) #show indices of min/max of topography
    #print('topo max:\t',   t.max(), np.unravel_index(t.argmax(), t.shape))
    #print('model min:\t',  zmin) #show min/max of entire model
    #print('model max:\t',  zmax)

    #Get z indices of land surface:
    #ind = sol.grid.topography._find_indices().copy()  #this has the wrong dimensions! do not use! it is stretched to rotate 90 degrees. instead:
    ind = (t - zmin) / dz              #calculate the cell index of each point in the dem array using the cell height (i.e. how many cells/layers up it is from the base)
    #ind = ind - 1                      #bump down by one to represent bottom of cell (this cuts lower than actual topo) 
    ind[ind==-1] = 0                   #zeros should stay zeros and not go negative
    ind = np.ceil(ind)                 #round up to nearest integer 
    ind = ind.astype(int)              #convert to integers for use as vertical indices
    #print('ind:\t', ind.shape, type(ind)) #debugging
    #print(np.unique(ind))

    #Crop off everything above the land surface:
    m = np.zeros((xres,yres))  #create array of zeroes of shape (xres,yres)
    gcrop = g.copy()           #make a copy bc need original as reference
    for x in range(xres):      #loop over x,y indices
        for y in range(yres):
            try:
                z = ind[x,y]            #get land surface elev at point (x,y)
                #m[x,y] = g[x,y,z]       #get lith value at land surface at point (x,y) - only if want to return a map
                gcrop[x,y,z:] = np.nan  #convert all lith values above land surface at point (x,y) to nans
            except:
                z = ind[y,x]
                gcrop[y,x,z:] = np.nan  #convert all lith values above land surface at point (x,y) to nans
    #print('map:\t', m.shape)       #for debugging
    #print('crop:\t', gcrop.shape)
    
    return gcrop



#############################################################################
def crop2raster(sol, mask_path, a=None, nanval=0):
    '''Crop the extent of the geologic model to the extent of an irregularly-shaped imported raster 
    with specified nan value indicating empty cells.
    
    sol:        gempy solution object containing lith_block and grids
    mask_path:  path to the raster file to use for cropping (can be bigger but not smaller than model extent)
    a:          optional: numpy array containing the gempy lithology to be cropped instead of sol.lith_block. 
                can supply an already-cropped array from mapping.crop2topo()
                must already have dimensions (xres,yres,zres)
    nanval:     value being used in raster file to indicate empty cells 
    
    returns:
    gcrop:      array of same dimensions as input, but with empty cells filled with np.nan'''
    
    #Get grid info
    xres = sol.grid.regular_grid.resolution[0]       #get x resolution
    yres = sol.grid.regular_grid.resolution[1]       #get y resolution
    if len(sol.grid.regular_grid.resolution)==3:     #if 3D grid
        zres = sol.grid.regular_grid.resolution[2]   #get z resolution
       
    #Get array to export:
    try:                                     #if array already provided, use array
        g = a.copy()                         #make a copy to avoid messing up original
    except:                                  #otherwise, use lith_block
        g = sol.lith_block.copy()            #make a copy to avoid messing up original
        g = np.reshape(g, (xres,yres,zres))  #reshape block to 3D
    g = g.round()                            #round to integers
        
    #Import & format raster to crop to:
    mask = gdal.Open(mask_path)     #import raster file
    mask = mask.ReadAsArray()       #read file as array (output will have an array for each color channel)
    mask = mask.astype(float)       #convert integer array to float array
    mask[mask==nanval] = np.nan     #replace zeros with NaNs (have to convert array to float first)
    if mask.shape != g.shape:
        #print('mask shape', mask.shape, 'does not match array shape', g.shape)
        mask = np.rot90(mask)       #rotate & flip to be same orientation as lith block
        mask = np.fliplr(mask)                 
        mask = np.flipud(mask)
        mask = mask[-xres:,0:yres]  #if raster larger than model xy extent, trim to size
        #plt.imshow(g[:,:,0])
        #plt.imshow(mask)
    
    #Crop input array:
    if g.shape == mask.shape:           #if arrays have same extent
        g[np.isnan(mask)] = np.nan      #crop lith to active cells in imported raster
    else:                                   #otherwise, assume lith is 3D
        g[np.isnan(mask),:] = np.nan    #crop over all z cells
        #acrop = np.rot90(acrop)             #rotate 90 degrees to give dimensions of (xres,yres,zres)
    
    return g



########################################################################################################
def export2gslib(a, filename, grid):
    '''Exports a numpy array to a gslib file (GeoModeller-style format), using gempy grid objects to get the correct dimensions.
    
    Inputs: 
    a:        numpy array to be exported 
    filename: filename or path to save to
    grid:     gempy grid object to use for the dimensions of the array, obtained from geo_model.grid.gridnamehere
              OR a list [xres,yres] or [xres,yres,zres]
              A list of available grid names can be found with geo_model.grid.grid_types. 
              The active grids in this list can be seen with geo_model.grid.active_grids
    
    Output: 
    filename.gslib: a gslib file
    '''
    
    #Format data:
    a[np.isnan(a)] = 0              #assign zeros to nan values (for SKS)
    a = np.round(a)                 #round to integers
    a = a.astype(int)               #convert from floats to integers (so that gslib file will have integers)
    
    #Get grid info:
    try:
        xres = grid.resolution[0]                          #get x resolution
        yres = grid.resolution[1]                          #get y resolution
        if len(grid.resolution)==3:                        #if 3D grid
            zres = grid.resolution[2]                      #get z resolution
    except:
        xres = grid[0]                          #get x resolution
        yres = grid[1]                          #get y resolution
        if len(grid)==3:                        #if 3D grid
            zres = grid[2]                      #get z resolution
       
    #Format array shape:
    try:                                                #try to reshape in 3D
        a = np.reshape(a,(xres,yres,zres))              #reshape 1D array to a 3D array with correct dimensions
        a = np.reshape(a, xres*yres*zres, order='F')    #reshape 3D array back to 1D array using Fortran indexing 
    except:                                             #if 3D doesn't work, try 2D 
        a = np.reshape(a,(xres,yres))                   #reshape 1D array to a 2D array with correct dimensions
        a = np.reshape(a, xres*yres, order='F')         #reshape 2D array back to 1D array using Fortran indexing 
        
    #Export:    
    df = pd.DataFrame(a)                            #store array in a pandas dataframe
    header = pd.DataFrame(['GemPy output',1,'lith'])    #set gslib file header
    df = header.append(df)                              #attach header and data
    df.to_csv(filename, header=False, index=False)      #write a text file in gslib format


#############################################################################
def importgslib(filename, grid):
    '''Imports a gslib file (GeoModeller-style) into a numpy array with dimensions taken from a GemPy grid object.
    
    Inputs:
    filename:  the name or path to the file to be imported
    grid:     gempy grid object to use for the dimensions of the array, obtained from geo_model.grid.gridnamehere
              OR a list of dimensions [xres,yres] or [xres,yres,zres]
              A list of available grid names can be found with geo_model.grid.grid_types. 
              The active grids in this list can be seen with geo_model.grid.active_grids
              
    Returns:
    a:        2D or 3D numpy array representing the gslib data, with dimensions taken from the grid'''
    
    #Get grid info:
    try:
        xres = grid.resolution[0]                          #get x resolution
        yres = grid.resolution[1]                          #get y resolution
        if len(grid.resolution)==3:                        #if 3D grid
            zres = grid.resolution[2]                      #get z resolution
    except:
        xres = grid[0]                          #get x resolution
        yres = grid[1]                          #get y resolution
        if len(grid)==3:                        #if 3D grid
            zres = grid[2]                      #get z resolution
    
    #Get data:
    a = pd.read_csv(filename, skiprows=2, dtype=float)  #read in gslib file to pandas df without header rows, as floats
    a = a.values                                        #get an array of the values
    a[a==0] = np.nan                                    #replace zeros with NaNs (array must be float first)

    #Reshape based on specified grid:
    try:                                             #for 3D files (blocks)
        a = np.reshape(a,(xres,yres,zres),order='F') #reshape to xyz grid using Fortran ordering
        a = np.rot90(a)                              #rotate 90 degrees CCW
    except:                                          #for 2D files (maps)
        a = np.reshape(a,(xres,yres),order='F')      #reshape to xyz grid using Fortran ordering
        #a = np.rot90(a)                              #rotate 90 degrees CCW
    
    return a



#############################################################################
def export2vtk(sol, vtk_path, a=None):
    '''Export gempy lith_block array to VTK file (for viewing in e.g. ParaView).
    
    Inputs:
    sol:      gempy solutions object containing lith_block and grids 
              OR list of model dimensions [xmin,xmax,xres, ymin,ymax,yres, zmin,zmax,zres]
    vtk_path: filepath to save VTK file to (must not include any file extension)
    a:        optional: if desired, specify a different array to export (such as a cropped lith_block)
    
    Returns:
    VTK file representing the exported array'''
    
    #Get coordinate info from grid & create VTK cells info:
    try:
        xmin = sol.grid.regular_grid.extent[0]        #min coordinate value (left)
        xmax = sol.grid.regular_grid.extent[1]        #max coordinate value (right)
        xres = sol.grid.regular_grid.resolution[0]    #number of pixels
        ymin = sol.grid.regular_grid.extent[2]
        ymax = sol.grid.regular_grid.extent[3]
        yres = sol.grid.regular_grid.resolution[1]
        zmin = sol.grid.regular_grid.extent[4]
        zmax = sol.grid.regular_grid.extent[5]
        zres = sol.grid.regular_grid.resolution[2]
    except:
        xmin = sol[0]        #min coordinate value (left)
        xmax = sol[1]        #max coordinate value (right)
        xres = sol[2]        #number of pixels
        ymin = sol[3]
        ymax = sol[4]
        yres = sol[5]
        zmin = sol[6]
        zmax = sol[7]
        zres = sol[8]
        
    dx   = (xmax-xmin)/xres                       #pixel width
    xvals = np.arange(xmin,xmax+dx,dx)            #calculate x coordinate values of the boundaries between cells
    dy   = (ymax-ymin)/yres
    yvals = np.arange(ymin,ymax+dy,dy)
    dz   = (zmax-zmin)/zres
    zvals = np.arange(zmin,zmax+dz,dz)

    #Format array for export:
    try:                                     #if array already provided, use array
        g = a.copy()                         #make a copy to avoid messing up original
    except:                                  #otherwise, use lith_block
        g = sol.lith_block.copy()            #make a copy to avoid messing up original
        g = g.round()                        #round to integers
        g = np.reshape(g, (xres,yres,zres))  #reshape block to 3D
    
    #Debugging checks:
    #print('x:', xmin,xmax,xres,dx)
    #print('y:', ymin,ymax,yres,dy)
    #print('z:', zmin,zmax,zres,dz)
    #print('shape of array to export:', g.shape)
    #plt.imshow(g[:,:,0])
    
    pyevtk.hl.gridToVTK(vtk_path, xvals, yvals, zvals, cellData={'data': g}) #export to VTK
    
    
#############################################################################
def plot_map(a, geo_model, plot_data=False, ref_points=None, ax=None):
    '''Plot 2D geologic map (generated by cropping lith_block with topography).
    Inputs:
    a:          Data representing the geologic map. Can be a 1D or 2D numpy array, or the path to a gslib file.
                To get 1D array: a = sol.geological_map. To get 2D array: a = np.reshape(sol.geological_map, (yres,xres)) (may need to also flip up-down)
                To import a gslib file: a = mapping.importgslib(filename, grid)
    geo_model:  gempy model object containing grid and solutions
    plot_data:  Whether or not to display the input data on the map. Defaults to false.
    ref_points: Optional path to a csv file of reference points to be added to the map (could be springs, peaks, towns, etc.).
                CSV must have the following labeled columns: Name, X, Y, Z
    ax:         Optional axis object to plot to. If not specified, new axes will be created.
                
    Returns:
    f,ax:       The figure and axis objects being plotted on.'''
    
    #Get grid info
    xres = geo_model.grid.topography.resolution[0]       #get x resolution from topo grid
    yres = geo_model.grid.topography.resolution[1]       #get y resolution
    extent = geo_model.grid.topography.extent            #get exent
        
    #Load data:
    try:
        a = np.reshape(a, (xres,yres))                        #if 1D, try reshaping to 2D
        print('reshaping 1D to 2D')
    except: 
        pass
    try:                                                          #if gslib file, try to load
        a = importgslib(a, geo_model.grid.topography)     #import geologic map from gslib file
        a = np.flipud(a)                                          #flip map up-down to account for gempy v2 bug
    except:                             
        pass

    #Generate colormap:
    cmap = matplotlib.colors.ListedColormap(geo_model.surfaces.colors.colordict.values())     #set colormap
    norm = matplotlib.colors.Normalize(vmin=1, vmax=len(geo_model.surfaces.colors.colordict)) #set normalization

    #Create axes:
    if ax==None:
        fig = plt.figure(figsize=(10,10))           #create empty figure
        ax = fig.add_subplot(111)                   #add subplot axes
        
    #Plot map:
    ax.imshow(a, extent=extent, cmap=cmap, norm=norm)              #plot map with coordinate extent and custom colormap

    #Plot data:
    if plot_data==True:
        pts = geo_model.orientations.df                                                                             #get orientations data points
        ax.scatter(pts.X,pts.Y, s=20, edgecolors='k', linewidths=1, c=pts.id, cmap=cmap, norm=norm, marker='<')    #plot as points with custom colormap
        pts = geo_model.surface_points.df                                                                           #get interface data points
        plt.scatter(pts.X,pts.Y, s=20, edgecolors='k', linewidths=1, c=pts.id, cmap=cmap, norm=norm, marker='o')    #plot as points with custom colormap
    if ref_points:
        pts = pd.read_csv(ref_points)                           #get reference points from file
        ax.scatter(pts.X,pts.Y, s=20, c='k')                   #plot as points 
        for i,pt in pts.iterrows():
            ax.annotate(pt.Name, (pt.X,pt.Y), fontsize=15)     #add labels to reference points

    #Add legend:
    patches = [matplotlib.patches.Patch(color=geo_model.surfaces.colors.colordict[key], label=key) for key in geo_model.surfaces.colors.colordict.keys()]  #loop to create a list of color patches & labels
    plt.legend(handles=patches, bbox_to_anchor=(1.35,1))  #add legend using color patches, at specific position



#############################################################################
def compare_geology(geo_model, geo_true, colordic=None, ref_points=None):
    
    '''Plots the GemPy-generated geologic map next to the true geologic map for comparison. 
    
    geo_model:    gempy model object containing grid and solutions
                  OR gempy solutions object containing geologic map
                  OR gslib file of lithologic unit values at the land surface, of dimensions (yres,xres)
    geo_true:     raster file of the actual geologic map, of same dimensions as geo_model map
    colordic:     optional dictionary of strings indicating which colors to map to which lithologic unit name (must be in order from youngest to oldest)
                  if not provided, uses gempy colormap. if unavailable, uses matplotlib colormap
    ref_points:   optional path to a csv file of reference points to be added to the map (could be springs, peaks, towns, etc.).
                  CSV must have the following labeled columns: Name, X, Y, Z
    '''
    
    #Import data & reformat:
    #Get true geologic map:
    trueRaster = gdal.Open(geo_true)   #raster must be rectangular tif with zeros as NaN values and formation values corresponding to GemPy fm values 
    true = trueRaster.ReadAsArray()    #convert raster to a numpy array (defaults to integers)
    true = true.astype(float)    #convert integer array to float array (to be able to use NaNs)
    true[true==0] = np.nan       #replace zeros with NaNs (have to convert array to float first)
    
    #Get grid info:
    try:     #from model object
        xres   = geo_model.grid.regular_grid.resolution[0]
        yres   = geo_model.grid.regular_grid.resolution[1]
        extent = geo_model.grid.regular_grid.extent[:4]
    except:  #from file
        #xres   = true.shape[0]
        #yres   = true.shape[1]
        xmin, dx, xskew, ymax, yskew, dx = trueRaster.GetGeoTransform() #get resolution and coordinate info (for some reason the order of skew and pixel size is flipped for y axis?!)
        xres = trueRaster.RasterXSize       #number of columns (aka number of cells along x axis)
        yres = trueRaster.RasterYSize       #number of rows (aka number of cells along y axis)
        xmax = xmin + (xres * abs(dx))      #lower right x coord = upper left x coord + (width of raster cells in x direction * number of raster cells in x direction)
        ymin = ymax + (yres * abs(dy))
        extent = [xmin, xmax, ymin, ymax]   #group into extent
        
    #Get modeled geologic map:
    try:        #from gempy geo_model 
        model = geo_model.solutions.geological_map              #get gempy geol map
        model = np.reshape(model, (yres, xres))                 #reshape to 2D
    except:  
        try:    #from gempy solutions
            model = geo_model.geological_map                        #get gempy geol map
            model = np.reshape(model, (yres, xres))                 #reshape to 2D
        except: #from file
            model = pd.read_csv(geo_model, skiprows=2, dtype=float) #read in gslib file to floats df, without header rows
            model = model.values                                    #get an array of the values
            model[model==0] = np.nan                                #replace zeros with NaNs (array must be float first)
            model = np.reshape(model,(yres,xres),order='F')         #reshape to xyz grid using Fortran ordering        
        model = np.flipud(model)                                    #flip to orient north

    #Plot:
    #Create figure:
    f = plt.figure(figsize=(20,20))                 #create figure
    
    #Generate colormap:
    try:    #set manually
        cmap = matplotlib.colors.ListedColormap(colordic.values())  #set colormap from dictionary provided
    except: 
        try: #set from gempy
            colordic = geo_model.surfaces.colors.colordict              #get gempy color dictionary
            cmap = matplotlib.colors.ListedColormap(colordic.values())  #set colormap from gempy
        except: #set from matplotlib
            cmap = 'viridis'                        #use one of matplotlib's default colormaps

    #Plot:
    plt.subplot(121)                                #create subplot in position 1 (1 row, 2 col, position 1)
    plt.imshow(model, extent=extent, cmap=cmap)     #plot model using colormap
    plt.title('model')
    if ref_points:
        pts = pd.read_csv(ref_points)                           #get reference points from file
        plt.scatter(pts.X,pts.Y, s=20, c='k')                   #plot as points 
        for i,pt in pts.iterrows():
            plt.annotate(pt.Name, (pt.X,pt.Y), fontsize=15)     #add labels to reference points
    plt.subplot(122)                                            #create 2nd subplot
    plt.imshow(true, extent=extent, cmap=cmap)                  #plot true geologic map
    plt.title('data')
    if ref_points:
        plt.scatter(pts.X,pts.Y, s=20, c='k')                   #plot as points 
        for i,pt in pts.iterrows():
            plt.annotate(pt.Name, (pt.X,pt.Y), fontsize=15)     #add labels

    #Add legend:
    try:   #if colordic provided
        patches = [matplotlib.patches.Patch(color=colordic[key], label=key) for key in colordic.keys()]  #loop to create a list of color patches & labels    
        plt.legend(handles=patches, bbox_to_anchor=(1.35,1))  #add legend using color patches, at specific position
    except: #otherwise don't plot legend
        pass
    return f


#############################################################################
def plotXsection(startpoints, endpoints, names, geo_model, lith, surflith, vscale=1, unitnames = None):
    '''Plots an approximate cross-section between the two specified points, using an elevation-cropped array (cells above the land surface should have nan values). Does not work well for N-S or nearly N-S xsections - use gempy built-in for that.
    startpoints: [[x1,y1],[x2,y2],...] float or [[col1,row1],[col1,row2],...] integer array of coordinates of starting points A
    endpoints:   [[x1,y1],[x2,y2],...] float or [[col1,row1],[col1,row2],...] integer array of coordinates of ending points B
    names:       [['A','B'],['C','D'],...] string array of names for starting and ending points
    grid_info:   [xmin,xmax,xres,dx,ymin,ymax,yres,dy,zmin,zmax,zres,dz] model grid info (can get these using importDEM() function if model grid is same as DEM grid)
    lith:        elevation-cropped array of lithologic unit indices of dimensions (nrow,ncol,nlay), i.e. (yres,xres,zres). Can use uncropped array, but will plot above land surface.
    surflith:    array of lithologic unit values at the land surface, of dimensions (yres,xres)
    vscale:      vertical exaggeration factor (y/x, defaults to 1)
    colors:      dictionary of color & unit names to use (in order from youngest to oldest) OR 'gempy' to use default gempy colormap
    '''
    
    #Get coordinate info from grid & create VTK cells info:
    xmin = geo_model.grid.regular_grid.extent[0]        #min coordinate value (left)
    xmax = geo_model.grid.regular_grid.extent[1]        #max coordinate value (right)
    xres = geo_model.grid.regular_grid.resolution[0]    #number of pixels
    dx   = (xmax-xmin)/xres                       #pixel width

    ymin = geo_model.grid.regular_grid.extent[2]
    ymax = geo_model.grid.regular_grid.extent[3]
    yres = geo_model.grid.regular_grid.resolution[1]
    dy   = (ymax-ymin)/yres

    zmin = geo_model.grid.regular_grid.extent[4]
    zmax = geo_model.grid.regular_grid.extent[5]
    zres = geo_model.grid.regular_grid.resolution[2]
    dz   = (zmax-zmin)/zres
    
    #Generate colormap:
    cmap = matplotlib.colors.ListedColormap(geo_model.surfaces.colors.colordict.values())     #set colormap
    norm = matplotlib.colors.Normalize(vmin=1, vmax=len(geo_model.surfaces.colors.colordict)) #set normalization
    
    #Plot geologic map once for reference:
    f1,ax1 = plt.subplots(1,1,figsize=(10,10))                  #create empty figure
    plt.imshow(surflith, cmap=cmap, norm=norm)                  #plot geology (normalized to gempy color range)
        
    f2,ax2 = plt.subplots(len(startpoints),1,figsize=(15,20))   #create figure and axes objects for subplots (one per xsection)
    
    for i in range(len(startpoints)):   #loop over number of sections
        #Get starting coordinates:
        xA = startpoints[i][0]   #get starting x coordinate
        yA = startpoints[i][1]   #get starting y coordinate
        xB = endpoints[i][0]     #get ending x coordinate
        yB = endpoints[i][1]     #get ending y coordinate

        #Calculate corresponding row,col
        if type(xA) != int:                     #if coordinates are NOT integers (i.e. not row,col numbers), convert them
            colA = (xA - xmin)//dx              #col:x calculate column index  c = (x1-x0)/dx 
            rowA = yres - ((yA - ymin)//dy)     #row:y calculate row index     r = ymax - (y1-y0)/dy
            colB = (xB - xmin)//dx                 
            rowB = yres - ((yB - ymin)//dy) 
        else:                                  #if coordinates are already in row,col format
            colA = xA
            rowA = yA
            colB = xB
            rowB = yB

        #Calculate line equation between points A and B:
        m = (rowB - rowA) / (colB - colA)   #calculate slope     m = (y2-y1)/(x2-x1)
        b = -m*colA + rowA                  #calculate intercept b = m*x1 + y1 (slope is neg here bc y axis is flipped)
        
        #Calculate true distance (not # of cells) between points A and B:
        #distance = ((xB-xA)**2 + (yB-yA)**2)**.5
        #xsizes.append(distance*.001)

        #Get xy indices for cells intersected by the x-sec line, then get z values for those xy points:
        xvals = np.arange(colA,colB)    #generate array of x values between the two points
        xvals = xvals.astype(int)       #convert to integer
        yvals = m*xvals + b             #calculate corresponding y values  y = mx + b 
        yvals = yvals.astype(int)       #convert to integers to be able to slice

        #xsec = lith[yvals,xvals,:].T    #select x-sec to plot and transpose to make it plot horizontally
        xsec = lith[xvals,yvals,:].T    #select x-sec to plot and transpose to make it plot horizontally
        
        #Plotting:
        #Add xsection lines to geologic map:
        plt.figure(f1.number)                           #make the map the active figure
        plt.plot([colA,colB],[rowA,rowB],'k')           #plot x-sec location line
        plt.annotate(names[i][0],xy=(colA,rowA),xytext=(colA-4,rowA+4)) #annotate start point
        plt.annotate(names[i][1],xy=(colB,rowB),xytext=(colB+1,rowB-1)) #annotate start point
        plt.ylim(bottom=yres, top=0) 
        plt.xlim(left=0, right=xres)

        #Plot cross-sections in a new figure:
        #Set and get correct subplot axes:
        if len(startpoints) == 1:               #check if there are more than 1 subplots (for indexing purposes)
            plt.sca(ax2)                        #make current subplot axes active (automatically makes fig active too)
            cax = plt.gca()                     #get current axes object
        else:          
            plt.sca(ax2[i])                     
            cax = plt.gca()
        cax.imshow(xsec, origin="lower", cmap=cmap, norm=norm)   #plot (with down=lower z indices)
        #cax.imshow(xsec, origin="lower", cmap=cmap)   #plot (with down=lower z indices)
        cax.set_aspect(vscale*dz/dx)                             #apply vertical exaggeration
        cax.set_ylim(bottom=0, top=zres)                         #set y limit to zres
        cax.set_title(names[i][0]+names[i][1])
        cax.set_anchor('W')                                      #align left (West)

        #Set ticks to accurately reflect elevation (masl):
        locs = cax.get_yticks()                               #get tick locations
        nlabels = len(cax.get_yticklabels())                  #get number of initial ticks 
        labels = np.linspace(zmin, zmax, nlabels)             #generate list of tick labels
        ticks = cax.set(yticks=locs,yticklabels=labels)       #set tick locations and labels
        
    return f1,ax1,f2,ax2


    
##############################################################################################
def xyz2rowcollay(X, Y, Z, grid, flip=False):

    '''Converts between X,Y,Z coordinates and row,col,lay coordinates.
    Inputs:
    X: list of x coordinates to convert to columns (or vice-versa)
    Y: list of y coordinates to convert to rows (or vice-versa)
    Z: list of z coordinates to convert to layers (or vice-versa)
    grid: gempy grid object containing resolution info etc., or list [xres,yres,zres,xmin,ymin,zmin,xmax,ymax,zmax]
    flip: if False (default), converts X,Y to row,col. If True, converts row,col to X,Y'''
    
    #Get grid info:
    try:                            #if gempy grid object given
        xres = grid.resolution[0]   #resolution (number of pixels)
        yres = grid.resolution[1]       
        zres = grid.resolution[2]       
        xmin = grid.extent[0]       #min coordinate value (left/bottom)
        ymin = grid.extent[2]
        zmin = grid.extent[4]
        xmax = grid.extent[1]       #max coordinate value (right/top)
        ymax = grid.extent[3]
        zmax = grid.extent[5]
    except:                         #if list given
        xres = grid[0]  #resolution (number of pixels)
        yres = grid[1]       
        zres = grid[2]       
        xmin = grid[3]  #min coordinate value (left/bottom)
        ymin = grid[4]
        zmin = grid[5]
        xmax = grid[6]  #max coordinate value (right/top)
        ymax = grid[7]
        zmax = grid[8]
    dx   = (xmax-xmin)/xres     #pixel width
    dy   = (ymax-ymin)/yres
    dz   = (zmax-zmin)/zres
    
    #Convert X,Y,Z to row,col,lay:
    if flip == False:
        cols = []
        for x in X:
            x   = x - 1             #move down by 1 so that pts on the upper boundary of a cell get counted as part of that cell & not the one above
            col = (x - xmin)//dx    #col:x calculate column index (0 = west/left) c = (x1-x0)/dx 
            if col < 0:             #if index assigned is below the lowest possible index (0)
                col = 0             #move it to be at the lowest index (0) (this only applies in edge cases)
            cols.append(col)

        rows = []
        for y in Y:
            y   = y - 1
            row = (y - ymin)//dy    #row:y calculate row index (0 = north/top)    r = (y1-y0)/dy
            if row < 0:
                row = 0                     
            rows.append(row)

        lays = []
        for z in Z:
            z   = z - 1
            lay = (z - zmin)//dz    #lay:z calculate layer index (0 = bottom)     l = (z1-z0)//dz
            if lay < 0:
                lay = 0
            lays.append(lay)
            
        return rows, cols, lays
    
    #Convert row,col,lay to X,Y,Z:
    else:
        cols = X  #object in X input position is list of columns
        X=[]      #create empty list to store new x values
        for col in cols:
            col = col+1                  #move over by 1 to avoid cutoffs
            x = xmin + (col*dx) - (0.5*dx)  #x: calculate x coordinates for center of cell (x = x0 + col*dx - 0.5dx)
            X.append(x)                  #store x coordinate

        rows = Y  #object in Y input position is list of rows
        Y=[]      #create empty list to store new y values
        for row in rows:
            row = row+1                  #move up by 1 to avoid cutoffs
            y = ymax - (row*dy) + (0.5*dy)  #y: calculate y coordinates for center of cell (y = y0 + row*dy - 0.5dy)
            Y.append(y)                  #store y coordinate
            
        lays = Z  #object in Z input position is list of layers
        Z=[]      #create empty list to store new z values
        for lay in lays:
            lay = lay+1                  #move up by 1 to avoid cutoffs
            z = zmin + (lay*dz) - (0.5*dz)  #z: calculate z coordinates for center of cell (z = z0 + lay*dz - 0.5dz)
            Z.append(z)                  #store y coordinate
            
            return X, Y, Z