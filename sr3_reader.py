"""
A reader for Computer Modelling Group Ltd. (CMG) SR3 output files.

Copyright 2023 Nikolai Andrianov, nia@geus.dk

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS
OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import h5py
import os, math
from datetime import datetime
import numpy as np
import pandas as pd
import time


class RawHDF:
    """
    A container for the raw data entries in the SR3 file.
    """

    class Grid:
        """
        A container for the grid structure.
        """

        class Cells:
            """
            A container for the cells structure.
            """

            def __init__(self):
                # Initialize empty structures for the cells' attributes
                self.centroids = np.array([])
                self.volumes = np.array([])
                self.corners = np.array([])
                self.inactive = np.array([])

        def __init__(self):
            # Initialize empty structures for cells
            self.cells = RawHDF.Grid.Cells()
            self.cart_dims = []
            self.grid_dim = 0
            self.n_cells = None
            self.n_active_cells = None

    def __init__(self, file, timeout):
        # Initialize empty structures for the datasets and their names
        self.data = {}
        self.names = []
        # Keep the SR3 file name
        self.file = file
        self.start_time = time.time()
        self.timeout = timeout
        # Time instants
        self.times = pd.DataFrame()
        # Component names
        self.comp_name = []
        # Names of spatial properties
        self.sp_name = []
        # Short and long descriptions of acronyms
        self.acronym_desc_short = {}
        self.acronym_desc_long = {}
        self.grid = RawHDF.Grid()
        self.units = {}

    def __call__(self, name, h5obj):
        # Only h5py datasets have dtype attribute, so we can search on this
        if hasattr(h5obj, 'dtype') and not name in self.names:
            self.names += [name]
            dataset = np.array(self.file[name][:])
            self.data[name] = dataset
            # Check if the file read does not take more than ...
            curr_time = time.time()
            dt_read = round(curr_time - self.start_time)
            if dt_read >= self.timeout:
                # Make sure that TimeSeries/WELLS/Origins is read (to be used in get_wells_timeseries)
                wells_datasets = ['TimeSeries/WELLS/Origins', 'TimeSeries/WELLS/Variables',
                                  'TimeSeries/WELLS/Data', 'TimeSeries/WELLS/Timesteps']
                for name in wells_datasets:
                    if not name in self.names:
                        self.names += [name]
                        dataset = np.array(self.file[name][:])
                        self.data[name] = dataset

                print('Warning: read %s datasets before the timeout of %s sec is exceeded..' %
                      (dt_read, len(self.names)))
                print('Increase the timeout in read_SR3 to get more datasets..')
                return 0
        # No return so that the visit function is recursive


class Completion:
    """
    A completion is associated with the (i,j,k) index of the connected grid cell, and with the timeseries data
    """

    def __init__(self, ijk):
        # Index of the connected grid cell
        self.ijk = ijk
        # Time series data
        self.data = pd.DataFrame()


def read_SR3(infile, timeout=60):
    """
    Reads the specified SR3 file and returns an instance of the RawHDF class with all the datasets from infile,
    unless it takes more than timeout seconds. If it takes more then timeout seconds to read the datasets,
    reading the SR3 file terminates and the returned instance of RawHDF contains only a part of the datasets.
    """
    try:
        sr3_file = h5py.File(infile, 'r')
        sr3 = RawHDF(sr3_file, timeout)
    except:
        raise IOError('Cannot open ' + infile + ' as a HDF file..')

    # Loop through all objects inside the hdf5 file and store the datasets in sr3
    sr3_file.visititems(sr3)

    # Extract report time instants from the master time table, which has entries <index, time_in_days, YYYYMMDD>
    assert ('General/MasterTimeTable' in sr3.data)
    td = []
    date = []
    for t in sr3.data['General/MasterTimeTable']:
        td.append(t[1])
        day = math.floor(t[2])
        frac = t[2] % 1
        hour = math.floor(frac * 24)
        frac = (frac * 24) % 1
        minute = math.floor(frac * 60)
        frac = (frac * 60) % 1
        second = math.floor(frac * 60)
        date_str = str(day) + ' ' + str(hour) + ':' + str(minute) + ':' + str(second)
        dt = datetime.strptime(date_str, "%Y%m%d %H:%M:%S")
        date.append(dt)

    sr3.times['Date'] = date
    sr3.times['Days'] = td

    # Augment the raw SR3 data with the component names, which are byte strings of the format (b'NAME',)
    for name in sr3.data['General/ComponentTable']:
        name_str = str(name)
        sr3.comp_name.append(name_str.split('\'')[1])

    # Augment the raw SR3 data with the names of spatial properties
    for key in sr3.data:
        if 'SpatialProperties/000000/' in key:
            sr3.sp_name.append(key.split('SpatialProperties/000000/')[1])

    # Augment the raw SR3 data with the description of acronyms
    for name in sr3.data['General/NameRecordTable']:
        acronym = str(name[0]).split('\'')[1]
        desc_short = str(name[1]).split('\'')[1]
        desc_long = str(name[2]).split('\'')[1]

        sr3.acronym_desc_short[acronym] = desc_short
        sr3.acronym_desc_long[acronym] = desc_long

    # Populate the grid structure with the properties
    get_grid(sr3)

    # Set up the dict with units
    get_units(sr3)

    return sr3


def search_acronym(sr3, str):
    """
    Prints the SR3 acronym(s), which either contain str, or whose description contains str.
    """
    print('\nSearching for the occurences of "' + str + '" among the descriptions of CMG acronyms..')
    for key in sr3.acronym_desc_long:
        if str.lower() in key.lower() or str.lower() in sr3.acronym_desc_long[key].lower():
            print('  ' + key + ': ' + sr3.acronym_desc_long[key])
    print('Done.\n')


def get_sector_timeseries(sr3):
    """
    Returns the timeseries for the available sectors.
    Currently, only for a single sector.
    """

    # Extract variable names from the byte strings of the format (b'NAME',)
    var_name = []
    for name in sr3.data['TimeSeries/SECTORS/Variables']:
        name_str = str(name)
        var_name.append(name_str.split('\'')[1])

    # Replace the indices in variable names with the corresponding components' names
    comp_fields = ['AQUCOMMOL', 'MINERAMOL', 'MINERACHG']
    for n, name in enumerate(var_name):
        cf = [c for c in comp_fields if c in name]
        if cf:
            # Get the index of the component
            ind = name[name.find('(') + 1:name.find(')')]
            try:
                i = int(ind)
            except:
                raise ValueError('Cannot convert the index ' + ind + ' for the component ' + name)

            # Replace the index with the component name
            var_name[n] = re.sub('\(.*?\)', '_' + sr3.comp_name[i - 1], name)

    # Sector can be FIELD, ..
    sector = sr3.data['TimeSeries/SECTORS/Origins']
    if len(sector) != 1:
        raise ValueError('More than one sector in the SR3 file.. Not implemented..')

    # Read the data for a single sector
    data = sr3.data['TimeSeries/SECTORS/Data']
    data = data.reshape((data.shape[0], data.shape[1]))

    # Augment data with time in days
    assert (len(sr3.times['Days']) == data.shape[0])
    data = np.c_[sr3.times['Days'], data]

    # Return a dataframe, indexed with date
    ts = pd.DataFrame(data=data, index=sr3.times['Date'], columns=['Days'] + var_name)
    return ts


def get_wells_timeseries(sr3):
    """
    Returns the timeseries for the available wells.
    """

    well_names = []
    for name in sr3.data['TimeSeries/WELLS/Origins']:
        name_str = str(name)
        well_names.append(name_str.split('\'')[1])

    # Extract variable names from the byte strings of the format (b'NAME',)
    var_name = []
    for name in sr3.data['TimeSeries/WELLS/Variables']:
        name_str = str(name)
        var_name.append(name_str.split('\'')[1])

    # Don't replace the indices in variable names with the corresponding components' names (for the moment).

    # Read the data for all wells
    data = sr3.data['TimeSeries/WELLS/Data']

    # Get the time instants when the well data is available
    ts_ind = sr3.data['TimeSeries/WELLS/Timesteps']
    t_sp = sr3.times['Days'].iloc[ts_ind]

    # Transform the data as dataframes for separate wells
    wells = {}
    for n, wn in enumerate(well_names):

        # Augment data with time in days
        assert (len(t_sp) == data.shape[0])
        well_data = np.c_[t_sp, data[:, :, n]]

        # Return a dataframe, indexed with date
        wells[wn] = pd.DataFrame(data=well_data, index=t_sp, columns=['Days'] + var_name)

    return wells


def get_layers_timeseries(sr3):
    """
    Returns the timeseries for the available layers.
    """

    # Make sure that LAYERS datasets are available
    layers_datasets = ['TimeSeries/LAYERS/Origins', 'TimeSeries/LAYERS/Variables',
                      'TimeSeries/LAYERS/Data', 'TimeSeries/LAYERS/Timesteps']
    for name in layers_datasets:
        if not name in sr3.names:
            print(name + ' is not available..')
            return None

    layers_names = []
    for name in sr3.data['TimeSeries/LAYERS/Origins']:
        name_str = str(name)
        layers_names.append(name_str.split('\'')[1])

    # Extract variable names from the byte strings of the format (b'NAME',)
    var_name = []
    for name in sr3.data['TimeSeries/LAYERS/Variables']:
        name_str = str(name)
        var_name.append(name_str.split('\'')[1])

    # Don't replace the indices in variable names with the corresponding components' names (for the moment).

    # Read the data for all wells
    data = sr3.data['TimeSeries/LAYERS/Data']

    # Get the time instants when the well data is available
    ts_ind = sr3.data['TimeSeries/LAYERS/Timesteps']
    t_sp = sr3.times['Days'].iloc[ts_ind]

    # Transform the data as dataframes for separate wells
    layers = {}
    for n, wn in enumerate(layers_names):

        # Augment data with time in days
        assert (len(t_sp) == data.shape[0])
        layers_data = np.c_[t_sp, data[:, :, n]]

        # Return a dataframe, indexed with date
        layers[wn] = pd.DataFrame(data=layers_data, index=t_sp, columns=['Days'] + var_name)

    return layers


def get_completions_timeseries(sr3):
    """
    Returns the timeseries for the available completions.

    Apparently completions in SR3 are represented as LAYERS with the names in the format <well_name>{ic,jc,kc}, where
    ic, jc, and kc are the indices of the grid block, connected to the well_name.
    """

    # Make sure that LAYERS datasets are available
    layers_datasets = ['TimeSeries/LAYERS/Origins', 'TimeSeries/LAYERS/Variables',
                      'TimeSeries/LAYERS/Data', 'TimeSeries/LAYERS/Timesteps']
    for name in layers_datasets:
        if not name in sr3.names:
            print(name + ' is not available..')
            return None

    layers_names = []
    for name in sr3.data['TimeSeries/LAYERS/Origins']:
        name_str = str(name)
        layers_names.append(name_str.split('\'')[1])

    # Extract variable names from the byte strings of the format (b'NAME',)
    # Don't replace the indices in variable names with the corresponding components' names (for the moment).
    var_name = []
    for name in sr3.data['TimeSeries/LAYERS/Variables']:
        name_str = str(name)
        var_name.append(name_str.split('\'')[1])

    # Get the time instants when the well data is available
    ts_ind = sr3.data['TimeSeries/LAYERS/Timesteps']
    t_sp = sr3.times['Days'].iloc[ts_ind]

    # Read the data for all completions
    data = sr3.data['TimeSeries/LAYERS/Data']
    assert (len(t_sp) == data.shape[0])

    # Extract well names and completion indices from layers_names (which have the format format <well_name>{ic,jc,kc})
    # and create the completions structure
    completions = {}
    for n, ln in enumerate(layers_names):
        try:
            # Parse layers_names
            wn = ln.split('{')[0]
            comp = ln.split('{')[1]
            comp = comp.split('}')[0]
            comp = comp.split(',')
            ind = [int(c) for c in comp]
        except:
            raise ValueError('Cannot extract wells and completion indices from ' + ln)

        # Write the well completions
        if wn not in completions:
            completions[wn] = []

        # Add the information on the connected block to the current completion
        completions[wn].append(Completion(ind))

        # Add the timeseries data, augmented with time in days
        layers_data = np.c_[t_sp, data[:, :, n]]

        # Return a dataframe, indexed with date
        completions[wn][-1].data = pd.DataFrame(data=layers_data, index=t_sp, columns=['Days'] + var_name)

    return completions


def get_wells_completions(sr3):
    """
    Returns a dict the indices of completed cells for the available wells.
    """

    well_names = []
    for name in sr3.data['TimeSeries/WELLS/Origins']:
        name_str = str(name)
        well_names.append(name_str.split('\'')[1])

    comp = {wn: [] for wn in well_names}
    for data in sr3.data['TimeSeries/LAYERS/LayerTable']:
        well_name = str(data[3]).split('\'')[1]
        comp_cell = str(data[2]).split('\'')[1]
        # Get to the 0-based indexing
        comp_cell_ind = np.array(comp_cell.split(','), dtype=int) - 1
        comp[well_name].append(comp_cell_ind)

    return comp


def list_spatial_properties(sr3):
    """
    Returns a list of spatial properties' names.
    """

    # Get a list of spatial properties
    prop_name = []
    for key in sr3.data:
        if 'SpatialProperties/000000/' in key:
            prop_name.append(key.split('SpatialProperties/000000/')[1])

    return prop_name


def list_data_categories(sr3, level, exclude_numeric):
    """
    Returns a list of data categories among the data names.

    The data names are strings, separated by slash "/", and category of n-th level is defined as substrings of data
    names between the (n-1)-th and n-th slash.
     E.g., the 0-th level categories are "General', 'Restart', 'SpaialProperties' etc.

     Exclude the numerical fields by setting exclude_numeric=True.
    """

    # Get a list of unique categories
    cat_name = []
    for key in sr3.data:
        token = key.split('/')
        if level < len(token):
            cat = token[level]
            if exclude_numeric:
                try:
                    num = int(cat)
                except:
                    cat_name.append(cat)
            else:
                cat_name.append(cat)

    # Keep the unique names only
    cat_name = set(cat_name)
    cat_name = list(cat_name)

    return cat_name


def list_wells_properties(sr3):
    """
    Returns a list of wells properties' names.
    """

    # Get a list of spatial properties
    prop_name = []
    for key in sr3.data:
        if 'WELLS' in key:
            val = sr3.data[key]
            prop_name.append(key)

    return prop_name


def list_timeseries_properties(sr3):
    """
    Returns a list of time series' names.
    """

    # Get a list of spatial properties
    prop_name = []
    for key in sr3.data:
        if 'TimeSeries' in key:
            val = sr3.data[key]
            prop_name.append(key)

    return prop_name


def get_spatial_properties(sr3, sel_names, activeonly=True, verbose=True):
    """
    Returns the tuple (sp_ind, sp), where sp is a dict of the time-dependent spatial properties with names, specified
    in sel_names, and sp_ind is a list of indices of corresponding time instants in sr3.times.

    If activeonly==True, the requested spatial parameters are defined on active cells only (e.g. pressure or saturation),
    otherwise the requested spatial parameters are defined for all cells (e.g. cell volumes).
    """

    if type(sel_names) is not list:
        if type(sel_names) is str:
            sel_names = [sel_names]
        else:
            print('Usage: get_spatial_properties(sr3, <list of CVM acronyms>')
            return {}, []

    # Get a list of time instants when the spatial properties are available
    sp = {}
    sp_ind = []
    for name in sel_names:
        prop = []
        for key in sr3.data:
            if 'SpatialProperties' in key:
                # Extract the time index and and the key variable name from the key of the form 'SpatialProperties/000000/VISG'
                key_parts = key.split('/')
                ind_time = key_parts[1]
                try:
                    key_var = key.split('SpatialProperties/' + ind_time + '/')[-1]
                except:
                    key_var = ''
                if name == key_var:
                    try:
                        i = int(ind_time)
                    except:
                        raise ValueError('Cannot convert the index ' + ind_time + ' for the component ' + name)
                    if i not in sp_ind:
                        sp_ind.append(i)

                    # Accumulate the spatial distribution at the current time step
                    prop.append(sr3.data[key])

        # Convert the time-dependent spatial distribution to a 3D numpy array
        if not prop:
            if verbose:
                print('Error: ' + name + ' not found among the SR3 spatial parameters! Empty array returned..')
            tmp = prop
        else:
            # If the dimension of the requested property is equal to the number of active cells, reshape the spatial property
            # to match active cells; the inactive cells are assigned with zero.
            # Otherwise return with the property as is.
            tmp = np.array(prop)
            ntimes = tmp.shape[0]
            n = tmp.shape[1]

            # The number of active cells is defined in get_grid()
            if sr3.grid.n_active_cells:
                if n == sr3.grid.n_active_cells:
                    ina = np.tile(sr3.grid.cells.inactive, ntimes)
                    tmp = np.zeros(sr3.grid.n_cells * ntimes)
                    tmp[~ina] = np.array(prop).flatten()
                    tmp = np.reshape(tmp, (ntimes, sr3.grid.n_cells))

        sp[name] = tmp

    return sp_ind, sp


def get_spatial_indexed(sr3, prefix, comps):
    """
    Returns the tuple (sp_ind, sp), where sp is a list of the time-dependent indexed properties (such as MOLALITY or
    MINERAL) of the specified components, and sp_ind is a list of indices of corresponding time instants in sr3.times.
    """

    # Get a list of internal names of the desired components
    sel_names = []
    for i, c in enumerate(comps):
        try:
            i = sr3.comp_name.index(c) + 1
            # ic.append(i)
            sel_names.append(prefix + '(' + str(i) + ')')
        except ValueError:
            print(c + ' is not found among available components...')

    (sp_ind, sp) = get_spatial_properties(sr3, sel_names)

    # Change the keys from prefix(i) to component names
    assert len(sp) == len(comps), 'The dimensions of the indexed spatial properties do not match the number of components!'
    sp_comps = {}
    for c, key in zip(comps, sp):
        sp_comps[c] = sp[key]

    return sp_ind, sp_comps


def list_molalities_minerals(sr3):
    """
    Returns the lists of spatially distributed molalities and minerals.
    """

    # Get the GEM indices between the brackets in MOLALITY(1) etc
    ind_mol = [int(name[name.find('(') + 1:name.find(')')]) for name in sr3.sp_name if 'MOLALITY' in name]
    molalities = [sr3.comp_name[i - 1] for i in ind_mol]

    # Get the GEM indices between the brackets in MINERAL(1) etc
    ind_min = [int(name[name.find('(') + 1:name.find(')')]) for name in sr3.sp_name if 'MINERAL' in name]
    minerals = [sr3.comp_name[i - 1] for i in ind_min]

    return molalities, minerals


def get_grid(sr3):
    """
    Assigns the grid structure with cells' volumes and centroids, and sets the grid dimensions.

    Not clear:
        * How the depths from SR3 are related to cell tops (the end values do not correspond to DTOP in the input)
        * What are the volumes (dimensions 2 times the cells' dimensions - due to dual porosity?)

    """

    # Set up the boolean array of inactive cells
    (sp_ind, sp) = get_spatial_properties(sr3, 'GRID/ICSTPS', activeonly=False)
    inactind = sp['GRID/ICSTPS'] == 0
    sr3.grid.cells.inactive = inactind[0]
    sr3.grid.n_active_cells = np.count_nonzero(~sr3.grid.cells.inactive)
    sr3.grid.n_cells = len(sr3.grid.cells.inactive)

    # Apparently GRID/NODES are only present in SR3 for corner-point grid geometry
    (sp_ind, sp) = get_spatial_properties(sr3, 'GRID/NODES', verbose=False)
    if len(sp['GRID/NODES']) > 0:
        print('Corner-point grid format detected, the grid geometry might not be visualized properly...')

    (sp_ind, sp) = get_spatial_properties(sr3, ['GRID/BLOCKDEPTH', 'GRID/BLOCKSIZE', # 'GRID/BVOL',
                                                'GRID/IGNTFR', 'GRID/IGNTGT', 'GRID/IGNTID', 'GRID/IGNTJD',
                                                'GRID/IGNTKD', 'GRID/IGNTNC', 'GRID/IGNTNS', 'GRID/IGNTZA'],
                                                 activeonly=False)

    # Assert that the grid properties are provided at t=0 only (can be different in geomechanics simulations)
    if len(sp_ind) != 1:
        print('The grid properties are not provided at several time instants, '
              'but the grid is constructed using the data at t=0 only!')

    assert sp_ind[0] == 0, 'The grid properties are not provided at t=0 only!'



    depths = np.array(sp['GRID/BLOCKDEPTH'][0])
    size = np.array(sp['GRID/BLOCKSIZE'][0])

    # Array dimensons do not match at least in the case of dual porosity (gmgmc080.sr3)
    #volumes = np.array(sp['GRID/BVOL'][0])

    # Get the block dimensions
    dx = size[0::3]
    dy = size[1::3]
    dz = size[2::3]

    # Get the number of Cartesian grid blocks
    Ni = int(sp['GRID/IGNTID'][0])
    Nj = int(sp['GRID/IGNTJD'][0])
    Nk = int(sp['GRID/IGNTKD'][0])
    sr3.grid.cart_dims = (Ni, Nj, Nk)
    # sr3.grid.n_cells = Ni * Nj * Nk
    sr3.grid.grid_dim = sum(1 for n in sr3.grid.cart_dims if n > 1)

    # Calculate the centroids and volumes
    dx = np.reshape(dx, sr3.grid.cart_dims, order='F')
    dy = np.reshape(dy, sr3.grid.cart_dims, order='F')
    dz = np.reshape(dz, sr3.grid.cart_dims, order='F')

    # The depths are related to cell centers
    d = np.reshape(depths, sr3.grid.cart_dims, order='F')
    depthsf = d.flatten(order='F')

    # Cell centroids
    xc = np.cumsum(dx, axis=0) - 0.5 * dx
    yc = np.cumsum(dy, axis=1) - 0.5 * dy
    zc = np.cumsum(dz, axis=2) - 0.5 * dz

    # Flatten so that i changes fastest, k - slowest
    xcf = xc.flatten(order='F')
    ycf = yc.flatten(order='F')
    zcf = zc.flatten(order='F')

    # Assign the sr3.grid.cells fields
    sr3.grid.cells.centroids = np.array([xcf, ycf, depthsf])

    volumes = dx * dy * dz
    sr3.grid.cells.volumes = volumes.flatten(order='F')

    # Get the cells' corners coordinates excluding the 0-th slices in x-, y, and z-directions
    x = np.cumsum(dx, axis=0)
    y = np.cumsum(dy, axis=1)
    z = np.cumsum(dz, axis=2)

    # Prepare the zeros arrays with proper dimensions
    xx = np.zeros((Ni + 1, Nj + 1, Nk + 1))
    yy = np.zeros((Ni + 1, Nj + 1, Nk + 1))
    zz = np.zeros((Ni + 1, Nj + 1, Nk + 1))

    # Immerse the available coordinates into the arrays of proper dimensions
    xx[1:, 1:, 1:] = x
    yy[1:, 1:, 1:] = y
    zz[1:, 1:, 1:] = z

    # Set up the coordinates for the 0-th slices in x-, y, and z-directions
    xx[:, 0, :] = xx[:, 1, :]
    xx[:, :, 0] = xx[:, :, 1]

    yy[0, :, :] = yy[1, :, :]
    yy[:, :, 0] = yy[:, :, 1]

    zz[0, :, :] = zz[1, :, :]
    zz[:, 0, :] = zz[:, 1, :]

    # Shift the depths to be negative so that abs(depths) grows downwards
    bottom = d[0, 0, 0] + dz[0, 0, 0] / 2
    zz = zz - bottom

    # Duplicate all coordinates ...
    xr = xx.copy()
    yr = yy.copy()
    zr = zz.copy()
    for n in range(3):
        xr = np.repeat(xr, 2, axis=n)
        yr = np.repeat(yr, 2, axis=n)
        zr = np.repeat(zr, 2, axis=n)

    # ... except the ones from the edges to get the coordinates of 8 corners per grid block
    for n in range(3):
        xr = np.delete(xr, [0, -1], axis=n)
        yr = np.delete(yr, [0, -1], axis=n)
        zr = np.delete(zr, [0, -1], axis=n)

    # Flatten so that i changes fastest, k - slowest
    xf = xr.flatten(order='F')
    yf = yr.flatten(order='F')
    zf = zr.flatten(order='F')

    # Grid the blocks' corners coordinates
    corners = np.array([xf, yf, zf])
    sr3.grid.cells.corners = np.transpose(corners)


def get_units(sr3):
    """
    Assigns the units structure with the values from SR3.

    The UnitsTable in SR3 has entries of the form (1, b'Time', b'day', b'day'), (2, b'Temperature', b'K', b'F'), etc.
    It seems that the units of the data are the 2nd element of the corresponding tuples in the 0-based indexing, and the
    3rd element is the unit in the system, selected in the GEM input file (SI or FIELD).

    The implementation below returns the internal CMG units.
    """

    for element in sr3.data['General/UnitsTable']:
        # Convert the byte strings of the format b'Time'
        name = str(element[1])
        unit = str(element[2])
        name = name.split('\'')[1]
        unit = unit.split('\'')[1]
        sr3.units[name] = unit
