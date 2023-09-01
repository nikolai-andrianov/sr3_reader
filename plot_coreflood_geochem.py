"""
Visualization of CMG GEM results for a 1D coreflooding experiment with geochemical reactions using the SR3 reader.

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

from sr3_reader import *
from matplotlib import pyplot
import numpy as np

# Read CMG results file
infile = 'coreflood_geochem.sr3'
sr3 = read_SR3(infile)

# Get the timeseries data from the wells (which are used to set up the boundary conditions)
wells_ts = get_wells_timeseries(sr3)

# Get the available molalities and minerals
molalities, minerals = list_molalities_minerals(sr3)
list_mol = ','.join(molalities)
list_min = ','.join(minerals)
print('Available molalities: ' + list_mol)
print('Available minerals: ' + list_min)

# The grid centers from a vertically (z-direction) located sample
x = sr3.grid.cells.centroids[2, :]
# Get the distance from the top
x = x - x[0]

# Convert x to cm
assert sr3.units['Length'] == 'm'
x = x * 100

# Get the time-dependent spatially-distributed porosity, permeability at start,
# and the water permeability resistance factor (RFW)
(sp_ind, sp) = get_spatial_properties(sr3, ['POROS', 'PERMK', 'RFW'])

# The time steps, at which the spatially-ditributed data is available
t_sp = sr3.times['Days'].iloc[sp_ind]

# Plot vs x at selected time instants
i_plot = np.arange(0, len(t_sp), len(t_sp) // 8)

poro = np.mean(sp['POROS'], axis=1)
perm0 = np.mean(sp['PERMK'], axis=1)
rfw = np.mean(sp['RFW'], axis=1)
perm = np.divide(perm0, rfw)
# Convert permeability to mD
assert sr3.units['Permeability'] == 'darcy'
perm *= 1000
props = [poro, perm]

# Plotting the average porosity & permeability vs time
labels = ['Porosity', 'Permeability (mD)']
figt, axst = pyplot.subplots(len(props), 1, num=1)
for n, (p, label) in enumerate(zip(props, labels)):
    axst[n].plot(t_sp, p, label=label)
    axst[n].legend()
    axst[n].xaxis.set_tick_params(labelbottom=False)

    # Plot the injection periods on the right y-axes
    ax2 = axst[n].twinx()
    ax2.fill_between(wells_ts['CO2-Injector']['Days'],
                     1 - wells_ts['CO2-Injector']['WELLSTATE'], color='red',
                     alpha=0.1)
    ax2.set_yticks([])
    ax2.xaxis.set_tick_params(labelbottom=False)

axst[-1].xaxis.set_tick_params(labelbottom=True)
axst[-1].set_xlabel('Time (' + sr3.units['Time'] + ')')
pyplot.show(block=False)

# Get the molalities
(sp_ind, sp) = get_spatial_indexed(sr3, 'MOLALITY', molalities)

# Plot the ratios of molalities to the molality of Cl- in the cell (inlet is the last cell #30, outlet is the first cell #1)
cell = 1
if 'Cl-' in molalities:
    n = 0
    figtr, axstr = pyplot.subplots(len(molalities) - 1, 1, num=10)
    for name in molalities:
        if name != 'Cl-':
            axstr[n].plot(t_sp, np.divide(sp[name][:, cell], sp['Cl-'][:, cell]), label=name + '/Cl-')
            axstr[n].legend()
            axstr[n].set_xlim(min(t_sp), max(t_sp))
            axstr[n].xaxis.set_tick_params(labelbottom=False)

            # Plot the injection periods on the right y-axes
            ax2 = axstr[n].twinx()
            ax2.fill_between(wells_ts['CO2-Injector']['Days'],
                             1 - wells_ts['CO2-Injector']['WELLSTATE'], color='red',
                             alpha=0.1)
            ax2.set_yticks([])
            ax2.set_xlim(min(t_sp), max(t_sp))
            ax2.xaxis.set_tick_params(labelbottom=False)
            
            n = n + 1

    axstr[-1].xaxis.set_tick_params(labelbottom=True)
    axstr[-1].set_xlabel('Time (' + sr3.units['Time'] + ')')
    axstr[0].set_title('Ratios of molalities at the outlet')
    figtr.set_figwidth(10)
    figtr.set_figheight(10)
    pyplot.show(block=False)


# Get the total mineral changes in moles
if minerals:
    (sp_ind, sp) = get_spatial_indexed(sr3, 'MINERAL', minerals)

    # Fix the minerals' names for plotting
    min_name = {m: m for m in minerals}
    for m in min_name:
        if min_name[m] == 'K-fe_fel':
            min_name[m] = 'Feldspar'
        if min_name[m] == 'Glaconit':
            min_name[m] = 'Glauconite'

    # Plotting the mean value for minerals vs time
    figt, axst = pyplot.subplots(len(minerals), 1, num=33)
    for n, name in enumerate(minerals):
        # Mean value for the whole core
        m = np.mean(sp[name], axis=1)
        axst[n].plot(t_sp, m, label=min_name[name])
        axst[n].legend()
        axst[n].set_xlim(min(t_sp), max(t_sp))
        axst[n].xaxis.set_tick_params(labelbottom=False)

        # Plot the injection periods on the right y-axes
        ax2 = axst[n].twinx()
        ax2.fill_between(wells_ts['CO2-Injector']['Days'],
                         1 - wells_ts['CO2-Injector']['WELLSTATE'], color='red',
                         alpha=0.1)
        ax2.set_yticks([])
        ax2.set_xlim(min(t_sp), max(t_sp))
        ax2.xaxis.set_tick_params(labelbottom=False)

    axst[-1].xaxis.set_tick_params(labelbottom=True)
    axst[-1].set_xlabel('Time (' + sr3.units['Time'] + ')')
    axst[0].set_title('Mineral change (mol/kgw)')
    pyplot.show(block=True)
