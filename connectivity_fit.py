"""
1. integrate layers
    map depths to layers
2. fit the data to what?
3. plot the connectivity

# Still don't know what these are:
data['roi']
# 0.00000e+00, 5.35920e-06, 9.28200e-06, 3.16764e-06, 4.53390e-06
data['roe']
# 2.1900e-07, 4.5240e-05, 4.7600e-05, 3.4358e-05, 5.2705e-05

"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy.optimize import curve_fit
import os
import distributions as pdfs

# Densities as given in the paper for layer 1, 23, 4, 5, 6
# Given in paragraph after equation (1)
CELL_DENSITY_E = [220, 45_240, 47_600, 34_360, 52_710]
CELL_DENSITY_I = [7_080, 12_760, 11_900, 7_540, 10_790]

# Helper dictionaries to map list index to actual layer and back
# Layer 1 is omitted from analysis, since no data were provided
IDX_TO_LAYER = {0: 23, 1: 4, 2: 5, 3: 6}
LAYER_TO_IDX = {IDX_TO_LAYER[key]: key for key in IDX_TO_LAYER.keys()}

# Data file
DATA_FILE = './data/Cat_Np.mat'

# DATA PROCESSING

def load_data():
    """Loads and preprocesses the data from the data file.
    """

    raw_data = sio.loadmat(DATA_FILE)  # loads matlab data

    data = dict()

    data['layers'] = raw_data['htz'].squeeze()  # borders of layers
    # 81.55475117,  587.07801261,  922.19052739, 1169.97051974, 1491.72
    data['z'] = np.squeeze(raw_data['z_int'])  # depth of index
    data['r'] = raw_data['Rxy_int'].squeeze().astype(float)  # lateral position

    data['layer_idx'] = dict()  # layer : indices values within the layer
    for i in range(len(data['layers'])-1):
        layer_name = IDX_TO_LAYER[i]
        data['layer_idx'][layer_name] = np.where(
                                        (data['z'] >= data['layers'][i]) 
                                         & (data['z'] < data['layers'][i+1]))[0]

    data['n'] = dict()  # expected number of potential synapses
    data['p'] = dict()  # probability of potential connectivity
    for conn, data_key in zip([data['n'], data['p']], ['Np', 'P']):
        for combination in ['ee', 'ei', 'ie', 'ii']:
            key = f'{data_key}_int_{combination}'
            conn[combination] = raw_data[key]

    return data

def split_layers(data: dict, layer_idx) -> dict:
    """Splits a dictionary of connectivity based on projections
    
    Returns:
        dict: A dictionary containing the connectivity data split by layers.
            e.g. new_data['ee']['4-4'] contains array of projections L4ExcL4Exc
            
    """
    new_data = dict()
    for key, value in data.items():
        new_data[key] = dict()
        for pre_l, pre_idx in layer_idx.items():
            for post_l, post_idx in layer_idx.items():
                conn_key = f'{pre_l}-{post_l}'
                new_data[key][conn_key] = value[pre_idx][:,post_idx]
    return new_data

def sum_layers(data, layer_idx):
    summed_data = dict()
    for key in data.keys():  # 'ee', 'ei', 'ie', 'ii'
        summed_data[key] = dict()

        for pre_l in layer_idx.keys():
            for post_l in layer_idx.keys():
                conn_key = f'{pre_l}-{post_l}'

                pre_idx = layer_idx[pre_l]
                post_idx = layer_idx[post_l]
                proj_data = data[key][pre_idx][:,post_idx]
                
                nan_values = np.isnan(proj_data)
                if not nan_values.all():
                    projections = proj_data.shape[0] * proj_data.shape[1]
                    invalid_projs = nan_values.sum(axis=(0,1))[0]

                    valid_projs = projections - invalid_projs
                    # valid_projs = (proj_data.size - nan_values.sum())/51.

                    # to sum properly, replace NaNs with 0
                    proj_data[nan_values] = 0
                    
                    # summation of data
                    summed_data[key][conn_key] = proj_data.sum(axis=(0,1)) / valid_projs
                else:
                    summed_data[key][conn_key] = np.zeros_like(proj_data[0,0])
    return summed_data

# HELPERS

def integrate(array, axis):
    # Very primitive, just to make code readable
    return array.sum(axis=axis)

def nan_const_along_axis(data, axis=2):
    """Test if NaNs are constant along a given axis.

    Args:
        data (dict): A dictionary containing the data as numpy arrays.
        axis (int, optional): The axis along which to test. Defaults to 2.
    Returns:
        bool: True if NaNs are constant along the axis.
    """
    const = []
    for key, value in data.items():
        nan_mask = np.isnan(value)
        nan_slice = nan_mask.take(0, axis=axis)
        c = all(np.all(nan_slice==nan_mask.take(i, axis=axis)) for i in range(value.shape[axis]))
        const.append(c)
    return all(const)

def data_available(data: dict) -> dict:
    """Calculates the fraction of data present in each projection

    Fraction is calculated as the count of NaNs divided by the total number of elements.

    Args:
        data (dict): A dictionary containing the data.

    Returns:
        dict: A dictionary containing the connectivity count for each connection.
    """
    data_available = dict()
    for key, sub_data in data.items():
        data_available[key] = dict()
        for conn_key, conn_data in sub_data.items():
            frac = 1 - np.isnan(conn_data).sum()/conn_data.size
            data_available[key][conn_key] = round(frac,2)
    return data_available

# REPLICATING PAPER FIGURES

def projection_heatmap(data, r_data, z_data, z_point, layers, fig_name):
    symmetric_data = np.concatenate((data[:,::-1], data[:,1:]), axis=-1)
    x_data = np.concatenate((-r_data[::-1], r_data[1:]), axis=-1)
    y_data = z_data.astype(float)

    fig, axis = plt.subplots()
    
    cs =  axis.contourf(x_data, y_data, symmetric_data)
    fig.colorbar(cs)
    axis.plot(0, z_point, 'ro')
    axis.invert_yaxis()
    axis.set_xlabel('r')
    axis.set_ylabel('z')
    for layer in layers:
        axis.axhline(y=layer, color='black', linestyle='--')
    fig.savefig(f'./plots/{fig_name}.png')
    plt.close(fig)

def incoming_heatmap(connectivity_data, z_post, z_data, r_data, layers, 
                     fig_name):
    z_idx = np.where(z_data > z_post)[0][0]
    projection_heatmap(connectivity_data[:,z_idx], r_data, z_data, z_post, 
                       layers, fig_name)

def outgoing_heatmap(connectivity_data, z_pre, z_data, r_data, layers, 
                     fig_name):
    z_idx = np.where(z_data > z_pre)[0][0]
    projection_heatmap(connectivity_data[z_idx], r_data, z_data, z_pre, layers, 
                       fig_name)

def matrix_heatmap(data, displacement,  z_data, r_data, layers, fig_name):
    r_idx = np.where(r_data > displacement)[0][0]

    fig, axis = plt.subplots()
    cs =  axis.contourf(z_data, z_data, data[:,:,r_idx])
    fig.colorbar(cs)
    axis.invert_yaxis()
    axis.set_xlabel('z_post')
    axis.set_ylabel('z_pre')
    for layer in layers:
        axis.axhline(y=layer, color='black', linestyle='--')
        axis.axvline(x=layer, color='black', linestyle='--')
    fig.savefig(f'./plots/{fig_name}.png')
    plt.close(fig)

def connections_in(data, r_data, layer_idx, pre_type : str='e'):
    data = data.copy()
    data[np.isnan(data)] = 0
    density = np.zeros(data.shape[0])
    if pre_type == 'e':
        cell_density = CELL_DENSITY_E
    elif pre_type == 'i':
        cell_density = CELL_DENSITY_I
    else:
        raise ValueError('pre_type must be either "e" or "i"')
    for layer, idx in layer_idx.items():
        density[idx] = cell_density[LAYER_TO_IDX[layer]+1]
    density = density[:, np.newaxis, np.newaxis]
    rho = r_data[np.newaxis, np.newaxis]
    return integrate(data*density*rho*2*np.pi, axis=(0,2))

def connections_out(data, r_data, layer_idx, post_type : str='e'):
    data = data.copy()
    data[np.isnan(data)] = 0
    density = np.zeros(data.shape[1])
    if post_type == 'e':
        cell_density = CELL_DENSITY_E
    elif post_type == 'i':
        cell_density = CELL_DENSITY_I
    else:
        raise ValueError('post_type must be either "e" or "i"')
    for layer, idx in layer_idx.items():
        density[idx] = cell_density[LAYER_TO_IDX[layer]+1]
    density = density[np.newaxis, :, np.newaxis]
    rho = r_data[np.newaxis, np.newaxis]
    return integrate(data*density*rho*2*np.pi, axis=(1,2))

def plot_fig_9a(convergence, divergence, z_data, layers, fig_name):
    fig, axis = plt.subplots()
    axis.plot(z_data, convergence, 'b', label='convergence')
    axis.plot(z_data, divergence, 'r', label='divergence')
    for layer in layers:
        axis.axvline(x=layer, color='black', linestyle='--')
    axis.legend()
    fig.savefig(f'./plots/{fig_name}.png')
    plt.close(fig)

# FITTING

def fit_projection(projection_data, x_data, distributions, inits):
    """Fits a projection_data to a given list of distributions
    """
    fits = []
    for pdf, init in zip(distributions, inits):
        try:
            fit_pars, fit_cov = curve_fit(pdf, x_data, projection_data, p0=init)
            fits.append(fit_pars)
        except RuntimeError:
            print(f"Fit failed for {pdf.__name__}")
            fits.append(None)
    return fits

def plot_fit(axis, proj_data, x_data, fits, distributions):
    axis.plot(x_data, proj_data, label='data')
    for pdf, fit_pars in zip(distributions, fits):
        if fit_pars is not None:
            label = f'{pdf.__name__}'
            axis.plot(x_data, pdf(x_data, *fit_pars), label=label)
    axis.legend()

def make_fit(proj_data, x_data, proj_name, distributions, inits, file_name, 
             axis=None):
    if np.all(proj_data==0):  # skip empty data
        return []

    fits = fit_projection(proj_data, x_data, distributions, inits)
    data_row = f"{proj_name} | "
    with open(f'./{file_name}', 'a') as f:
        for pdf, fit_pars in zip(distributions, fits):
            f.write(data_row + pdf.__name__ + ' | '
                    + ' | '.join([str(x) for x in fit_pars])
                    + '\n')        

    if axis is None:
        fig, axis = plt.subplots()
        fig.suptitle(proj_name)
    else:
        fig = None
    plot_fit(axis, proj_data, x_data, fits, distributions)
    if fig is not None:
        fig.savefig(f'./plots/{file_name}_{proj_name}.png')
        plt.close(fig)

    return fits

def make_bulk_fit(proj_dict, x_data, distributions, inits, file_name):
    """
    example p_summed['ee']

    """
    fig, axes = plt.subplots(4, 4, figsize=(12, 12))
    for projection, projection_data in proj_dict.items():

        # prepares fig for table plot
        row, col = [LAYER_TO_IDX[int(x)] for x in projection.split('-')]
        ax = axes[row, col]
        if row == 0:
            ax.set_title(f'Layer {IDX_TO_LAYER[col]}')
        if col == 0:
            ax.set_ylabel(f'Layer {IDX_TO_LAYER[row]}')

        fits = make_fit(projection_data, x_data, projection, distributions, inits, file_name)

        if np.all(projection_data==0):  # skip empty data
            ax.axis('off')
            # ax.set_visible(False)
        else:
            plot_fit(ax, projection_data, x_data, fits, distributions)
    fig.savefig(f'./plots/table/{file_name}.png')
    plt.close(fig)


def main():
    # TODO: Quality of fit?
    # TODO: how to fit the normalization?
    # Or how is the probability normalized?
    # TODO: relationship between N and P
    # TODO: Fit the data for the model!

    # Shape fits the pictures but the values are slightly off, why?
    # eg divergence, convergence seem off by 0.5
    # 
    # suggest number of conns is also a bit off
    # 
    # z_pre = np.where(data['z'] == 405)[0][0]  # 40
    # z_post = np.where(data['z'] == 725)[0][0]  # 72
    # rho = np.where(data['r'] == 70)[0][0]  # 7
    # data['n']['ee'][z_pre, z_post, rho]  # 1.438085

    #####################################
    # LOADING DATA AND REPLICATING FIGS #
    #####################################

    data = load_data()

    # Plots Fig 5A
    outgoing_heatmap(data['n']['ee'], 650, data['z'], data['r'], 
                     data['layers'], 'fig_5A')

    # Plots Fig 5D
    outgoing_heatmap(data['n']['ee'], 400, data['z'], data['r'], 
                     data['layers'], 'fig_5D')

    # Plots Fig 10A
    matrix_heatmap(data['n']['ee'], 0, data['z'], data['r'], 
                   data['layers'], 'fig_10A')

    # Plots Fig 9A
    convergence = connections_in(data['n']['ee'], data['r'], data['layer_idx'],
                                 pre_type='e')
    divergence = connections_out(data['n']['ee'], data['r'], data['layer_idx'], 
                                 post_type='e')
    plot_fig_9a(convergence, divergence, data['z'], data['layers'], 'fig_9A')

    #######################################
    # THIS PART IS FOR EXPLORING THE DATA #
    #######################################

    # NOT TESTED!
    # nan_const_along_axis(data['n'])  # True
    # nan_const_along_axis(data['p'])  # True

    n_splitted = split_layers(data['n'], data['layer_idx'])
    p_splitted = split_layers(data['p'], data['layer_idx'])

    # n_available = data_available(n_splitted)
    # for key, value in n_available['ii'].items():
    #     print(f'{key}: {value}')

    ###########################################################
    # THIS PART IS FOR FITTING THE PROJECTIONS BETWEEN LAYERS #
    ###########################################################

    n_summed = sum_layers(data['n'], data['layer_idx'])
    p_summed = sum_layers(data['p'], data['layer_idx'])

    # Set up the fitting
    distributions = [pdfs.hyperbolic_pdf]
    inits = [[0.01, 180.]]

    # To run all the fitting run this:
    make_bulk_fit(p_summed['ee'], data['r'], distributions, inits,
                  'fit_pars_test.txt')

    # To fit a single projection run this:
    make_fit(p_summed['ee']['6-4'], data['r'], 'L6ExcL4Exc', distributions, 
             inits,'fit_pars_test.txt')



if __name__ == '__main__':
    main()