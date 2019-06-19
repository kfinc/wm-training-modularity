import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def matrix_networks_plot(M, network_colors, dpi = 300, colorbar = False, group = None, ses = None, suffix = None, out_dir = None):
    """Creates and saves matrixplot with networks color labels. """
    
    small = 15
    medium = 15
    bigger = 15
    
    plt.rc('font', size=small)          # controls default text sizes
    plt.rc('axes', titlesize=small)     # fontsize of the axes title
    plt.rc('axes', linewidth=2.2)
    plt.rc('axes', labelsize=medium)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
    plt.rc('legend', fontsize=small)    # legend fontsize
    plt.rc('figure', titlesize=bigger)  # fontsize of the figure title
    plt.rc('lines', linewidth=2.2, color='gray')
    
    g = sns.clustermap(M,
                   cmap="RdBu_r",
                   row_cluster=False,
                   col_cluster=False,
                   row_colors=network_colors,
                   col_colors=network_colors,
                   linewidths=0,
                   yticklabels=False,
                   xticklabels=False,
                   vmax = 0.8)

    # Adjust the postion of the main colorbar for the heatmap
    g.cax.set_position([.97, .2, .03, .45])
    g.fig.suptitle(f'{group}: {ses}', size = 20)
    #g.ax_heatmap.set_title(f'{group}: {ses}', size = 15)
    g.cax.set_visible(colorbar)
    
    if out_dir == None:
        "Figure not saved"
    else:
    
        if suffix != None:
            g.savefig(f'{out_dir}{group}_{ses}_{suffux}.pdf', dpi=dpi)
        else:
            g.savefig(f'{out_dir}{group}_{ses}.pdf', dpi=dpi)


def swarm_box_plot(x, y, hue, data):
    plt.style.use('seaborn-white')
    plt.rcParams['font.family'] = 'Helvetica'

    plt.figure(figsize = (8, 6))

    ax = sns.swarmplot(x = x, y = y, hue = hue, data = data, dodge = True, alpha = 0.8, size = 8)
    ax = sns.boxplot(x = x, y = y, hue = hue, data = data, dodge = True,
            showcaps = False, boxprops = {'facecolor':'None'},
            showfliers = False)
    plt.xticks(np.arange(4), ('1', '2', '3', '4'))
    ax.set(xlabel='Scan')
    ax.tick_params(axis='both', color = 'black', length = 5, width = 2)
    
    return ax


def swarm_box_plot_integ(x, hue, data, net1, net2):
    
    plt.style.use('seaborn-white')
    plt.rcParams['font.family'] = 'Helvetica'

    small = 25
    medium = 25
    bigger = 25

    plt.rc('font', size=small)          # controls default text sizes
    plt.rc('axes', titlesize=small)     # fontsize of the axes title
    plt.rc('axes', linewidth=2.2)
    plt.rc('axes', labelsize=medium)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
    plt.rc('legend', fontsize=small)    # legend fontsize
    plt.rc('figure', titlesize=bigger)  # fontsize of the figure title
    plt.rc('lines', linewidth=2.2, color='gray')
    
    plt.figure(figsize = (8, 6))

    ax = sns.swarmplot(x = x, y = f'{net1}_integration', hue = hue, data = data[data['Network'] == net2], dodge = True, alpha = 0.8, size = 8)
    ax = sns.boxplot(x = x, y = f'{net1}_integration', hue = hue, data = data[data['Network'] == net2], dodge = True,
            showcaps = False, boxprops = {'facecolor':'None'},
            showfliers = False)
    plt.xticks(np.arange(4), ('Naive', 'Early', 'Middle', 'Late'))
    ax.set(xlabel='Scan')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.tick_params(axis='both', color = 'black', length = 5, width = 2)

    
    return ax


def matrix_networks_slope_plot(M, network_colors, dpi=300, colorbar=False, group=None, suffix=None, out_dir=None, vmin = -0.08, vmax = 0.08):
    """Creates and saves matrixplot with networks color labels. """
    
    small = 15
    medium = 15
    bigger = 15

    plt.rc('font', size=small)          # controls default text sizes
    plt.rc('axes', titlesize=small)     # fontsize of the axes title
    plt.rc('axes', linewidth=2.2)
    plt.rc('axes', labelsize=medium)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
    plt.rc('legend', fontsize=small)    # legend fontsize
    plt.rc('figure', titlesize=bigger)  # fontsize of the figure title
    plt.rc('lines', linewidth=2.2, color='gray')
    g = sns.clustermap(M,
                   cmap="RdBu_r",
                   row_cluster=False,
                   col_cluster=False,
                   row_colors=network_colors,
                   col_colors=network_colors,
                   linewidths=0,
                   yticklabels=False,
                   xticklabels=False, vmin = vmin, vmax = vmax)

    # Adjust the postion of the main colorbar for the heatmap
    g.cax.set_position([.97, .2, .03, .45])
    g.fig.suptitle(f'{group}', size=20)
    #g.ax_heatmap.set_title(f'{group}: {ses}', size = 15)
    g.cax.set_visible(colorbar)
    
    if out_dir == None:
        "Figure not saved"
    else:
    
        if suffix != None:
            g.savefig(f'{out_dir}{group}_{ses}_{suffux}.pdf', dpi=dpi)
        else:
            g.savefig(f'{out_dir}{group}_{ses}.pdf', dpi=dpi)