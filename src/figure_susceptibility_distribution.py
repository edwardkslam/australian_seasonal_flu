#!/usr/bin/env python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy as np
import plotting_style as ps
import model_parameters as mp

from plotting_functions import get_position, add_bounding_subplot
from flu_final_size_model import get_dist_susses, compare_final_sizes, calc_R0_from_Reff
from susceptibility_models import pick_sus_func

def plot_susceptibility(
        n_classes=5,
        cmap=plt.cm.Blues,
        axis=None,
        susceptibility_model="linear",
        homotypic_protection=1,
        broadly_neutralizing_sus=None,
        **kwargs):
    if axis is None:
        fig, axis = plt.subplots()

    max_var = n_classes - 2 - (broadly_neutralizing_sus is not None)
    variants = list(range(max_var))
    escape_factors = [0.1, 0.2, 0.3]
    color_vals = np.linspace(0.4, 0.8, len(escape_factors))
        
    for e, col_val in zip(escape_factors, color_vals):
        sus_func = pick_sus_func(susceptibility_model,
                                 e,
                                 homotypic_protection=homotypic_protection)

        susses = [0] + [sus_func(var) for var in variants]
        axis.plot([-1] + variants, susses, "-o",
                  color=cmap(col_val),
                  label=e,
                  **kwargs)
        axis.axvline(max_var - 0.5,
                     linestyle="dashed",
                     color='k',
                     lw=2)
        axis.plot([max_var], [1],
                  "o",
                  color=cmap(col_val),
                  **kwargs)

        if broadly_neutralizing_sus is not None:
            axis.plot([max_var + 1], [broadly_neutralizing_sus],
                      "o",
                      color=cmap(col_val),
                      **kwargs)
            
    axis.legend(frameon=True,
                framealpha=1,
                fancybox=True,
                title="escape",
                loc='upper left')
    axis.set_xticks(list(range(-1, n_classes)))
    axis.set_xticklabels(["R"] +
                         list(range(n_classes - 3)) +
                         ["S", "B"])


def plot_distribution(
        distribution,
        axis=None,
        cmap=plt.cm.Blues):
    if axis is None:
        fig, axis = plt.subplots()

    # hackily reorder classes for plotting
    sus_classes = list(distribution[2:-1])
    sus_dist = ([distribution[0]] + sus_classes +
                [distribution[1], distribution[-1]])
    
    xlabs = (["R"] + [str(num) for num in range(0, len(sus_dist) - 3)] +
             ["S", "B"])
    print(len(sus_dist))
    print(len(xlabs))
    axis.bar(xlabs,
             sus_dist,
             color=cmap(0.8))


def plot_differences_by_R(
        distribution,
        cmap=plt.cm.Blues,
        susceptibility_model="linear",
        escapes=[0.1, 0.2, 0.3],
        homotypic_protection=1,
        axes=None,
        plot_wt=False,
        broadly_neutralizing_sus=mp.broadly_neutralizing_sus,
        lw=None,
        max_R=2,
        by="R0",
        **kwargs):
    
    if axes is None:
        fig, axes = plt.subplots(2, 1,
                                 sharex=True,
                                 sharey=False)

    Rs = np.linspace(1, max_R, 50)
    color_vals = np.linspace(0.4, 0.8, len(escapes))
    
    for e, c_val in zip(escapes, color_vals):
        if by == "Reff":
            mut_susses = get_dist_susses(
                -1,
                distribution,
                e,
                susceptibility_model=susceptibility_model,
                homotypic_protection=homotypic_protection,
                broadly_neutralizing_sus=broadly_neutralizing_sus)
            
            R0s = np.array([calc_R0_from_Reff(mut_susses,
                                              distribution,
                                              R)
                            for R in Rs])
        elif by == "R0":
            R0s = Rs
        else:
            raise ValueError(" 'by' must be either 'Reff' or 'R0'")
        
        vals = [compare_final_sizes(
            distribution,
            e,
            R0,
            homotypic_protection=homotypic_protection,
            susceptibility_model=susceptibility_model,
            broadly_neutralizing_sus=broadly_neutralizing_sus)
                for R0 in R0s]

        diffs = [val[3] for val in vals]
        wt_final_sizes = [np.sum(val[0]) for val in vals]
        mut_final_sizes = [np.sum(val[1]) for val in vals]
        
        axes[1].plot(Rs, diffs,
                     color=cmap(c_val),
                     label=e,
                     lw=lw,
                     **kwargs)
        
        axes[0].plot(Rs, mut_final_sizes,
                     color=cmap(c_val),
                     label=e,
                     lw=lw,
                     **kwargs)
        if plot_wt:
            axes[0].plot(Rs, wt_final_sizes,
                         color=cmap(c_val),
                         linestyle="dashed",
                         lw=lw/2,
                         **kwargs)

    axes[1].grid(b=True)
    axes[0].grid(b=True)
    axes[0].set_xlim(left=1)
    axes[0].set_ylim(bottom=0)
    axes[1].set_xlim(left=1)
    axes[1].set_ylim(bottom=0)

    axes[0].legend(frameon=True,
                   fancybox=True,
                   framealpha=1,
                   title="escape")

def make_susceptibility_distribution_figure(
        homotypic_protection=mp.default_homotypic_protection,
        broadly_neutralizing_sus=mp.broadly_neutralizing_sus,
        example_dists=mp.example_sus_distributions,
        output_path="../out/figure_susceptibility_distribution.pdf"):

    # do input checking 
    class_totals = np.array([sum(dist) for dist in example_dists])
    if not np.all(np.abs(class_totals - 1) < 1e-7):
        raise ValueError("all susceptibility "
                         "distributions must "
                         "sum to 1\n\n"
                         "Here we have the "
                         "following sums: {}"
                         "\n\n".format(class_totals))
                         
    class_nums = np.array([len(dist) for dist in example_dists])
    if not np.all(class_nums == class_nums[0]):
        raise ValueError("all susceptibility "
                         "distributions must "
                         "have the same number "
                         "of host classes\n\n"
                         "Here we have the "
                         "following counts: {}"
                         "\n\n".format(class_nums))
    
    # model params
    max_R = 2
    max_classes = class_nums[0]

    # figure params / setup
    width = 10
    height = 15
    figsize=(width, height)
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(4, 4)
    n_rows = 4
    n_cols = 4
    all_dists = add_bounding_subplot(fig, gs[0,:])
    all_susses = add_bounding_subplot(fig, gs[1,:])
    all_final_sizes = add_bounding_subplot(fig, gs[2,:])
    all_diffs = add_bounding_subplot(fig, gs[3,:])
    
    dists = [fig.add_subplot(n_rows, n_cols, get_position(1, 1, n_cols))]
    dists += [fig.add_subplot(n_rows, n_cols,
                              get_position(1, k, n_cols),
                              sharex=dists[0],
                              sharey=dists[0])
              for k in range(2, n_cols+1)]
    
    susses = [fig.add_subplot(n_rows, n_cols, get_position(2, 1, n_cols))]
    susses += [fig.add_subplot(n_rows, n_cols,
                               get_position(2, k, n_cols),
                              sharex=susses[0],
                              sharey=susses[0])
               for k in range(2, n_cols+1)]

    final_sizes = [fig.add_subplot(n_rows, n_cols,
                                   get_position(3, 1, n_cols))]
    final_sizes += [fig.add_subplot(n_rows, n_cols,
                              get_position(3, k, n_cols),
                              sharex=final_sizes[0],
                              sharey=final_sizes[0])
                    for k in range(2, n_cols+1)]
    
    diffs = [fig.add_subplot(n_rows, n_cols, get_position(4, 1, n_cols))]
    diffs += [fig.add_subplot(n_rows, n_cols,
                              get_position(4, k, n_cols),
                              sharex=diffs[0],
                              sharey=diffs[0])
              for k in range(2, n_cols+1)]

    upper_axes = susses + dists
    lower_axes = final_sizes + diffs
    all_axes = susses + dists + final_sizes + diffs

    cmaps = [plt.cm.Blues, plt.cm.Blues, plt.cm.Greens, plt.cm.Greens]
    sus_models = ["linear", "multiplicative"] * 2

    for i, distribution in enumerate(example_dists):
        print(i)
        print(distribution)
        plot_susceptibility(
            n_classes=max_classes,
            axis=susses[i],
            cmap=cmaps[i],
            susceptibility_model=sus_models[i],
            homotypic_protection=homotypic_protection,
            broadly_neutralizing_sus=broadly_neutralizing_sus)
        plot_distribution(distribution,
                          axis=dists[i],
                          cmap=cmaps[i])
        plot_differences_by_R(distribution,
                              axes=[final_sizes[i],
                                    diffs[i]],
                              cmap=cmaps[i],
                              susceptibility_model=sus_models[i],
                              homotypic_protection=homotypic_protection,
                              lw=ps.lw,
                              by="Reff",
                              max_R=max_R,
                              plot_wt=True)

    # styling
    susses[0].set_xlim(left=-1.5, right=max_classes - 1)
    dists[0].set_xlim(left=-0.5, right=max_classes)
    susses[0].set_ylim(bottom=-0.05, top=1.05)

    stripped_num_fmt = mpl.ticker.FormatStrFormatter('%g')
    diffs[0].set_ylim(bottom=0, top=0.5)

    for axis in lower_axes:
        axis.set_xlim(left=1, right=max_R)
        axis.set_xticks(np.arange(1, max_R + 0.25, 0.25))
        axis.xaxis.set_major_formatter(stripped_num_fmt)

    for axis in all_axes:
        axis.grid(b=True)
        axis.yaxis.set_major_formatter(stripped_num_fmt)
        
    all_susses.set_ylabel("susceptibility to cluster 0")
    all_dists.set_xlabel("most recent cluster seen",
                         fontsize="xx-large")
    all_susses.set_xlabel("most recent cluster seen",
                         fontsize="xx-large")
    all_dists.set_ylabel("fraction of population")
    
    all_final_sizes.set_ylabel("attack rate")
    all_final_sizes.set_xlabel("mutant effective reproduction number ($\mathcal{R}_e$)",
                               fontsize="xx-large")
    all_diffs.set_xlabel("mutant effective reproduction number ($\mathcal{R}_e$)",
                         fontsize="xx-large")
    all_diffs.set_ylabel("$\Delta$ attack rate\n(new -- old)")

    fig.tight_layout(h_pad=0.1)
    fig.savefig(output_path)

if __name__ == "__main__":
    make_susceptibility_distribution_figure(
        output_path=sys.argv[1])
