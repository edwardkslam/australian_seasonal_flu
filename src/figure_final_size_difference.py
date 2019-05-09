#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import plotting_style as ps
import model_parameters as mp
import sys

from flu_final_size_model import flu_final_size, compare_final_sizes



def get_sus_shared(pct_naive,
                   pct_wt_imm,
                   mutant_escape_from_wt,
                   desired_avg_sus):
    numer = pct_naive + mutant_escape_from_wt * pct_wt_imm - desired_avg_sus 
    denom = pct_naive + pct_wt_imm - 1
    result = numer / denom
    if result < 0 or result > 1:
        print("no choice of shared susceptibility"
              "can achieve desired average")
        return None
    else:
        return result

def plot_differences_by_wt_imm(escape_factors=[0.1, 0.2, 0.3],
                               shared_sus=1,
                               imm_sus=0,
                               R0=1.8,
                               cmap=plt.cm.Greens,
                               axes=None,
                               **kwargs):
    if axes is None:
        fig, axes = plt.subplots(3, 1, sharex=True, sharey=True)

    color_vals = np.linspace(0.4, 0.8, len(escape_factors))

    wt_imms = np.linspace(0, 1, 100)
    distributions = [[wt_imm, 1 - wt_imm]
                     for wt_imm in wt_imms]
    for k_escape, escape_factor in enumerate(escape_factors):
        sus_0 = np.array([imm_sus + escape_factor, shared_sus])
        sus_1 = np.array([imm_sus, shared_sus])

        final_size_0s = np.array(
            [np.sum(compare_final_sizes(distribution,
                                        escape_factor,
                                        R0,
                                        sus_dist_0=sus_0,
                                        sus_dist_1=sus_1)[0])
             for distribution in distributions])
        final_size_1s = np.array(
            [np.sum(compare_final_sizes(distribution,
                                        escape_factor,
                                        R0,
                                        sus_dist_0=sus_0,
                                        sus_dist_1=sus_1)[1])
             for distribution in distributions])
        
        axes[0].plot(wt_imms, final_size_0s,
                     label="{:.0f}\%".format(escape_factor * 100),
                     color=cmap(color_vals[k_escape]),
                     **kwargs)
        axes[1].plot(wt_imms, final_size_1s,
                     label="{:.0f}\%".format(escape_factor * 100),
                     color=cmap(color_vals[k_escape]),
                     **kwargs)
        axes[2].plot(wt_imms, final_size_0s - final_size_1s,
                     label="{:.0f}\%".format(escape_factor * 100),
                     color=cmap(color_vals[k_escape]),
                     **kwargs)
    axes[0].legend(title="escape",
                   frameon=True,
                   fancybox=True,
                   framealpha=1)


def plot_final_size_difference(shared_sus=mp.shared_sus,
                               R0=mp.R0,
                               output_path="../out/final_size_diff.pdf",
                               **kwargs):
    escapes = np.array(mp.final_size_diff_escape_factors)
    imm_sus = mp.final_size_diff_immune_susceptibilities
    cmaps = [plt.cm.Greens, plt.cm.Blues, plt.cm.Reds]

    figsize=(9, 9)
    
    fig, axes = plt.subplots(3, 3, sharex=True, sharey=True,
                             figsize=figsize)

    for i, imm, cmap in zip(range(3), imm_sus, cmaps):
        ax_slice = axes[:, i]
        plot_differences_by_wt_imm(escape_factors=escapes,
                                   shared_sus=shared_sus,
                                   imm_sus=imm,
                                   R0=R0,
                                   cmap=cmap,
                                   axes=ax_slice,
                                   **kwargs)
    xticks = [0, 0.2, 0.4, 0.6, 0.8, 1]
    for ax in axes.flatten():
        ax.grid(b=True)
        ax.set_ylim(bottom=0, top=0.35)
        ax.set_xlim(left=0, right=1)
        ax.set_xticks(xticks)
        ax.label_outer()
        
    axes[2][1].set_xlabel("percentage of strongly immune hosts")
    axes[0][0].set_ylabel("new variant\nattack rate ")
    axes[1][0].set_ylabel("old variant\nattack rate")
    axes[2][0].set_ylabel("$\Delta$ attack rate\n(new - old)")
    axes[0][0].set_title("strongly immune hosts are\n{:.0f}% susceptible"
                       "".format(imm_sus[0]*100))
    axes[0][1].set_title("strongly immune hosts are\n{:.0f}% susceptible"
                       "".format(imm_sus[1]*100))
    axes[0][2].set_title("strongly immune hosts are\n{:.0f}% susceptible"
                       "".format(imm_sus[2]*100))
    fig.savefig(output_path)


if __name__ == "__main__":
    plot_final_size_difference(output_path=sys.argv[1],
                               lw=5)
