import numpy as np
import matplotlib.pyplot as plt
import plotting_style as ps

def add_bounding_subplot(figure, position=None):
    if position is None:
        position = 111
    ax = figure.add_subplot(position)
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w',
                   top=False,
                   bottom=False,
                   left=False,
                   right=False)
    return ax



def get_position(row,
                 col,
                 n_cols):
    return (row - 1) * n_cols + col

def setup_multipanel(fig, plot_positions,
                     gridspec=None,
                     letters=True,
                     letter_loc=None):

    if gridspec is None:
        n_rows = np.max([val["row"] for val in plot_positions.values()])
        n_cols = np.max([val["col"] for val in plot_positions.values()])
        
    plots = {}

    letters = ["a", "b", "c", "d", "e", "f",
               "g", "h", "i", "j", "k", "l",
               "m", "n", "o", "p", "q", "r",
               "s", "t", "u", "v", "w", "x",
               "y", "z", "aa", "bb", "cc"]

    iteration = 0
    for plot_identity, position_dict in plot_positions.items():

        sharex = position_dict["sharex"]
        sharey = position_dict["sharey"]

        if sharex is not None:
            plot_sharex = plots[sharex]
        else:
            plot_sharex = sharex
            
        if sharey is not None:
            plot_sharey = plots[sharey]
        else:
            plot_sharey = sharey

        if gridspec is None:
            position_id = get_position(position_dict["row"],
                                       position_dict["col"],
                                       n_cols)
            ax = fig.add_subplot(n_rows,
                                 n_cols,
                                 position_id,
                                 sharex=plot_sharex,
                                 sharey=plot_sharey)
        else:
            grid_pos = position_dict["grid_position"]
            ax = fig.add_subplot(gridspec[grid_pos],
                                 sharex=plot_sharex,
                                 sharey=plot_sharey)
        if letters:
            letter_locations = position_dict.get("letter_loc", letter_loc) 
            if letter_locations is None:
                letter_x, letter_y = ps.letter_loc
            else:
                letter_x, letter_y = letter_locations
            ax.text(letter_x, letter_y, letters[iteration],
                    transform=ax.transAxes,
                    fontsize=ps.letter_size,
                    fontweight='bold', va='top')
        
        plots[plot_identity] = ax
        iteration += 1

    return plots
