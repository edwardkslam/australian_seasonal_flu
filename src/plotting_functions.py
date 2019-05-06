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

