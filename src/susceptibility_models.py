#!/usr/bin/env python3

# filename: susceptibility_models.py
# author: Dylan Morris <dhmorris@princeton.edu>
#
# description: implements susceptibility models
# as a function of antigenic distance

import numpy as np
from scipy.optimize import fsolve

def linear_susceptibility(distance, escape, homotypic_protection=1):
    return min(1, (1 - homotypic_protection + distance * escape))

def multiplicative_susceptibility(distance,
                                  escape,
                                  homotypic_protection=1):
    return 1 - (((1 - escape)**distance) * homotypic_protection)

def pick_sus_func(susceptibility_model, escape,
                  homotypic_protection=1):
    sus_funcs = {
        "multiplicative": multiplicative_susceptibility,
        "linear": linear_susceptibility
    }
    def sus_func(dist):
        return sus_funcs[susceptibility_model](dist, escape, homotypic_protection=homotypic_protection)
        
    return sus_func
