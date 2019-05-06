#!/usr/bin/env python3

import numpy as np

final_size_diff_escape_factors = [0.1, 0.2, 0.3]
final_size_diff_immune_susceptibilities = [0.0, 0.1, 0.2]
shared_sus = .6
R0 = 2

example_sus_distributions = [
    np.array([0.0, 0.1, 0.2, 0.2, 0.2, 0.2, 0.1, 0]),
    np.array([0.0, 0.1, 0.2, 0.2, 0.2, 0.2, 0.1, 0]),
    np.array([0.1, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]),
    np.array([0.1, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    ]


broadly_neutralizing_sus = 0.5
default_homotypic_protection = 1
