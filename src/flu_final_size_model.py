import numpy as np
from final_size import calc_n_class_S_final
from susceptibility_models import pick_sus_func

def flu_final_size(susceptibilities,
                  class_distribution,
                  R0,
                  init_inf=None):
    if init_inf is None:
        init_inf = np.zeros_like(class_distribution)
        init_inf[0] = 0.00001
    susceptibilities = np.array(susceptibilities)
    class_distribution = np.array(class_distribution)
    n_classes = susceptibilities.size
    imm_mat = np.diag(susceptibilities)
    R0_mat = R0 * np.ones_like(imm_mat)
    beta_mat = imm_mat @ R0_mat
    class_distribution = class_distribution / np.sum(class_distribution)
    s_finals = calc_n_class_S_final(beta_mat,
                                    class_distribution,
                                    init_inf)
    r_finals = class_distribution - s_finals
    return (s_finals, r_finals)


def calc_Reff_simple(susceptibilities,
                     class_distribution,
                     R0):
    susceptibilities = np.array(susceptibilities)
    class_distribution = np.array(class_distribution) / np.sum(class_distribution)
    return np.sum(R0 * susceptibilities * class_distribution)

def get_dist_susses(
        variant,
        distribution,
        escape_factor,
        susceptibility_model="linear",
        homotypic_protection=1,
        broadly_neutralizing_sus=1):
    """
    For the paraticular parametrization of 
    susceptiblity classes used in this model, namely:
    [removed, fully naive, memory to wt (0th cluster), 
    memory to 1st cluster ... memory to kth previous cluster,
    has broadly neutralizing antibodies]
    return the susceptibilities to variant 0 (wild type)
    variant -1 (mutant), or further variants (later mutants)
    """
    sus_func = pick_sus_func(susceptibility_model,
                             escape_factor,
                             homotypic_protection=homotypic_protection)

    distribution = np.array(distribution)
    n_non_B_cell_classes = 3 # removed, naive, broadly neutralizing
    n_B_cell_classes = distribution.size - n_non_B_cell_classes

    return np.array(
        [0, 1] +
        [sus_func(abs(variant - k)) for k in range(n_B_cell_classes)] +
        [broadly_neutralizing_sus])


def calc_R0_from_Reff(susceptibilities,
                      class_distribution,
                      desired_Reff):
    susceptibilities = np.array(susceptibilities)
    class_distribution = np.array(class_distribution) / np.sum(class_distribution)

    potential = np.sum(susceptibilities * class_distribution)
    R0 = desired_Reff / potential
    return R0

def compare_final_sizes(distribution, escape_factor, R0,
                        susceptibility_model="linear",
                        homotypic_protection=1,
                        sus_dist_0=None,
                        sus_dist_1=None,
                        broadly_neutralizing_sus=1):
    
    if sus_dist_0 is None:
        sus_dist_0 = get_dist_susses(
            0,
            distribution,
            escape_factor,
            susceptibility_model=susceptibility_model,
            homotypic_protection=homotypic_protection,
            broadly_neutralizing_sus=broadly_neutralizing_sus)
    if sus_dist_1 is None:
        sus_dist_1 = get_dist_susses(
            -1,
            distribution,
            escape_factor,
            susceptibility_model=susceptibility_model,
            homotypic_protection=homotypic_protection,
            broadly_neutralizing_sus=broadly_neutralizing_sus)


    #print(sus_dist_0, sus_dist_1)
    _, final_size_0 = flu_final_size(sus_dist_0,
                                     distribution,
                                     R0)
    _, final_size_1 = flu_final_size(sus_dist_1,
                                     distribution,
                                     R0)
    diffs = final_size_1 - final_size_0
    total_diff = np.sum(diffs)
    return (final_size_0, final_size_1, diffs, total_diff)
