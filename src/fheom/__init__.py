"""
fHEOM - Factorized Hierarchical Equations of Motion

A novel numerical method for handling spatially correlated environmental baths
in HEOM simulations through low-rank decomposition.

Key Innovation:
    Reduces HEOM hierarchy size from (Nk+1)^N to (Nk+1)^r modes
    through low-rank factorization of spatial correlation matrices.

Author: Roberto Ignacio Henriquez-Perozo
License: MIT
"""

from .fheom import (
    FactorizationResult,
    spatial_correlation_matrix,
    factorize_correlation_matrix,
    construct_factorized_bath_operators,
    get_factorized_bath,
    compare_full_vs_factorized,
    fheom,
    factorize,
    spatial_corr,
)

from .heom_utils import (
    run_heom_simulation,
    compute_t2,
    compute_t2_fft,
    cm_to_angular_freq,
)

__version__ = "1.0.0"
__author__ = "Roberto Ignacio Henriquez-Perozo"
__email__ = "roberto@henriquezperozo.com"

__all__ = [
    'FactorizationResult',
    'spatial_correlation_matrix',
    'factorize_correlation_matrix', 
    'construct_factorized_bath_operators',
    'get_factorized_bath',
    'compare_full_vs_factorized',
    'fheom',
    'factorize',
    'spatial_corr',
    'run_heom_simulation',
    'compute_t2',
    'compute_t2_fft',
    'cm_to_angular_freq',
]
