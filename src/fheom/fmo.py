"""
FMO Complex (Fenna-Matthews-Olson) - Photosynthetic Energy Transfer

Experimentally validated quantum coherence in biological photosynthesis.
7-site excitonic model based on Adolphs & Renger (2006) structure.

References:
    - Engel et al., Nature 446, 782-786 (2007) - Experimental evidence
    - Ishizaki & Fleming, PNAS 106, 17255 (2009) - HEOM simulation
    - Adolphs & Renger, Biophys. J. 91, 2778 (2006) - Hamiltonian
"""

import numpy as np
from qutip import Qobj
import qutip as qt
from .heom_utils import cm_to_angular_freq


def build_fmo_hamiltonian():
    """
    Builds the 7-site FMO Hamiltonian using the parameters from Adolphs & Renger (2006).
    This is the primary, cited Hamiltonian for the paper's validation.

    Source: Adolphs, J., & Renger, T. (2006). How proteins trigger excitation
    energy transfer in the FMO complex of green sulfur bacteria. Biophysical
    Journal, 91(8), 2778–2797. DOI: 10.1529/biophysj.105.079483
    
    Returns Hamiltonian in rad/s (compatible with time in seconds).
    """
    # Convert cm^-1 to rad/s: ω = 2πcν where c = 2.998e10 cm/s
    cm_to_rad_s = 2 * np.pi * 2.998e10  # rad/s per cm^-1

    # Site energies in cm^-1, from Table 1 in Adolphs & Renger (2006)
    site_energies_cm = np.array([
        12410, 12530, 12210, 12320, 12480, 12620, 12440
    ])

    # Coupling matrix in cm^-1, from Table 1 in Adolphs & Renger (2006)
    couplings_cm = np.array([
        [0,    -87.7,   5.5,  -5.9,   6.7, -13.7,  -9.9],
        [-87.7,    0,  30.8,   8.2,   0.7,  11.8,   4.3],
        [  5.5, 30.8,     0, -53.5,  -2.2,  -9.6,   6.0],
        [ -5.9,  8.2, -53.5,     0, -70.7, -17.0, -63.3],
        [  6.7,  0.7,  -2.2, -70.7,     0,  81.1,  -1.3],
        [-13.7, 11.8,  -9.6, -17.0,  81.1,     0,  39.7],
        [ -9.9,  4.3,   6.0, -63.3,  -1.3,  39.7,     0]
    ])

    hamiltonian_cm = np.diag(site_energies_cm) + couplings_cm
    hamiltonian_rad_s = hamiltonian_cm * cm_to_rad_s
    hamiltonian = qt.Qobj(hamiltonian_rad_s)
    return hamiltonian


def bath_operator() -> qt.Qobj:
    """
    Bath coupling operator for FMO.
    
    Couples to local site populations (dephasing noise).
    Sum over all site projection operators.
    
    Returns:
        7x7 QuTiP Qobj for bath coupling
    """
    Q = qt.Qobj(np.zeros((7, 7), dtype=complex))
    
    # Couple to each site population
    for i in range(7):
        Q += qt.basis(7, i) * qt.basis(7, i).dag()
    
    return Q


def default_initial_state() -> qt.Qobj:
    """
    Default initial state for FMO.
    
    Site 0 (BChl 1) excited, all others ground.
    Represents initial excitation from baseplate.
    
    Returns:
        7x7 density matrix
    """
    return qt.ket2dm(qt.basis(7, 0))


def observable_site0() -> qt.Qobj:
    """Observable for site 0 population."""
    return qt.basis(7, 0) * qt.basis(7, 0).dag()


def get_default_params() -> dict:
    """
    Get default FMO simulation parameters.
    
    Based on Ishizaki & Fleming (2009) "fast bath" model.
    
    Returns:
        Dictionary of standard parameters
    """
    return {
        'lam_cm': 35.0,      # Reorganization energy (cm^-1)
        'gamma_cm': 667.0,   # Bath cutoff (cm^-1), τ_c ≈ 8 fs
        'T_K': 77.0,         # Liquid nitrogen temperature (K)
        'Nk': 5,             # Hierarchy depth
        't_max_fs': 500.0    # Simulation time (fs)
    }


def get_experimental_T2() -> tuple:
    """
    Get experimental T2 range from literature.
    
    Returns:
        (T2_low, T2_high) in femtoseconds
        
    References:
        - Engel et al., Nature 446, 782 (2007): ~60-80 fs at 77K
    """
    return (60.0, 80.0)


def get_site_coordinates():
    """Returns the 3D coordinates of the FMO sites in Angstroms."""
    return np.array([
        [28.3, 18.3, 19.3],
        [32.4, 29.5, 13.6],
        [35.1, 21.3, 5.1],
        [22.4, 25.5, 2.7],
        [16.8, 16.3, 10.5],
        [11.9, 25.1, 14.5],
        [20.4, 32.7, 21.3]
    ])

def get_correlated_bath_q_ops(correlation_length=10.0):
    """
    Constructs bath coupling operators for a spatially correlated environment.
    The correlation is assumed to decay exponentially with distance.

    Args:
        correlation_length (float): The characteristic length (in Angstroms)
                                    over which environmental correlations decay.

    Returns:
        A list of QuTiP Qobj operators, one for each correlated bath mode.
    """
    coords = get_site_coordinates()
    num_sites = coords.shape[0]
    
    # Calculate pairwise distance matrix
    dist_matrix = np.linalg.norm(coords[:, np.newaxis, :] - coords, axis=2)
    
    # Construct spatial correlation matrix
    correlation_matrix = np.exp(-dist_matrix / correlation_length)
    
    # Find eigenvectors of the correlation matrix
    eigenvalues, eigenvectors = np.linalg.eigh(correlation_matrix)
    
    # The eigenvectors are the collective modes. Create a Qobj for each.
    q_ops = []
    for i in range(num_sites):
        op = qt.Qobj(np.zeros((num_sites, num_sites), dtype=float))
        for j in range(num_sites):
            # The operator for mode i is a sum of projectors weighted by the eigenvector components
            op += eigenvectors[j, i] * qt.basis(num_sites, j) * qt.basis(num_sites, j).dag()
        q_ops.append(op)
        
    return q_ops


def build_fmo_hamiltonian_ishizaki():
    """
    Builds the 7-site FMO Hamiltonian using the parameters from Ishizaki & Fleming (2009).
    This is a widely respected alternative to the Adolphs & Renger (2006) model.

    Source: Ishizaki, A., & Fleming, G. R. (2009). Unified treatment of quantum
    coherent and incoherent hopping dynamics in electronic energy transfer:
    Reduced hierarchy equation approach. The Journal of Chemical Physics, 130(23), 234111.
    DOI: 10.1063/1.3155372
    """
    cm_to_rad_ps = 2 * np.pi * 29979245800 * 1e-12

    # Site energies in cm^-1
    site_energies_cm = np.array([
        12410, 12530, 12210, 12320, 12480, 12620, 12440
    ])

    # Add the standard 210 cm^-1 offset from the paper
    site_energies_cm = site_energies_cm - 210

    # Coupling matrix in cm^-1
    couplings_cm = np.array([
        [0,    -87.7,   5.5,  -5.9,   6.7, -13.7,  -9.9],
        [-87.7,    0,  30.8,   8.2,   0.7,  11.8,   4.3],
        [  5.5, 30.8,     0,  -8.2,   0.9,  -7.0,  -3.3],
        [ -5.9,  8.2,  -8.2,     0, -62.3,  -9.6,   6.0],
        [  6.7,  0.7,   0.9, -62.3,     0,  81.1,  -1.3],
        [-13.7, 11.8,  -7.0,  -9.6,  81.1,     0,  -2.2],
        [ -9.9,  4.3,  -3.3,   6.0,  -1.3,  -2.2,     0]
    ])

    hamiltonian_cm = np.diag(site_energies_cm) + couplings_cm
    hamiltonian = qt.Qobj(hamiltonian_cm * cm_to_rad_ps)
    return hamiltonian


def build_fmo_hamiltonian_original_code():
    """
    Builds the 7-site FMO Hamiltonian using the parameters
    originally present in the codebase. These differ slightly from
    the cited Adolphs & Renger (2006) paper.
    """
    cm_to_rad_ps = 2 * np.pi * 29979245800 * 1e-12

    # Site energies in cm^-1
    site_energies_cm = np.array([12410, 12530, 12210, 12320, 12480, 12630, 12440])

    # Coupling matrix in cm^-1
    couplings_cm = np.array([
        [0, -87.7, 5.5, -5.9, 6.7, -13.7, -9.9],
        [-87.7, 0, 30.8, 8.2, 0.7, 11.8, 4.3],
        [5.5, 30.8, 0, -53.5, -2.2, -9.6, 6.0],
        [-5.9, 8.2, -53.5, 0, -70.7, -17.0, -63.3],
        [6.7, 0.7, -2.2, -70.7, 0, 81.1, -1.3],
        [-13.7, 11.8, -9.6, -17.0, 81.1, 0, -2.2],
        [-9.9, 4.3, 6.0, -63.3, -1.3, -2.2, 0]
    ])
    hamiltonian_cm = np.diag(site_energies_cm) + couplings_cm
    hamiltonian = Qobj(hamiltonian_cm * cm_to_rad_ps)
    return hamiltonian

