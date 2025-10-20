"""
Factorized HEOM (fHEOM) — Low-rank decomposition of correlated baths.

This module implements a numerical method for handling spatially correlated
environmental baths in HEOM simulations. By decomposing the correlation matrix
into a low-rank representation, we reduce the number of effective bath modes from
N_sites to r ≪ N_sites, reducing computational scaling from exponential to polynomial.

Key Innovation:
    For a correlation matrix C (N×N), we compute:
        C ≈ L L^T  (Cholesky or eigendecomposition)
    
    Then construct r effective bath operators:
        Q_eff[k] = sum_i L[i,k] * Q_site[i]
    
    Result: HEOM hierarchy size (Nk+1)^r instead of (Nk+1)^N.

Reference:
    Low-rank factorization method for correlated baths in HEOM simulations
    of photosynthetic quantum coherence systems. See README.md for details.
"""
import numpy as np
import qutip as qt
from typing import List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class FactorizationResult:
    """Result of correlation matrix factorization."""
    q_eff: List[qt.Qobj]  # Effective bath operators
    weights: np.ndarray    # Rescaling factors for each mode
    rank: int              # Number of effective modes retained
    reconstruction_error: float  # Relative Frobenius norm error
    eigenvalues: np.ndarray  # All eigenvalues (for diagnostics)
    total_variance: float    # Sum of all eigenvalues
    explained_variance: float  # Fraction of variance retained


def spatial_correlation_matrix(
    coordinates: np.ndarray,
    correlation_length: float = 10.0,
    kernel: str = 'exponential'
) -> np.ndarray:
    """
    Construct spatial correlation matrix from site coordinates.
    
    Args:
        coordinates: N×3 array of site positions (Angstroms)
        correlation_length: Decay length scale (Angstroms)
        kernel: 'exponential', 'gaussian', or 'power_law'
    
    Returns:
        C: N×N symmetric correlation matrix (positive semidefinite)
    
    Examples:
        >>> coords = np.array([[0,0,0], [5,0,0], [10,0,0]])  # 3 sites in line
        >>> C = spatial_correlation_matrix(coords, correlation_length=5.0)
        >>> C[0,0]  # Self-correlation
        1.0
        >>> 0.3 < C[0,1] < 0.4  # Nearest neighbor
        True
        >>> C[0,2] < C[0,1]  # Further site weaker
        True
    """
    N = len(coordinates)
    C = np.zeros((N, N))
    
    for i in range(N):
        for j in range(N):
            r_ij = np.linalg.norm(coordinates[i] - coordinates[j])
            
            if kernel == 'exponential':
                C[i, j] = np.exp(-r_ij / correlation_length)
            elif kernel == 'gaussian':
                C[i, j] = np.exp(-(r_ij / correlation_length)**2)
            elif kernel == 'power_law':
                C[i, j] = 1.0 / (1.0 + (r_ij / correlation_length)**2)
            else:
                raise ValueError(f"Unknown kernel: {kernel}")
    
    # Ensure symmetry (numerical precision)
    C = 0.5 * (C + C.T)
    
    return C


def factorize_correlation_matrix(
    corr_matrix: np.ndarray,
    rank: Optional[int] = None,
    variance_threshold: float = 0.99,
    method: str = 'eigendecomposition'
) -> Tuple[np.ndarray, np.ndarray, dict]:
    """
    Factorize correlation matrix into low-rank form.
    
    Args:
        corr_matrix: N×N symmetric positive semidefinite correlation matrix
        rank: Number of modes to retain (None = auto from variance_threshold)
        variance_threshold: Retain modes capturing this fraction of variance
        method: 'eigendecomposition' or 'cholesky'
    
    Returns:
        L: N×r factorization matrix (C ≈ L L^T)
        eigenvalues: r eigenvalues (for weighting bath parameters)
        info: dict with diagnostics (rank, error, explained_variance)
    
    Notes:
        - Eigendecomposition: Optimal in L2 sense, always works
        - Cholesky: Faster, but requires C strictly positive definite
    """
    N = corr_matrix.shape[0]
    
    if method == 'eigendecomposition':
        # Eigen-decomposition (stable for semidefinite matrices)
        eigvals, eigvecs = np.linalg.eigh(corr_matrix)
        
        # Sort descending
        idx = np.argsort(eigvals)[::-1]
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]
        
        # Filter out negative eigenvalues (numerical noise)
        eigvals = np.maximum(eigvals, 0.0)
        
        # Auto-select rank
        total_variance = np.sum(eigvals)
        if rank is None:
            cumsum = np.cumsum(eigvals) / total_variance
            rank = int(np.searchsorted(cumsum, variance_threshold)) + 1
            rank = min(rank, N)  # Cap at full rank
        
        # Construct low-rank factorization
        L = eigvecs[:, :rank] @ np.diag(np.sqrt(eigvals[:rank]))
        retained_eigvals = eigvals[:rank]
        
    elif method == 'cholesky':
        # Cholesky decomposition (faster but less robust)
        try:
            L_full = np.linalg.cholesky(corr_matrix)  # C = L L^T
            
            # SVD to get rank reduction
            U, s, Vh = np.linalg.svd(L_full, full_matrices=False)
            
            total_variance = np.sum(s**2)
            if rank is None:
                cumsum = np.cumsum(s**2) / total_variance
                rank = int(np.searchsorted(cumsum, variance_threshold)) + 1
            
            L = U[:, :rank] @ np.diag(s[:rank])
            retained_eigvals = s[:rank]**2
            eigvals = s**2  # For diagnostics
            
        except np.linalg.LinAlgError:
            # Fall back to eigendecomposition
            return factorize_correlation_matrix(
                corr_matrix, rank, variance_threshold, method='eigendecomposition'
            )
    else:
        raise ValueError(f"Unknown method: {method}")
    
    # Compute reconstruction error
    C_approx = L @ L.T
    error = np.linalg.norm(corr_matrix - C_approx, 'fro') / np.linalg.norm(corr_matrix, 'fro')
    
    explained_variance = np.sum(retained_eigvals) / total_variance
    
    info = {
        'rank': rank,
        'reconstruction_error': error,
        'explained_variance': explained_variance,
        'total_variance': total_variance,
        'eigenvalues': eigvals,
    }
    
    return L, retained_eigvals, info


def construct_factorized_bath_operators(
    q_site_ops: List[qt.Qobj],
    factorization_matrix: np.ndarray,
    normalize: bool = True
) -> List[qt.Qobj]:
    """
    Construct effective bath operators from factorization.
    
    Args:
        q_site_ops: List of N site bath operators [Q_0, Q_1, ..., Q_{N-1}]
        factorization_matrix: N×r matrix from factorize_correlation_matrix
        normalize: Whether to normalize each effective operator
    
    Returns:
        q_eff: List of r effective bath operators
    
    Math:
        Q_eff[k] = sum_{i=0}^{N-1} L[i,k] * Q_site[i]
    
    Example:
        >>> # 3 sites, rank-2 factorization
        >>> q_sites = [qt.sigmaz() for _ in range(3)]
        >>> L = np.array([[0.8, 0.1], [0.6, 0.3], [0.4, 0.5]])
        >>> q_eff = construct_factorized_bath_operators(q_sites, L)
        >>> len(q_eff)
        2
    """
    N, rank = factorization_matrix.shape
    
    if len(q_site_ops) != N:
        raise ValueError(f"Expected {N} site operators, got {len(q_site_ops)}")
    
    q_eff = []
    
    for k in range(rank):
        # Linear combination of site operators
        op_k = sum(
            factorization_matrix[i, k] * q_site_ops[i]
            for i in range(N)
        )
        
        if normalize:
            # Normalize to preserve bath strength
            norm = op_k.norm()
            if norm > 1e-12:
                op_k = op_k / norm
        
        q_eff.append(op_k)
    
    return q_eff


def get_factorized_bath(
    q_site_ops: List[qt.Qobj],
    coordinates: np.ndarray,
    rank: Optional[int] = None,
    correlation_length: float = 10.0,
    variance_threshold: float = 0.99,
    kernel: str = 'exponential',
    method: str = 'eigendecomposition'
) -> FactorizationResult:
    """
    One-shot function: Construct factorized bath from site operators and coordinates.
    
    Args:
        q_site_ops: List of N site bath operators
        coordinates: N×3 array of site positions (Angstroms)
        rank: Number of effective modes (None = auto)
        correlation_length: Spatial decay length (Angstroms)
        variance_threshold: Fraction of variance to retain (if rank=None)
        kernel: Correlation kernel ('exponential', 'gaussian', 'power_law')
        method: Factorization method ('eigendecomposition', 'cholesky')
    
    Returns:
        FactorizationResult with q_eff, weights, and diagnostics
    
    Example:
        >>> from fheom.fmo import get_site_coordinates, bath_operator
        >>> coords = get_site_coordinates()  # FMO 7-site
        >>> q_sites = [bath_operator() for _ in range(7)]  # Use same operator for each site
        >>> result = get_factorized_bath(q_sites, coords, rank=3)
        >>> len(result.q_eff)
        3
        >>> result.explained_variance > 0.95
        True
    """
    N = len(q_site_ops)
    
    # 1. Construct correlation matrix
    C = spatial_correlation_matrix(coordinates, correlation_length, kernel)
    
    # 2. Factorize
    L, eigvals, info = factorize_correlation_matrix(C, rank, variance_threshold, method)
    
    # 3. Construct effective operators
    q_eff = construct_factorized_bath_operators(q_site_ops, L, normalize=False)
    
    # 4. Weights for bath parameter rescaling
    # Each mode gets a weight proportional to sqrt(eigenvalue)
    weights = np.sqrt(eigvals)
    
    return FactorizationResult(
        q_eff=q_eff,
        weights=weights,
        rank=info['rank'],
        reconstruction_error=info['reconstruction_error'],
        eigenvalues=info['eigenvalues'],
        total_variance=info['total_variance'],
        explained_variance=info['explained_variance']
    )


def compare_full_vs_factorized(
    H: qt.Qobj,
    q_site_ops: List[qt.Qobj],
    coordinates: np.ndarray,
    lam_cm: float,
    gamma_cm: float,
    T_K: float,
    Nk: int,
    tlist: np.ndarray,
    rho0: qt.Qobj,
    e_ops: List[qt.Qobj],
    rank_values: List[int] = [2, 3, 4],
    correlation_length: float = 10.0
) -> dict:
    """
    Benchmark: Compare full correlated bath vs fHEOM at different ranks.
    
    Args:
        H, q_site_ops, lam_cm, gamma_cm, T_K, Nk, tlist, rho0, e_ops: Standard HEOM params
        coordinates: N×3 site positions
        rank_values: List of ranks to test
        correlation_length: Spatial correlation decay length
    
    Returns:
        results: dict with keys:
            'full': Full simulation result
            'factorized': dict mapping rank → (result, error, speedup, memory_reduction)
    
    This function is used to generate Figure 7 in the paper (fHEOM performance).
    """
    from fheom.heom_utils import run_heom_simulation
    from fheom.fmo import get_correlated_bath_q_ops
    import time
    
    results = {}
    
    # 1. Run full correlated bath (reference)
    print("Running full correlated bath (reference)...")
    t0 = time.time()
    
    q_full = get_correlated_bath_q_ops()  # Full N-mode correlation
    
    try:
        res_full = run_heom_simulation(
            H, q_full, lam_cm, gamma_cm, T_K, Nk, tlist, rho0, e_ops=e_ops
        )
        runtime_full = time.time() - t0
        trace_full = res_full.expect[0]  # Site 0 population
        
        results['full'] = {
            'result': res_full,
            'runtime': runtime_full,
            'trace': trace_full,
            'n_modes': len(q_full),
        }
        print(f"  ✓ Full: {len(q_full)} modes, {runtime_full:.1f}s")
        
    except Exception as e:
        print(f"  ✗ Full simulation failed: {e}")
        results['full'] = None
        trace_full = None
    
    # 2. Run factorized versions
    results['factorized'] = {}
    
    for rank in rank_values:
        print(f"\nRunning fHEOM with rank={rank}...")
        t0 = time.time()
        
        # Get factorized bath
        fheom = get_factorized_bath(
            q_site_ops, coordinates, rank=rank,
            correlation_length=correlation_length
        )
        
        try:
            res_fact = run_heom_simulation(
                H, fheom.q_eff, lam_cm, gamma_cm, T_K, Nk, tlist, rho0, e_ops=e_ops
            )
            runtime_fact = time.time() - t0
            trace_fact = res_fact.expect[0]
            
            # Compute error vs full
            if trace_full is not None:
                error = np.linalg.norm(trace_fact - trace_full) / np.linalg.norm(trace_full)
            else:
                error = np.nan
            
            # Speedup and memory estimates
            speedup = runtime_full / runtime_fact if runtime_full else np.nan
            
            # Memory scales as (Nk+1)^n_modes approximately
            n_full = len(q_full) if results['full'] else 7
            memory_reduction = (Nk + 1)**n_full / (Nk + 1)**rank
            
            results['factorized'][rank] = {
                'result': res_fact,
                'runtime': runtime_fact,
                'trace': trace_fact,
                'error': error,
                'speedup': speedup,
                'memory_reduction': memory_reduction,
                'explained_variance': fheom.explained_variance,
                'reconstruction_error': fheom.reconstruction_error,
            }
            
            print(f"  ✓ Rank {rank}: {runtime_fact:.1f}s, "
                  f"error={error:.2%}, speedup={speedup:.1f}×")
            
        except Exception as e:
            print(f"  ✗ Rank {rank} failed: {e}")
            results['factorized'][rank] = None
    
    return results


# Convenience aliases
fheom = get_factorized_bath
factorize = factorize_correlation_matrix
spatial_corr = spatial_correlation_matrix

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
]





