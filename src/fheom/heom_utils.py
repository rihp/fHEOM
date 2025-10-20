"""
HEOM Simulation Utilities - Shared API for Quantum Biology Models

Provides unified interface for running HEOM simulations across all biological models.
Based on QuTiP's hierarchical equations of motion solver.
"""

import numpy as np
import qutip as qt
from typing import List, Optional, Tuple
import time
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.fft import fft, fftfreq


def cm_to_angular_freq(cm_inv: float) -> float:
    """Convert wavenumber (cm^-1) to angular frequency (rad/s)."""
    c = 2.99792458e10  # cm/s
    return 2 * np.pi * c * cm_inv


def run_heom_simulation(
    H: qt.Qobj,
    Q: qt.Qobj | list[qt.Qobj],
    lam_cm: float,
    gamma_cm: float,
    T_K: float,
    Nk: int,
    tlist: np.ndarray,
    rho0: qt.Qobj,
    e_ops: Optional[List[qt.Qobj]] = None,
    options: Optional[dict] = None
) -> qt.solver.Result:
    """
    Run HEOM simulation with Drude-Lorentz spectral density.
    
    Args:
        H: System Hamiltonian
        Q: Bath coupling operator (which system operator couples to bath)
        lam_cm: Reorganization energy (cm^-1)
        gamma_cm: Bath cutoff frequency (cm^-1)
        T_K: Temperature (Kelvin)
        Nk: Hierarchy depth
        tlist: Time points (seconds)
        rho0: Initial density matrix
        e_ops: List of expectation value operators
        
    Returns:
        QuTiP Result object with .expect, .states, etc.
        
    References:
        - Ishizaki & Fleming, J. Chem. Phys. 130, 234111 (2009)
        - Tanimura, J. Phys. Soc. Jpn. 75, 082001 (2006)
    """
    try:
        from qutip.solver.heom import HEOMSolver, DrudeLorentzBath
    except ImportError:
        raise ImportError(
            "HEOM solver not available. "
            "Requires QuTiP >= 5.0 with qutip.solver.heom module."
        )
    
    # Convert bath parameters to angular frequency
    lam_rad = cm_to_angular_freq(lam_cm)
    gamma_rad = cm_to_angular_freq(gamma_cm)
    
    # Convert temperature to energy units (rad/s)
    kB = 1.380649e-23  # J/K
    hbar = 1.054571817e-34  # J⋅s
    T_energy = kB * T_K / hbar  # rad/s
    
    # Create the bath or list of baths
    if isinstance(Q, list):
        baths = []
        for q_op in Q:
            baths.append(DrudeLorentzBath(
                q_op, 
                lam=lam_rad, 
                gamma=gamma_rad, 
                T=T_energy, 
                Nk=Nk
            ))
        bath_arg = baths
    else:
        bath_arg = DrudeLorentzBath(
            Q, 
            lam=lam_rad, 
            gamma=gamma_rad, 
            T=T_energy, 
            Nk=Nk
        )
    
    # Create HEOM solver
    solver = HEOMSolver(H, bath_arg, max_depth=Nk, options=options)
    
    # Run simulation
    result = solver.run(rho0, tlist, e_ops=e_ops)
    
    return result


def compute_t2_fft(trace: np.ndarray, t_fs: np.ndarray) -> Tuple[float, Optional[dict]]:
    """
    Calculate T2 coherence time by analyzing the frequency spectrum (FFT) of the trace.
    This method is more robust than curve fitting for noisy or complex signals.

    1. Performs a Fast Fourier Transform (FFT) on the population trace.
    2. Finds the dominant frequency peak corresponding to quantum beating.
    3. Calculates the Full Width at Half Maximum (FWHM) of this peak.
    4. T2 is calculated from the FWHM (T2 = 1 / (pi * FWHM)).
    """
    N = len(trace)
    if N < 2:
        return np.nan, None

    dt = (t_fs[1] - t_fs[0]) * 1e-15  # Time step in seconds
    
    # Detrend the signal and apply a Hann window for cleaner FFT
    signal = trace - np.mean(trace)
    window = np.hanning(N)
    signal_windowed = signal * window
    
    # FFT
    yf = fft(signal_windowed)
    xf = fftfreq(N, dt)[:N//2]
    
    # Power spectrum
    power = 2.0/N * np.abs(yf[0:N//2])
    
    # Find the main peak
    try:
        peaks, _ = find_peaks(power, height=0.01)
        if not len(peaks):
            return np.nan, None
        
        main_peak_idx = peaks[np.argmax(power[peaks])]
        
        # FWHM
        peak_power = power[main_peak_idx]
        half_max = peak_power / 2.0
        
        # Find where the power drops to half max
        left_idx = np.where(power[:main_peak_idx] < half_max)[0]
        right_idx = np.where(power[main_peak_idx:] < half_max)[0]
        
        if not len(left_idx) or not len(right_idx):
            return np.nan, None
            
        fwhm_left = xf[left_idx[-1]]
        fwhm_right = xf[main_peak_idx + right_idx[0]]
        fwhm_hz = fwhm_right - fwhm_left

        if fwhm_hz <= 0:
            return np.nan, None
            
        # T2 from FWHM of a Lorentzian peak
        t2_sec = 1 / (np.pi * fwhm_hz)
        t2_fs = t2_sec * 1e15

        analysis_data = {
            "peak_freq_hz": xf[main_peak_idx],
            "fwhm_hz": fwhm_hz
        }
        
        return t2_fs, analysis_data
        
    except Exception:
        return np.nan, None


def gaussian_decay_envelope(t, T2, a, omega, phi, c):
    """Model for fitting the Gaussian decay envelope of population oscillations."""
    return a * np.exp(-(t / T2)**2) * np.cos(omega * t + phi) + c


def compute_t2(trace: np.ndarray, t_fs: np.ndarray) -> float:
    """
    Calculate T2 coherence time. This now uses the robust FFT method.
    """
    return compute_t2_fft(trace, t_fs)


def compute_t2_threshold(trace: np.ndarray, t_fs: np.ndarray, threshold: float = 1/np.e) -> float:
    """Fallback: Calculate T2 coherence time (1/e decay)."""
    # Normalize the trace by its initial value
    normalized_trace = trace / trace[0]
    below = np.where(normalized_trace < threshold)[0]
    if len(below) > 0:
        return float(t_fs[below[0]])
    return np.nan


def save_result(result_dict: dict, path: str):
    """
    Save simulation result to JSON.
    
    Args:
        result_dict: Dictionary of results
        path: Output file path
    """
    import json
    from pathlib import Path
    
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(path, 'w') as f:
        json.dump(result_dict, f, indent=2)
    
    print(f"[✓] Saved {path.name} ({time.strftime('%H:%M:%S')})")


def print_result_summary(result_dict: dict):
    """Print formatted result summary."""
    print("\n" + "="*60)
    print(f"Model: {result_dict.get('model', 'Unknown')}")
    print("="*60)
    
    for key, value in result_dict.items():
        if key in ['model', 'time_fs', 'trace']:
            continue
        if isinstance(value, float):
            if 'runtime' in key:
                print(f"  {key}: {value:.2f}s")
            elif 'T2' in key or 'fs' in key:
                if np.isnan(value):
                    print(f"  {key}: >simulation_length")
                else:
                    print(f"  {key}: {value:.1f} fs")
            elif value < 0.01:
                print(f"  {key}: {value:.6f}")
            else:
                print(f"  {key}: {value:.3f}")
        else:
            print(f"  {key}: {value}")
    
    print("="*60)

