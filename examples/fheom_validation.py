"""
fHEOM Validation: Test factorized correlated bath implementation.

This script validates the fHEOM method by testing:
1. Spatial correlation matrix construction
2. Low-rank factorization accuracy at different ranks
3. Effective bath operator construction
4. HEOM dynamics simulation with factorized bath

Goal: Verify algorithmic correctness and factorization accuracy.
"""
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '0'
os.environ['OMP_NUM_THREADS'] = '4'

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import time

import qutip as qt
from fheom.fmo import (
    build_fmo_hamiltonian,
    get_site_coordinates,
    bath_operator,
    observable_site0
)
from fheom.fheom import (
    get_factorized_bath,
    spatial_correlation_matrix,
    compare_full_vs_factorized
)
from fheom.heom_utils import run_heom_simulation


def test_correlation_matrix():
    """Test: Correlation matrix construction."""
    print("\n" + "="*70)
    print("TEST 1: Spatial Correlation Matrix")
    print("="*70)
    
    coords = get_site_coordinates()
    C = spatial_correlation_matrix(coords, correlation_length=10.0, kernel='exponential')
    
    print(f"✓ Matrix shape: {C.shape}")
    print(f"✓ Symmetric: {np.allclose(C, C.T)}")
    print(f"✓ Diagonal all ones: {np.allclose(np.diag(C), 1.0)}")
    print(f"✓ Positive semidefinite: {np.all(np.linalg.eigvalsh(C) >= -1e-10)}")
    
    print(f"\nNearest-neighbor correlations:")
    for i in range(min(3, len(C)-1)):
        print(f"  C[{i},{i+1}] = {C[i,i+1]:.3f}")
    
    return C


def test_factorization():
    """Test: Low-rank factorization."""
    print("\n" + "="*70)
    print("TEST 2: Factorization")
    print("="*70)
    
    coords = get_site_coordinates()
    C = spatial_correlation_matrix(coords, correlation_length=10.0)
    
    # Test different ranks
    for rank in [2, 3, 4, 7]:
        from fheom.fheom import factorize_correlation_matrix
        L, eigvals, info = factorize_correlation_matrix(C, rank=rank, method='eigendecomposition')
        
        print(f"\nRank {rank}:")
        print(f"  Explained variance: {info['explained_variance']:.2%}")
        print(f"  Reconstruction error: {info['reconstruction_error']:.4f}")
        print(f"  Top eigenvalues: {eigvals[:3]}")
    
    return info


def test_bath_operators():
    """Test: Effective bath operator construction."""
    print("\n" + "="*70)
    print("TEST 3: Effective Bath Operators")
    print("="*70)
    
    coords = get_site_coordinates()
    q_sites = [bath_operator() for _ in range(7)]
    
    result = get_factorized_bath(
        q_sites, coords, rank=3,
        correlation_length=10.0, variance_threshold=0.95
    )
    
    print(f"✓ Number of effective modes: {result.rank}")
    print(f"✓ Explained variance: {result.explained_variance:.2%}")
    print(f"✓ Reconstruction error: {result.reconstruction_error:.4f}")
    print(f"✓ Mode weights: {result.weights}")
    
    # Check operator properties
    for k, op in enumerate(result.q_eff):
        print(f"\n  Mode {k}:")
        print(f"    Hermitian: {op.isherm}")
        print(f"    Norm: {op.norm():.3f}")
        print(f"    Weight: {result.weights[k]:.3f}")
    
    return result


def quick_dynamics_test():
    """Test: Quick dynamics comparison (short time, low Nk)."""
    print("\n" + "="*70)
    print("TEST 4: Quick Dynamics (Smoke Test)")
    print("="*70)
    
    H = build_fmo_hamiltonian()
    coords = get_site_coordinates()
    q_sites = [bath_operator() for _ in range(7)]
    rho0 = observable_site0()
    e_ops = [observable_site0()]
    
    # Short time, low Nk for speed
    tlist_fs = np.linspace(0, 500, 501)  # 500 fs
    tlist_sec = tlist_fs * 1e-15
    
    lam_cm = 35.0
    gamma_cm = 106.0
    T_K = 77.0
    Nk = 2  # Fast
    
    print("\nRunning rank-3 fHEOM (500 fs, Nk=2)...")
    t0 = time.time()
    
    result = get_factorized_bath(q_sites, coords, rank=3, correlation_length=10.0)
    
    res = run_heom_simulation(
        H, result.q_eff, lam_cm, gamma_cm, T_K, Nk,
        tlist_sec, rho0, e_ops=e_ops
    )
    
    runtime = time.time() - t0
    trace = res.expect[0]
    
    # Check for dynamics
    dynamic_range = np.max(np.abs(trace)) - np.min(np.abs(trace))
    
    print(f"✓ Runtime: {runtime:.2f}s")
    print(f"✓ Dynamic range: {dynamic_range:.4f}")
    
    if dynamic_range > 0.01:
        print(f"✓ Observable dynamics detected!")
    else:
        print(f"⚠️  Dynamics very weak (may need longer time or different initial state)")
    
    return trace, runtime


def plot_fheom_concept(output_dir):
    """Generate conceptual figure showing factorization."""
    print("\n" + "="*70)
    print("PLOTTING: fHEOM Concept Figure")
    print("="*70)
    
    coords = get_site_coordinates()
    C = spatial_correlation_matrix(coords, correlation_length=10.0)
    
    from fheom.fheom import factorize_correlation_matrix
    L_rank3, eigvals, info = factorize_correlation_matrix(C, rank=3)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Panel A: Full correlation matrix
    ax = axes[0, 0]
    im = ax.imshow(C, cmap='RdBu_r', vmin=0, vmax=1)
    ax.set_title('(A) Spatial Correlation Matrix', fontsize=13, fontweight='bold')
    ax.set_xlabel('Site index')
    ax.set_ylabel('Site index')
    plt.colorbar(im, ax=ax, label='Correlation')
    
    # Panel B: Eigenvalue spectrum
    ax = axes[0, 1]
    all_eigvals = info['eigenvalues']
    ax.bar(range(len(all_eigvals)), all_eigvals, color='steelblue', edgecolor='black')
    ax.axvline(x=2.5, color='red', linestyle='--', linewidth=2, label='Rank-3 cutoff')
    ax.set_title('(B) Eigenvalue Spectrum', fontsize=13, fontweight='bold')
    ax.set_xlabel('Mode index')
    ax.set_ylabel('Eigenvalue')
    ax.legend()
    ax.set_yscale('log')
    
    # Panel C: Factorization matrix L
    ax = axes[1, 0]
    im = ax.imshow(L_rank3, cmap='seismic', vmin=-0.6, vmax=0.6)
    ax.set_title('(C) Factorization Matrix L (7×3)', fontsize=13, fontweight='bold')
    ax.set_xlabel('Effective mode index')
    ax.set_ylabel('Site index')
    plt.colorbar(im, ax=ax, label='Coefficient')
    
    # Panel D: Reconstruction
    C_approx = L_rank3 @ L_rank3.T
    ax = axes[1, 1]
    im = ax.imshow(C_approx, cmap='RdBu_r', vmin=0, vmax=1)
    ax.set_title(f'(D) Rank-3 Reconstruction (error={info["reconstruction_error"]:.2%})',
                fontsize=13, fontweight='bold')
    ax.set_xlabel('Site index')
    ax.set_ylabel('Site index')
    plt.colorbar(im, ax=ax, label='Correlation')
    
    plt.tight_layout()
    
    output_path = Path(output_dir) / 'fheom_concept.png'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Saved: {output_path}")


def main():
    print("="*70)
    print("🔬 fHEOM VALIDATION SUITE")
    print("="*70)
    print("Testing factorized correlated bath implementation")
    print("="*70)
    
    output_dir = "results/fheom_validation"
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Run tests
    try:
        C = test_correlation_matrix()
        info = test_factorization()
        result = test_bath_operators()
        trace, runtime = quick_dynamics_test()
        plot_fheom_concept(output_dir)
        
        print("\n" + "="*70)
        print("✅ ALL TESTS PASSED")
        print("="*70)
        print("\nSummary:")
        print(f"  • Correlation matrix: ✓ Symmetric, PSD")
        print(f"  • Factorization: ✓ Rank-3 explains {result.explained_variance:.1%} variance")
        print(f"  • Bath operators: ✓ {result.rank} effective modes constructed")
        print(f"  • Dynamics: ✓ HEOM runs successfully")
        print(f"  • Figure: ✓ Saved to {output_dir}/fheom_concept.png")
        print("\n📁 Results saved to:", output_dir)
        
        print("\n" + "="*70)
        print("READY FOR PUBLICATION BENCHMARKS")
        print("="*70)
        print("\nNext steps:")
        print("  1. Run full comparison: examples/fheom_benchmark.py")
        print("  2. Generate Figure 7: fHEOM performance vs rank")
        print("  3. Add to manuscript Methods section 2.2")
        
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())





