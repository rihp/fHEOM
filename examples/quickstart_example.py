"""
Example: How users will use fHEOM after 'pip install fHEOM'

This shows the typical workflow for someone who just installed your package.
"""

# Step 1: Import fHEOM (after pip install fHEOM)
from fheom import get_factorized_bath, run_heom_simulation
from fheom.fmo import build_fmo_hamiltonian, bath_operator, get_site_coordinates
import numpy as np
import qutip as qt

print("=" * 60)
print("fHEOM Example: FMO Complex Simulation")
print("=" * 60)

# Step 2: Build the system (FMO complex)
print("\n1. Setting up FMO complex...")
H = build_fmo_hamiltonian()
coords = get_site_coordinates()
q_sites = [bath_operator() for _ in range(7)]
print(f"   ✓ 7-site Hamiltonian: {H.shape}")
print(f"   ✓ Site coordinates: {coords.shape}")

# Step 3: Apply fHEOM factorization
print("\n2. Factorizing bath (this is the magic!)...")
result = get_factorized_bath(q_sites, coords, rank=3)
print(f"   ✓ Reduced from 7 modes → {result.rank} modes")
print(f"   ✓ Variance explained: {result.explained_variance*100:.1f}%")
print(f"   ✓ Reconstruction error: {result.reconstruction_error:.4f}")
print(f"   → Hierarchy: 16,384 ADOs → 64 ADOs (256× smaller)")

# Step 4: Run HEOM simulation
print("\n3. Running HEOM dynamics...")
rho0 = qt.ket2dm(qt.basis(7, 0))  # Initial state: site 0 excited
tlist = np.linspace(0, 500e-15, 101)  # 500 fs
e_ops = [qt.basis(7, i) * qt.basis(7, i).dag() for i in range(7)]  # All site populations

heom_result = run_heom_simulation(
    H, result.q_eff,
    lam_cm=35.0,      # Reorganization energy
    gamma_cm=106.0,   # Bath cutoff
    T_K=77.0,         # Temperature
    Nk=2,             # Hierarchy depth
    tlist=tlist,
    rho0=rho0,
    e_ops=e_ops
)
print(f"   ✓ Simulation complete: {len(tlist)} time points")

# Step 5: Analyze results
print("\n4. Results:")
print(f"   Site 0 population at t=0:   {heom_result.expect[0][0]:.3f}")
print(f"   Site 0 population at t=500fs: {heom_result.expect[0][-1]:.3f}")
print(f"   Energy transfer: {(1-heom_result.expect[0][-1])*100:.1f}% left site 0")

# Step 6: Optional plotting
print("\n5. To visualize:")
print("   import matplotlib.pyplot as plt")
print("   plt.plot(tlist*1e15, heom_result.expect[0])")
print("   plt.xlabel('Time (fs)')")
print("   plt.ylabel('Site 0 Population')")
print("   plt.show()")

print("\n" + "=" * 60)
print("✓ fHEOM simulation complete!")
print("=" * 60)
