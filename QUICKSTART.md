# fHEOM Quick Start Guide

## 🚀 Get Running in 2 Minutes

### Step 1: Install
```bash
cd fHEOM
pip install -e .
```

### Step 2: Import and Use
```python
from fheom import get_factorized_bath
from fheom.fmo import build_fmo_hamiltonian, bath_operator, get_site_coordinates

# Setup FMO complex
H = build_fmo_hamiltonian()
coords = get_site_coordinates()
q_sites = [bath_operator() for _ in range(7)]

# Get factorized bath (rank-3 instead of 7 modes!)
result = get_factorized_bath(q_sites, coords, rank=3)

print(f"✓ Reduced to {result.rank} effective modes")
print(f"✓ Variance explained: {result.explained_variance:.1%}")
print(f"✓ Error: {result.reconstruction_error:.4f}")
```

### Step 3: Run Full Validation
```bash
python examples/fheom_validation.py
```
Generates visualization at: `results/fheom_validation/fheom_concept.png`

## 📚 Core API Reference

### Main Function: `get_factorized_bath()`
**One-shot factorization of spatial correlation**

```python
from fheom import get_factorized_bath

result = get_factorized_bath(
    q_site_ops,              # List of N bath operators
    coordinates,             # N×3 array of site positions
    rank=3,                  # Number of effective modes (auto if None)
    correlation_length=10.0, # Spatial decay (Angstroms)
    variance_threshold=0.99, # Variance to retain
    kernel='exponential'     # 'exponential', 'gaussian', 'power_law'
)

# Access results
result.q_eff                # List of r effective bath operators
result.rank                 # Rank used
result.explained_variance   # Fraction of variance retained (0-1)
result.reconstruction_error # L2 error in correlation matrix
result.weights              # Mode weights for bath scaling
```

### Supporting Functions

```python
# Build correlation matrix from coordinates
C = spatial_correlation_matrix(
    coordinates,
    correlation_length=10.0,
    kernel='exponential'
)

# Factorize: C ≈ LL^T
L, eigvals, info = factorize_correlation_matrix(
    C,
    rank=3,
    variance_threshold=0.99,
    method='eigendecomposition'  # or 'cholesky'
)

# Construct effective operators
q_eff = construct_factorized_bath_operators(
    q_site_ops,
    L,
    normalize=True
)

# Run HEOM simulation
result = run_heom_simulation(
    H,                      # System Hamiltonian
    q_eff,                  # Bath operators
    lam_cm=35.0,           # Reorganization energy (cm⁻¹)
    gamma_cm=106.0,        # Bath cutoff (cm⁻¹)
    T_K=77.0,              # Temperature (K)
    Nk=2,                  # Hierarchy depth
    tlist=times,           # Time points
    rho0=initial_state,    # Initial density matrix
    e_ops=[observable]     # Expectation value operators
)
```

## 🎯 Key Features

| Feature | Benefit |
|---------|---------|
| **Low-rank factorization** | Reduce modes from N to r ≪ N |
| **Spatial correlation** | Capture realistic bath structure |
| **Automatic rank selection** | Variance-threshold or manual |
| **Multiple kernels** | Exponential, Gaussian, power-law |
| **HEOM integration** | Works seamlessly with QuTiP |
| **CPU-based** | No GPU required |

## 📊 Computational Complexity

**FMO Complex (N=7 sites, Nk=3 hierarchy depth):**

| Method | Modes | ADO Count | Scaling |
|--------|-------|-----------|---------|
| Full HEOM | 7 | 16,384 | (Nk+1)^N |
| fHEOM rank-2 | 2 | 16 | (Nk+1)^r |
| fHEOM rank-3 | 3 | 64 | (Nk+1)^r |
| fHEOM rank-4 | 4 | 256 | (Nk+1)^r |

**Note**: Runtime and memory depend on solver implementation and hardware.

## 🔬 Validation

The repository includes comprehensive tests:

```bash
# Run all validation tests
python examples/fheom_validation.py

# Tests included:
# 1. Correlation matrix construction
# 2. Low-rank factorization accuracy
# 3. Effective bath operator properties
# 4. HEOM dynamics simulation
# 5. Concept visualization
```

## 🐍 System Requirements

- Python 3.8+
- QuTiP 5.2+ (HEOM solver)
- NumPy, SciPy (CPU-based operations)
- Matplotlib (visualization)
- No GPU or CUDA required

## 📖 Documentation

- **README.md** - Full algorithm description and references
- **MANIFEST.md** - What's included in the package
- **VERIFY.md** - Verification checklist
- **REPLICATION_SUMMARY.md** - Implementation details

## 🤔 Common Questions

**Q: How do I choose the rank?**  
A: Use `rank=None` for automatic selection based on `variance_threshold` (default 0.99), or specify manually (e.g., `rank=3`).

**Q: What coordinates format?**  
A: N×3 NumPy array in Angstroms. Example for FMO:
```python
coords = np.array([
    [28.3, 18.3, 19.3],
    [32.4, 29.5, 13.6],
    ...
])
```

**Q: Which kernel should I use?**  
A: 'exponential' (default) is standard. Try 'gaussian' for sharper decay or 'power_law' for long-range effects.

**Q: Can I use my own system?**  
A: Yes! Provide your Hamiltonian, site coordinates, and bath operators. See FMO implementation for reference.

## ⚡ Quick Examples

### Example 1: Get Effective Modes
```python
from fheom import get_factorized_bath
from fheom.fmo import *

H = build_fmo_hamiltonian()
coords = get_site_coordinates()
q_sites = [bath_operator() for _ in range(7)]

# Rank-3 factorization
result = get_factorized_bath(q_sites, coords, rank=3)
print(f"Modes: {len(result.q_eff)}, Variance: {result.explained_variance:.1%}")
```

### Example 2: Compare Ranks
```python
for rank in [2, 3, 4]:
    result = get_factorized_bath(q_sites, coords, rank=rank)
    print(f"Rank {rank}: {result.explained_variance:.1%} variance")
```

### Example 3: Full HEOM Simulation
```python
from fheom import run_heom_simulation
import numpy as np

# Setup
H = build_fmo_hamiltonian()
rho0 = observable_site0()  # Initial state
e_ops = [observable_site0()]  # Observable
tlist = np.linspace(0, 500e-15, 501)  # 500 fs

# Get factorized bath
result = get_factorized_bath(q_sites, coords, rank=3)

# Run simulation
heom_result = run_heom_simulation(
    H, result.q_eff,
    lam_cm=35.0, gamma_cm=106.0, T_K=77.0, Nk=2,
    tlist=tlist, rho0=rho0, e_ops=e_ops
)

# Plot
import matplotlib.pyplot as plt
plt.plot(tlist*1e15, heom_result.expect[0])
plt.xlabel('Time (fs)')
plt.ylabel('Site 0 population')
plt.show()
```

---

**Ready to get started? Run:**
```bash
cd fHEOM
python examples/fheom_validation.py
```

For more details, see **README.md** 📖
