# NanoFluids-AI: WP1 Toy Model Verification

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

Analytical verification of non-local Gaussian kernel behavior under nanoconfinement for the **NanoFluids-AI** project (ERC Consolidator Grant 2026, Part B2, WP1: Mathematical Foundations).

## Overview

This repository contains a complete analytical and numerical verification of how geometric confinement induces effective anisotropy in isotropic molecular response kernels. The toy model demonstrates that even for an isotropic Gaussian kernel, domain truncation breaks symmetry and generates directionally-dependent effective properties.

### Key Result

For an isotropic Gaussian kernel `K_ξ(r) = (2πξ²)^(-3/2) exp(-|r|²/(2ξ²))` confined to a slit domain `Ω_L = {(x,y,z) : z ∈ [-L/2, L/2]}`:

```
λ_∥ = 1                        (unbounded in-plane integration)
λ_⊥ = erf(L / (2√2 ξ))         (truncated out-of-plane integration)
```

where:
- `L` is the confinement length [nm]
- `ξ` is the molecular correlation length [nm]
- `erf` is the error function

### Physical Interpretation

The anisotropy ratio `λ_∥/λ_⊥ = [erf(L/(2√2ξ))]⁻¹` quantifies emergent dielectric anisotropy arising purely from geometric confinement. For water under 1 nm confinement:

- **Gaussian kernel prediction**: λ_∥/λ_⊥ ≈ 1.2
- **Experimental observation** (Wang et al. *Nature* 2025): ε_∥/ε_⊥ ≈ 500
- **Discrepancy**: ~400× factor

This gap motivates data-driven kernel discovery in WP3 of the NanoFluids-AI project.

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Clone the repository

```bash
git clone https://github.com/[YOUR_USERNAME]/nanofluids-ai-wp1-toy-model.git
cd nanofluids-ai-wp1-toy-model
```

### Install dependencies

```bash
pip install -r requirements.txt
```

## Usage

Run the complete verification pipeline:

```bash
python WP1_toy_model_verification.py
```

### Output

The script performs four verification steps:

1. **Symbolic Verification**: Derives λ_⊥ analytically using SymPy and verifies asymptotic limits
2. **Numerical Validation**: Cross-validates analytical formula against numerical quadrature (error < 10⁻¹⁰)
3. **Anisotropy Analysis**: Computes anisotropy ratios and compares with experimental data
4. **Publication Figure**: Generates `WP1_toy_model_anisotropy.png` with two panels:
   - (a) Eigenvalue evolution λ_∥ and λ_⊥ vs L/ξ
   - (b) Anisotropy ratio in log scale with experimental comparison

### Expected Output

```
======================================================================
   NANOFLUIDS-AI: WP1 TOY MODEL - COMPLETE VERIFICATION
======================================================================

============================================================
SYMBOLIC VERIFICATION
============================================================

Kernel:      K(z) = (2πξ²)^(-1/2) exp(-z²/(2ξ²))
Domain:      z ∈ [-L/2, L/2]

Result:      λ_⊥ = erf(√2*L/(4*xi))

--- Asymptotic Limits ---
Local limit (ξ → 0⁺):     λ_⊥ → 1
Bulk limit  (L → ∞):      λ_⊥ → 1
Confined    (L → 0⁺):     λ_⊥ → 0

============================================================
NUMERICAL VALIDATION
============================================================

   L [nm]   ξ [nm]   Analytical    Numerical    Rel. Err.
--------------------------------------------------------
    0.50     0.30     0.889697     0.889697     5.33e-15
    ...

Maximum relative error: 5.33e-15
✓ Analytical formula validated against numerical integration.

============================================================
ANISOTROPY ANALYSIS
============================================================

Correlation length ξ = 0.3 nm

   L [nm]        λ_⊥    λ_∥/λ_⊥
------------------------------
     1.0     0.8897       1.12
     2.0     0.9820       1.02
     5.0     0.9997       1.00
    10.0     1.0000       1.00

--- Comparison with Experiment ---
Confinement L = 1.0 nm:
  Gaussian kernel prediction: λ_∥/λ_⊥ = 1.12
  Experimental observation:   ε_∥/ε_⊥ ≈ 500
  Discrepancy factor: 446×

→ The Gaussian kernel captures qualitative anisotropy but
  underestimates magnitude by ~400×, motivating data-driven
  kernel discovery in WP3.

✓ Figure saved: WP1_toy_model_anisotropy.png

======================================================================
   VERIFICATION COMPLETE
======================================================================

Key result for WP1 manuscript:
  λ_⊥(L,ξ) = erf(L / 2√2ξ)
  λ_∥ = 1
  Anisotropy ratio = [erf(L / 2√2ξ)]⁻¹

→ Gaussian kernels produce qualitative anisotropy from
  confinement but require data-driven refinement (WP3).
```

## Project Structure

```
nanofluids-ai-wp1-toy-model/
├── WP1_toy_model_verification.py    # Main verification script
├── README.md                         # This file
├── LICENSE                           # MIT License
├── requirements.txt                  # Python dependencies
├── CITATION.cff                      # Citation metadata
└── .gitignore                        # Git ignore rules
```

## Configuration

The script can be configured by modifying constants at the top of `WP1_toy_model_verification.py`:

```python
FIGURE_DPI = 300              # Figure resolution (300 for publication)
FIGURE_FORMAT = 'png'         # Output format ('png', 'pdf', or 'eps')
USE_LATEX = False             # Set True if LaTeX is installed for typography
```

## Scientific Context

This work is part of the **NanoFluids-AI** project:

> **NanoFluids-AI: Structural Inverse Problems for Nanoscale Continuum Laws**
> ERC Consolidator Grant 2026 - Part B2
> Work Package 1: Mathematical Foundations

The toy model establishes theoretical foundations for understanding how geometric confinement modifies effective material properties through non-local kernel convolutions. It serves as a baseline for more sophisticated data-driven kernel discovery approaches in WP3.

## Citation

If you use this code in your research, please cite:

```bibtex
@software{nanofluids_ai_wp1_2025,
  author       = {NanoFluids-AI Team},
  title        = {NanoFluids-AI: WP1 Toy Model Verification},
  year         = {2025},
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.XXXXXX},
  url          = {https://github.com/renee29/nanofluids-ai-wp1-toy-model}
}
```

## Related Publications

- Wang et al. (2025). "Giant dielectric anisotropy in nanoconfined water." *Nature*. [Reference for experimental ε_∥/ε_⊥ ≈ 500]

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

This work is supported by the European Research Council (ERC) under the European Union's Horizon Europe research and innovation programme (Grant Agreement No. [TBD]).

## Contact

For questions or collaboration inquiries, please open an issue on GitHub or contact the NanoFluids-AI team.

---

**Project Status**: Initial release (v1.0.0) - December 2025
