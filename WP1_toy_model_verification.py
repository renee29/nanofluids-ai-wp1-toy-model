#!/usr/bin/env python3
"""
================================================================================
NANOFLUIDS-AI: WP1 Preliminary Result
Analytical Verification of the Non-local Kernel Toy Model
================================================================================

This script verifies the asymptotic behaviour of Gaussian convolution kernels
under nanoconfinement, as described in:

    NanoFluids-AI: Structural Inverse Problems for Nanoscale Continuum Laws
    ERC Consolidator Grant 2026 - Part B2, WP1 (Mathematical Foundations)

Key Result:
    For an isotropic Gaussian kernel K_ξ confined to a slit domain Ω_L,
    domain-induced symmetry breaking generates effective anisotropy:
    
        λ_∥ = 1  (unbounded in-plane integration)
        λ_⊥ = erf(L / (2√2 ξ))  (truncated out-of-plane integration)

Physical Interpretation:
    The ratio λ_∥/λ_⊥ quantifies the emergent dielectric anisotropy arising
    purely from geometric confinement, even for an isotropic molecular response.

Author: NanoFluids-AI Team
License: MIT
Repository: https://github.com/[TBD]/nanofluids-ai
================================================================================
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate as scipy_integrate
from scipy.special import erf
from sympy import symbols, exp, integrate, sqrt, pi, erf as sym_erf, simplify, limit, oo

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
FIGURE_DPI = 300
FIGURE_FORMAT = 'png'  # For publication: 'pdf' or 'eps'
USE_LATEX = False      # Set True if LaTeX fully installed


def setup_matplotlib():
    """Configure matplotlib for publication-quality figures."""
    if USE_LATEX:
        plt.rcParams.update({
            'text.usetex': True,
            'font.family': 'serif',
            'font.serif': ['Computer Modern Roman'],
            'font.size': 10,
            'axes.labelsize': 11,
            'axes.titlesize': 11,
            'legend.fontsize': 9,
            'xtick.labelsize': 9,
            'ytick.labelsize': 9,
        })
    plt.rcParams.update({
        'figure.figsize': (3.5, 2.8),  # Single-column width
        'figure.dpi': FIGURE_DPI,
        'axes.linewidth': 0.8,
        'lines.linewidth': 1.2,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
    })


# -----------------------------------------------------------------------------
# Symbolic Verification (SymPy)
# -----------------------------------------------------------------------------
def verify_symbolic():
    """
    Analytically compute λ_⊥ using symbolic integration.
    
    Returns
    -------
    lambda_perp_expr : sympy expression
        Symbolic expression for λ_⊥(L, ξ)
    """
    print("=" * 60)
    print("SYMBOLIC VERIFICATION")
    print("=" * 60)
    
    # Define symbolic variables
    z, L, xi = symbols('z L xi', real=True, positive=True)
    
    # Gaussian kernel (1D marginal in z-direction, normalised)
    # Full 3D kernel: K(r) = (2πξ²)^(-3/2) exp(-|r|²/(2ξ²))
    # Marginal in z:  K(z) = (2πξ²)^(-1/2) exp(-z²/(2ξ²))
    kernel_z = (1 / (sqrt(2 * pi) * xi)) * exp(-z**2 / (2 * xi**2))
    
    # Integration bounds: slit domain z ∈ [-L/2, L/2]
    z_min, z_max = -L/2, L/2
    
    # Compute λ_⊥ = ∫_{-L/2}^{L/2} K(z) dz
    lambda_perp = integrate(kernel_z, (z, z_min, z_max))
    lambda_perp = simplify(lambda_perp)
    
    print(f"\nKernel:      K(z) = (2πξ²)^(-1/2) exp(-z²/(2ξ²))")
    print(f"Domain:      z ∈ [-L/2, L/2]")
    print(f"\nResult:      λ_⊥ = {lambda_perp}")
    
    # Verify asymptotic limits
    print("\n--- Asymptotic Limits ---")
    
    # Limit 1: ξ → 0 (local limit, kernel → δ)
    # Should recover λ_⊥ → 1 (full integral of delta)
    limit_local = limit(lambda_perp, xi, 0, '+')
    print(f"Local limit (ξ → 0⁺):     λ_⊥ → {limit_local}")
    
    # Limit 2: L → ∞ (bulk limit, no confinement)
    # Should recover λ_⊥ → 1 (full integral of Gaussian)
    limit_bulk = limit(lambda_perp, L, oo)
    print(f"Bulk limit  (L → ∞):      λ_⊥ → {limit_bulk}")
    
    # Limit 3: L → 0 (extreme confinement)
    # Should recover λ_⊥ → 0 (no space to integrate)
    limit_confined = limit(lambda_perp, L, 0, '+')
    print(f"Confined    (L → 0⁺):     λ_⊥ → {limit_confined}")
    
    return lambda_perp


# -----------------------------------------------------------------------------
# Numerical Validation (SciPy)
# -----------------------------------------------------------------------------
def lambda_perp_analytical(L, xi):
    """
    Analytical formula for λ_⊥.
    
    Parameters
    ----------
    L : float or array
        Confinement length [nm]
    xi : float
        Correlation length [nm]
    
    Returns
    -------
    lambda_perp : float or array
        Out-of-plane eigenvalue
    """
    return erf(L / (2 * np.sqrt(2) * xi))


def lambda_perp_numerical(L, xi, n_points=1000):
    """
    Numerical integration for validation.
    
    Parameters
    ----------
    L : float
        Confinement length [nm]
    xi : float
        Correlation length [nm]
    n_points : int
        Number of quadrature points
    
    Returns
    -------
    result : float
        Numerically integrated λ_⊥
    """
    def kernel(z):
        return (1 / (np.sqrt(2 * np.pi) * xi)) * np.exp(-z**2 / (2 * xi**2))
    
    result, _ = scipy_integrate.quad(kernel, -L/2, L/2)
    return result


def validate_numerical():
    """
    Cross-validate analytical formula against numerical quadrature.
    
    Returns
    -------
    max_error : float
        Maximum relative error between analytical and numerical
    """
    print("\n" + "=" * 60)
    print("NUMERICAL VALIDATION")
    print("=" * 60)
    
    # Test grid
    L_values = np.linspace(0.5, 10, 20)  # nm
    xi_values = [0.3, 0.5, 1.0]          # nm
    
    max_error = 0
    print(f"\n{'L [nm]':>8} {'ξ [nm]':>8} {'Analytical':>12} {'Numerical':>12} {'Rel. Err.':>12}")
    print("-" * 56)
    
    for xi in xi_values:
        for L in L_values[::5]:  # Sample every 5th point for display
            ana = lambda_perp_analytical(L, xi)
            num = lambda_perp_numerical(L, xi)
            rel_err = abs(ana - num) / num if num > 1e-10 else 0
            max_error = max(max_error, rel_err)
            print(f"{L:8.2f} {xi:8.2f} {ana:12.6f} {num:12.6f} {rel_err:12.2e}")
    
    print(f"\nMaximum relative error: {max_error:.2e}")
    assert max_error < 1e-10, "Numerical validation failed!"
    print("✓ Analytical formula validated against numerical integration.")
    
    return max_error


# -----------------------------------------------------------------------------
# Anisotropy Analysis
# -----------------------------------------------------------------------------
def compute_anisotropy_ratio(L, xi):
    """
    Compute the anisotropy ratio λ_∥/λ_⊥.
    
    Since λ_∥ = 1 (unbounded in-plane), the ratio is simply 1/λ_⊥.
    
    Parameters
    ----------
    L : float or array
        Confinement length [nm]
    xi : float
        Correlation length [nm]
    
    Returns
    -------
    ratio : float or array
        Anisotropy ratio λ_∥/λ_⊥
    """
    lambda_perp = lambda_perp_analytical(L, xi)
    # Avoid division by zero for very small L
    lambda_perp = np.maximum(lambda_perp, 1e-10)
    return 1.0 / lambda_perp


def analyze_anisotropy():
    """
    Analyze anisotropy ratio and compare with experimental data.
    """
    print("\n" + "=" * 60)
    print("ANISOTROPY ANALYSIS")
    print("=" * 60)
    
    # Representative parameters
    xi = 0.3  # nm (molecular correlation length for water)
    L_values = [1.0, 2.0, 5.0, 10.0]  # nm
    
    print(f"\nCorrelation length ξ = {xi} nm")
    print(f"\n{'L [nm]':>8} {'λ_⊥':>10} {'λ_∥/λ_⊥':>10}")
    print("-" * 30)
    
    for L in L_values:
        lp = lambda_perp_analytical(L, xi)
        ratio = compute_anisotropy_ratio(L, xi)
        print(f"{L:8.1f} {lp:10.4f} {ratio:10.2f}")
    
    # Comparison with experiment (Wang et al. 2025)
    print("\n--- Comparison with Experiment ---")
    L_exp = 1.0  # nm (bilayer water)
    ratio_exp = 500  # ε_∥/ε_⊥ from Wang et al. Nature 2025
    ratio_theory = compute_anisotropy_ratio(L_exp, xi)
    
    print(f"Confinement L = {L_exp} nm:")
    print(f"  Gaussian kernel prediction: λ_∥/λ_⊥ = {ratio_theory:.2f}")
    print(f"  Experimental observation:   ε_∥/ε_⊥ ≈ {ratio_exp}")
    print(f"  Discrepancy factor: {ratio_exp / ratio_theory:.0f}×")
    print("\n→ The Gaussian kernel captures qualitative anisotropy but")
    print("  underestimates magnitude by ~400×, motivating data-driven")
    print("  kernel discovery in WP3.")


# -----------------------------------------------------------------------------
# Publication Figure
# -----------------------------------------------------------------------------
def generate_figure(save=True):
    """
    Generate publication-quality figure showing λ_⊥(L/ξ) and anisotropy.
    
    Parameters
    ----------
    save : bool
        If True, save figure to file
    
    Returns
    -------
    fig : matplotlib Figure
    """
    setup_matplotlib()
    
    # Dimensionless confinement parameter
    L_over_xi = np.linspace(0.1, 10, 200)
    
    # Compute eigenvalues (for ξ = 1, so L/ξ = L)
    lambda_perp = lambda_perp_analytical(L_over_xi, xi=1.0)
    lambda_para = np.ones_like(L_over_xi)
    ratio = lambda_para / lambda_perp
    
    # Create figure with two panels
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 2.8))
    
    # Panel (a): Eigenvalues
    ax1.plot(L_over_xi, lambda_para, 'b-', label=r'$\lambda_\parallel = 1$')
    ax1.plot(L_over_xi, lambda_perp, 'r-', label=r'$\lambda_\perp = \mathrm{erf}(L/2\sqrt{2}\xi)$')
    ax1.axhline(1, color='gray', linestyle=':', linewidth=0.5)
    ax1.set_xlabel(r'$L/\xi$')
    ax1.set_ylabel(r'Eigenvalue $\lambda$')
    ax1.set_xlim(0, 10)
    ax1.set_ylim(0, 1.1)
    ax1.legend(loc='right', fontsize=8, framealpha=0.9)
    ax1.text(-0.25, 1.05, '(a)', transform=ax1.transAxes, 
             fontsize=11, fontweight='bold', verticalalignment='top')
    
    # Panel (b): Anisotropy ratio
    ax2.semilogy(L_over_xi, ratio, 'k-', linewidth=1.5)
    ax2.axhline(500, color='orange', linestyle='--', linewidth=1, 
                label=r'Expt. $\varepsilon_\parallel/\varepsilon_\perp \approx 500$')
    ax2.axvline(1/0.3, color='green', linestyle=':', linewidth=1,
                label=r'$L = 1$ nm, $\xi = 0.3$ nm')
    ax2.set_xlabel(r'$L/\xi$')
    ax2.set_ylabel(r'Anisotropy ratio $\lambda_\parallel/\lambda_\perp$')
    ax2.set_xlim(0, 10)
    ax2.set_ylim(1, 1000)
    #ax2.legend(loc='upper right', fontsize=8)
    ax2.legend(loc='center right', fontsize=8, framealpha=0.9) 
    #ax2.legend(loc='center right', bbox_to_anchor=(0.98, 0.6), fontsize=9)
    ax2.text(-0.25, 1.05, '(b)', transform=ax2.transAxes,
             fontsize=11, fontweight='bold', verticalalignment='top')
    
    plt.tight_layout()
    
    if save:
        filename = f'WP1_toy_model_anisotropy.{FIGURE_FORMAT}'
        fig.savefig(filename, dpi=FIGURE_DPI, format=FIGURE_FORMAT)
        print(f"\n✓ Figure saved: {filename}")
    
    return fig


# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------
def main():
    """Run complete verification pipeline."""
    print("\n" + "=" * 70)
    print("   NANOFLUIDS-AI: WP1 TOY MODEL - COMPLETE VERIFICATION")
    print("=" * 70)
    
    # Step 1: Symbolic derivation
    lambda_expr = verify_symbolic()
    
    # Step 2: Numerical cross-validation
    validate_numerical()
    
    # Step 3: Physical interpretation
    analyze_anisotropy()
    
    # Step 4: Generate publication figure
    fig = generate_figure(save=True)
    
    print("\n" + "=" * 70)
    print("   VERIFICATION COMPLETE")
    print("=" * 70)
    print("\nKey result for WP1 manuscript:")
    print("  λ_⊥(L,ξ) = erf(L / 2√2ξ)")
    print("  λ_∥ = 1")
    print("  Anisotropy ratio = [erf(L / 2√2ξ)]⁻¹")
    print("\n→ Gaussian kernels produce qualitative anisotropy from")
    print("  confinement but require data-driven refinement (WP3).")
    
    return fig


if __name__ == "__main__":
    fig = main()
    plt.show()
