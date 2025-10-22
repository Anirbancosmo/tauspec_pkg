#!/usr/bin/env python3
"""Example usage of TauSpec package.

This script demonstrates how to:
1. Set up configurations
2. Calculate ionization fraction power spectra
3. Calculate optical depth power spectra
4. Save and plot results
"""

import numpy as np
import matplotlib.pyplot as plt
from tauspec import TauSpecCalculator
from tauspec.config import CosmologyConfig, ReionizationConfig, BubbleConfig

# Configure plotting
plt.style.use('seaborn-v0_8-darkgrid' if 'seaborn-v0_8-darkgrid' in plt.style.available else 'default')


def main():
    """Run example calculations."""
    
    print("=" * 70)
    print("TauSpec Example: Calculating Optical Depth Power Spectra")
    print("=" * 70)
    
    # 1. Set up configurations
    print("\n1. Setting up configurations...")
    
    cosmo_config = CosmologyConfig(
        h=0.69,
        omega_m=0.2877,
        omega_b=0.04757,
        omega_lambda=0.721,
        n_s=0.9667,
        sigma8=0.80
    )
    
    reion_config = ReionizationConfig(
        model='tanh',
        z_re=7.5,
        delta_z=2.0
    )
    
    bubble_config = BubbleConfig(
        Rb=5.0,
        sigma_lnr=np.log(2),
        b=6.0,
        T0=2e4
    )
    
    print("   ✓ Cosmology: h={}, Ωm={}, Ωb={}".format(
        cosmo_config.h, cosmo_config.omega_m, cosmo_config.omega_b))
    print("   ✓ Reionization: z_re={}, Δz={}".format(
        reion_config.z_re, reion_config.delta_z))
    print("   ✓ Bubbles: Rb={} Mpc, b={}".format(
        bubble_config.Rb, bubble_config.b))
    
    # 2. Create calculator
    print("\n2. Initializing calculator...")
    calc = TauSpecCalculator(cosmo_config, reion_config, bubble_config)
    print("   ✓ Calculator initialized")
    
    # 3. Calculate ionization fraction history
    print("\n3. Calculating ionization fraction history...")
    z_xe = np.linspace(0, 20, 100)
    xe = calc.xe(z_xe)
    
    print("   ✓ xe(z=7.5) = {:.4f}".format(calc.xe(7.5)))
    print("   ✓ xe(z=10) = {:.4f}".format(calc.xe(10.0)))
    
    # 4. Calculate ionization fraction power spectrum
    print("\n4. Calculating ionization fraction power spectrum...")
    z_pxe = 7.5
    k = np.logspace(-2, 1, 30)
    pxe = calc.p_xexe(z_pxe, k)
    
    print("   ✓ P_xe calculated at z={} for {} k values".format(z_pxe, len(k)))
    print("   ✓ P_xe(k=0.1) = {:.4e} (Mpc/h)^3".format(
        calc.p_xexe(z_pxe, 0.1)))
    
    # 5. Calculate optical depth power spectrum
    print("\n5. Calculating optical depth power spectrum...")
    ell = np.logspace(np.log10(2), np.log10(5000), 15)
    cl_tau = calc.cltau(ell)
    
    print("   ✓ C_l^tau calculated for {} multipoles".format(len(ell)))
    print("   ✓ C_l^tau(l=100) = {:.4e}".format(
        np.interp(100, ell, cl_tau)))
    
    # 6. Save outputs
    print("\n6. Saving outputs...")
    calc.save_output('xe', filename='xe_example.dat')
    calc.save_output('pxe', z=z_pxe, k=k, filename='pxe_example.dat')
    calc.save_output('cltau', ell=ell, filename='cltau_example.dat')
    print("   ✓ Outputs saved to xe_example.dat, pxe_example.dat, cltau_example.dat")
    
    # 7. Create plots
    print("\n7. Creating plots...")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Ionization fraction history
    ax = axes[0, 0]
    ax.plot(z_xe, xe, 'b-', linewidth=2)
    ax.axhline(y=0.5, color='r', linestyle='--', alpha=0.5, label='xe=0.5')
    ax.axvline(x=reion_config.z_re, color='r', linestyle='--', 
               alpha=0.5, label=f'z_re={reion_config.z_re}')
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('Ionization Fraction $x_e$')
    ax.set_title('Reionization History')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Power spectrum P_xe(k)
    ax = axes[0, 1]
    delta_xe = k ** 3 * pxe / (2 * np.pi ** 2)
    ax.loglog(k, delta_xe, 'b-', linewidth=2)
    ax.set_xlabel('Wavenumber k [h/Mpc]')
    ax.set_ylabel(r'$\Delta_{x_e}^2(k) = k^3 P_{x_e}(k)/(2\pi^2)$')
    ax.set_title(f'Ionization Fraction Power Spectrum (z={z_pxe})')
    ax.grid(True, alpha=0.3, which='both')
    
    # Plot 3: C_l^tau
    ax = axes[1, 0]
    ax.loglog(ell, ell * (ell + 1) * cl_tau / (2 * np.pi), 'b-', linewidth=2)
    ax.set_xlabel(r'Multipole $\ell$')
    ax.set_ylabel(r'$\ell(\ell+1)C_\ell^{\tau\tau}/(2\pi)$')
    ax.set_title('Optical Depth Power Spectrum')
    ax.grid(True, alpha=0.3, which='both')
    
    # Plot 4: Information summary
    ax = axes[1, 1]
    ax.axis('off')
    
    info_text = f"""
    TauSpec Calculation Summary
    
    Cosmology:
    • h = {cosmo_config.h}
    • Ωm = {cosmo_config.omega_m}
    • Ωb = {cosmo_config.omega_b}
    • σ8 = {cosmo_config.sigma8}
    
    Reionization:
    • Model: {reion_config.model}
    • z_re = {reion_config.z_re}
    • Δz = {reion_config.delta_z}
    
    Bubble Properties:
    • Rb = {bubble_config.Rb} Mpc
    • σ_lnR = {bubble_config.sigma_lnr:.3f}
    • Bias b = {bubble_config.b}
    • Temperature T0 = {bubble_config.T0:.0e} K
    """
    
    ax.text(0.1, 0.5, info_text, fontsize=10, family='monospace',
            verticalalignment='center', transform=ax.transAxes)
    
    plt.tight_layout()
    plt.savefig('tauspec_example.png', dpi=150, bbox_inches='tight')
    print("   ✓ Plot saved to tauspec_example.png")
    
    print("\n" + "=" * 70)
    print("Example completed successfully!")
    print("=" * 70)
    
    # Display plot
    plt.show()


if __name__ == '__main__':
    main()
