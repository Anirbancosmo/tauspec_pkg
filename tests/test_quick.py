#!/usr/bin/env python3
"""Quick test of TauSpec functionality."""

import numpy as np
from tauspec import TauSpecCalculator
from tauspec.config import CosmologyConfig, ReionizationConfig, BubbleConfig
import time

print("Testing TauSpec Package")
print("=" * 60)

# Create configurations
print("\n1. Creating configurations...")
cosmo_config = CosmologyConfig(h=0.69, omega_m=0.2877)
reion_config = ReionizationConfig(z_re=7.5, delta_z=2.0)
bubble_config = BubbleConfig(Rb=5.0, b=6.0)
print("   ✓ Configurations created")

# Initialize calculator
print("\n2. Initializing calculator...")
start = time.time()
calc = TauSpecCalculator(cosmo_config, reion_config, bubble_config)
print(f"   ✓ Calculator initialized in {time.time()-start:.3f} seconds")

# Test ionization fraction
print("\n3. Testing ionization fraction xe(z)...")
z_test = np.array([5.0, 7.5, 10.0, 15.0])
xe_test = calc.xe(z_test)
for z, xe in zip(z_test, xe_test):
    print(f"   xe(z={z:.1f}) = {xe:.4f}")

# Test power spectrum calculation
print("\n4. Testing power spectrum P_xe(k, z)...")
start = time.time()
k = np.logspace(-1, 0, 10)
pxe = calc.p_xexe(7.5, k)
print(f"   ✓ P_xe calculated for {len(k)} k values in {time.time()-start:.3f} seconds")
print(f"   P_xe(k=0.1, z=7.5) = {calc.p_xexe(7.5, 0.1):.4e} (Mpc/h)^3")

# Test C_l calculation
print("\n5. Testing C_l^tau calculation...")
start = time.time()
ell = np.array([2, 10, 100, 1000])
cl_tau = calc.cltau(ell)
elapsed = time.time() - start
print(f"   ✓ C_l^tau calculated for {len(ell)} multipoles in {elapsed:.3f} seconds")
for l, cl in zip(ell, cl_tau):
    print(f"   C_l^tau(l={l}) = {cl:.4e}")

# Test file output
print("\n6. Testing file output...")
calc.save_output('xe', filename='test_xe.dat')
calc.save_output('pxe', z=7.5, k=np.logspace(-2, 1, 20), filename='test_pxe.dat')
print("   ✓ Output files created successfully")

print("\n" + "=" * 60)
print("All tests passed! ✓")
print("=" * 60)
