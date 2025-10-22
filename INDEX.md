# TauSpec Package - Getting Started Guide

Welcome to TauSpec 0.2.0! This guide will help you get started with the modernized and optimized optical depth power spectrum calculator.

## 📚 Documentation Files

Start here to understand the package:

1. **[SUMMARY.md](SUMMARY.md)** - Complete overview of the package ⭐ START HERE
2. **[README.md](README.md)** - Main documentation and API reference
3. **[MIGRATION.md](MIGRATION.md)** - Guide for users of the old code
4. **[PERFORMANCE.md](PERFORMANCE.md)** - Performance improvements and benchmarks

## 🚀 Quick Start (5 minutes)

### 1. Install
```bash
cd tauspec_pkg
pip install -e .
```

### 2. Run Example
```bash
python test_quick.py
```

### 3. Try It Yourself
```python
import numpy as np
from tauspec import TauSpecCalculator
from tauspec.config import CosmologyConfig, ReionizationConfig, BubbleConfig

# Configure
cosmo = CosmologyConfig(h=0.69, omega_m=0.2877)
reion = ReionizationConfig(z_re=7.5, delta_z=2.0)
bubble = BubbleConfig(Rb=5.0, b=6.0)

# Calculate
calc = TauSpecCalculator(cosmo, reion, bubble)
pxe = calc.p_xexe(7.5, np.logspace(-2, 1, 50))
```

## 📂 Directory Structure

```
tauspec_pkg/
├── 📖 Documentation
│   ├── README.md           - Main documentation
│   ├── SUMMARY.md          - Complete overview
│   ├── MIGRATION.md        - Migration guide
│   ├── PERFORMANCE.md      - Performance details
│   └── INDEX.md            - This file
│
├── 📦 Package Source
│   └── tauspec/
│       ├── __init__.py     - Package interface
│       ├── calculator.py   - Main calculations
│       ├── cosmology.py    - Astropy cosmology
│       ├── reionization.py - Reionization models
│       └── config.py       - Configuration classes
│
├── 🧪 Tests
│   └── tests/
│       └── test_config.py  - Configuration tests
│
├── 📝 Examples
│   ├── examples/
│   │   ├── basic_usage.py      - Simple example with plots
│   │   └── tutorial.ipynb      - Comprehensive notebook
│   ├── test_quick.py           - Quick functionality test
│   └── compare_performance.py  - Performance demo
│
└── ⚙️  Configuration
    ├── setup.py            - Package installation
    ├── pyproject.toml      - Modern packaging
    ├── requirements.txt    - Dependencies
    ├── LICENSE             - MIT License
    └── .gitignore          - Git configuration
```

## 🎯 Common Tasks

### Calculate Reionization History
```python
z = np.linspace(5, 15, 100)
xe = calc.xe(z)
```

### Calculate Power Spectrum at Redshift
```python
k = np.logspace(-2, 1, 50)
pxe = calc.p_xexe(z=7.5, k=k)
```

### Calculate Optical Depth Power Spectrum
```python
ell = np.logspace(np.log10(2), np.log10(5000), 20)
cl_tau = calc.cltau(ell)
```

### Save Results
```python
calc.save_output('pxe', z=7.5, k=k, filename='my_pxe.dat')
calc.save_output('cltau', ell=ell, filename='my_cltau.dat')
calc.save_output('xe', filename='my_xe.dat')
```

### Use Custom Reionization History
```python
# Create file with z, xe columns
# Then:
reion = ReionizationConfig(model='custom', custom_file='my_xe_history.dat')
calc = TauSpecCalculator(cosmo, reion, bubble)
```

## 📊 Examples to Run

### 1. Quick Test (1 minute)
```bash
python test_quick.py
```
Tests basic functionality and shows performance.

### 2. Basic Usage (2 minutes)
```bash
python examples/basic_usage.py
```
Creates plots and saves outputs.

### 3. Performance Demo (2 minutes)
```bash
python compare_performance.py
```
Demonstrates speed improvements and features.

### 4. Full Tutorial (10 minutes)
```bash
jupyter notebook examples/tutorial.ipynb
```
Comprehensive walkthrough with explanations.

## 🔧 Configuration Options

### Cosmology
- `h` - Hubble parameter (0 < h < 1)
- `omega_m` - Total matter density
- `omega_b` - Baryon density
- `omega_lambda` - Dark energy density
- `n_s` - Scalar spectral index
- `sigma8` - Fluctuation amplitude

### Reionization
- `model` - 'tanh' or 'custom'
- `z_re` - Reionization redshift
- `delta_z` - Duration
- `custom_file` - Path to xe(z) file (for custom model)

### Bubbles
- `Rb` - Characteristic size (Mpc)
- `sigma_lnr` - Distribution width
- `b` - Linear bias
- `T0` - Temperature (K)

## 🎨 Plotting Examples

### Ionization Fraction
```python
import matplotlib.pyplot as plt

z = np.linspace(0, 20, 200)
xe = calc.xe(z)

plt.plot(z, xe)
plt.xlabel('Redshift z')
plt.ylabel('Ionization Fraction $x_e$')
plt.show()
```

### Power Spectrum
```python
k = np.logspace(-2, 1, 50)
pxe = calc.p_xexe(7.5, k)

plt.loglog(k, k**3 * pxe / (2*np.pi**2))
plt.xlabel('k [h/Mpc]')
plt.ylabel(r'$\Delta_{x_e}^2(k)$')
plt.show()
```

### Optical Depth Spectrum
```python
ell = np.logspace(np.log10(2), np.log10(5000), 20)
cl_tau = calc.cltau(ell)

plt.loglog(ell, ell*(ell+1)*cl_tau/(2*np.pi))
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1)C_\ell^{\tau\tau}/(2\pi)$')
plt.show()
```

## 🐛 Troubleshooting

### Installation Issues
```bash
# Make sure pip is up to date
pip install --upgrade pip

# Install with verbose output
pip install -e . -v
```

### Import Errors
```python
# Check installation
import tauspec
print(tauspec.__version__)  # Should print "0.2.0"
```

### Performance Issues
- Reduce `n_chi_samples` in `IntegrationConfig` for faster (less accurate) results
- Use fewer k or ell values for quick estimates
- Reuse calculator objects instead of creating new ones

### Integration Warnings
The warnings about integration subdivisions are normal for some parameter ranges.
They indicate the integrand is difficult but results are still reliable.

## 📖 Learning Path

1. **Beginner** (30 min)
   - Read SUMMARY.md
   - Run test_quick.py
   - Try basic example from Quick Start

2. **Intermediate** (1 hour)
   - Read README.md
   - Run examples/basic_usage.py
   - Modify configurations and explore

3. **Advanced** (2+ hours)
   - Work through tutorial.ipynb
   - Read PERFORMANCE.md
   - Try custom reionization history
   - Explore source code in tauspec/

## 🤝 Support

- **Documentation**: Start with SUMMARY.md
- **Examples**: See examples/ directory
- **API Reference**: See docstrings in source code
- **Migration**: See MIGRATION.md for upgrading from old code

## 📦 Dependencies

All dependencies are modern, well-maintained packages:
- `numpy >= 1.20.0` - Array operations
- `scipy >= 1.7.0` - Integration and interpolation
- `astropy >= 5.0.0` - Cosmology calculations

## ⚡ Key Features

- ✅ **Fast**: 2-5x faster than original
- ✅ **Modern**: Uses Astropy instead of deprecated cosmolopy
- ✅ **Reliable**: Type hints and validation
- ✅ **Flexible**: Easy to configure
- ✅ **Documented**: Comprehensive docs and examples
- ✅ **Tested**: Unit tests included
- ✅ **Installable**: Proper Python package

## 🎓 Scientific Usage

This package calculates optical depth power spectra during the Epoch of Reionization
based on the patchy reionization model. Key physics:

- **Ionization fraction**: xe(z) with tanh or custom model
- **Bubble model**: Log-normal size distribution
- **Power spectrum**: P_xe(k,z) with linear bias
- **Optical depth**: C_l^tau from line-of-sight integration

For methodology, see:
- Bharadwaj & Ali (2005): https://arxiv.org/abs/astro-ph/0511141

## 📝 Citation

If you use this code in research, please cite:

```bibtex
@software{tauspec,
  author = {Roy, Anirban},
  title = {TauSpec: Optical Depth Power Spectra Calculator},
  year = {2025},
  version = {0.2.0},
  url = {https://github.com/yourusername/tauspec}
}
```

## 🚀 Next Steps

1. Read [SUMMARY.md](SUMMARY.md) for complete overview
2. Run `python test_quick.py` to test installation
3. Try `python examples/basic_usage.py` for a real example
4. Explore `tutorial.ipynb` for comprehensive guide
5. Read [MIGRATION.md](MIGRATION.md) if updating from old code

## 📄 License

MIT License - Free to use, modify, and distribute.

---

**Ready to calculate some optical depth power spectra? Start with test_quick.py!** 🎉
