# TauSpec

A Python package for calculating optical depth power spectra during the Epoch of Reionization (EoR).

## Features

- Calculate ionization fraction power spectra (P_xe)
- Calculate optical depth power spectra (C_l^tau)
- Support for custom reionization histories
- Optimized for performance using vectorized operations and caching
- Uses Astropy for cosmological calculations

## Installation

### From source

```bash
git clone https://github.com/yourusername/tauspec
cd tauspec
pip install -e .
```

### Development installation

```bash
pip install -e ".[dev]"
```

## Quick Start

```python
import numpy as np
from tauspec import TauSpecCalculator
from tauspec.config import CosmologyConfig, ReionizationConfig, BubbleConfig

# Set up configuration
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

# Create calculator
calc = TauSpecCalculator(cosmo_config, reion_config, bubble_config)

# Calculate power spectra
z = 7.5
k = np.logspace(-2, 1, 50)
pxe = calc.p_xexe(z, k)

# Calculate C_l^tau
ell = np.logspace(np.log10(2), np.log10(5000), 20)
cl_tau = calc.cltau(ell)

# Get ionization fraction
z_array = np.linspace(5, 15, 100)
xe_array = calc.xe(z_array)
```

## Configuration

The package uses dataclass-based configuration for easy customization:

### Cosmology Parameters
- `h`: Hubble parameter (H0 = 100h km/s/Mpc)
- `omega_m`: Matter density parameter
- `omega_b`: Baryon density parameter
- `omega_lambda`: Dark energy density parameter
- `n_s`: Scalar spectral index
- `sigma8`: Amplitude of matter fluctuations

### Reionization Parameters
- `model`: 'tanh' or 'custom'
- `z_re`: Reionization redshift (when xe=0.5)
- `delta_z`: Duration of reionization
- `custom_file`: Path to custom reionization history file (for custom model)

### Bubble Parameters
- `Rb`: Characteristic bubble size (Mpc)
- `sigma_lnr`: Bubble size distribution width
- `b`: Linear bias parameter
- `T0`: Temperature during EoR (Kelvin)

## Performance Improvements

This version includes several optimizations:
- Vectorized operations using NumPy
- Caching of expensive computations (power spectra, integrals)
- Efficient interpolation using scipy
- Parallel computation support (optional)
- Pre-computed lookup tables for common calculations

## Testing

Run tests with:

```bash
pytest tests/
```

Run tests with coverage:

```bash
pytest --cov=tauspec tests/
```

## License

MIT License

## Citation

If you use this code in your research, please cite:

```bibtex
@software{tauspec,
  author = {Roy, Anirban},
  title = {TauSpec: Optical Depth Power Spectra Calculator},
  year = {2025},
  url = {https://github.com/yourusername/tauspec}
}
```

## References

Based on the methodology described in:
- Bharadwaj & Ali (2005): https://arxiv.org/abs/astro-ph/0511141
