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
git clone https://github.com/Anirbancosmo/tauspec_pkg
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

```
@article{Roy:2018gcv,
    author = "Roy, A. and Lapi, A. and Spergel, D. and Baccigalupi, C.",
    title = "{Observing patchy reionization with future CMB polarization experiments}",
    eprint = "1801.02393",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1088/1475-7516/2018/05/014",
    journal = "JCAP",
    volume = "05",
    pages = "014",
    year = "2018"
}
```

and 

```
@article{Namikawa:2021zhh,
    author = "Namikawa, Toshiya and Roy, Anirban and Sherwin, Blake D. and Battaglia, Nicholas and Spergel, David N.",
    title = "{Constraining reionization with the first measurement of the cross-correlation between the CMB optical-depth fluctuations and the Compton y-map}",
    eprint = "2102.00975",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1103/PhysRevD.104.063514",
    journal = "Phys. Rev. D",
    volume = "104",
    number = "6",
    pages = "063514",
    year = "2021"
}
```




## References

- Roy et al. (2018): https://arxiv.org/abs/1801.02393
- Namikawa, Roy, et al. (2021): https://arxiv.org/abs/2102.00975
