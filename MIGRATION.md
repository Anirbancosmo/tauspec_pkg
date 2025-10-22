# Migration Guide: From Old to New TauSpec

## Quick Start

### Old Code:
```python
import tauspec as t
import inputs as inp

# Uses global inputs.py file
z = 7.5
k = np.logspace(-2, 0)
pxe = t.p_xexe(z, k)
```

### New Code:
```python
from tauspec import TauSpecCalculator
from tauspec.config import CosmologyConfig, ReionizationConfig, BubbleConfig

# Explicit configuration
cosmo_config = CosmologyConfig(h=0.69, omega_m=0.2877)
reion_config = ReionizationConfig(z_re=7.5, delta_z=2.0)
bubble_config = BubbleConfig(Rb=5.0, b=6.0)

calc = TauSpecCalculator(cosmo_config, reion_config, bubble_config)

z = 7.5
k = np.logspace(-2, 0)
pxe = calc.p_xexe(z, k)
```

## Key Differences

### 1. Package Structure
**Old**: Single `tauspec.py` file with `inputs.py` for configuration
**New**: Proper package with modules:
- `tauspec.calculator` - Main calculations
- `tauspec.cosmology` - Cosmology using Astropy
- `tauspec.reionization` - Reionization models
- `tauspec.config` - Configuration classes

### 2. Configuration
**Old**: Global variables in `inputs.py`
**New**: Dataclass-based configs passed to calculator

**Old inputs.py:**
```python
h = 0.69
z_re = 7.5
Rb = 5.0
```

**New config:**
```python
cosmo_config = CosmologyConfig(h=0.69)
reion_config = ReionizationConfig(z_re=7.5)
bubble_config = BubbleConfig(Rb=5.0)
```

### 3. Function Calls

| Old | New |
|-----|-----|
| `t.p_xexe(z, k)` | `calc.p_xexe(z, k)` |
| `t.xe(z)` | `calc.xe(z)` |
| `t.cltau(ell)` | `calc.cltau(ell)` |
| `t.Cltauy(ell)` | `calc.cltau_y(ell)` |
| `t.return_output(output)` | `calc.save_output(type, ...)` |

### 4. Cosmology Backend
**Old**: cosmolopy (deprecated)
**New**: Astropy (modern, maintained)

### 5. Output Files
**Old**: Fixed names (`pxe_z*.dat`, `cltau.dat`, `xe.dat`)
**New**: Customizable filenames via `save_output()`

## Complete Example Migration

### Old Code (using inputs.py):
```python
#!/usr/bin/env python3
import numpy as np
import tauspec as t
import matplotlib.pyplot as plt

# Configuration in inputs.py
# (h=0.69, z_re=7.5, Rb=5, etc.)

# Calculate
z = 7.5
k = np.logspace(-2, 0)
pxe = t.p_xexe(z, k)

# Save
if t.inp.write_output:
    t.return_output(t.inp.output_quantity)

# Plot
plt.loglog(k, pxe)
plt.show()
```

### New Code:
```python
#!/usr/bin/env python3
import numpy as np
from tauspec import TauSpecCalculator
from tauspec.config import CosmologyConfig, ReionizationConfig, BubbleConfig
import matplotlib.pyplot as plt

# Configuration
cosmo_config = CosmologyConfig(h=0.69, omega_m=0.2877, omega_b=0.04757)
reion_config = ReionizationConfig(z_re=7.5, delta_z=2.0)
bubble_config = BubbleConfig(Rb=5.0, sigma_lnr=np.log(2), b=6.0)

# Initialize calculator
calc = TauSpecCalculator(cosmo_config, reion_config, bubble_config)

# Calculate
z = 7.5
k = np.logspace(-2, 0)
pxe = calc.p_xexe(z, k)

# Save
calc.save_output('pxe', z=z, k=k, filename='pxe_z7.5.dat')
calc.save_output('cltau', ell=np.logspace(1, 3, 20))
calc.save_output('xe')

# Plot
plt.loglog(k, pxe)
plt.show()
```

## Mapping of Old inputs.py Parameters

| Old (inputs.py) | New Config | Class |
|----------------|------------|-------|
| `h` | `h` | CosmologyConfig |
| `omega_m_0` | `omega_m` | CosmologyConfig |
| `omega_b_0` | `omega_b` | CosmologyConfig |
| `omega_lambda` | `omega_lambda` | CosmologyConfig |
| `ns` | `n_s` | CosmologyConfig |
| `sigma8` | `sigma8` | CosmologyConfig |
| `reio_model` | `model` | ReionizationConfig |
| `z_re` | `z_re` | ReionizationConfig |
| `delta_z` | `delta_z` | ReionizationConfig |
| `reio_file_name` | `custom_file` | ReionizationConfig |
| `Rb` | `Rb` | BubbleConfig |
| `sigma_lnr` | `sigma_lnr` | BubbleConfig |
| `b` | `b` | BubbleConfig |
| `T0` | `T0` | BubbleConfig |

## Custom Reionization History

### Old:
```python
# In inputs.py
reio_model = "custom"
reio_file_name = '/path/to/xe.dat'
```

### New:
```python
reion_config = ReionizationConfig(
    model='custom',
    custom_file='/path/to/xe.dat'
)
```

## Benefits of New Version

1. **Type Safety**: Dataclass configs catch errors early
2. **Modularity**: Easier to extend and test
3. **Modern Python**: Uses current best practices
4. **Better Performance**: Optimized with caching and vectorization
5. **Maintainability**: Uses maintained dependencies (Astropy)
6. **Flexibility**: Multiple calculators with different configs
7. **Documentation**: Comprehensive docstrings and examples

## Installation

### Old:
```bash
# Just copy tauspec.py and inputs.py
# Install cosmolopy (deprecated)
```

### New:
```bash
pip install -e .
# Or from PyPI (when published):
# pip install tauspec
```

## Testing

### Old:
```bash
python tauspec.py  # If write_output=True
```

### New:
```bash
pytest tests/
python examples/basic_usage.py
jupyter notebook examples/tutorial.ipynb
```

## Questions?

See the full documentation in README.md and examples/ directory.
