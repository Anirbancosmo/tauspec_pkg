# Changelog

## [0.2.0] - 2025-10-22

### Major Changes
- Complete package restructure from single-file script to proper Python package
- Replaced cosmolopy (deprecated) with Astropy (modern, maintained)
- Performance optimizations: 2-5x faster overall

### Added
- Proper package structure with `tauspec` module
- Configuration system using dataclasses (`CosmologyConfig`, `ReionizationConfig`, `BubbleConfig`)
- Type hints throughout the codebase
- Comprehensive test suite in `tests/`
- Multiple usage examples in `examples/`
- Full documentation (README, MIGRATION, PERFORMANCE, INDEX, SUMMARY, QUICKREF)
- `pyproject.toml` for modern Python packaging
- `setup.py` for pip installation
- MIT License
- `.gitignore` and project configuration files

### Changed
- **Cosmology backend**: cosmolopy → Astropy
  - Pre-computed interpolation tables for z(chi), a(chi)
  - More accurate distance calculations
  - Better maintained and supported

- **API structure**: 
  - Old: `import tauspec as t; t.p_xexe(z, k)`
  - New: `calc = TauSpecCalculator(...); calc.p_xexe(z, k)`
  
- **Configuration**:
  - Old: Global `inputs.py` file
  - New: Type-safe dataclass configs passed to calculator

- **Function names**:
  - `Cltauy()` → `cltau_y()` (lowercase for consistency)
  - `return_output()` → `save_output()` (clearer name)

### Improved
- **Performance**:
  - Vectorized all array operations (3-5x faster)
  - LRU caching for power spectrum calculations
  - Pre-computed bubble distribution integrals
  - Optimized integration with Simpson's rule
  
- **Code quality**:
  - Modular design with separate files for concerns
  - Comprehensive docstrings
  - Error handling and validation
  - Clear variable names
  - PEP 8 compliant

- **Usability**:
  - Easy installation with pip
  - Clear error messages
  - Flexible output filenames
  - Support for both scalar and array inputs
  - Better documentation

### Documentation
- Added comprehensive README.md
- Added MIGRATION.md for users of old code
- Added PERFORMANCE.md with benchmarks
- Added INDEX.md as navigation guide
- Added SUMMARY.md with complete overview
- Added QUICKREF.md as quick reference
- Added inline docstrings throughout
- Added tutorial.ipynb Jupyter notebook
- Added basic_usage.py example script

### Testing
- Added test_config.py for configuration validation
- Added test_quick.py for functionality testing
- Added compare_performance.py for benchmarking
- Ready for pytest integration

### Dependencies
- **Removed**: cosmolopy (deprecated)
- **Added**: astropy >= 5.0.0
- **Updated**: numpy >= 1.20.0, scipy >= 1.7.0

### File Structure
```
Old:
tauspec.py
inputs.py
get_pk.py

New:
tauspec_pkg/
├── tauspec/
│   ├── __init__.py
│   ├── calculator.py
│   ├── cosmology.py
│   ├── reionization.py
│   └── config.py
├── tests/
├── examples/
├── Documentation
└── Configuration files
```

### Breaking Changes
- Requires explicit configuration objects instead of global inputs
- Different import structure: `from tauspec import TauSpecCalculator`
- Calculator must be instantiated before use
- Some function names changed for consistency
- Requires Astropy instead of cosmolopy

### Migration Path
See MIGRATION.md for detailed migration instructions from version 0.1.0.

---

## [0.1.0] - Original Version

Initial implementation with:
- Single tauspec.py file
- inputs.py for configuration
- cosmolopy for cosmology calculations
- Basic functionality for P_xe and C_l calculations
