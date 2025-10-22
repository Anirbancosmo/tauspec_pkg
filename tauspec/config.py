"""Configuration classes for TauSpec package."""

from dataclasses import dataclass, field
from typing import Optional
import numpy as np


@dataclass
class CosmologyConfig:
    """Cosmological parameters configuration.
    
    Attributes:
        h: Hubble parameter (H0 = 100h km/s/Mpc)
        omega_m: Total matter density parameter
        omega_b: Baryon density parameter
        omega_lambda: Dark energy density parameter
        omega_k: Curvature density parameter
        n_s: Scalar spectral index
        sigma8: Amplitude of matter fluctuations at 8 Mpc/h
    """
    h: float = 0.69
    omega_m: float = 0.2877
    omega_b: float = 0.04757
    omega_lambda: float = 0.721
    omega_k: float = 0.0
    n_s: float = 0.9667
    sigma8: float = 0.80
    
    def __post_init__(self):
        """Validate cosmological parameters."""
        if not 0 < self.h < 1:
            raise ValueError("h must be between 0 and 1")
        if not 0 < self.omega_m < 1:
            raise ValueError("omega_m must be between 0 and 1")
        if not 0 < self.omega_b < self.omega_m:
            raise ValueError("omega_b must be positive and less than omega_m")


@dataclass
class ReionizationConfig:
    """Reionization history configuration.
    
    Attributes:
        model: Type of reionization model ('tanh' or 'custom')
        z_re: Reionization redshift (when xe=0.5)
        delta_z: Duration of reionization
        custom_file: Path to custom reionization history file (z, xe columns)
    """
    model: str = "tanh"
    z_re: float = 7.5
    delta_z: float = 2.0
    custom_file: Optional[str] = None
    
    def __post_init__(self):
        """Validate reionization parameters."""
        if self.model not in ['tanh', 'custom']:
            raise ValueError("model must be 'tanh' or 'custom'")
        if self.model == 'custom' and self.custom_file is None:
            raise ValueError("custom_file must be provided for custom model")
        if self.z_re < 0:
            raise ValueError("z_re must be non-negative")
        if self.delta_z <= 0:
            raise ValueError("delta_z must be positive")


@dataclass
class BubbleConfig:
    """Bubble size distribution and bias configuration.
    
    Attributes:
        Rb: Characteristic bubble size (Mpc)
        sigma_lnr: Width of log-normal bubble size distribution
        b: Linear bias parameter
        T0: Temperature during EoR (Kelvin)
        R_int_min: Minimum radius for integration (Mpc)
        R_int_max: Maximum radius for integration (Mpc)
    """
    Rb: float = 5.0
    sigma_lnr: float = np.log(2)
    b: float = 6.0
    T0: float = 2e4
    R_int_min: float = 0.1
    R_int_max: float = 500.0
    
    def __post_init__(self):
        """Validate bubble parameters."""
        if self.Rb <= 0:
            raise ValueError("Rb must be positive")
        if self.sigma_lnr <= 0:
            raise ValueError("sigma_lnr must be positive")
        if self.T0 <= 0:
            raise ValueError("T0 must be positive")
        if self.R_int_min >= self.R_int_max:
            raise ValueError("R_int_min must be less than R_int_max")


@dataclass
class IntegrationConfig:
    """Integration parameters configuration.
    
    Attributes:
        z_min: Minimum redshift for integration
        z_max: Maximum redshift for integration
        n_chi_samples: Number of comoving distance samples
        n_k_samples: Number of wavenumber samples for interpolation
    """
    z_min: float = 5.0
    z_max: float = 20.0
    n_chi_samples: int = 50
    n_k_samples: int = 100
    
    def __post_init__(self):
        """Validate integration parameters."""
        if self.z_min >= self.z_max:
            raise ValueError("z_min must be less than z_max")
        if self.n_chi_samples < 10:
            raise ValueError("n_chi_samples must be at least 10")
