"""Reionization history models."""

import numpy as np
from scipy.interpolate import interp1d
from .config import ReionizationConfig


class ReionizationHistory:
    """Reionization history models.
    
    Provides methods to calculate ionization fraction as a function of redshift
    using different models (tanh or custom).
    """
    
    def __init__(self, config: ReionizationConfig):
        """Initialize reionization history.
        
        Args:
            config: ReionizationConfig object with reionization parameters
        """
        self.config = config
        
        if config.model == 'custom':
            self._load_custom_history()
    
    def _load_custom_history(self):
        """Load custom reionization history from file."""
        data = np.loadtxt(self.config.custom_file)
        z_file = data[:, 0]
        xe_file = data[:, 1]
        
        # Create interpolator
        self._xe_interp = interp1d(z_file, xe_file, kind='cubic',
                                   bounds_error=False, 
                                   fill_value=(xe_file[-1], xe_file[0]))
    
    def xe(self, z):
        """Calculate ionization fraction at given redshift.
        
        Args:
            z: Redshift (scalar or array)
            
        Returns:
            Ionization fraction xe(z)
        """
        z = np.atleast_1d(z)
        
        if self.config.model == 'tanh':
            result = self._xe_tanh(z)
        elif self.config.model == 'custom':
            result = self._xe_interp(z)
        else:
            raise ValueError(f"Unknown model: {self.config.model}")
        
        return result if z.size > 1 else float(result)
    
    def _xe_tanh(self, z):
        """Calculate ionization fraction using tanh model.
        
        Args:
            z: Redshift (array)
            
        Returns:
            Ionization fraction xe(z)
        """
        z_re = self.config.z_re
        delta_z = self.config.delta_z
        
        # Following convention from original code
        delta_y = 1.5 * np.sqrt(1 + z_re) * delta_z
        y = (1 + z) ** (3.0 / 2.0)
        y_re = (1 + z_re) ** (3.0 / 2.0)
        
        xe = 0.5 * (1 + np.tanh((y_re - y) / delta_y))
        
        return xe
    
    def tau(self, z_start, z_end=0.0, n_points=1000):
        """Calculate optical depth between two redshifts.
        
        Args:
            z_start: Starting redshift
            z_end: Ending redshift (default: 0)
            n_points: Number of integration points
            
        Returns:
            Optical depth tau
        """
        # This would require cosmology object
        # Placeholder for future implementation
        raise NotImplementedError("tau calculation requires Cosmology object")
