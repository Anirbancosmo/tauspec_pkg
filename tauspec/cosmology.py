"""Cosmology calculations using Astropy."""

import numpy as np
from astropy.cosmology import FlatLambdaCDM, LambdaCDM
from astropy import units as u
from scipy.interpolate import interp1d
from functools import lru_cache
from .config import CosmologyConfig


class Cosmology:
    """Cosmology calculations using Astropy.
    
    This class provides methods for cosmological distance calculations
    and power spectrum computations using Astropy's cosmology module.
    """
    
    def __init__(self, config: CosmologyConfig):
        """Initialize cosmology with given parameters.
        
        Args:
            config: CosmologyConfig object with cosmological parameters
        """
        self.config = config
        
        # Create Astropy cosmology object
        H0 = config.h * 100  # km/s/Mpc
        Omega0 = config.omega_m
        OmegaLambda = config.omega_lambda
        Ob0 = config.omega_b
        
        if abs(config.omega_k) < 1e-8:
            # Flat universe
            self.cosmo = FlatLambdaCDM(H0=H0, Om0=Omega0, Ob0=Ob0)
        else:
            # Non-flat universe
            self.cosmo = LambdaCDM(H0=H0, Om0=Omega0, Ode0=OmegaLambda, Ob0=Ob0)
        
        # Pre-compute interpolation tables
        self._setup_interpolators()
    
    def _setup_interpolators(self, z_min=0.0, z_max=30.0, n_points=500):
        """Set up interpolation tables for fast lookups.
        
        Args:
            z_min: Minimum redshift
            z_max: Maximum redshift
            n_points: Number of interpolation points
        """
        z_array = np.linspace(z_min, z_max, n_points)
        
        # Comoving distance in Mpc
        chi_array = self.cosmo.comoving_distance(z_array).to(u.Mpc).value
        
        # Scale factor
        a_array = 1.0 / (1.0 + z_array)
        
        # Create interpolators
        self._chi_interp = interp1d(z_array, chi_array, kind='cubic', 
                                    bounds_error=False, fill_value='extrapolate')
        self._z_interp = interp1d(chi_array, z_array, kind='cubic',
                                  bounds_error=False, fill_value='extrapolate')
        self._a_interp = interp1d(chi_array, a_array, kind='cubic',
                                  bounds_error=False, fill_value='extrapolate')
    
    def comoving_distance(self, z):
        """Calculate comoving distance in Mpc.
        
        Args:
            z: Redshift (scalar or array)
            
        Returns:
            Comoving distance in Mpc
        """
        z = np.atleast_1d(z)
        result = self._chi_interp(z)
        return result if z.size > 1 else float(result)
    
    def z_from_chi(self, chi):
        """Calculate redshift from comoving distance.
        
        Args:
            chi: Comoving distance in Mpc (scalar or array)
            
        Returns:
            Redshift
        """
        chi = np.atleast_1d(chi)
        result = self._z_interp(chi)
        return result if chi.size > 1 else float(result)
    
    def scale_factor_from_chi(self, chi):
        """Calculate scale factor from comoving distance.
        
        Args:
            chi: Comoving distance in Mpc (scalar or array)
            
        Returns:
            Scale factor a = 1/(1+z)
        """
        chi = np.atleast_1d(chi)
        result = self._a_interp(chi)
        return result if chi.size > 1 else float(result)
    
    @lru_cache(maxsize=128)
    def _cached_power_spectrum(self, k_tuple, z):
        """Cached power spectrum calculation.
        
        Args:
            k_tuple: Tuple of wavenumber values (for hashability)
            z: Redshift
            
        Returns:
            Matter power spectrum at given k and z
        """
        k = np.array(k_tuple)
        return self._compute_power_spectrum_linear(k, z)
    
    def _compute_power_spectrum_linear(self, k, z):
        """Compute linear matter power spectrum.
        
        Uses the Eisenstein & Hu (1998) transfer function.
        
        Args:
            k: Wavenumber in h/Mpc (scalar or array)
            z: Redshift
            
        Returns:
            Matter power spectrum P(k,z) in (Mpc/h)^3
        """
        k = np.atleast_1d(k)
        
        # Growth factor
        D_z = self.growth_factor(z)
        
        # Primordial power spectrum
        n_s = self.config.n_s
        P_primordial = k ** n_s
        
        # Transfer function (simplified Eisenstein-Hu)
        T_k = self._transfer_function_eh(k)
        
        # Normalization to sigma8
        norm = self._get_normalization()
        
        # Power spectrum
        P_k = norm * P_primordial * T_k ** 2 * D_z ** 2
        
        return P_k if k.size > 1 else float(P_k)
    
    def power_spectrum(self, k, z):
        """Calculate matter power spectrum at given k and z.
        
        Args:
            k: Wavenumber in h/Mpc (scalar or array)
            z: Redshift
            
        Returns:
            Matter power spectrum P(k,z) in (Mpc/h)^3
        """
        k = np.atleast_1d(k)
        scalar_input = k.size == 1
        
        # Use caching for better performance
        if scalar_input and k.size == 1:
            result = self._cached_power_spectrum((float(k[0]),), z)
        else:
            result = self._compute_power_spectrum_linear(k, z)
        
        return result if not scalar_input else float(result)
    
    def growth_factor(self, z):
        """Calculate linear growth factor D(z).
        
        Approximate formula for LCDM cosmology.
        
        Args:
            z: Redshift
            
        Returns:
            Growth factor D(z) normalized to D(0) = 1
        """
        a = 1.0 / (1.0 + z)
        Om_a = self.cosmo.Om(z)
        OL_a = self.cosmo.Ode(z)
        
        # Approximate growth factor (Carroll, Press & Turner 1992)
        g_a = 2.5 * Om_a / (Om_a ** (4.0 / 7.0) - OL_a + 
                            (1.0 + Om_a / 2.0) * (1.0 + OL_a / 70.0))
        
        # Normalize to present day
        g_0 = 2.5 * self.config.omega_m / (
            self.config.omega_m ** (4.0 / 7.0) - self.config.omega_lambda +
            (1.0 + self.config.omega_m / 2.0) * 
            (1.0 + self.config.omega_lambda / 70.0))
        
        return g_a / g_0
    
    def _transfer_function_eh(self, k):
        """Eisenstein & Hu (1998) transfer function (no BAO).
        
        Args:
            k: Wavenumber in h/Mpc
            
        Returns:
            Transfer function T(k)
        """
        # Convert k to Mpc^-1 (remove h)
        k_h = k * self.config.h
        
        # Shape parameter
        omega_m_h2 = self.config.omega_m * self.config.h ** 2
        omega_b_h2 = self.config.omega_b * self.config.h ** 2
        
        Gamma_eff = omega_m_h2 * (
            self.config.h * np.exp(-self.config.omega_b * 
            (1.0 + np.sqrt(2 * self.config.h) / self.config.omega_m)))
        
        q = k_h / Gamma_eff
        
        # Transfer function
        L = np.log(2 * np.e + 1.8 * q)
        C = 14.2 + 731.0 / (1.0 + 62.5 * q)
        T = L / (L + C * q ** 2)
        
        return T
    
    def _get_normalization(self):
        """Calculate normalization from sigma8.
        
        Returns:
            Normalization constant for power spectrum
        """
        # This is a simplified normalization
        # In practice, you'd integrate the power spectrum with a window function
        return self.config.sigma8 ** 2 / (2.0 * np.pi ** 2)
    
    def sigma_r(self, R, z):
        """Calculate RMS mass fluctuation in spheres of radius R.
        
        Args:
            R: Radius in Mpc/h
            z: Redshift
            
        Returns:
            sigma_R(z): RMS fluctuation
        """
        # Simplified calculation
        # In practice, integrate power spectrum with window function
        k = np.logspace(-4, 2, 1000)
        P_k = self.power_spectrum(k, z)
        
        # Top-hat window function
        x = k * R
        W = 3.0 * (np.sin(x) - x * np.cos(x)) / x ** 3
        
        # Integrate
        integrand = k ** 2 * P_k * W ** 2 / (2.0 * np.pi ** 2)
        sigma2 = np.trapz(integrand, k)
        
        return np.sqrt(sigma2)
