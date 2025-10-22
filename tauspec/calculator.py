"""Main calculator for optical depth power spectra.

This module contains the TauSpecCalculator class which provides methods
for calculating ionization fraction power spectra and optical depth
power spectra during reionization.
"""

import numpy as np
from scipy.integrate import quad, simpson
from scipy.interpolate import interp1d
from functools import lru_cache
from typing import Union, Tuple

from .config import CosmologyConfig, ReionizationConfig, BubbleConfig, IntegrationConfig
from .cosmology import Cosmology
from .reionization import ReionizationHistory


# Physical constants
C_LIGHT_M = 3e8  # Speed of light in m/s
K_B_OVER_MEC2 = 1.689e-10  # k_B/(m_e c^2) in K^{-1}
N_P0 = 0.22 * (3.086e+22) ** 3  # Proton number density (Mpc^-3)
SIGMA_T = 6.652e-29 / (3.086e+22) ** 2  # Thomson cross section (Mpc^2)
PI = np.pi


class TauSpecCalculator:
    """Calculator for optical depth power spectra.
    
    This class provides methods to calculate:
    - Ionization fraction power spectra P_xe(k, z)
    - Optical depth power spectra C_l^tau
    - Temperature-weighted optical depth power spectra C_l^tau_y
    
    The calculations are optimized using:
    - Vectorized NumPy operations
    - Caching of expensive computations
    - Pre-computed interpolation tables
    """
    
    def __init__(
        self,
        cosmo_config: CosmologyConfig,
        reion_config: ReionizationConfig,
        bubble_config: BubbleConfig,
        integration_config: IntegrationConfig = None
    ):
        """Initialize the calculator.
        
        Args:
            cosmo_config: Cosmology configuration
            reion_config: Reionization configuration
            bubble_config: Bubble size distribution configuration
            integration_config: Integration parameters (optional)
        """
        self.cosmo_config = cosmo_config
        self.reion_config = reion_config
        self.bubble_config = bubble_config
        self.integration_config = integration_config or IntegrationConfig()
        
        # Initialize cosmology and reionization
        self.cosmo = Cosmology(cosmo_config)
        self.reion = ReionizationHistory(reion_config)
        
        # Pre-compute some values
        self._sigma_t_np0_sq = (SIGMA_T * N_P0) ** 2
        
        # Cache for bubble distribution functions
        self._V_av_cache = None
        self._I_cache = {}
        self._F_cache = {}
    
    def xe(self, z: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Calculate ionization fraction.
        
        Args:
            z: Redshift (scalar or array)
            
        Returns:
            Ionization fraction xe(z)
        """
        return self.reion.xe(z)
    
    def te(self, z: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Calculate electron temperature.
        
        Args:
            z: Redshift (scalar or array)
            
        Returns:
            Electron temperature in Kelvin
        """
        # Constant temperature model
        z = np.atleast_1d(z)
        result = np.full_like(z, self.bubble_config.T0, dtype=float)
        return result if z.size > 1 else float(result)
    
    def _P_bubble(self, R: np.ndarray) -> np.ndarray:
        """Bubble size distribution function.
        
        Log-normal distribution:
        P(R) = 1/(R sqrt(2π σ_lnR^2)) exp(-(ln(R/Rb))^2/(2σ_lnR^2))
        
        Args:
            R: Bubble radius in Mpc (array)
            
        Returns:
            Probability density P(R)
        """
        Rb = self.bubble_config.Rb
        sigma_lnr = self.bubble_config.sigma_lnr
        
        norm = 1.0 / (R * np.sqrt(2 * PI * sigma_lnr ** 2))
        exponent = -(np.log(R / Rb) ** 2) / (2 * sigma_lnr ** 2)
        
        return norm * np.exp(exponent)
    
    def _W_tophat(self, R: np.ndarray, k: np.ndarray) -> np.ndarray:
        """Top-hat window function in Fourier space.
        
        W(R,k) = 3/(kR)^3 [sin(kR) - kR cos(kR)]
        
        Args:
            R: Radius in Mpc (array)
            k: Wavenumber in h/Mpc (array)
            
        Returns:
            Window function W(R,k)
        """
        kR = np.outer(k, R) if R.ndim == 1 and k.ndim == 1 else k * R
        
        # Avoid division by zero
        with np.errstate(divide='ignore', invalid='ignore'):
            W = 3.0 * (np.sin(kR) - kR * np.cos(kR)) / kR ** 3
            # Handle kR -> 0 limit: W -> 1
            W = np.where(np.abs(kR) < 1e-6, 1.0, W)
        
        return W
    
    def _V_bubble(self, R: np.ndarray) -> np.ndarray:
        """Bubble volume.
        
        Args:
            R: Bubble radius in Mpc (array)
            
        Returns:
            Volume in Mpc^3
        """
        return (4.0 / 3.0) * PI * R ** 3
    
    def _V_average(self) -> float:
        """Calculate average bubble volume.
        
        <V> = ∫ P(R) V(R) dR
        
        Returns:
            Average volume in Mpc^3
        """
        if self._V_av_cache is not None:
            return self._V_av_cache
        
        R_min = self.bubble_config.R_int_min
        R_max = self.bubble_config.R_int_max
        
        def integrand(R):
            return self._P_bubble(np.array([R]))[0] * self._V_bubble(np.array([R]))[0]
        
        result, _ = quad(integrand, R_min, R_max, limit=100)
        self._V_av_cache = result
        
        return result
    
    def _I_function(self, k: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Calculate I(k) function.
        
        I(k) = b/<V> ∫ P(R) V(R) W(R,k) dR
        
        Args:
            k: Wavenumber in h/Mpc (scalar or array)
            
        Returns:
            I(k) values
        """
        k = np.atleast_1d(k)
        scalar_input = k.size == 1
        
        # Check cache
        result = np.zeros_like(k, dtype=float)
        for i, k_val in enumerate(k):
            k_key = round(k_val, 8)  # Round for cache key
            if k_key in self._I_cache:
                result[i] = self._I_cache[k_key]
            else:
                # Compute
                R_min = self.bubble_config.R_int_min
                R_max = self.bubble_config.R_int_max
                
                def integrand(R):
                    P = self._P_bubble(np.array([R]))[0]
                    V = self._V_bubble(np.array([R]))[0]
                    W = self._W_tophat(np.array([R]), np.array([k_val]))[0]
                    return P * V * W
                
                integral, _ = quad(integrand, R_min, R_max, limit=100)
                value = self.bubble_config.b * integral / self._V_average()
                
                self._I_cache[k_key] = value
                result[i] = value
        
        return result if not scalar_input else float(result[0])
    
    def _F_function(self, k: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Calculate F(k) function.
        
        F(k) = 1/<V> ∫ P(R) V(R)^2 W(R,k)^2 dR
        
        Args:
            k: Wavenumber in h/Mpc (scalar or array)
            
        Returns:
            F(k) values
        """
        k = np.atleast_1d(k)
        scalar_input = k.size == 1
        
        # Check cache
        result = np.zeros_like(k, dtype=float)
        for i, k_val in enumerate(k):
            k_key = round(k_val, 8)
            if k_key in self._F_cache:
                result[i] = self._F_cache[k_key]
            else:
                # Compute
                R_min = self.bubble_config.R_int_min
                R_max = self.bubble_config.R_int_max
                
                def integrand(R):
                    P = self._P_bubble(np.array([R]))[0]
                    V = self._V_bubble(np.array([R]))[0]
                    W = self._W_tophat(np.array([R]), np.array([k_val]))[0]
                    return P * V ** 2 * W ** 2
                
                integral, _ = quad(integrand, R_min, R_max, limit=100)
                value = integral / self._V_average()
                
                self._F_cache[k_key] = value
                result[i] = value
        
        return result if not scalar_input else float(result[0])
    
    def _G_function(self, k: Union[float, np.ndarray], z: float) -> Union[float, np.ndarray]:
        """Calculate G(k,z) function.
        
        Non-linear correction using approximation from Bharadwaj & Ali (2005).
        
        Args:
            k: Wavenumber in h/Mpc (scalar or array)
            z: Redshift
            
        Returns:
            G(k,z) values
        """
        k = np.atleast_1d(k)
        
        P_k = self.cosmo.power_spectrum(k, z)
        sigma_r_sq = self.cosmo.sigma_r(self.bubble_config.Rb, z) ** 2
        V_av = self._V_average()
        
        V_sigma = V_av * sigma_r_sq
        
        result = P_k * V_sigma / np.sqrt(P_k ** 2 + V_sigma ** 2)
        
        return result if k.size > 1 else float(result)
    
    def p_xexe(
        self,
        z: float,
        k: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """Calculate ionization fraction power spectrum.
        
        P_xe(k,z) = xe(1-xe)[F(k) + G(k,z)] + 
                    [(1-xe)ln(1-xe)I(k) - xe]^2 P(k,z)
        
        Args:
            z: Redshift
            k: Wavenumber in h/Mpc (scalar or array)
            
        Returns:
            Power spectrum P_xe(k,z) in (Mpc/h)^3
        """
        k = np.atleast_1d(k)
        scalar_input = k.size == 1
        
        # Get ionization fraction
        xe = self.xe(z)
        
        # Get functions
        I_k = self._I_function(k)
        F_k = self._F_function(k)
        G_k = self._G_function(k, z)
        P_k = self.cosmo.power_spectrum(k, z)
        
        # Calculate power spectrum
        term1 = xe * (1 - xe) * (F_k + G_k)
        term2_factor = (1 - xe) * np.log(np.maximum(1 - xe, 1e-10)) * I_k - xe
        term2 = term2_factor ** 2 * P_k
        
        result = term1 + term2
        
        return result if not scalar_input else float(result)
    
    def cltau(
        self,
        ell: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """Calculate optical depth power spectrum.
        
        C_l^tau = ∫ (σ_T n_p0)^2 / (χ^2 a^4) P_xe(k=l/χ, z(χ)) dχ
        
        Args:
            ell: Multipole (scalar or array)
            
        Returns:
            Optical depth power spectrum C_l^tau
        """
        ell = np.atleast_1d(ell)
        scalar_input = ell.size == 1
        
        # Set up integration grid
        z_min = self.integration_config.z_min
        z_max = self.integration_config.z_max
        n_samples = self.integration_config.n_chi_samples
        
        chi_min = self.cosmo.comoving_distance(z_min)
        chi_max = self.cosmo.comoving_distance(z_max)
        
        # Use log spacing for better sampling
        chi_array = np.logspace(
            np.log10(chi_min + 1),
            np.log10(chi_max - 1),
            n_samples
        )
        
        # Create grid for vectorized computation
        chi_grid = chi_array[:, np.newaxis]  # Shape: (n_chi, 1)
        ell_grid = ell[np.newaxis, :]  # Shape: (1, n_ell)
        
        # Calculate k = ell/chi for each combination
        k_grid = ell_grid / chi_grid  # Shape: (n_chi, n_ell)
        
        # Get redshifts and scale factors
        z_array = self.cosmo.z_from_chi(chi_array)
        a_array = self.cosmo.scale_factor_from_chi(chi_array)
        
        # Calculate P_xe for each (chi, ell) pair
        integrand = np.zeros_like(k_grid)
        for i, (z_val, a_val, chi_val) in enumerate(zip(z_array, a_array, chi_array)):
            k_values = k_grid[i, :]
            P_xe = self.p_xexe(z_val, k_values)
            integrand[i, :] = (self._sigma_t_np0_sq / (chi_val ** 2 * a_val ** 4)) * P_xe
        
        # Integrate using Simpson's rule
        result = simpson(integrand, x=chi_array, axis=0)
        
        return result if not scalar_input else float(result)
    
    def cltau_y(
        self,
        ell: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """Calculate temperature-weighted optical depth power spectrum.
        
        C_l^tau_y = ∫ k_B/(m_e c^2) (σ_T n_p0)^2 T_e / (χ^2 a^4) P_xe(k=l/χ, z(χ)) dχ
        
        Args:
            ell: Multipole (scalar or array)
            
        Returns:
            Temperature-weighted optical depth power spectrum C_l^tau_y
        """
        ell = np.atleast_1d(ell)
        scalar_input = ell.size == 1
        
        # Set up integration grid
        z_min = self.integration_config.z_min
        z_max = self.integration_config.z_max
        n_samples = self.integration_config.n_chi_samples
        
        chi_min = self.cosmo.comoving_distance(z_min)
        chi_max = self.cosmo.comoving_distance(z_max)
        
        chi_array = np.logspace(
            np.log10(chi_min + 1),
            np.log10(chi_max - 1),
            n_samples
        )
        
        # Create grid
        chi_grid = chi_array[:, np.newaxis]
        ell_grid = ell[np.newaxis, :]
        k_grid = ell_grid / chi_grid
        
        # Get arrays
        z_array = self.cosmo.z_from_chi(chi_array)
        a_array = self.cosmo.scale_factor_from_chi(chi_array)
        te_array = self.te(z_array)
        
        # Calculate integrand
        integrand = np.zeros_like(k_grid)
        for i, (z_val, a_val, chi_val, te_val) in enumerate(
            zip(z_array, a_array, chi_array, te_array)
        ):
            k_values = k_grid[i, :]
            P_xe = self.p_xexe(z_val, k_values)
            integrand[i, :] = (
                K_B_OVER_MEC2 * self._sigma_t_np0_sq * te_val / 
                (chi_val ** 2 * a_val ** 4)
            ) * P_xe
        
        # Integrate
        result = simpson(integrand, x=chi_array, axis=0)
        
        return result if not scalar_input else float(result)
    
    def save_output(
        self,
        output_type: str,
        z: float = None,
        k: np.ndarray = None,
        ell: np.ndarray = None,
        filename: str = None
    ):
        """Save calculation output to file.
        
        Args:
            output_type: Type of output ('pxe', 'cltau', 'xe')
            z: Redshift for power spectrum (for 'pxe')
            k: Wavenumber array (for 'pxe')
            ell: Multipole array (for 'cltau')
            filename: Output filename (optional)
        """
        if output_type == 'pxe':
            if z is None:
                z = self.reion_config.z_re
            if k is None:
                k = np.logspace(-2, 1, 50)
            
            pxe = self.p_xexe(z, k)
            
            if filename is None:
                filename = f"pxe_z{z:.2f}.dat"
            
            np.savetxt(
                filename,
                np.column_stack([k, pxe]),
                header="k [h/Mpc]\t\tP_xe [(Mpc/h)^3]",
                fmt='%e'
            )
            
        elif output_type == 'cltau':
            if ell is None:
                ell = np.logspace(np.log10(2), np.log10(5000), 20)
            
            cl = self.cltau(ell)
            
            if filename is None:
                filename = "cltau.dat"
            
            np.savetxt(
                filename,
                np.column_stack([ell, cl]),
                header="ell\t\tC_l^tau",
                fmt='%e'
            )
            
        elif output_type == 'xe':
            z_array = np.linspace(0, 20, 100)
            xe_array = self.xe(z_array)
            
            if filename is None:
                filename = "xe.dat"
            
            np.savetxt(
                filename,
                np.column_stack([z_array, xe_array]),
                header="z\t\txe",
                fmt='%e'
            )
