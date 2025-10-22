"""Tests for configuration classes."""

import pytest
import numpy as np
from tauspec.config import (
    CosmologyConfig,
    ReionizationConfig,
    BubbleConfig,
    IntegrationConfig
)


class TestCosmologyConfig:
    """Tests for CosmologyConfig."""
    
    def test_default_values(self):
        """Test default cosmology configuration."""
        config = CosmologyConfig()
        assert config.h == 0.69
        assert config.omega_m == 0.2877
        assert config.omega_b == 0.04757
    
    def test_custom_values(self):
        """Test custom cosmology configuration."""
        config = CosmologyConfig(h=0.70, omega_m=0.30)
        assert config.h == 0.70
        assert config.omega_m == 0.30
    
    def test_validation(self):
        """Test parameter validation."""
        with pytest.raises(ValueError):
            CosmologyConfig(h=1.5)  # Invalid h
        
        with pytest.raises(ValueError):
            CosmologyConfig(omega_m=1.5)  # Invalid omega_m
        
        with pytest.raises(ValueError):
            CosmologyConfig(omega_b=0.5, omega_m=0.3)  # omega_b > omega_m


class TestReionizationConfig:
    """Tests for ReionizationConfig."""
    
    def test_tanh_model(self):
        """Test tanh model configuration."""
        config = ReionizationConfig(model='tanh', z_re=7.5, delta_z=2.0)
        assert config.model == 'tanh'
        assert config.z_re == 7.5
    
    def test_custom_model_without_file(self):
        """Test that custom model requires file."""
        with pytest.raises(ValueError):
            ReionizationConfig(model='custom')
    
    def test_invalid_model(self):
        """Test invalid model name."""
        with pytest.raises(ValueError):
            ReionizationConfig(model='invalid')


class TestBubbleConfig:
    """Tests for BubbleConfig."""
    
    def test_default_values(self):
        """Test default bubble configuration."""
        config = BubbleConfig()
        assert config.Rb == 5.0
        assert config.sigma_lnr == np.log(2)
        assert config.b == 6.0
    
    def test_validation(self):
        """Test parameter validation."""
        with pytest.raises(ValueError):
            BubbleConfig(Rb=-1)  # Negative Rb
        
        with pytest.raises(ValueError):
            BubbleConfig(R_int_min=100, R_int_max=10)  # min > max


class TestIntegrationConfig:
    """Tests for IntegrationConfig."""
    
    def test_default_values(self):
        """Test default integration configuration."""
        config = IntegrationConfig()
        assert config.z_min == 5.0
        assert config.z_max == 20.0
    
    def test_validation(self):
        """Test parameter validation."""
        with pytest.raises(ValueError):
            IntegrationConfig(z_min=15, z_max=10)  # min > max
        
        with pytest.raises(ValueError):
            IntegrationConfig(n_chi_samples=5)  # Too few samples
