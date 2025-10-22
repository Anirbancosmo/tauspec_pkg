"""TauSpec: Optical Depth Power Spectra Calculator.

A package for calculating optical depth power spectra during the Epoch of Reionization.
"""

__version__ = "0.2.0"
__author__ = "Anirban Roy"

from .calculator import TauSpecCalculator
from .config import CosmologyConfig, ReionizationConfig, BubbleConfig
from .cosmology import Cosmology
from .reionization import ReionizationHistory

__all__ = [
    "TauSpecCalculator",
    "CosmologyConfig",
    "ReionizationConfig",
    "BubbleConfig",
    "Cosmology",
    "ReionizationHistory",
]
