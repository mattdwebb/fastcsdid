"""Group-time difference-in-differences for stationary repeated cross sections."""
from .core import fit, IdentificationError, influence_function, rif

__all__ = ["fit", "IdentificationError", "influence_function", "rif"]
