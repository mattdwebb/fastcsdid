"""Native CSDID reference implementation. No R or Stata runtime dependency."""
from .core import fit, IdentificationError, influence_function, rif

__all__ = ["fit", "IdentificationError", "influence_function", "rif"]
