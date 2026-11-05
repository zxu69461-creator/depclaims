"""depclaims: score published cancer-cell growth claims against DepMap screens of the same cell line."""
from .core import Atlas, Score, build, scale, THRESHOLD

__all__ = ['Atlas', 'Score', 'build', 'scale', 'THRESHOLD']
__version__ = '0.1.0'
