from importlib.metadata import version as _version

__all__ = ["__version__"]
__version__ = _version("flowsig")

from . import preprocessing as pp
from . import tools as tl
from . import plotting as pl
