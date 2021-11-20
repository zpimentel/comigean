from pkg_resources import get_distribution, DistributionNotFound

__version__ = None

try:
    __version__ = get_distribution('comigean').version
except DistributionNotFound:
    pass
