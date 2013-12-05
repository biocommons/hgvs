import pkg_resources
import warnings

try:
    __version__ = pkg_resources.get_distribution(__package__).version
except pkg_resources.DistributionNotFound as e:
    warnings.warn("can't get __version__ because %s package isn't installed" % __package__,Warning)
    __version__ = None
