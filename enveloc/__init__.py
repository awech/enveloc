from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("enveloc")
except PackageNotFoundError:
    __version__ = "unknown"
