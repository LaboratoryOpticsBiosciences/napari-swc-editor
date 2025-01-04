try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"
from ._reader import napari_get_reader, read_swc
from ._sample_data import make_empty_sample, make_sample_data
from ._widget import SWCEditorWidget
from ._writer import write_multiple, write_single_image

__all__ = (
    "napari_get_reader",
    "read_swc",
    "write_single_image",
    "write_multiple",
    "make_sample_data",
    "make_empty_sample",
    "SWCEditorWidget",
)
