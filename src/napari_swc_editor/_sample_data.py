"""
This module is an example of a barebones sample data provider for napari.

It implements the "sample data" specification.
see: https://napari.org/stable/plugins/guides.html?#sample-data

Replace code below according to your needs.
"""

from __future__ import annotations

import tempfile
from urllib.request import urlretrieve

from ._reader import reader_function


def make_sample_data():
    # Download data from neuromorpho.org
    url = "https://raw.githubusercontent.com/LaboratoryOpticsBiosciences/napari-swc-reader/refs/heads/main/data/204-2-6nj.CNG.swc"
    urlretrieve(url, "204-2-6nj.CNG.swc")
    result = reader_function("204-2-6nj.CNG.swc")
    return result


def make_empty_sample():
    content = "# SWC created using napari-swc-editor\n"
    # create temp file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".swc", delete=False
    ) as f:
        f.write(content)
        result = reader_function(f.name)
        return result
