"""
This module is an example of a barebones numpy reader plugin for napari.

It implements the Reader specification, but your plugin may choose to
implement multiple readers or even other plugin contributions. see:
https://napari.org/stable/plugins/guides.html?#readers
"""

import napari

from .bindings import bind_layers_with_events
from .swc_io import parse_data_from_swc_file


def napari_get_reader(path):
    """A basic implementation of a Reader contribution.

    Parameters
    ----------
    path : str or list of str
        Path to file, or list of paths.

    Returns
    -------
    function or None
        If the path is a recognized format, return a function that accepts the
        same path or list of paths, and returns a list of layer data tuples.
    """
    if isinstance(path, list):
        # reader plugins may be handed single path, or a list of paths.
        # if it is a list, it is assumed to be an image stack...
        # so we are only going to look at the first file.
        path = path[0]

    # if we know we cannot read the file, we immediately return None.
    if not path.endswith(".swc"):
        return None

    # otherwise we return the *function* that can read ``path``.
    return reader_function


def reader_function(path):
    """Take a path or list of paths and return a list of LayerData tuples.
    The first layer will be a points layer, the second layer will be a shapes layer.

    Parameters
    ----------
    path : str or list of str
        Path to file, or list of paths.

    Returns
    -------
    layer_data : list of tuples
        A list of LayerData tuples. Each path will return two tuples, one for points and one for shapes.
    """
    # handle both a string and a list of strings
    paths = [path] if isinstance(path, str) else path

    viewer = napari.current_viewer()

    if viewer is None:
        viewer = napari.Viewer()

    for _path in paths:
        with open(_path) as f:
            file_content = f.read()

        nodes, radius, lines = parse_data_from_swc_file(file_content)

        shape_layer = viewer.add_shapes(
            lines, shape_type="line", edge_width=radius
        )

        add_kwargs_points = {
            "n_dimensional": True,
            "size": radius,
            "metadata": {
                "raw_swc": file_content,
                "shape_layer": shape_layer,
            },
        }

        point_layer = viewer.add_points(nodes, **add_kwargs_points)

        bind_layers_with_events(point_layer, shape_layer)

    return [point_layer, shape_layer]