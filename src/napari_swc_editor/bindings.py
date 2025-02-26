import napari
import numpy as np
from napari.layers.base._base_constants import ActionType

from .swc_io import (
    add_edge,
    add_points,
    create_line_data_from_swc_data,
    get_branch_from_node,
    get_index_from_treenode_id,
    get_treenode_id_from_index,
    move_points,
    parse_data_from_swc_file,
    parse_swc_content,
    remove_edge,
    remove_points,
    sort_edge_indices,
    structure_id_to_symbol,
    symbol_to_structure_id,
    update_point_properties,
)


def add_napari_layers_from_swc_content(
    file_content: str, viewer: napari.Viewer
):
    """Create layers from a swc file

    Parameters
    ----------
    file_content : swc_content
        Content of the swc file
        Must have the following columns:
            - treenode_id
            - structure_id
            - x
            - y
            - z
            - r
            - parent_treenode_id
    viewer : napari.Viewer

    Returns
    -------
    layers : list of tuples
        List of layers to be added to the napari viewer
    """

    points, radius, lines, structure = parse_data_from_swc_file(file_content)

    structure_symbol = structure_id_to_symbol(structure)

    shape_layer = viewer.add_shapes(
        lines, shape_type="line", edge_width=radius, ndim=3
    )

    add_kwargs_points = {
        "n_dimensional": True,
        "size": radius,
        "blending": "additive",
        "symbol": structure_symbol,
        "face_color": "transparent",
        "metadata": {
            "raw_swc": file_content,
            "shape_layer": shape_layer,
            "Ctrl_activated": False,  # link input from keyboard
            "widget_link_activated": False,  # link input from widget
            # delay the update of the vector layer used for faster loading of the swc file
        },
    }

    point_layer = viewer.add_points(points, **add_kwargs_points)

    # only update the vector layer when the 3D viewer is active for faster loading
    viewer.dims.events.ndisplay.connect(
        lambda e: update_edges(point_layer) if e.value == 3 else None
    )
    # when the shape layer is visible, update the edges
    shape_layer.events.visible.connect(lambda _: update_edges(point_layer))

    bind_layers_with_events(point_layer, shape_layer)

    return [point_layer, shape_layer]


def update_edges(point_layer, lines=None, edge_width=None):
    """Update the edges of the shape layer based on the swc content

    Parameters
    ----------
    point_layer : napari.layers.Points
        Points layer
    """

    # if the shape layer is not visible, do not update the edges
    # this is for faster loading of the swc file
    if point_layer.metadata["shape_layer"].visible is False:
        return

    if lines is None or edge_width is None:
        raw_swc = point_layer.metadata["raw_swc"]
        df = parse_swc_content(raw_swc)

        lines, edge_width = create_line_data_from_swc_data(df)

    # when updating the shape layer directly, the previous data
    # is not removed correctly. So we remove it first
    point_layer.metadata["shape_layer"].data = []
    point_layer.metadata["shape_layer"].add_lines(lines, edge_width=edge_width)


def bind_layers_with_events(point_layer, shape_layer):
    """Bind the events of the point layer with the swc content

    Parameters
    ----------
    point_layer : napari.layers.Points
        Points layer
    shape_layer : napari.layers.Shapes
        Shapes layer
    """

    # check if "swc_data" is already in the metadata
    if "swc_data" not in point_layer.metadata:
        # if not, parse the swc content
        point_layer.metadata["swc_data"] = parse_swc_content(
            point_layer.metadata["raw_swc"]
        )

    point_layer.events.data.connect(event_add_points)
    point_layer.events.data.connect(event_move_points)
    point_layer.events.data.connect(event_remove_points)
    point_layer.events.size.connect(event_update_point_properties)
    point_layer.events.symbol.connect(event_update_point_properties)

    point_layer.bind_key("l")(event_add_edge)
    point_layer.bind_key("u")(event_remove_edge)
    point_layer.bind_key("b", select_branch)
    point_layer.bind_key("Shift+b", select_upstream_branch)

    point_layer.bind_key("Control", linked_point)

    point_layer.metadata["shape_layer"] = shape_layer

    # useful to ajust in 2D the third dimension
    point_layer.bind_key("up", on_arrow_up, overwrite=True)
    point_layer.bind_key("down", on_arrow_down, overwrite=True)


def on_arrow_up(layer):
    """Move the selected points up by 1 unit in the z direction"""
    selection_indices = list(layer.selected_data)
    layer.data[selection_indices] += np.array([1, 0, 0])
    layer.refresh()
    layer.events.data(
        value=layer.data,
        action=ActionType.CHANGED,
        data_indices=tuple(selection_indices),
        vertex_indices=((),),
    )


def on_arrow_down(layer):
    """Move the selected points down by 1 unit in the z direction"""
    selection_indices = list(layer.selected_data)
    layer.data[selection_indices] -= np.array([1, 0, 0])
    layer.refresh()
    layer.events.data(
        value=layer.data,
        action=ActionType.CHANGED,
        data_indices=tuple(selection_indices),
        vertex_indices=((),),
    )


def select_upstream_branch(layer):
    """Select the upstream branch of the selected point

    Parameters
    ----------
    layer : napari.layers.Points
        Points layer
    treenode_id : int
        Index of the selected point
    """

    select_branch(layer, upstream=True, downstream=False)


def select_branch(layer, upstream=True, downstream=True):
    """Select the branch of the selected point

    Parameters
    ----------
    layer : napari.layers.Points
        Points layer
    treenode_id : int
        Index of the selected point
    """

    raw_swc = layer.metadata["raw_swc"]
    df = parse_swc_content(raw_swc)

    selected_branches = []
    for selected in list(layer.selected_data):
        indices = get_treenode_id_from_index([selected], df)[0]
        branch_indices = get_branch_from_node(
            indices, df, upstream, downstream
        )
        selected_branches.extend(branch_indices)

    selected_branches = np.unique(selected_branches)
    indices = get_index_from_treenode_id(selected_branches, df)
    layer.selected_data = indices


def linked_point(layer):
    """Activate the Ctrl key for the next event.
    Used to link points together"""

    layer.metadata["Ctrl_activated"] = True
    yield
    layer.metadata["Ctrl_activated"] = False


def event_add_points(event):

    if event.action == "added":
        raw_swc = event.source.metadata["raw_swc"]

        df = parse_swc_content(raw_swc)

        # indices must be deduced because the event.data_indices has only (-1,) see
        # https://github.com/napari/napari/issues/7507
        indices = [
            -i - 1 for i in reversed(range(len(event.source.data) - len(df)))
        ]
        new_pos = event.source.data[indices]
        new_radius = event.source.size[indices]
        new_structure = event.source.symbol[indices]
        new_parents = -1

        # if shift is activated, the add the new edges from previous selected point
        if (
            event.source.metadata["Ctrl_activated"]
            or event.source.metadata["widget_link_activated"]
        ) and len(event.source.selected_data) > 0:

            previous_selected = list(event.source.selected_data)[-1]
            new_parents = get_treenode_id_from_index([previous_selected], df)[
                0
            ]

        structure_id = symbol_to_structure_id(
            [structure.value for structure in new_structure]
        )

        new_swc, df = add_points(
            raw_swc, new_pos, new_radius, structure_id, new_parents, df
        )

        event.source.metadata["raw_swc"] = new_swc
        event.source.metadata["swc_data"] = df

        if new_parents != -1:

            new_lines, new_r = create_line_data_from_swc_data(df)
            update_edges(event.source, new_lines, new_r)


def event_move_points(event):

    if event.action == "changed":
        raw_swc = event.source.metadata["raw_swc"]
        df = parse_swc_content(raw_swc)

        new_pos = event.source.data[list(event.data_indices)]
        indices = get_treenode_id_from_index(list(event.data_indices), df)

        new_swc, new_lines, df = move_points(raw_swc, indices, new_pos, df)

        event.source.metadata["raw_swc"] = new_swc
        update_edges(event.source, new_lines)
        event.source.metadata["swc_data"] = df


def event_remove_points(event):

    if event.action == "removed":
        raw_swc = event.source.metadata["raw_swc"]
        df = parse_swc_content(raw_swc)

        indices = get_treenode_id_from_index(list(event.data_indices), df)

        new_swc, new_lines, new_r, df = remove_points(raw_swc, indices, df)

        event.source.metadata["raw_swc"] = new_swc
        update_edges(event.source, new_lines, new_r)
        event.source.metadata["swc_data"] = df


def event_update_point_properties(event):
    """Update the properties (`size` -> `radius` and `symbol` -> `structure_id`)
    of a point to the swc content point

    Parameters
    ----------
    layer : napari.layers.Points
        Points layer
    """

    raw_swc = event.source.metadata["raw_swc"]
    df = parse_swc_content(raw_swc)

    indices = get_treenode_id_from_index(list(event.source.selected_data), df)

    structure_object = event.source.symbol[list(event.source.selected_data)]
    structure_id = symbol_to_structure_id(
        [structure.value for structure in structure_object]
    )

    properties = {
        "r": event.source.size[list(event.source.selected_data)],
        "structure_id": structure_id,
    }

    new_swc, new_lines, new_r, df = update_point_properties(
        raw_swc,
        indices,
        properties,
        df,
    )

    event.source.metadata["raw_swc"] = new_swc
    update_edges(event.source, new_lines, new_r)
    event.source.metadata["swc_data"] = df


def event_add_edge(
    layer, treenode_id=None, parent_treenode_id=None, sort=True
):
    """Add an edge between two selected points

    Parameters
    ----------
    layer : napari.layers.Points
        Points layer
    sort : bool, optional
        If True, the indices will be sorted so soma are linked to other points,
        by default True
    """

    raw_swc = layer.metadata["raw_swc"]
    df = parse_swc_content(raw_swc)

    if treenode_id is None:
        treenode_id = get_treenode_id_from_index(list(layer.selected_data), df)

    if parent_treenode_id is None:
        # if no parent is given, the parent will be the following selected points
        if sort:
            treenode_id = sort_edge_indices(raw_swc, treenode_id, df)
        parent_treenode_id = treenode_id[:-1].astype(int)
        treenode_id = treenode_id[1:].astype(int)

    new_swc, new_lines, new_r, df = add_edge(
        raw_swc, treenode_id, parent_treenode_id, df
    )

    layer.metadata["raw_swc"] = new_swc
    update_edges(layer, new_lines, new_r)
    layer.metadata["swc_data"] = df


def event_remove_edge(layer, sort=True):
    """Add an edge between two selected points

    Parameters
    ----------
    layer : napari.layers.Points
        Points layer
    sort : bool, optional
        If True, the indices will be sorted so soma are linked to other points,
        by default True
    """

    raw_swc = layer.metadata["raw_swc"]
    df = parse_swc_content(raw_swc)

    indices = get_treenode_id_from_index(list(layer.selected_data), df)

    if sort:
        indices = sort_edge_indices(raw_swc, indices, df)
    new_swc, new_lines, new_r, df = remove_edge(raw_swc, indices, df)

    layer.metadata["raw_swc"] = new_swc
    update_edges(layer, new_lines, new_r)
    layer.metadata["swc_data"] = df
