import napari

from .swc_io import (
    add_edge,
    add_points,
    create_line_data_from_swc_data,
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
        lines, shape_type="line", edge_width=radius
    )

    add_kwargs_points = {
        "n_dimensional": True,
        "size": radius,
        "blending": "additive",
        "symbol": structure_symbol,
        "metadata": {
            "raw_swc": file_content,
            "shape_layer": shape_layer,
            "Ctrl_activated": False,
        },
    }

    point_layer = viewer.add_points(points, **add_kwargs_points)

    bind_layers_with_events(point_layer, shape_layer)

    return [point_layer, shape_layer]


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
    point_layer.bind_key("Shift-l")(event_add_edge_wo_sort)
    point_layer.bind_key("u")(event_remove_edge)

    point_layer.metadata["shape_layer"] = shape_layer


def event_add_points(event):

    if event.action == "added":
        raw_swc = event.source.metadata["raw_swc"]

        df = parse_swc_content(raw_swc)
        new_pos = event.source.data[list(event.data_indices)]
        new_radius = event.source.size[list(event.data_indices)]
        new_structure = event.source.symbol[list(event.data_indices)]
        new_parents = -1

        # if shift is activated, the add the new edges from previous selected point
        if (
            event.source.metadata["Ctrl_activated"]
            and len(event.source.selected_data) > 0
        ):

            previous_selected = list(event.source.selected_data)[-1]
            new_parents = get_treenode_id_from_index([previous_selected], df)[
                0
            ]

        new_swc, df = add_points(
            raw_swc, new_pos, new_radius, new_structure, new_parents, df
        )

        event.source.metadata["raw_swc"] = new_swc
        event.source.metadata["swc_data"] = df

        if new_parents != -1:
            new_lines, new_r = create_line_data_from_swc_data(df)
            event.source.metadata["shape_layer"].data = []
            event.source.metadata["shape_layer"].add_lines(
                new_lines, edge_width=new_r
            )


def event_move_points(event):

    if event.action == "changed":
        raw_swc = event.source.metadata["raw_swc"]
        df = parse_swc_content(raw_swc)

        new_pos = event.source.data[list(event.data_indices)]
        indices = get_treenode_id_from_index(list(event.data_indices), df)

        new_swc, new_lines, df = move_points(raw_swc, indices, new_pos, df)

        event.source.metadata["raw_swc"] = new_swc
        event.source.metadata["shape_layer"].data = new_lines
        # event.source.metadata["shape_layer"].edge_width = df
        event.source.metadata["swc_data"] = df


def event_remove_points(event):

    if event.action == "removed":
        raw_swc = event.source.metadata["raw_swc"]
        df = parse_swc_content(raw_swc)

        indices = get_treenode_id_from_index(list(event.data_indices), df)

        new_swc, new_lines, new_r, df = remove_points(raw_swc, indices, df)

        event.source.metadata["raw_swc"] = new_swc
        event.source.metadata["shape_layer"].data = []
        event.source.metadata["shape_layer"].add_lines(
            new_lines, edge_width=new_r
        )
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
    # when updating the shape layer directly, the previous data
    # is not removed correctly. So we remove it first
    event.source.metadata["shape_layer"].data = []
    event.source.metadata["shape_layer"].add_lines(new_lines, edge_width=new_r)
    event.source.metadata["swc_data"] = df


def event_add_edge_wo_sort(layer):
    """Add an edge between two selected points without sorting the indices

    Parameters
    ----------
    layer : napari.layers.Points
        Points layer
    """
    ### This function is not used in the current version of the plugin
    ### Because napari does not support history of the selected points
    event_add_edge(layer, sort=False)


def event_add_edge(layer, indices=None, sort=True):
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

    if indices is None:
        indices = get_treenode_id_from_index(list(layer.selected_data), df)

    if sort:
        indices = sort_edge_indices(raw_swc, indices, df)
    new_swc, new_lines, new_r, df = add_edge(raw_swc, indices, df)

    layer.metadata["raw_swc"] = new_swc
    # when updating the shape layer directly, the previous data
    # is not removed correctly. So we remove it first
    layer.metadata["shape_layer"].data = []
    layer.metadata["shape_layer"].add_lines(new_lines, edge_width=new_r)
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
    # when updating the shape layer directly, the previous data
    # is not removed correctly. So we remove it first
    layer.metadata["shape_layer"].data = []
    layer.metadata["shape_layer"].add_lines(new_lines, edge_width=new_r)
    layer.metadata["swc_data"] = df
