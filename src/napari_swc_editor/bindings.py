from .swc_io import (
    add_edge,
    add_points,
    get_treenode_id_from_index,
    move_points,
    parse_swc_content,
    remove_edge,
    remove_points,
    sort_edge_indices,
    symbol_to_structure_id,
    update_node_properties,
)


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
    point_layer.events.size.connect(event_update_node_properties)
    point_layer.events.symbol.connect(event_update_node_properties)

    point_layer.bind_key("l")(event_add_edge)
    point_layer.bind_key("Shift-l")(event_add_edge_wo_sort)
    # point_layer.bind_key('u')(event_add_edge)

    point_layer.metadata["shape_layer"] = shape_layer


def event_add_points(event):

    if event.action == "added":
        raw_swc = event.source.metadata["raw_swc"]

        df = parse_swc_content(raw_swc)
        new_pos = event.source.data[list(event.data_indices)]
        new_radius = event.source.size[list(event.data_indices)]
        new_structure = event.source.symbol[list(event.data_indices)]

        new_swc = add_points(raw_swc, new_pos, new_radius, new_structure, df)

        event.source.metadata["raw_swc"] = new_swc
        event.source.metadata["swc_data"] = df


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


def event_update_node_properties(event):
    """Update the properties (`size` -> `radius` and `symbol` -> `structure_id`)
    of a point to the swc content node

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

    new_swc, new_lines, new_r, df = update_node_properties(
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
    event_add_edge(layer, sort=False)


def event_add_edge(layer, sort=True):
    """Add an edge between two selected points

    Parameters
    ----------
    layer : napari.layers.Points
        Points layer
    sort : bool, optional
        If True, the indices will be sorted so soma are linked to other nodes,
        by default True
    """

    raw_swc = layer.metadata["raw_swc"]
    df = parse_swc_content(raw_swc)

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
        If True, the indices will be sorted so soma are linked to other nodes,
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
