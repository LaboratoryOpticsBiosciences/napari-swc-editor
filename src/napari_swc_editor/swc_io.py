import io

import numpy as np
import pandas as pd

SWC_SYMBOL = {
    0: "clobber",  # undefined
    1: "star",  # soma
    2: "disc",  # axon
    3: "triangle_down",  # basal dendrite
    4: "triangle_up",  # apical dendrite
}


def parse_swc_content(file_content):
    """Parse a swc file and return a dataframe with the data.
    Must have the following columns:
    - treenode_id
    - structure_id
    - x
    - y
    - z
    - r
    - parent_treenode_id

    Parameters
    ----------
    file_content : swc_content
        Content of the swc file

    Returns
    -------
    df : pd.DataFrame
        Dataframe with the data extracted from the swc file
    """

    df = pd.read_csv(
        io.StringIO(file_content),
        sep=r"\s+",  # separator is any whitespace
        comment="#",
        # set columns names according to SWC format
        # http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html
        names=[
            "treenode_id",
            "structure_id",
            "x",
            "y",
            "z",
            "r",
            "parent_treenode_id",
        ],
        index_col=0,
    )

    return df


def structure_id_to_symbol(structure_ids, swc_structure=SWC_SYMBOL):
    """Convert list structure id to a names

    Parameters
    ----------
    structure_id : list
        List of structure ids
    swc_structure : dict
        Dictionary with the structure id as key and the name as value
        Default is SWC_STRUCTURE

    Returns
    -------
    name : str
        Name of the structure
    """

    names = [
        swc_structure.get(structure_id, "x") for structure_id in structure_ids
    ]

    return names


def symbol_to_structure_id(symbol, swc_structure=SWC_SYMBOL):
    """Convert a symbol to a structure id

    Parameters
    ----------
    symbol : str
        Symbol of the structure
    swc_structure : dict
        Dictionary with the structure id as key and the name as value
        Default is SWC_STRUCTURE

    Returns
    -------
    structure_id : int
        Id of the structure
    """

    # invert the dictionary
    swc_structure = {v: k for k, v in swc_structure.items()}
    structure_id = [swc_structure.get(s, 0) for s in symbol]

    return structure_id


def parse_data_from_swc_file(file_content):
    """Create layers from a swc file

    Parameters
    ----------
    file_content : swc_content
        Content of the swc file

    Returns
    -------
    points : np.ndarray
        All positions of the points
    radius : np.ndarray
        Radius of the points
    lines : np.ndarray
        All lines connecting the points
    structure : np.ndarray
        Structure of the points
    """

    df = parse_swc_content(file_content)

    points, radius, structure = create_point_data_from_swc_data(df)
    lines, _ = create_line_data_from_swc_data(df)

    return points, radius, lines, structure


def create_point_data_from_swc_data(df):
    """Take a dataframe extracted from a swc and create point data

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe extracted from a swc file. Should have the following columns:
        - x: x coordinate of the node
        - y: y coordinate of the node
        - z: z coordinate of the node
        - structure_id: id of the structure
        - r: radius of the node

    Returns
    -------
    points : np.ndarray
        All positions of the points. The coordinates are in the napari order (z, y, x)
    radius : np.ndarray
        Radius of the points
    structure : np.ndarray
        Structure of the points
    """

    radius = df["r"].values

    structure = df["structure_id"].values

    # for each node create a point
    points = df[["z", "y", "x"]].values

    return points, radius, structure


def create_line_data_from_swc_data(df):
    """Take a dataframe extracted from a swc and create line data

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe extracted from a swc file. Should have the following columns:
        - x: x coordinate of the node
        - y: y coordinate of the node
        - z: z coordinate of the node
        - r: radius of the node
        - parent_treenode_id: id of the parent node

    Returns
    -------
    lines : np.ndarray
        All lines connecting the points
    radius : np.ndarray
        Radius of the lines
    """

    # for each nodes create a point
    points = df[["z", "y", "x"]].values

    # for each edge create a line
    edges = df["parent_treenode_id"].values

    # remove all soma nodes
    points = points[edges != -1]
    edges = edges[edges != -1]

    # for each id in edges, get the corresponding node according to its index
    prev_point = df.loc[edges, ["z", "y", "x"]].values

    lines = np.array([points, prev_point])
    lines = np.moveaxis(lines, 0, 1)

    radius = df.loc[edges, "r"].values

    return lines, radius


def write_swc_content(df, swc_content=None):
    """Write a dataframe to a swc file content

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with the data to be written to the swc file
        Should contain:
            - treenode_id: id of the node
            - structure_id: id of the structure
            - x: x coordinate of the node
            - y: y coordinate of the node
            - z: z coordinate of the node
            - r: radius of the node
            - parent_treenode_id: id of the parent node
    swc_content : str
        Original content of the swc file. If provided, the header lines will be kept.

    Returns
    -------
    new_swc_content : str
        New content of the swc file
    """
    # get header lines starting with #
    header_lines = [
        line for line in swc_content.split("\n") if line.startswith("#")
    ]
    # create a new swc content
    new_swc_content = "\n".join(header_lines) + "\n"

    if df.size > 0:
        result_string = "\n".join(
            df.reset_index().astype(str).apply(" ".join, axis=1)
        )
        new_swc_content += result_string
    return new_swc_content


def add_points(
    swc_content,
    new_positions,
    new_radius,
    structure_id=0,
    parent_treenode_id=-1,
    swc_df=None,
):
    """Add a point to the swc content

    Parameters
    ----------
    swc_content : swc_content
        Content of the swc file
    new_positions : np.ndarray
        New positions to be added in napari order (z, y, x)
    new_radius : np.ndarray
        Radius of the new positions
    structure_id : np.ndarray
        Structure of the new positions, see SWC_SYMBOL
        Default is 0 (undefined)
    parent_treenode_id : np.ndarray
        Parent of the new positions
        Default is -1 (no parent)
    swc_df : pd.DataFrame
        Dataframe extracted from the swc file. Should have the following columns:
        - x: x coordinate of the node
        - y: y coordinate of the node
        - z: z coordinate of the node
        - r: radius of the node
        - parent_treenode_id: id of the parent node

    Returns
    -------
    new_swc_content : swc_content
        New content of the swc file
    """

    if new_positions.ndim == 1:
        new_positions = new_positions[np.newaxis]

    if swc_df is None:
        swc_df = parse_swc_content(swc_content)

    # change napari position order to swc order
    new_points = pd.DataFrame(new_positions, columns=["z", "y", "x"])
    new_points["r"] = new_radius
    new_points["structure_id"] = structure_id
    new_points["parent_treenode_id"] = parent_treenode_id

    # order columns to respect swc format
    new_points = new_points[
        ["structure_id", "x", "y", "z", "r", "parent_treenode_id"]
    ]

    if swc_df.size > 0:
        previous_max = swc_df.index.max()
        max_index = np.array(
            previous_max + np.arange(1, len(new_points) + 1)
        ).astype(int)
        new_points.index = max_index

        new_df = pd.concat([swc_df, new_points])
    else:
        new_df = new_points.copy()

    new_df.index.name = "treenode_id"
    new_swc_content = write_swc_content(new_df, swc_content)
    return new_swc_content, new_df


def move_points(swc_content, index, new_positions, swc_df=None):
    """Move a point in the swc content

    Parameters
    ----------
    swc_content : swc_content
        Content of the swc file
    index : int
        Index of the point to be moved
    new_positions : np.ndarray
        New positions of the point in napari order (z, y, x)
    swc_df : pd.DataFrame
        Dataframe extracted from the swc file. Should have the following columns:
        - x: x coordinate of the node
        - y: y coordinate of the node
        - z: z coordinate of the node
        - r: radius of the node
        - parent_treenode_id: id of the parent node

    Returns
    -------
    new_swc_content : swc_content
        New content of the swc file
    moved_lines : np.ndarray
        New lines connecting the points
    swc_df : pd.DataFrame
        Dataframe extracted from the new swc file
    """

    if swc_df is None:
        swc_df = parse_swc_content(swc_content)

    swc_df.loc[index, ["z", "y", "x"]] = new_positions

    new_swc_content = write_swc_content(swc_df, swc_content)
    moved_lines, _ = create_line_data_from_swc_data(swc_df)

    return new_swc_content, moved_lines, swc_df


def remove_points(swc_content, indices, swc_df=None):
    """Delete points in the swc content

    Parameters
    ----------
    swc_content : swc_content
        Content of the swc file
    indices : list of int
        Indices of the points to be deleted
    swc_df : pd.DataFrame
        Dataframe extracted from the swc file. Should have the following columns:
        - x: x coordinate of the node
        - y: y coordinate of the node
        - z: z coordinate of the node
        - r: radius of the node
        - parent_treenode_id: id of the parent node

    Returns
    -------
    new_swc_content : swc_content
        New content of the sw
    moved_lines : np.ndarray
        New lines connecting the points
    new_r : np.ndarray
        New radius of the lines
    swc_df : pd.DataFrame
        Dataframe extracted from the swc file
    """

    if swc_df is None:
        swc_df = parse_swc_content(swc_content)

    swc_df = swc_df.drop(indices)

    mask = swc_df["parent_treenode_id"].isin(indices)
    swc_df.loc[mask, "parent_treenode_id"] = -1

    new_swc_content = write_swc_content(swc_df, swc_content)

    moved_lines, new_r = create_line_data_from_swc_data(swc_df)

    return new_swc_content, moved_lines, new_r, swc_df


def get_treenode_id_from_index(iloc, df):
    """Get the treenode_id from the iloc index

    Parameters
    ----------
    iloc : int or list of int
        Index of the row in the dataframe
    df : pd.DataFrame
        Dataframe extracted from a swc file. Should have the following columns:
        - treenode_id as index
        - parent_treenode_id: id of the parent node

    Returns
    -------
    indices : np.ndarray
        Treenode_id of the selected index
    """

    if isinstance(iloc, int):
        iloc = [iloc]

    indices = df.iloc[iloc].index.values

    return indices


def add_edge(swc_content, indices, swc_df=None):
    """Add an edge between two or more indices in order

    Parameters
    ----------
    swc_content : str
        Content of the swc file
    indices : list of int
        Indices of the points to be connected at least two indices are needed
    swc_df : pd.DataFrame, optional
        Dataframe extracted from a swc file. Should have the following columns:
        - treenode_id as index
        - parent_treenode_id: id of the parent node

    Returns
    -------
    new_swc_content : str
        New content of the swc file
    new_lines : np.ndarray
        New lines connecting the nodes
    new_r : np.ndarray
        New radius of the lines
    swc_df : pd.DataFrame
        Dataframe extracted from the swc file
    """

    assert len(indices) >= 2, "At least two indices are needed to create edges"

    if swc_df is None:
        swc_df = parse_swc_content(swc_content)

    for i in range(1, len(indices)):
        swc_df.loc[indices[i], "parent_treenode_id"] = indices[i - 1]

    new_lines, new_r = create_line_data_from_swc_data(swc_df)

    new_swc_content = write_swc_content(swc_df, swc_content)

    return new_swc_content, new_lines, new_r, swc_df


def remove_edge(swc_content, indices, swc_df=None):
    """Remove an edge between two or more indices in order

    Parameters
    ----------
    swc_content : str
        Content of the swc file
    indices : list of int
        Indices of the points with edges to be removed at least one indices
        are needed
    swc_df : pd.DataFrame, optional
        Dataframe extracted from a swc file. Should have the following columns:
        - treenode_id as index
        - parent_treenode_id: id of the parent node

    Returns
    -------
    new_swc_content : str
        New content of the swc file
    new_lines : np.ndarray
        New lines connecting the points
    new_r : np.ndarray
        New radius of the lines
    swc_df : pd.DataFrame
        Dataframe extracted from the swc file
    """

    assert len(indices) >= 1, "At least one indices are needed to remove edges"

    if swc_df is None:
        swc_df = parse_swc_content(swc_content)

    for i in range(1, len(indices)):
        swc_df.loc[indices[i], "parent_treenode_id"] = -1

    new_lines, new_r = create_line_data_from_swc_data(swc_df)

    new_swc_content = write_swc_content(swc_df, swc_content)

    return new_swc_content, new_lines, new_r, swc_df


def sort_edge_indices(swc_content, indices, swc_df=None):
    """Sort the indices of the edges
    With only two indices:
        - if one is the parent of the other, the parent should be the first
        - if one is the soma, the soma should be the last
    else:
        keep as it is

    Parameters
    ----------
    swc_content : str
        Content of the swc file
    indices : list of int
        Indices of the points with edges to be removed at least one indices
        are needed
    swc_df : pd.DataFrame, optional
        Dataframe extracted from a swc file. Should have the following columns:
        - treenode_id as index
        - parent_treenode_id: id of the parent node

    Returns
    -------
    sorted_indices : np.ndarray
        Sorted indices
    """

    assert len(indices) >= 2, "At least two indices are needed to remove edges"

    if swc_df is None:
        swc_df = parse_swc_content(swc_content)

    # check if one indices has a parent_treenode_id in the list
    parent_ids = swc_df.loc[indices, "parent_treenode_id"].values

    soma_id = indices[parent_ids == -1]
    non_soma_id = indices[parent_ids != -1]

    new_indices = []

    # ideal case when only two indices with one is a soma, then the soma should be the last
    if len(indices) == 2 and len(soma_id) == 1:
        first_node = non_soma_id[0]
        new_indices.append(first_node)

        second_node = soma_id[0]
        new_indices.append(second_node)

        return np.array(new_indices)

    # ideal case when only two indices with one is the parent of the other
    if len(indices) == 2 and (
        indices[0] in parent_ids or indices[1] in parent_ids
    ):
        if indices[0] in parent_ids:
            first_node = indices[0]
            second_node = indices[1]
        else:
            first_node = indices[1]
            second_node = indices[0]

        new_indices.append(first_node)
        new_indices.append(second_node)

        return np.array(new_indices)

    # else, we need to sort the indices
    sorted_indices = np.sort(indices)

    return sorted_indices


def update_point_properties(swc_content, indices, new_properties, swc_df=None):
    """Update the properties of the nodes

    Parameters
    ----------
    swc_content : str
        Content of the swc file
    indices : list of int
        Indices of the points to be updated
    new_properties : dict
        Properties to be updated. Such as: `r` and `structure_id`
    swc_df : pd.DataFrame, optional
        Dataframe extracted from a swc file. Should have the following columns:
        - treenode_id as index
        - parent_treenode_id: id of the parent node

    Returns
    -------
    new_swc_content : str
        New content of the swc file
    swc_df : pd.DataFrame
        Dataframe extracted from the swc file
    """

    if swc_df is None:
        swc_df = parse_swc_content(swc_content)

    for key, value in new_properties.items():
        swc_df.loc[indices, key] = value

    new_lines, new_r = create_line_data_from_swc_data(swc_df)

    new_swc_content = write_swc_content(swc_df, swc_content)

    return new_swc_content, new_lines, new_r, swc_df
