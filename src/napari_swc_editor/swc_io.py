import io

import numpy as np
import pandas as pd


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


def parse_data_from_swc_file(file_content):
    """Create layers from a swc file

    Parameters
    ----------
    file_content : swc_content
        Content of the swc file

    Returns
    -------
    nodes : np.ndarray
        All positions of the nodes
    radius : np.ndarray
        Radius of the nodes
    lines : np.ndarray
        All lines connecting the nodes
    """

    df = parse_swc_content(file_content)

    nodes, radius = create_point_data_from_swc_data(df)
    lines, _ = create_line_data_from_swc_data(df)

    return nodes, radius, lines


def create_point_data_from_swc_data(df):
    """Take a dataframe extracted from a swc and create point data

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe extracted from a swc file. Should have the following columns:
        - x: x coordinate of the node
        - y: y coordinate of the node
        - z: z coordinate of the node
        - r: radius of the node

    Returns
    -------
    nodes : np.ndarray
        All positions of the nodes
    radius : np.ndarray
        Radius of the nodes
    """

    radius = df["r"].values

    # for each node create a point
    nodes = df[["x", "y", "z"]].values

    return nodes, radius


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
        All lines connecting the nodes
    radius : np.ndarray
        Radius of the lines
    """

    # for each node create a point
    nodes = df[["x", "y", "z"]].values

    # for each edge create a line
    edges = df["parent_treenode_id"].values

    # remove all soma nodes
    nodes = nodes[edges != -1]
    edges = edges[edges != -1]

    # for each id in edges, get the corresponding node according to its index
    prev_node = df.loc[edges, ["x", "y", "z"]].values

    lines = np.array([nodes, prev_node])
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


def add_points(swc_content, new_positions, new_radius, swc_df=None):
    """Add a point to the swc content

    Parameters
    ----------
    swc_content : swc_content
        Content of the swc file
    new_positions : np.ndarray
        New positions to be added
    new_radius : np.ndarray
        Radius of the new positions
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

    new_nodes = pd.DataFrame(new_positions, columns=["x", "y", "z"])
    new_nodes["r"] = new_radius
    new_nodes["structure_id"] = 0  # 0 is undefined structure
    new_nodes["parent_treenode_id"] = -1

    if swc_df.size > 0:
        previous_max = swc_df.index.max()
        max_index = np.array(
            previous_max + np.arange(1, len(new_nodes) + 1)
        ).astype(int)
        new_nodes.index = max_index

        new_df = pd.concat([swc_df, new_nodes])
    else:
        new_df = new_nodes

    new_df.index.name = "treenode_id"
    new_swc_content = write_swc_content(new_df, swc_content)
    return new_swc_content


def move_points(swc_content, index, new_positions, swc_df=None):
    """Move a point in the swc content

    Parameters
    ----------
    swc_content : swc_content
        Content of the swc file
    index : int
        Index of the point to be moved
    new_positions : np.ndarray
        New positions of the point
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
        New lines connecting the nodes
    swc_df : pd.DataFrame
        Dataframe extracted from the new swc file
    """

    if swc_df is None:
        swc_df = parse_swc_content(swc_content)

    swc_df.loc[index, ["x", "y", "z"]] = new_positions

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
        New lines connecting the nodes
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
        New lines connecting the nodes
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

    return new_swc_content, new_lines, new_r


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