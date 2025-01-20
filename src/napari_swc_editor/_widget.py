"""
This module contains four napari widgets declared in
different ways:

- a pure Python function flagged with `autogenerate: true`
    in the plugin manifest. Type annotations are used by
    magicgui to generate widgets for each parameter. Best
    suited for simple processing tasks - usually taking
    in and/or returning a layer.
- a `magic_factory` decorated function. The `magic_factory`
    decorator allows us to customize aspects of the resulting
    GUI, including the widgets associated with each parameter.
    Best used when you have a very simple processing task,
    but want some control over the autogenerated widgets. If you
    find yourself needing to define lots of nested functions to achieve
    your functionality, maybe look at the `Container` widget!
- a `magicgui.widgets.Container` subclass. This provides lots
    of flexibility and customization options while still supporting
    `magicgui` widgets and convenience methods for creating widgets
    from type annotations. If you want to customize your widgets and
    connect callbacks, this is the best widget option for you.
- a `QWidget` subclass. This provides maximal flexibility but requires
    full specification of widget layouts, callbacks, events, etc.

References:
- Widget specification: https://napari.org/stable/plugins/guides.html?#widgets
- magicgui docs: https://pyapp-kit.github.io/magicgui/

Replace code below according to your needs.
"""

from typing import TYPE_CHECKING

from magicgui.widgets import Container, Table, create_widget

if TYPE_CHECKING:
    import napari


class SWCEditorWidget(Container):
    def __init__(self, viewer: "napari.viewer.Viewer"):
        super().__init__()
        self._viewer = viewer
        # use create_widget to generate widgets from type annotations
        self._point_layer_combo = create_widget(
            label="Points", annotation="napari.layers.Points"
        )
        self._link_previous_node_checkbox = create_widget(
            label="link previous node with new node (same as using CTRL+Click)",
            widget_type="CheckBox",
        )

        self._show_table_button = create_widget(
            label="Show swc table", annotation=bool, widget_type="PushButton"
        )

        self._get_layer_data()

        # connect your own callbacks
        self._point_layer_combo.changed.connect(self._get_layer_data)
        self._link_previous_node_checkbox.changed.connect(
            self._set_link_previous_node
        )

        self._show_table_button.changed.connect(self._set_table)

        # append into/extend the container with your widgets
        self.extend(
            [
                self._point_layer_combo,
                self._link_previous_node_checkbox,
                self._show_table_button,
            ]
        )

    def _set_table(self):
        if not hasattr(self, "_table") or self._table is None:
            table = Table(
                value=self._point_layer_combo.value.metadata[
                    "swc_data"
                ].reset_index()
            )
            table.native.setEditTriggers(table.native.NoEditTriggers)
            table.native.clicked.connect(self._table_clicked)
            table.native.doubleClicked.connect(self._table_double_clicked)

            self._point_layer_combo.value.selected_data.events._current.connect(
                self._change_current_cell
            )
            self._table = table
            self._viewer.window.add_dock_widget(
                self._table.native, area="right"
            )
        else:
            self._table.value = self._point_layer_combo.value.metadata[
                "swc_data"
            ].reset_index()

    def _change_current_cell(self, current_selected):
        if current_selected is not None:
            column = (
                self._table.native.currentColumn()
                if self._table.native.currentColumn() != -1
                else 0
            )
            self._table.native.setCurrentCell(current_selected, column)

    def _get_layer_data(self):
        layer = self._point_layer_combo.value
        if layer is None:
            self._link_previous_node_checkbox.value = False
            self._link_previous_node_checkbox.enabled = False
            self._show_table_button.enabled = False
        else:
            self._link_previous_node_checkbox.value = layer.metadata.get(
                "link_previous_node", False
            )
            self._link_previous_node_checkbox.enabled = True
            self._show_table_button.enabled = True
            layer.events.data.connect(self._set_table)

    def _set_link_previous_node(self, value):
        layer = self._point_layer_combo.value
        if layer is None:
            return
        layer.metadata["widget_link_activated"] = value

    def _table_clicked(self, event):
        row = event.row()
        x = self._table["x"][row]
        y = self._table["y"][row]
        z = self._table["z"][row]

        # get the current viewer
        viewer = self._viewer

        current_step = list(viewer.dims.current_step)
        current_step[-1] = row
        viewer.dims.current_step = current_step
        # use negative index to get the order of the dimensions even if len(ndims) > 3
        current_step[-3] = int(z)
        current_step[-2] = int(y)
        current_step[-1] = int(x)

        viewer.dims.point = current_step

        # set camera to the clicked point
        viewer.camera.center = [
            0,
            current_step[viewer.dims.order[1]],
            current_step[viewer.dims.order[2]],
        ]

        self._point_layer_combo.value.selected_data = [row]

    def _table_double_clicked(self, event):
        # same as clicking on the table but zoom
        self._table_clicked(event)
        self._viewer.camera.zoom = 80
