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

from magicgui.widgets import Container, create_widget

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
            annotation=float,
            widget_type="CheckBox",
        )

        self._get_layer_data()

        # connect your own callbacks
        self._link_previous_node_checkbox.changed.connect(
            self._set_link_previous_node
        )
        self._link_previous_node_checkbox.changed.connect(self._get_layer_data)

        # append into/extend the container with your widgets
        self.extend(
            [
                self._point_layer_combo,
                self._link_previous_node_checkbox,
            ]
        )

    def _get_layer_data(self):
        layer = self._point_layer_combo.value
        if layer is None:
            self._link_previous_node_checkbox.value = False
            self._link_previous_node_checkbox.enabled = False
        else:
            self._link_previous_node_checkbox.value = layer.metadata.get(
                "link_previous_node", False
            )
            self._link_previous_node_checkbox.enabled = True

    def _set_link_previous_node(self, value):
        layer = self._point_layer_combo.value
        if layer is None:
            return
        layer.metadata["link_previous_node"] = value
