#
# pypeid_modes.py -- keyboard modes for manipulating PypeIt plots in Ginga
#
#
from ginga.gw.PlotView import PlotViewBase
from ginga.misc import Bunch
from ginga.modes.mode_base import Mode
from ginga.modes.plot2d import Plot2DMode


class Spec1DMode(Plot2DMode):
    """
    Spec1DMode enables special bindings for viewing PypeIt spec1d files.

    Enter the mode by
    -----------------
    * Space, then "1"

    Exit the mode by
    ----------------
    * Esc

    Default bindings in mode
    ------------------------
    * Shift + left click : set pan position
    * middle click : set pan position
    * scroll : zoom in/out
    * ctrl + scroll : zoom in/out X axis only
    * shift + scroll : zoom in/out Y axis only
    * alt + scroll : zoom in/out Y axis only
    * meta + scroll : zoom in/out at cursor
    * equals : zoom in one zoom level
    * ctrl + equals : zoom in X axis one zoom level
    * shift + equals (plus) : zoom in Y axis one zoom level
    * minus : zoom out one zoom level
    * ctrl + minus : zoom out X axis one zoom level
    * shift + minus (underscore): zoom out Y axis one zoom level
    * 9 : zoom out maintaining cursor position
    * ctrl + 9 : zoom out X axis maintaining cursor position
    * shift + 9 (left paren): zoom out Y axis maintaining cursor position
    * 0 : zoom in maintaining cursor position
    * ctrl + 0 : zoom in X axis maintaining cursor position
    * shift + 0 (right paren): zoom in Y axis maintaining cursor position
    * backquote : fit plot to window
    * 1 : fit x axis only to window
    * 2 : fit y axis only to window
    * left arrow : pan left
    * right arrow : pan right
    * up arrow : pan up
    * down arrow : pan down
    * k : set lower X range to X value at cursor
    * l : set upper X range to X value at cursor
    * K : set lower Y range to Y value at cursor
    * L : set upper Y range to Y value at cursor
    """

    # Needs to be set by reference viewer (via set_shell_ref) before any
    # channel viewers are created
    fv = None

    @classmethod
    def set_shell_ref(cls, fv):
        cls.fv = fv

    @classmethod
    def is_compatible_viewer(cls, viewer):
        return isinstance(viewer, PlotViewBase)

    def __init__(self, viewer, settings=None):
        super().__init__(viewer, settings=settings)

        self.actions = dict(
            dmod_spec1d=['__1', None, 'spec1d'],

            ms_showxy=['spec1d+nobtn'],
            ms_panset2d=['spec1d+middle', 'spec1d+shift+left'],

            sc_zoom2d=['spec1d+*+scroll'],

            kp_pan_left=['spec1d+*+left'],
            kp_pan_right=['spec1d+*+right'],
            kp_pan_up=['spec1d+*+up'],
            kp_pan_down=['spec1d+*+down'],

            kp_zoom_in=['spec1d++', 'spec1d+='],
            kp_zoom_out=['spec1d+-', 'spec1d+_'],
            kp_zoom_cursor_in=['spec1d+*+0', 'spec1d+*+)'],
            kp_zoom_cursor_out=['spec1d+*+9', 'spec1d+*+('],

            kp_zoom_fit=['spec1d+backquote'],
            kp_zoom_fit_x=['spec1d+1'],
            kp_zoom_fit_y=['spec1d+2'],

            kp_cut_x_lo=['spec1d+k'],
            kp_cut_x_hi=['spec1d+l'],
            kp_cut_y_lo=['spec1d+K'],
            kp_cut_y_hi=['spec1d+L'],

            plot_zoom_rate=1.2,
            plot_pan_pct=0.10,
        )

    def __str__(self):
        return 'spec1d'

    def start(self):
        pass

    def stop(self):
        pass

    #####  SCROLL ACTION CALLBACKS #####

    #####  MOUSE ACTION CALLBACKS #####

    def ms_showxy(self, viewer, event, data_x, data_y):
        """Motion event in the channel viewer window.  Show the pointing
        information under the cursor.
        """
        self.fv.showxy(viewer, data_x, data_y)
        return False
