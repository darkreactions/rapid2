"""
Interactive tool to draw mask on an image or image-like array.
Adapted from https://gist.github.com/tonysyu/3090704
"""

import os
import sys

import cv2
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.path import Path
import numpy as np
import pandas as pd
from shapely.geometry import Polygon as shapely_polygon
import seaborn as sns

from calibration_image_processor import (get_image_filenames_to_process, get_subsetted_image_filenames, Calibrator,
                                         RESET_CALIBRATION, DONE, ImageHandler, HALF_CM_NUM_PIXELS, DiffusionMeasurer,
                                         LIQUID_HEIGHT_CM, ELAPSED_HOURS, get_filetime_from_filename)


POLYGON_COLOR = 'b'
CRYSTAL_SIZE = 'crystal_visible_surface_area'


class CrystalMeasurer(ImageHandler):
    """An interactive polygon editor.
    Parameters
    ----------
    poly_xy : list of (float, float)
        List of (x, y) coordinates used as vertices of the polygon.
    max_ds : float
        Max pixel distance to count as a vertex hit.
    Key-bindings
    ------------
    't' : toggle vertex markers on and off.  When vertex markers are on,
          you can move them, delete them
    'd' : delete the vertex under point
    'i' : insert a vertex at point.  You must be within max_ds of the
          line connecting two existing vertices
    """

    def __init__(self, image_filepath, calibration, poly_xy=None, max_ds=10, polygon_color=POLYGON_COLOR):
        self.showverts = True
        self.image = self.load_tif(image_filepath)

        self.max_ds = max_ds
        self.fig, self.ax = plt.subplots(figsize=(12, 8))
        self.visible_image = self.ax.imshow(self.image)
        self.status = None
        self.calibration = calibration

        self.verts = None
        if poly_xy is None:
            poly_xy = self.default_vertices(self.ax)
        self.poly = Polygon(poly_xy, animated=True, fc=polygon_color, ec='none', alpha=0.0)

        self.ax.add_patch(self.poly)
        self.ax.set_clip_on(False)
        self.ax.set_title("Surround the crystal with the points. Click and drag a point to move it;\n"
                          "While hovering over a point, type 'i' to insert new point; 'd' to delete.\n"
                          "'t' toggles visibility of vertex points\n"
                          # "'r' resets calibration\n"
                          "'enter' moves to next task")

        x, y = zip(*self.poly.xy)
        self.line = plt.Line2D(x, y, color='none', marker='o', mfc='r', alpha=0.6, animated=True)
        self._update_line()
        self.ax.add_line(self.line)

        self.poly.add_callback(self.poly_changed)
        self._ind = None  # the active vert

        canvas = self.poly.figure.canvas
        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas = canvas
        self.background = None
        self.redraw_masked_canvas()
        for annotation in self.calibration.calibration_events.values():
            self.ax.hlines(annotation[1], xmin=annotation[0] - 30, xmax=annotation[0] + 30, color='green')
        plt.show()

    def get_mask(self):
        """Return image mask given by mask creator"""
        h, w, color = self.image.shape
        y, x = np.mgrid[:h, :w]
        points = np.transpose((x.ravel(), y.ravel()))
        path = Path(self.verts)
        mask = path.contains_points(points)
        return mask.reshape(h, w)

    def poly_changed(self, poly):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        # Artist.update_from(self.line, poly)
        self.line.set_visible(vis)  # don't use the poly visibility state

    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        self.canvas.blit(self.ax.bbox)

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        ignore = not self.showverts or event.inaxes is None or event.button != 1
        if ignore:
            return
        self._ind = self.get_ind_under_cursor(event)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        ignore = not self.showverts or event.button != 1
        if ignore:
            return
        self._ind = None
        self.redraw_masked_canvas()

    def key_press_callback(self, event):
        """whenever a key is pressed"""
        # if not event.inaxes:
        #     return
        if event.key == 't':
            # toggle visibility of vertices
            self.showverts = not self.showverts
            self.line.set_visible(self.showverts)
            if not self.showverts:
                self._ind = None
        elif event.key == 'd':
            # delete outline point for crystal
            ind = self.get_ind_under_cursor(event)
            if ind is None:
                return
            if ind == 0 or ind == self.last_vert_ind:
                print("Cannot delete root node")
                return
            self.poly.xy = [tup for i, tup in enumerate(self.poly.xy)
                            if i != ind]
            self._update_line()
        # elif event.key == 'r':
        #     # reset calibration of the scale
        #     self.status = RESET_CALIBRATION
        #     plt.close()
        elif event.key == 'enter':
            # todo: segmentation faults when I try to save image
            # plt.savefig(f'test{image_filename}.png')
            self.status = DONE
            plt.close()
            return
        elif event.key == 'i':
            # insert a new vertex
            xys = self.poly.get_transform().transform(self.poly.xy)
            p = event.x, event.y  # cursor coords
            for i in range(len(xys) - 1):
                s0 = xys[i]
                s1 = xys[i + 1]
                d = dist_point_to_segment(p, s0, s1)
                if d <= self.max_ds:
                    self.poly.xy = np.array(
                        list(self.poly.xy[:i + 1]) +
                        [(event.xdata, event.ydata)] +
                        list(self.poly.xy[i + 1:]))
                    self._update_line()
                    break
        self.redraw_masked_canvas()

    def redraw_masked_canvas(self):
        mask = self.get_mask()
        new_visible_image = self.image.copy()
        new_visible_image[~mask] = np.uint8(np.clip(new_visible_image[~mask] - 100., 0, 255))
        self.visible_image.set_data(new_visible_image)

        self.canvas.draw()

    def motion_notify_callback(self, event):
        'on mouse movement'
        ignore = (not self.showverts or event.inaxes is None or
                  event.button != 1 or self._ind is None)
        if ignore:
            return
        x, y = event.xdata, event.ydata

        if self._ind == 0 or self._ind == self.last_vert_ind:
            self.poly.xy[0] = x, y
            self.poly.xy[self.last_vert_ind] = x, y
        else:
            self.poly.xy[self._ind] = x, y
        self._update_line()

        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        self.canvas.blit(self.ax.bbox)

    def _update_line(self):
        # save verts because polygon gets deleted when figure is closed
        self.verts = self.poly.xy
        self.last_vert_ind = len(self.poly.xy) - 1
        self.line.set_data(zip(*self.poly.xy))

    def get_ind_under_cursor(self, event):
        'get the index of the vertex under cursor if within max_ds tolerance'
        # display coords
        xy = np.asarray(self.poly.xy)
        xyt = self.poly.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.sqrt((xt - event.x) ** 2 + (yt - event.y) ** 2)
        indseq = np.nonzero(np.equal(d, np.amin(d)))[0]
        ind = indseq[0]
        if d[ind] >= self.max_ds:
            ind = None
        return ind

    @staticmethod
    def default_vertices(ax):
        """Default to rectangle that has a quarter-width/height border."""
        xlims = ax.get_xlim()
        ylims = ax.get_ylim()
        w = np.diff(xlims)
        h = np.diff(ylims)
        x1, x2 = xlims + w // 4 * np.array([1, -1])
        y1, y2 = ylims + h // 4 * np.array([1, -1])
        return ((x1, y1), (x1, y2), (x2, y2), (x2, y1))


def calculate_area_of_polygon(vertices):
    try:
        polygon = shapely_polygon(vertices)
        return polygon.area
    except ValueError:
        # there are not enough points in the Polygon to calculate area
        return 0


if __name__ == '__main__':
    """
    Creates an edited crystal growth analysis output file that contains the height of the crystal nucleation point
    
    Usage:
    python add_crystal_height_to_existing_files.py analysis_output_filename image_folder
    I.e., 
    python add_crystal_height_to_existing_files.py {output csv from calculate_growth_of_crystal.py} {path to folder with images on which that analysis has been done}
    e.g.: python add_crystal_height_to_existing_files.py calculate_growth_of_crystal_MA_336_3_with_FAH.csv /Users/nick.leiby/repos/chaos-perovskite/exp/diffusion_image_processing/data/MA_336_3_with_FAH_CyclohexMethyl
    
    """

    analysis_output_filename = sys.argv[1]
    image_folder = sys.argv[2]

    # find filename with first crystal from output file from
    df = pd.read_csv(analysis_output_filename)
    df['crystal_vertices'] = None
    df['height_cm_of_mean_crystal_vertex_from_bottom'] = None
    filenames_with_visible_crystals = df['filename'].values

    first_image_file_with_crystal = filenames_with_visible_crystals[0]
    image_filepath = os.path.join(image_folder, first_image_file_with_crystal)
    # open the calibrator to label scale and label the vertices at the start of nucleation

    calibration = Calibrator(image_filename=image_filepath)
    previous_polygon_vertices = None

    mc = CrystalMeasurer(image_filepath, calibration=calibration)
    mean_crystal_vertices_location = np.mean(mc.verts, axis=0)[1]
    num_pixels_from_bottom = mean_crystal_vertices_location - calibration.solvent_tube_bottom
    crystal_height = abs(num_pixels_from_bottom/calibration.half_cm_pixel_distance)
    df.at[0, 'height_cm_of_mean_crystal_vertex_from_bottom'] = crystal_height
    df.at[0, 'crystal_vertices'] = mc.verts
    output_filename = "initial_crystal_height_" + analysis_output_filename
    df.to_csv(output_filename, index=False)
