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

def dist(x, y):
    """
    Return the distance between two points.
    """
    d = x - y
    return np.sqrt(np.dot(d, d))


def dist_point_to_segment(p, s0, s1):
    """
    Get the distance of a point to a segment.
      *p*, *s0*, *s1* are *xy* sequences
    This algorithm from
    http://geomalgorithms.com/a02-_lines.html
    """
    v = s1 - s0
    w = p - s0
    c1 = np.dot(w, v)
    if c1 <= 0:
        return dist(p, s0)
    c2 = np.dot(v, v)
    if c2 <= c1:
        return dist(p, s1)
    b = c1 / c2
    pb = s0 + b * v
    return dist(p, pb)


class CrystalMeasurer(ImageHandler):
    """An interactive polygon editor.
    Parameters
    ----------
    poly_xy : None or list of (float, float)
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
        self.ax.patch.set_facecolor('grey')
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
                          "'enter' moves to next task", color='red')

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


def binary_search_interaction(event):
    global has_crystal
    if event.key not in ['y', 'n']:
        return
    if event.key == 'y':
        has_crystal = True
        print('crystal seen')
        plt.close()
    if event.key == 'n':
        has_crystal = False
        print('no crystal')
        plt.close()


def binary_search_for_first_crystal_visible(image_folder, image_filenames):
    # binary search:
    current_ind = len(image_filenames) - 1
    first_ind_seen_with_crystal = None
    last_ind_seen_without_crystal = 0
    found_first_image_with_crystal = False
    global has_crystal
    has_crystal = False
    while found_first_image_with_crystal is False:
        print(
            f"previous image has_crystal {has_crystal}, current_ind {current_ind}, "
            f"last_ind_seen_without_crystal {last_ind_seen_without_crystal}, "
            f"first_ind_seen_with_crystal {first_ind_seen_with_crystal}")
        image_filename = image_filenames[current_ind]
        image = ImageHandler.load_tif(os.path.join(image_folder, image_filename))
        fig, ax = plt.subplots(figsize=(12, 18))
        ax.patch.set_facecolor('grey')
        ax.set_title("Searching for first appearance of crystal\n Hit 'y' if there is crystal, 'n' if not", color='red')
        binary_searcher = fig.canvas.mpl_connect('key_press_event', binary_search_interaction)
        ax.imshow(image)
        plt.show()
        if has_crystal:
            first_ind_seen_with_crystal = current_ind
        elif not has_crystal:
            last_ind_seen_without_crystal = current_ind
        if first_ind_seen_with_crystal is None:
            print("No crystal seen in the last image file- no crystal found")
            break
        if first_ind_seen_with_crystal - last_ind_seen_without_crystal == 1:
            found_first_image_with_crystal = True
            print(f"First crystal found in file {image_filenames[current_ind]}, ind {current_ind}")
            return current_ind
        current_ind = (first_ind_seen_with_crystal + last_ind_seen_without_crystal) // 2


def plot_change_in_levels(observations, output_filename):
    fig, axes = plt.subplots(2, 1, figsize=(8, 6))
    sns.lineplot(x=ELAPSED_HOURS, y=LIQUID_HEIGHT_CM, hue='measurement_name', data=observations, ax=axes[0], marker='o')
    sns.lineplot(x=ELAPSED_HOURS, y=CRYSTAL_SIZE, data=observations, ax=axes[1], marker='o')
    experiment = os.path.splitext(output_filename)[0]
    axes[0].set_title("Liquid height in tube (cm) over time")
    axes[1].set_title("Crystal Size (cm^2) over time")
    fig.tight_layout()
    fig.savefig(f"Liquid height and crystal size over time for {experiment}")



if __name__ == '__main__':
    """
    Does a binary search of the images to find the first witha  crystal, Calibrates the image to a scale,
     tracks the meniscus levels, and outlines the crystal to measure its size
    Usage:
    python calculate_growth_of_crystal.py {path to directory of images} {output csv filename for data} {sampling frequency (every number of files/minutes) to consider for crystal growth}
    e.g.:
    python calculate_growth_of_crystal.py "data/MA_336_1_No FAH" MA_336_1_crystal_size.csv 30
    """



    image_folder = sys.argv[1]
    csv_output_filename = sys.argv[2]
    subset_input_images = int(sys.argv[3])

    # image_folder = "data/MA_336_1_No FAH"
    # csv_output_filename = "MA_336_1_crystal_size.csv"
    # subset_input_images = 30
    image_filenames, _ = get_image_filenames_to_process(image_folder)
    # times are in the filename, so we save the first filename to figure the start time of the experiment
    first_filename_of_experiment = image_filenames[0]

    first_crystal_found_filename_ind = binary_search_for_first_crystal_visible(image_folder, image_filenames)
    # first_crystal_found_filename_ind = 2477

    image_filenames_with_crystal = image_filenames[first_crystal_found_filename_ind:]
    print("These image files should have crystal present:", image_filenames_with_crystal)
    # take one image with crystal per subset_input_images minutes
    subsetted_image_filenames_with_crystal = get_subsetted_image_filenames(image_filenames_with_crystal,
                                                                           subset_input_images)
    print(f"measuring crystal size for  {len(subsetted_image_filenames_with_crystal)} files")

    previous_polygon_vertices = None
    crystal_size_data = []
    calibration = Calibrator(image_filename=os.path.join(image_folder, subsetted_image_filenames_with_crystal[0]))

    for image_filename in subsetted_image_filenames_with_crystal:
        print(f'Processing {image_filename}')
        image_filepath = os.path.join(image_folder, image_filename)
        dm = DiffusionMeasurer(image_filepath, calibration)
        while dm.status != DONE:
            # generate new calibration inputs
            if dm.status == RESET_CALIBRATION:
                calibration = Calibrator(image_filename=os.path.join(image_folder, image_filename))
            dm = DiffusionMeasurer(image_filepath, calibration)

        mc = CrystalMeasurer(image_filepath, calibration=calibration, poly_xy=previous_polygon_vertices)
        # while mc.status != DONE:
        #     # generate new calibration inputs
        #     if mc.status == RESET_CALIBRATION:
        #         calibration = Calibrator(image_filename=os.path.join(image_folder, image_filename))
        #     mc = CrystalMeasurer(image_filepath, calibration=calibration, poly_xy=previous_polygon_vertices)

        previous_polygon_vertices = mc.verts
        area = calculate_area_of_polygon(mc.verts)
        crystal_data = {
            "filename": image_filename,
            "crystal_pixel_area": area,
            'pixels_in_square_cm': (2 * mc.calibration.half_cm_pixel_distance) ** 2,
            'crystal_vertices': mc.verts,
            "y_location_of_mean_crystal_verts": np.mean(mc.verts, axis=0)[1]
        }
        for x in dm.diffusion_data:
            x.update(crystal_data)
            crystal_size_data.append(x)
        print(f"{image_filename} Area of polygon: {area}")
    crystal_size_data_df = pd.DataFrame(crystal_size_data)
    crystal_size_data_df[CRYSTAL_SIZE] = crystal_size_data_df["crystal_pixel_area"] / crystal_size_data_df['pixels_in_square_cm']
    crystal_size_data_df.to_csv(csv_output_filename)

    all_data = pd.read_csv(csv_output_filename)
    all_data['liquid_height_pixels'] = all_data['y'] - all_data['tube_bottom']
    height_pixels_of_mean_crystal_vertex_from_bottom = abs(all_data['y_location_of_mean_crystal_verts'] - all_data['tube_bottom'])
    all_data['height_cm_of_mean_crystal_vertex_from_bottom'] = height_pixels_of_mean_crystal_vertex_from_bottom / all_data[HALF_CM_NUM_PIXELS]
    all_data['file_time'] = all_data.filename.apply(get_filetime_from_filename)
    experiment_start_time = get_filetime_from_filename(first_filename_of_experiment)
    all_data['experiment_start_time'] = experiment_start_time
    all_data[LIQUID_HEIGHT_CM] = abs(all_data['liquid_height_pixels'] / (2 * all_data[HALF_CM_NUM_PIXELS]))
    all_data[ELAPSED_HOURS] = all_data.apply(lambda row: (row['file_time'] - experiment_start_time).total_seconds() / 3600, axis=1)

    all_data.to_csv(csv_output_filename, index=False)

    all_data = pd.read_csv(csv_output_filename)
    plot_change_in_levels(all_data, csv_output_filename)
