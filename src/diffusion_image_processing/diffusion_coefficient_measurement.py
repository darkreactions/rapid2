"""
Copyright [2020] [Two Six Labs, LLC]

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from calculate_growth_of_crystal import ImageHandler
from calibration_image_processor import (
    get_filetime_from_filename,
    get_subsetted_image_filenames,
    get_image_filenames_to_process,
)


def get_title(title):
    return f"{title}\nMove on to next image with 'enter'.\nReset the clicks on this image by hitting 'x'.\nRecalibrate the ruler scale with 'r'. Finish processing images with 'd'"


import os
import sys

import cv2
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.path import Path
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point
import seaborn as sns

SCALE_MEASURES = ["12_cm", "13_cm", "14_cm", "15_cm"]
RESET_CALIBRATION = "reset_calibration"
DONE = "done"
CM_NUM_PIXELS = "cm_num_pixels"
# LIQUID_HEIGHT_CM = 'liquid height (cm)'
ELAPSED_HOURS = "elapsed hours"
FINISH_PROCESSING = "FINISH_PROCESSING"

def dist(x, y):
    """
    Return the distance between two points.
    """
    d = np.array(x) - np.array(y)
    return np.sqrt(np.dot(d, d))

class Calibrator(ImageHandler):
    def __init__(self, image_filename):
        self.image_filename = image_filename
        self.calibration_events = {}
        self.scale_measures = SCALE_MEASURES
        self.events_df = None
        self.nclicks = 0
        self.status = None

        self.fig, self.ax = plt.subplots(figsize=(12, 18))
        self.scale_getter = self.fig.canvas.mpl_connect(
            "button_press_event", self.get_scale
        )
        self.image_resetter = self.fig.canvas.mpl_connect(
            "key_press_event", self.reset_image
        )

        plt.title(get_title(f"click {self.scale_measures[self.nclicks]} line"))
        image = self.load_tif(self.image_filename)
        self.visible_image = self.ax.imshow(image)
        plt.show()
        self.cm_pixel_distance = self.get_cm_from_scale()
        # self.antisolvent_tube_bottom, self.solvent_tube_bottom = self.get_bottom_of_tube()

    def get_scale(self, event):

        if self.nclicks < len(self.scale_measures):
            ann = self.ax.hlines(
                event.ydata, xmin=event.xdata - 30, xmax=event.xdata + 30, color="green"
            )

            self.calibration_events[self.scale_measures[self.nclicks]] = (
                event.xdata,
                event.ydata,
            )
        self.nclicks += 1
        try:
            plt.title(get_title(f"click {self.scale_measures[self.nclicks]} line"))
        except IndexError:
            plt.title(get_title(f"Hit 'enter' to save calibration and continue"))

        plt.draw()

    def reset_image(self, event):
        if event.key not in ["r", "enter"]:
            return
        if event.key == "r":
            self.status = None
            plt.close()

        if event.key == "enter":
            if self.nclicks >= len(self.scale_measures):
                self.status = DONE
                plt.close()

    def get_cm_from_scale(self):
        """
        given a df of observations, and the calibration number corresponding to the image, return the number of pixels associated with 1 cm
        """
        # get calibration points corresponding to scale, sorted by the number on the scale
        measurement_names = sorted(
            [k for k in self.calibration_events.keys() if k.endswith("cm")],
            key=lambda k: float(k.split("_")[0]),
        )
        # get y value for scale measures
        measurements = np.array(
            [self.calibration_events[k][1] for k in measurement_names]
        )
        cm = abs(np.mean(np.diff(measurements)))
        return cm

    #
    # def get_bottom_of_tube(self):
    #     """
    #     :param calibration_number: the index of a calibration definition
    #     :return: a dict valued with the pixel location of the bottom of the tubes in the image
    #     """
    #     return self.calibration_events['bottom left tube'][1], self.calibration_events['bottom right tube'][1]


class DiffusionMeasurer(ImageHandler):
    def __init__(self, image_filename, calibration_object):
        """

        :param image_folder:
        :param output_filename:
        :param out_folder:
        :param subset_input_images: if specified, samples an image every subset_input_images (sorted by time) from a folder with many images
        """
        self.image_filename = image_filename
        self.fig, self.ax = plt.subplots(figsize=(12, 18))
        self.meniscus_getter = self.fig.canvas.mpl_connect(
            "button_press_event", self.meniscus_click
        )
        self.image_resetter = self.fig.canvas.mpl_connect(
            "key_press_event", self.reset_image
        )

        self.calibration = calibration_object
        self.clicks_data = []
        self.diffusion_diffraction_measures = []
        self.nclicks = 0

        self.first_crystal_index = None
        self.status = None

        image = self.load_tif(self.image_filename)
        self.max_y, self.max_x, _ = image.shape
        self.visible_image = self.ax.imshow(image)
        for annotation in self.calibration.calibration_events.values():
            self.ax.hlines(
                annotation[1],
                xmin=annotation[0] - 30,
                xmax=annotation[0] + 30,
                color="green",
            )
        plt.title(get_title("click straight laser line bottom left"))
        self.laser_line_points = []
        self.half_max_points = []
        self.base_width_points = []

        plt.show()

    def meniscus_click(self, event):
        # this is awful, sorry to future self, just not worth optimizing
        if self.nclicks == 0:
            self.clicks_data.append(
                {
                    "filename": self.image_filename,
                    "x": event.xdata,
                    "y": event.ydata,
                    "measurement_name": "bottom_left_line",
                    CM_NUM_PIXELS: self.calibration.cm_pixel_distance,
                }
            )
            # directions for next action
            plt.title(get_title("click straight laser line top right"))
            self.laser_line_points.append((event.xdata, event.ydata))
        elif self.nclicks == 1:
            self.clicks_data.append(
                {
                    "filename": self.image_filename,
                    "x": event.xdata,
                    "y": event.ydata,
                    "measurement_name": "top_right_line",
                    CM_NUM_PIXELS: self.calibration.cm_pixel_distance,
                }
            )
            # directions for next action
            plt.title(get_title("Click top of curve"))
            self.laser_line_points.append((event.xdata, event.ydata))
            # annotate plot with line for straight laser
            x, y = zip(*self.laser_line_points)
            self.ax.plot(x, y, "r-")
        elif self.nclicks == 2:
            self.clicks_data.append(
                {
                    "filename": self.image_filename,
                    "x": event.xdata,
                    "y": event.ydata,
                    "measurement_name": "top_of_curve",
                    CM_NUM_PIXELS: self.calibration.cm_pixel_distance,
                }
            )
            curve_point = (event.xdata, event.ydata)

            x, y = zip(*self.laser_line_points)
            from sklearn.linear_model import LinearRegression
            lm = LinearRegression().fit(np.array(x).reshape(-1, 1), y)
            curve_peak_laser_line_intersect_prediction = lm.predict(
                np.array(curve_point[0]).reshape(1, -1))

            bottom_limit_point = (event.xdata, self.max_y)

            vertical_line_points = [curve_point, bottom_limit_point]
            x, y = zip(*vertical_line_points)
            self.ax.plot(x, y, "g-")

           # plt.title(get_title("Hit enter to save and move to next image"))
            plt.title(get_title("Click on the intersection of the blue half max line with the laser line"))
            vertical_point_of_intersection_coords =  [curve_point[0], curve_peak_laser_line_intersect_prediction[0]]
            vertical_point_of_intersection = Point(*vertical_point_of_intersection_coords)
            self.vertical_point_of_intersection = vertical_point_of_intersection
            self.clicks_data.append(
                {
                    "filename": self.image_filename,
                    "x": vertical_point_of_intersection.x,
                    "y": vertical_point_of_intersection.y,
                    "measurement_name": "vertical_point_of_intersection",
                    CM_NUM_PIXELS: self.calibration.cm_pixel_distance,
                }
            )
            self.ax.plot(
                vertical_point_of_intersection.x, vertical_point_of_intersection.y, "go"
            )

            vertical_distance = dist(vertical_point_of_intersection, curve_point)
            self.vertical_distance = vertical_distance

            self.diffusion_diffraction_measures.append({"filename": self.image_filename,
                                             "measurement_name": "vertical_distance",
                                             "distance_pixels": vertical_distance,
                                             "distance_cm": vertical_distance / self.calibration.cm_pixel_distance
                                                          })

            half_max_center = np.array(vertical_point_of_intersection_coords) - np.array([0, vertical_distance/2])
            self.ax.hlines(half_max_center[1], xmin=half_max_center[0] - 100, xmax=half_max_center[0] + 100, color='blue')

        # get half max measures
        elif self.nclicks == 3:
            self.clicks_data.append(
                {
                    "filename": self.image_filename,
                    "x": event.xdata,
                    "y": event.ydata,
                    "measurement_name": "half_max_left_point",
                    CM_NUM_PIXELS: self.calibration.cm_pixel_distance,
                }
            )
            self.ax.plot(event.xdata, event.ydata, "bo")
            self.half_max_points.append((event.xdata, event.ydata))

        elif self.nclicks == 4:
            self.clicks_data.append(
                {
                    "filename": self.image_filename,
                    "x": event.xdata,
                    "y": event.ydata,
                    "measurement_name": "half_max_right_point",
                    CM_NUM_PIXELS: self.calibration.cm_pixel_distance,
                }
            )
            self.ax.plot(event.xdata, event.ydata, "bo")
            self.half_max_points.append((event.xdata, event.ydata))

            plt.title(get_title("Click laser line at left and right base of curve to measure width"))
            half_max_width = dist(self.half_max_points[0], self.half_max_points[1])
            self.diffusion_diffraction_measures.append({"filename": self.image_filename,
                                                    "measurement_name": "half_max_width",
                                                    "distance_pixels": half_max_width,
                                                    "distance_cm": half_max_width / self.calibration.cm_pixel_distance
                                                    })
            # get half max measures
        elif self.nclicks == 5:
            self.clicks_data.append(
                {
                    "filename": self.image_filename,
                    "x": event.xdata,
                    "y": event.ydata,
                    "measurement_name": "base_width_left_point",
                    CM_NUM_PIXELS: self.calibration.cm_pixel_distance,
                }
            )
            self.ax.plot(event.xdata, event.ydata, "bo")
            self.base_width_points.append((event.xdata, event.ydata))

        elif self.nclicks == 6:
            self.clicks_data.append(
                {
                    "filename": self.image_filename,
                    "x": event.xdata,
                    "y": event.ydata,
                    "measurement_name": "base_width_right_point",
                    CM_NUM_PIXELS: self.calibration.cm_pixel_distance,
                }
            )
            self.ax.plot(event.xdata, event.ydata, "bo")
            self.base_width_points.append((event.xdata, event.ydata))

            plt.title(get_title("Hit enter to save and move to next image"))
            base_width = dist(self.base_width_points[0], self.base_width_points[1])
            self.diffusion_diffraction_measures.append({"filename": self.image_filename,
                                                        "measurement_name": "base_width_diagonal",
                                                        "distance_pixels": base_width,
                                                        "distance_cm": base_width / self.calibration.cm_pixel_distance
                                                        })
            # does not follow the laser line, just measures horizontal width
            base_width_horizontal = abs(self.base_width_points[0][0] - self.base_width_points[1][0])

            self.diffusion_diffraction_measures.append({"filename": self.image_filename,
                                                        "measurement_name": "base_width_horizontal",
                                                        "distance_pixels": base_width_horizontal,
                                                        "distance_cm": base_width_horizontal / self.calibration.cm_pixel_distance
                                                        })

        else:
            # plt.close()
            return
        self.nclicks += 1
        plt.draw()

    def reset_image(self, event):
        # something has gone wrong with user input- user asks to try again with same image
        if event.key not in ["x", "r", "enter", 'd']:
            return
        elif event.key == "x":
            print("resetting image clicks")
            self.status = None
            plt.close()
        elif event.key == "r":
            self.status = RESET_CALIBRATION
            plt.close()
        if event.key == "d":
            self.status = FINISH_PROCESSING
            plt.close()
        elif event.key == "enter":
            self.status = DONE
            plt.close()



def plot_laser_peak_height(observations, output_filename):
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.lineplot(x=ELAPSED_HOURS, y="distance_cm", hue='measurement_name', data=observations, ax=ax, marker='o')
    experiment = os.path.splitext(output_filename)[0]
    ax.set_title("Laser peak height (cm)")
    fig.tight_layout()
    fig.savefig(f"Laser peak height for {experiment}")




if __name__ == "__main__":
    """
    Usage 
    python diffusion_coefficient_measurement.py {path to images folder} {output filename} {frame sample rate}

    python diffusion_coefficient_measurement.py "data/Aug_5_DCM_DMSO_1" Aug_5_DCM_DMSO_1.csv 60
    """
    image_folder = sys.argv[1]
    csv_output_filename = sys.argv[2]
    subset_input_images = int(sys.argv[3])


    print(
        f"Location of images: {image_folder}, output file {csv_output_filename}, frame subset rate {subset_input_images}"
    )

    diffusion_measurer_data = []
    image_filenames, all_events = get_image_filenames_to_process(
        image_folder, subset_input_images
    )
    image_filename = image_filenames[0]
    calibration = Calibrator(image_filename=os.path.join(image_folder, image_filename))
    for image_filename in image_filenames:
        dm = DiffusionMeasurer(os.path.join(image_folder, image_filename), calibration)
        # reset was triggered or window was closed
        while dm.status not in [DONE, FINISH_PROCESSING]:
            # generate new calibration inputs
            if dm.status == RESET_CALIBRATION:
                calibration = Calibrator(
                    image_filename=os.path.join(image_folder, image_filename)
                )
            # if the plot
            dm = DiffusionMeasurer(
                os.path.join(image_folder, image_filename), calibration
            )
        # add the rows to the data output
        diffusion_measurer_data.append(pd.DataFrame(dm.diffusion_diffraction_measures))
        if dm.status == FINISH_PROCESSING:
            break
    all_data = pd.concat(diffusion_measurer_data).reset_index()

    all_data['file_time'] = all_data.filename.apply(get_filetime_from_filename)
    start_time = all_data['file_time'].sort_values()[0]
    all_data[ELAPSED_HOURS] = all_data.apply(lambda row: (row['file_time'] - start_time).total_seconds()/3600, axis=1)
    csv_output_filename = f"laser_diffusion_coefficient_measurement_{image_folder.split('/')[-1]}.csv"
    all_data.to_csv(csv_output_filename, index=False)
    #
    all_data = pd.read_csv(csv_output_filename)
    plot_laser_peak_height(all_data, csv_output_filename)
