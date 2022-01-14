import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cv2
import seaborn as sns

import os
import re
import datetime
import sys

SCALE_MEASURES = ["0_cm", "0.5_cm", "1_cm", "1.5_cm", "2_cm", 'bottom left tube', 'bottom right tube']
RESET_CALIBRATION = 'reset_calibration'
DONE = 'done'
HALF_CM_NUM_PIXELS = 'half_cm_num_pixels'
LIQUID_HEIGHT_CM = 'liquid height (cm)'
ELAPSED_HOURS = 'elapsed hours'


def get_title(title):
    return f"{title}\nReset the clicks on this image by hitting 'x'. Recalibrate the ruler scale with 'r'"


def get_filetime_from_filename(filepath):
    """
    :param filename:
    :return: datetime from string as formatted by the Haver
    """
    FILENAME_PATTERN = "(\d{8}-\d{6})-\S+.tif"
    TIME_STRING_FORMAT = "%Y%m%d-%H%M%S"
    filename = os.path.basename(filepath)
    match = re.match(FILENAME_PATTERN, filename)
    if match:
        time_string = match.groups()[0]
        file_time = datetime.datetime.strptime(time_string, TIME_STRING_FORMAT)
        return file_time


def get_image_filenames_to_process(image_folder, subset_input_images: int = 1):
    # # if we've partially processed the data, don't make us do it again
    files_already_processed = []
    existing_events_dfs = []
    # if os.path.exists(output_filename):
    #     print(f'calibration file {output_filename} found')
    #     existing_events_df = pd.read_csv(output_filename)
    #     existing_events_dfs = [existing_events_df]
    #     print(existing_events_df.head())
    #     files_already_processed = existing_events_df['filename'].unique()
    #     print(f"Already processed: {files_already_processed}")
    image_dir_list = os.listdir(image_folder)
    if subset_input_images:
        print(f"Subsetting images, using every {subset_input_images} images")
        image_dir_list = get_subsetted_image_filenames(image_dir_list, subset_input_images)
    image_filenames = [x for x in image_dir_list if
                       x.endswith('.tif') and x not in files_already_processed]
    print(f"Now processing {len(image_filenames)} {image_filenames}")
    return image_filenames, existing_events_dfs


def get_subsetted_image_filenames(image_dir_list, subset_input_images: int =1):
    # get timestamps from files
    file_times = []
    for filename in image_dir_list:
        file_time = get_filetime_from_filename(filename)
        if file_time:
            file_times.append((filename, file_time))
    # sort the files by creation timestamp and take every nth one for n=subset_input_images
    subsetted_files = sorted(file_times, key=lambda x: x[1])[::subset_input_images]
    return [x[0] for x in subsetted_files]


class ImageHandler:
    @staticmethod
    def load_tif(image_filepath):
        image = cv2.imread(image_filepath)
        # flip the color channels to match the tif input
        image = image[:, :, ::-1]
        return image


class Calibrator(ImageHandler):

    def __init__(self, image_filename):
        self.image_filename = image_filename
        self.calibration_events = {}
        self.scale_measures = SCALE_MEASURES
        self.events_df = None
        self.nclicks = 0
        self.status = None

        self.fig, self.ax = plt.subplots(figsize=(12, 18))
        self.scale_getter = self.fig.canvas.mpl_connect('button_press_event', self.get_scale)
        self.image_resetter = self.fig.canvas.mpl_connect('key_press_event', self.reset_image)

        plt.title(get_title(f'click {self.scale_measures[self.nclicks]} line'))
        image = self.load_tif(self.image_filename)
        self.visible_image = self.ax.imshow(image)
        plt.show()
        self.half_cm_pixel_distance = self.get_half_cm_from_scale()
        self.antisolvent_tube_bottom, self.solvent_tube_bottom = self.get_bottom_of_tube()

    def get_scale(self, event):

        if self.nclicks < len(self.scale_measures):
            ann = self.ax.hlines(event.ydata, xmin=event.xdata - 30, xmax=event.xdata + 30, color='green')

            self.calibration_events[self.scale_measures[self.nclicks]] = (event.xdata, event.ydata)
        self.nclicks += 1
        try:
            plt.title(get_title(f'click {self.scale_measures[self.nclicks]} line'))
        except IndexError:
            plt.title(get_title(f"Hit 'enter' to save calibration and continue"))

        plt.draw()

    def reset_image(self, event):
        if event.key not in ['r', 'enter']:
            return
        if event.key == 'r':
            self.status = None
            plt.close()
        if event.key == 'enter':
            if self.nclicks >= len(self.scale_measures):
                self.status = DONE
                plt.close()

    def get_half_cm_from_scale(self):
        """
        given a df of observations, and the calibration number corresponding to the image, return the number of pixels associated with 0.5 cm
        """
        # get calibration points corresponding to scale, sorted by the number on the scale
        measurement_names = sorted([k for k in self.calibration_events.keys() if k.endswith('cm')], key= lambda k: float(k.split('_')[0]))
        # get y value for scale measures
        measurements = np.array([self.calibration_events[k][1] for k in measurement_names])
        half_cm = abs(np.mean(np.diff(measurements)))
        return half_cm

    def get_bottom_of_tube(self):
        """
        :param calibration_number: the index of a calibration definition
        :return: a dict valued with the pixel location of the bottom of the tubes in the image
        """
        return self.calibration_events['bottom left tube'][1], self.calibration_events['bottom right tube'][1]


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
        self.meniscus_getter = self.fig.canvas.mpl_connect('button_press_event', self.meniscus_click)
        self.image_resetter = self.fig.canvas.mpl_connect('key_press_event', self.reset_image)

        self.calibration = calibration_object
        self.diffusion_data = []
        self.nclicks = 0

        self.first_crystal_index = None
        self.status = None

        image = self.load_tif(self.image_filename)
        self.visible_image = self.ax.imshow(image)
        for annotation in self.calibration.calibration_events.values():
            self.ax.hlines(annotation[1], xmin=annotation[0] - 30, xmax=annotation[0] + 30, color='green')
        plt.title(get_title('click left tube meniscus'))
        plt.show()

    def meniscus_click(self, event):
        if self.nclicks == 0:
            ann = self.ax.hlines(event.ydata, xmin=event.xdata - 30, xmax=event.xdata + 30, color='red')
            self.diffusion_data.append(
                {"filename": self.image_filename, 'x': event.xdata, 'y': event.ydata, 'measurement_name': 'antisolvent',
                 HALF_CM_NUM_PIXELS: self.calibration.half_cm_pixel_distance,
                 'tube_bottom': self.calibration.antisolvent_tube_bottom})
            # directions for next action
            plt.title(get_title('click right tube meniscus'))
        elif self.nclicks == 1:
            ann = self.ax.hlines(event.ydata, xmin=event.xdata - 30, xmax=event.xdata + 30, color='blue')
            self.diffusion_data.append(
                {"filename": self.image_filename, 'x': event.xdata, 'y': event.ydata, 'measurement_name': 'solvent',
                 HALF_CM_NUM_PIXELS: self.calibration.half_cm_pixel_distance,
                 'tube_bottom': self.calibration.solvent_tube_bottom})
            # directions for next action
            plt.title(get_title('Hit enter to save and move to next image'))
        else:
            plt.close()
            return
        self.nclicks += 1
        plt.draw()

    def reset_image(self, event):
        # something has gone wrong with user input- user asks to try again with same image
        if event.key not in ['x', 'r', 'enter']:
            return
        elif event.key == 'x':
            print('resetting image clicks')
            self.status = None
            plt.close()
        elif event.key == 'r':
            self.status = RESET_CALIBRATION
            plt.close()
        elif event.key == 'enter':
            self.status = DONE
            plt.close()


def plot_change_in_levels(observations, output_filename):
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.lineplot(x=ELAPSED_HOURS, y=LIQUID_HEIGHT_CM, hue='measurement_name', data=observations, ax=ax, marker='o')
    experiment = os.path.splitext(output_filename)[0]
    ax.set_title("Liquid height in tube (cm) over time")
    fig.savefig(f"Liquid height level over time for {experiment}")


if __name__ == '__main__':
    """
    Usage 
    python calibration_image_processor.py {path to images folder} {output filename} {frame sample rate}

    python calibration_image_processor.py "data/MA_336_1_No FAH" MA_336_1_NoFAH.csv 60
    """
    image_folder = sys.argv[1]
    csv_output_filename = sys.argv[2]
    subset_input_images = int(sys.argv[3])

    # image_folder = "data/MA_336_1_No FAH"
    # csv_output_filename = "MA_336_1_NoFAH.csv"
    # out_folder = 'out'
    # subset_input_images = 60


    print(
        f"Location of images: {image_folder}, output file {csv_output_filename}, frame subset rate {subset_input_images}")

    diffusion_measurer_data = []
    image_filenames, all_events = get_image_filenames_to_process(image_folder, subset_input_images)
    image_filename = image_filenames[0]
    # todo: if you try to reset the calibration on the first go it just exits
    calibration = Calibrator(image_filename=os.path.join(image_folder, image_filename))
    for image_filename in image_filenames:
        dm = DiffusionMeasurer(os.path.join(image_folder, image_filename), calibration)
        # reset was triggered or window was closed
        while dm.status != DONE:
            # generate new calibration inputs
            if dm.status == RESET_CALIBRATION:
                calibration = Calibrator(image_filename=os.path.join(image_folder, image_filename))
            # if the plot
            dm = DiffusionMeasurer(os.path.join(image_folder, image_filename), calibration)
        # add the rows to the data output
        diffusion_measurer_data.append(pd.DataFrame(dm.diffusion_data))

    # process data with calibration, output to graph and csv
    all_data = pd.concat(diffusion_measurer_data).reset_index()
    all_data['liquid_height_pixels'] = all_data['y'] - all_data['tube_bottom']
    all_data['file_time'] = all_data.filename.apply(get_filetime_from_filename)
    all_data[LIQUID_HEIGHT_CM] = abs(all_data['liquid_height_pixels'] / (2 * all_data[HALF_CM_NUM_PIXELS]))
    start_time = all_data['file_time'].sort_values()[0]
    all_data[ELAPSED_HOURS] = all_data.apply(lambda row: (row['file_time'] - start_time).total_seconds()/3600, axis=1)
    all_data.to_csv(csv_output_filename, index=False)

    all_data = pd.read_csv(csv_output_filename)
    plot_change_in_levels(all_data, csv_output_filename)
