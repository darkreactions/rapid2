import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


import os
import sys

SCALE_MEASURES = [
    "0_cm",
    "0.5_cm",
    "1_cm",
    "1.5_cm",
    "2_cm",
    "bottom left tube",
    "bottom right tube",
]
RESET_CALIBRATION = "reset_calibration"
DONE = "done"
HALF_CM_NUM_PIXELS = "half_cm_num_pixels"
LIQUID_HEIGHT_CM = "liquid height (cm)"
ELAPSED_HOURS = "elapsed hours"

from calibration_image_processor import (
    ImageHandler,
    DiffusionMeasurer,
    Calibrator,
    get_image_filenames_to_process,
    get_filetime_from_filename,
    plot_change_in_levels,
)
from calculate_growth_of_crystal import CrystalMeasurer


def binary_search_interaction(event):
    global has_crystal
    if event.key not in ["y", "n"]:
        return
    if event.key == "y":
        has_crystal = True
        print("crystal seen")
        plt.close()
    if event.key == "n":
        has_crystal = False
        print("no crystal")
        plt.close()
    if event.key == "enter":
        plt.close()

        raise KeyError("No crystal in this data set, abort")


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
            f"first_ind_seen_with_crystal {first_ind_seen_with_crystal}"
        )
        image_filename = image_filenames[current_ind]
        image = ImageHandler.load_tif(os.path.join(image_folder, image_filename))
        fig, ax = plt.subplots(figsize=(12, 18))
        ax.patch.set_facecolor('grey')

        ax.set_title(
            "Searching for first appearance of crystal\n Hit 'y' if there is crystal, 'n' if not",
            color='red'
        )
        binary_searcher = fig.canvas.mpl_connect(
            "key_press_event", binary_search_interaction
        )

        ax.imshow(image)
        plt.show()
        if has_crystal:
            first_ind_seen_with_crystal = current_ind
        elif not has_crystal:
            last_ind_seen_without_crystal = current_ind

        if first_ind_seen_with_crystal is None:
            # we start with the last frame in the while loop, so if there is no crystal we exit immediately
            print("No crystal seen in the last image file- no crystal found")
            break
        if first_ind_seen_with_crystal - last_ind_seen_without_crystal == 1:
            found_first_image_with_crystal = True
            print(
                f"First crystal found in file {image_filenames[current_ind]}, ind {first_ind_seen_with_crystal}"
            )
            return first_ind_seen_with_crystal
        current_ind = (first_ind_seen_with_crystal + last_ind_seen_without_crystal) // 2


if __name__ == "__main__":
    """
    Script to collect all liquid heights and crystal start time and height, not not crystal growth over time
    
    Usage 
    python get_liquid_height_and_crystal_start.py {path to images folder} {output filename} {frame sample rate}

    python get_liquid_height_and_crystal_start.py "data/MA_336_1_No FAH" MA_336_1_NoFAH.csv" 60
    """
    image_folder = sys.argv[1]
    csv_output_filename = sys.argv[2]
    subset_input_images = int(sys.argv[3])

    # image_folder = "data/MA_336_1_No FAH"
    # csv_output_filename = "MA_336_1_NoFAH.csv"
    # out_folder = 'out'
    # subset_input_images = 60

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
        while dm.status != DONE:
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
        diffusion_measurer_data.append(pd.DataFrame(dm.diffusion_data))

    # process data with calibration, output to graph and csv
    all_data = pd.concat(diffusion_measurer_data).reset_index(drop=True)
    all_data['crystal_present'] = 0
    all_data["y_location_of_mean_crystal_verts"] = None
    all_data['height_cm_of_mean_crystal_vertex_from_bottom'] = None


    # ask if this experiment has crystallization, and if so, find the start time and height of the crystal
    crystal_image_filenames, all_events = get_image_filenames_to_process(
        image_folder, subset_input_images=1
    )
    first_crystal_found_filename_ind = binary_search_for_first_crystal_visible(
        image_folder, crystal_image_filenames)
    mc = None
    if first_crystal_found_filename_ind:
        # skip this if htere is no crystal
        previous_polygon_vertices = None
        first_crystal_found_image_filename = crystal_image_filenames[first_crystal_found_filename_ind]
        crystal_image_filepath = os.path.join(image_folder, first_crystal_found_image_filename)
        dm = DiffusionMeasurer(crystal_image_filepath, calibration)

        mc = CrystalMeasurer(crystal_image_filepath, calibration=calibration,
                             poly_xy=previous_polygon_vertices)
        nucleation_rows = pd.DataFrame(dm.diffusion_data)
        nucleation_rows['crystal_present'] = 1
        nucleation_rows["y_location_of_mean_crystal_verts"] = np.mean(mc.verts, axis=0)[1]
        height_pixels_of_mean_crystal_vertex_from_bottom = abs(
            nucleation_rows['y_location_of_mean_crystal_verts'] - nucleation_rows['tube_bottom'])

        nucleation_rows['height_cm_of_mean_crystal_vertex_from_bottom'] = height_pixels_of_mean_crystal_vertex_from_bottom / (2*nucleation_rows[HALF_CM_NUM_PIXELS])
        all_data = pd.concat([all_data, nucleation_rows]).sort_values("filename").reset_index(drop=True)
        first_crystal_index = (all_data.crystal_present.values == 1).argmax()
        all_data.loc[all_data.index > first_crystal_index, 'crystal_present'] = 1
    all_data["liquid_height_pixels"] = all_data["y"] - all_data["tube_bottom"]
    all_data["file_time"] = all_data.filename.apply(get_filetime_from_filename)
    all_data[LIQUID_HEIGHT_CM] = abs(
        all_data["liquid_height_pixels"] / (2*all_data[HALF_CM_NUM_PIXELS])
    )
    start_time = all_data["file_time"].sort_values()[0]
    all_data[ELAPSED_HOURS] = all_data.apply(
        lambda row: (row["file_time"] - start_time).total_seconds() / 3600, axis=1
    )
    all_data["elapsed_minutes"] = all_data.apply(
        lambda row: (row["file_time"] - start_time).total_seconds() / 60, axis=1
    )

    all_data.to_csv(csv_output_filename, index=False)

    all_data = pd.read_csv(csv_output_filename)
    plot_change_in_levels(all_data, csv_output_filename)
