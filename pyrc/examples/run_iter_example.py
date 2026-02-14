import sys
import os
import subprocess
import glob
from time import time

from pyrc import utils
from pyrc.gain_calibration import GainCalibration
from pyrc.gain_calibration_validation import Validation
from pyrc.utils import delete_validated_files_recursive
from pyrc import (
    RadioCartographerConfig,
    RadioCartographerRunner,
    Channel,
    CalibrationMethod,
    CoordinateSystem,
    TimeShiftMode,
    Receiver,
    CentroidType,
)


def process_fits_file(executable_path: str, input_path: str, output_dir: str) -> None:
    """
    Run gain calibration and radio-cartographer on a single FITS file.

    All radio-cartographer outputs will be written into `output_dir`
    by temporarily changing the working directory.
    """
    start_time = time()
    basename = os.path.basename(input_path)
    stem, ext = os.path.splitext(basename)

    print(f"\n=== Processing {basename} ===")
    #validation

    v = Validation(input_path)
    validated_path = v.validate()

    if not os.path.exists(validated_path):
        print(f"  [WARN] Validated file not found: {validated_path}")
        print("        Using the input file as the validated path instead.")
        validated_path = input_path

    # Gain calibration
    t0 = time()
    g_c = GainCalibration(validated_path, 0, 0, None, None, None, None)
    pre_cal_delta_2, post_cal_delta_2 = g_c.calibrate_gain()

    g_c_1 = GainCalibration(validated_path, 0, 1, None, None, None, None)
    pre_cal_delta_1, post_cal_delta_1 = g_c_1.calibrate_gain()

    print("  Gain calibration deltas:")
    print(f"    precal 1:  {pre_cal_delta_1}")
    print(f"    precal 2:  {pre_cal_delta_2}")
    print(f"    postcal 1: {post_cal_delta_1}")
    print(f"    postcal 2: {post_cal_delta_2}")
    print("  Validation+calibration time:", round(time() - t0, 3), "seconds")

    config = RadioCartographerConfig(
        channel=Channel.COMPOSITE,
        receiver=Receiver.HI,
        calibration_method=CalibrationMethod.INTERPOLATED,
        coordinate_system=CoordinateSystem.EQUATORIAL,
        time_shift_mode=TimeShiftMode.AUTO,
        time_shift_value=0.0,
        generate_raw_map=False,
        min_freq=1350.0,
        max_freq=1750.0,
        bg_scale=6.0,
        rfi_scale=0.1,
        m10_plus_processing=False,
        weight_scale=0.333333,
        photometry_enabled=False,
        photo_inner_radius=1.25,
        photo_outer_radius=5.0,
        photo_centroid_type=CentroidType.BRIGHTEST,
        trim_size=0.0,
        lss_mapping=False,
        gain_delta_start_1=pre_cal_delta_1,
        gain_delta_end_1=post_cal_delta_1,
        gain_delta_start_2=pre_cal_delta_2,
        gain_delta_end_2=post_cal_delta_2,
    )

    runner = RadioCartographerRunner(executable_path)

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Run RC with cwd set to output_dir so all outputs land there
    original_cwd = os.getcwd()
    try:
        os.chdir(output_dir)
        print("  Starting radio-cartographer execution...")
        runner.run(config, input_path, check=True, capture_output=False)
        print("  Execution completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"  [ERROR] Command failed with exit code {e.returncode}")
        raise
    finally:
        os.chdir(original_cwd)

    print("  Total time for this file:", round(time() - start_time, 3), "seconds")
def find_fits_files(input_dir: str):
    """
    Recursively find all .fits files under input_dir, excluding *_validated.fits.
    """
    fits_paths = []
    for root, dirs, files in os.walk(input_dir):
        for name in files:
            lower = name.lower()
            if lower.endswith(".fits") and not lower.endswith("_validated.fits"):
                fits_paths.append(os.path.join(root, name))
    return sorted(fits_paths)

def run_batch(executable_path: str, input_dir: str, output_dir: str) -> None:
    """
    Iterate over all FITS files in input_dir and run processing,
    placing outputs into output_dir.
    """
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)

    print(f"Using executable: {executable_path}")
    print(f"Input directory:  {input_dir}")
    print(f"Output directory: {output_dir}")

    # --- New: sanity check on input directory ---
    if not os.path.isdir(input_dir):
        print(f"[ERROR] Input directory does NOT exist: {input_dir}")
        parent = os.path.dirname(input_dir) or "/"
        print(f"[DEBUG] Parent directory: {parent}")
        if os.path.isdir(parent):
            print("[DEBUG] Contents of parent directory:")
            for entry in os.listdir(parent):
                print("  ", entry)
        else:
            print("[DEBUG] Parent directory also does not exist.")
        return


    # --- New: use recursive search for .fits files ---
    fits_paths = find_fits_files(input_dir)

    if not fits_paths:
        print("No FITS files found in input directory (after recursive search).")
        return

    print(f"Found {len(fits_paths)} FITS files to process:")
    for p in fits_paths:
        print(" ", os.path.relpath(p, input_dir))

    for path in fits_paths:
        try:
            process_fits_file(executable_path, path, output_dir)
        except Exception as e:
            # Continue with other files even if one fails
            print(f"  [ERROR] Failed to process {os.path.basename(path)}: {e}")



if __name__ == "__main__":
    EXECUTABLE_PATH = "/skynet/radio-cartographer/radio-cartographer"
    INPUT_DIR = "testing/test_files/raw_daisy_200s"
    OUTPUT_DIR = "testing/test_files/processed_daisy_200s_rfi01"
    # ------------------------------------------------------------

    run_batch(EXECUTABLE_PATH, INPUT_DIR, OUTPUT_DIR)
    delete_validated_files_recursive(INPUT_DIR, dry_run=False)
