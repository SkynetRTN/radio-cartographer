import sys
import os
import argparse
import subprocess
from pyrc import utils
from time import time
from pyrc.gain_calibration import Gain_Calibration
from pyrc.gain_calibration_validation import Validation
# Add the project root to sys.path
# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from pyrc import (
    RadioCartographerConfig,
    RadioCartographerRunner,
    Channel,
    CalibrationMethod,
    CoordinateSystem,
    TimeShiftMode,
    Receiver,
    CentroidType
)

def run_example(executable_path="radio-cartographer"):
    """
    Runs the radio-cartographer on the 'orion' standard input file using
    the parameters requested for verification.
    """
    print(f"Using executable: {executable_path}")
    start_time = time()
    # Use the standard path provided in the request
    filename = "testing/test_files/raw_daisy_200s/0144858.fits"
    v = Validation(filename)
    validated_path = v.validate()
    g_c = Gain_Calibration(validated_path, 0, 0, None, None,
                        None, None)
    precaldelta2, postcaldelta2 = g_c.Gain_calibration()

    g_c_1 =     g_c = Gain_Calibration(validated_path, 0, 1, None, None,
                                       None, None)
    precaldelta1, postcaldelta1  = g_c.Gain_calibration()

    print(precaldelta1, precaldelta2, postcaldelta1, postcaldelta2)
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
        exclusion_bands=[],
        bg_scale=6.0,
        rfi_scale=0.35,
        m10_plus_processing=False,
        weight_scale=0.333333,
        photometry_enabled=False,
        photo_inner_radius=1.25,
        photo_outer_radius=5.0,
        photo_centroid_type=CentroidType.BRIGHTEST,
        trim_size=0.0,
        lss_mapping=False,
        gain_delta_start_1= precaldelta1,
        gain_delta_end_1=postcaldelta1,
        gain_delta_start_2=precaldelta2,
        gain_delta_end_2=postcaldelta2,
    )

    
    # Initialize runner

    runner = RadioCartographerRunner(executable_path)
    
    try:
        # Run with capture_output=False to stream output to stdout/stderr directly,
        # or True to capture and print after. Streaming is usually better for long jobs.
        # But here let's capture to show we can handle it.
        # Actually, let's stream it so the user sees progress if they run it.
        print("Starting execution...")
        runner.run(config, filename, check=True, capture_output=False)
        print("Execution completed successfully.")
        
    except subprocess.CalledProcessError as e:
        print(f"Error: Command failed with exit code {e.returncode}")
        sys.exit(e.returncode)
    except FileNotFoundError:
        print(f"Error: Executable '{executable_path}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Radio Cartographer Example")
    parser.add_argument("--bin", default="radio-cartographer", help="Path to the radio-cartographer binary")
    args = parser.parse_args()
    
    run_example(args.bin)
