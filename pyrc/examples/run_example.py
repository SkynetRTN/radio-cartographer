import sys
import argparse
import subprocess
import tomllib
from pathlib import Path

# Add the parent directory to sys.path to resolve the package
sys.path.append(str(Path(__file__).parent.parent))

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
from pyrc.gain_calibration import Gain_Calibration
from pyrc.gain_calibration_validation import Validation


def load_config(config_path="run_config.toml"):
    """Loads configuration from a TOML file."""
    with open(config_path, "rb") as f:
        return tomllib.load(f)


def run_example(executable_path, config_file="run_config.toml"):
    """
    Runs the radio-cartographer using parameters from a configuration file.
    """
    print(f"Using executable: {executable_path}")

    # Load configuration
    config_path = Path(__file__).parent / config_file
    if not config_path.exists():
        print(f"Error: Configuration file '{config_path}' not found.")
        sys.exit(1)

    print(f"Loading configuration from: {config_path}")
    toml_config = load_config(config_path)

    rc_config = toml_config.get("radio_cartographer", {})
    input_config = toml_config.get("input", {})
    filenames = input_config.get("filenames", [])

    if not filenames:
        print("Error: No input filenames specified in configuration.")
        sys.exit(1)

    # Validation
    # Validate the first file for gain cal
    v = Validation(filenames[0])
    validated_path = v.validate()

    # Gain Calibration on the first file
    # Note: Gain calibration logic remains hardcoded as per original example structure
    # unless we want to parameterize these steps too.
    g_c = Gain_Calibration(validated_path, 0, 0, None, None, None, None)
    precaldelta2, postcaldelta2 = g_c.Gain_calibration()

    g_c_1 = Gain_Calibration(validated_path, 0, 1, None, None, None, None)
    precaldelta1, postcaldelta1 = g_c_1.Gain_calibration()

    print(f"Gain deltas: {precaldelta1}, {precaldelta2}, {postcaldelta1}, {postcaldelta2}")

    # Map string enums to actual Enum values
    # Helper to map string to enum, default to provided default if not found or invalid
    def get_enum(enum_cls, key, default):
        val = rc_config.get(key)
        if val:
            try:
                return enum_cls[val]
            except KeyError:
                print(f"Warning: Invalid enum value '{val}' for {key}. Using default.")
        return default

    config = RadioCartographerConfig(
        channel=get_enum(Channel, "channel", Channel.COMPOSITE),
        receiver=get_enum(Receiver, "receiver", Receiver.HI),
        calibration_method=get_enum(CalibrationMethod, "calibration_method", CalibrationMethod.INTERPOLATED),
        coordinate_system=get_enum(CoordinateSystem, "coordinate_system", CoordinateSystem.EQUATORIAL),
        time_shift_mode=get_enum(TimeShiftMode, "time_shift_mode", TimeShiftMode.AUTO),
        time_shift_value=rc_config.get("time_shift_value", 0.0),
        generate_raw_map=rc_config.get("generate_raw_map", False),
        min_freq=rc_config.get("min_freq", 1350.0),
        max_freq=rc_config.get("max_freq", 1750.0),
        exclusion_bands=rc_config.get("exclusion_bands", []),
        bg_scale=rc_config.get("bg_scale", 6.0),
        rfi_scale=rc_config.get("rfi_scale", 0.7),
        m10_plus_processing=rc_config.get("m10_plus_processing", False),
        weight_scale=rc_config.get("weight_scale", 0.333333),
        photometry_enabled=rc_config.get("photometry_enabled", False),
        photo_inner_radius=rc_config.get("photo_inner_radius", 1.25),
        photo_outer_radius=rc_config.get("photo_outer_radius", 5.0),
        photo_centroid_type=get_enum(CentroidType, "photo_centroid_type", CentroidType.BRIGHTEST),
        trim_size=rc_config.get("trim_size", 0.0),
        lss_mapping=rc_config.get("lss_mapping", False),
        gain_delta_start_1=precaldelta1,
        gain_delta_end_1=postcaldelta1,
        gain_delta_start_2=precaldelta2,
        gain_delta_end_2=postcaldelta2,
        skip_time_shift=rc_config.get("skip_time_shift", False),
        skip_background_subtraction=rc_config.get("skip_background_subtraction", False),
        skip_rfi_removal=rc_config.get("skip_rfi_removal", False),
        skip_surface_modeling=rc_config.get("skip_surface_modeling", False),
    )

    # Initialize runner
    runner = RadioCartographerRunner(executable_path)

    try:
        print("Starting execution...")
        runner.run(config, filenames, check=True, capture_output=False)
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
    parser.add_argument("--bin", default="/skynet/radio-cartographer/radio-cartographer",
                        help="Path to the radio-cartographer binary")
    parser.add_argument("--config", default="run_config.toml",
                        help="Path to the configuration file (relative to script or absolute)")
    args = parser.parse_args()

    run_example(args.bin, args.config)
