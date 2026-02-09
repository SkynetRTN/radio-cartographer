import sys
import os
import subprocess
import time

# Add pyrc to path
sys.path.append("/skynet/radio-cartographer/radio-cartographer/pyrc")
from pyrc.parameters import RadioCartographerConfig, Channel, Receiver, CalibrationMethod, CoordinateSystem, TimeShiftMode, CentroidType

def run_test(output_suffix, skips=None):
    if skips is None:
        skips = {}

    config = RadioCartographerConfig()
    config.channel = Channel.COMPOSITE
    config.calibration_method = CalibrationMethod.INTERPOLATED
    config.coordinate_system = CoordinateSystem.EQUATORIAL
    config.trim_size = 0.0
    
    config.receiver = Receiver.HI
    config.min_freq = 1350.0
    config.max_freq = 1435.0
    config.exclusion_bands = []

    config.rfi_scale = 0.8
    config.weight_scale = 1.0
    config.bg_scale = 0.7
    
    config.photometry_enabled = False
    config.m10_plus_processing = False # Default off
    config.lss_mapping = False
    config.generate_raw_map = False
    
    config.time_shift_mode = TimeShiftMode.AUTO
    config.time_shift_value = 0.0

    # Apply skips
    if skips.get("skip_time_shift"):
        config.skip_time_shift = True
    if skips.get("skip_background_subtraction"):
        config.skip_background_subtraction = True
    if skips.get("skip_rfi_removal"):
        config.skip_rfi_removal = True
    if skips.get("skip_surface_modeling"):
        config.skip_surface_modeling = True
    
    # Generate args
    args = config.to_args()
    
    # Prepend filename
    input_file = "/skynet/test_standards_static/gbo20_L_lores_daisy_cyga.fits"
    exec_path = "/skynet/radio-cartographer/radio-cartographer"
    
    cmd = [exec_path, input_file] + args
    
    print(f"Running command: {' '.join(cmd)}")
    
    start_time = time.time()
    try:
        # Run process
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("Stdout:", result.stdout[-500:]) # Last 500 chars
    except subprocess.CalledProcessError as e:
        print("Error running command:")
        print(e.stderr)
        return False
        
    duration = time.time() - start_time
    print(f"Success! Duration: {duration:.2f}s")
    
    # Check for output file
    # Base name handling in logic: strip dir, strip ext, add _processed.fits
    basename = os.path.basename(input_file).rsplit('.', 1)[0]
    output_file = f"{basename}_processed.fits"
    
    if os.path.exists(output_file):
        print(f"Output file created: {output_file}")
    else:
        print("Output file NOT found.")
        
    return True

if __name__ == "__main__":
    print("Running baseline test (No Skips)...")
    if not run_test("baseline"):
        print("Baseline failed!")
        exit(1)

    print("\nRunning Skip TS...")
    run_test("skip_ts", {"skip_time_shift": True})

    print("\nRunning Skip BG...")
    run_test("skip_bg", {"skip_background_subtraction": True})
    
    print("\nRunning Skip RFI...")
    run_test("skip_rfi", {"skip_rfi_removal": True})
    
    print("\nRunning Skip Surface Modeling...")
    run_test("skip_surface", {"skip_surface_modeling": True})
    
    print("\nRunning Skip All...")
    run_test("skip_all", {
        "skip_time_shift": True,
        "skip_background_subtraction": True, 
        "skip_rfi_removal": True,
        "skip_surface_modeling": True
    })
