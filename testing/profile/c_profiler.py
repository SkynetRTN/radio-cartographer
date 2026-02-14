import os
import sys
import subprocess
import argparse
import shutil
import re
from pathlib import Path

# Add the repository root to sys.path to allow importing from pyrc
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
# structure is repo/pyrc/pyrc (package) and repo/pyrc/examples (scripts)
# We want to add repo/pyrc to sys.path so 'import pyrc' works
sys.path.append(str(REPO_ROOT / "pyrc"))
# We also want to import run_example, which is in repo/pyrc/examples
sys.path.append(str(REPO_ROOT / "pyrc" / "examples"))

import run_example


BUILD_DIR = REPO_ROOT / "build"
EXECUTABLE_NAME = "radio-cartographer"
EXECUTABLE_PATH = BUILD_DIR / EXECUTABLE_NAME
PROFILE_RESULTS_DIR = Path(__file__).parent / "results"
GMON_OUT = PROFILE_RESULTS_DIR / "gmon.out"

def clean_build():
    """Removes the build directory."""
    if BUILD_DIR.exists():
        print(f"Cleaning build directory: {BUILD_DIR}")
        shutil.rmtree(BUILD_DIR)

def configure_and_build(rebuild=False):
    """Configures CMake with profiling enabled and builds the project."""
    if rebuild:
        clean_build()

    BUILD_DIR.mkdir(exist_ok=True)
    
    # Configure
    cmd_cmake = [
        "cmake",
        "-DENABLE_PROFILING=ON",
        "-DBUILD_PYTHON_BINDINGS=OFF", 
        ".."
    ]
    print(f"Configuring CMake: {' '.join(cmd_cmake)}")
    subprocess.check_call(cmd_cmake, cwd=BUILD_DIR)

    # Build
    cmd_build = ["cmake", "--build", "."]
    print(f"Building project: {' '.join(cmd_build)}")
    subprocess.check_call(cmd_build, cwd=BUILD_DIR)

def run_profiling_session():
    """Runs the example script to generate profile data."""
    print("Running profiling session...")
    
    # Ensure previous gmon.out is gone from both locations
    cwd_gmon = Path.cwd() / "gmon.out"
    if GMON_OUT.exists():
        GMON_OUT.unlink()
    if cwd_gmon.exists():
        cwd_gmon.unlink()

    # We use the run_example function directly, passing our profiled executable
    # The config is assumed to be the default one in pyrc/examples/run_config.toml
    print(f"Using executable: {EXECUTABLE_PATH}")
    
    try:
        # run_example(executable_path, config_file="run_config.toml")
        # We need to make sure run_example can find run_config.toml
        # run_example uses Path(__file__).parent / config_file
        # so it should be fine if we just call it.
        run_example.run_example(str(EXECUTABLE_PATH))
    except SystemExit as e:
        if e.code != 0:
            print(f"Execution failed with code {e.code}")
            sys.exit(e.code)
    except Exception as e:
        print(f"Execution failed: {e}")
        sys.exit(1)
    
    # The C++ executable generates gmon.out in its working directory (CWD)
    # Move it to the results directory if found in CWD
    if cwd_gmon.exists():
        PROFILE_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        shutil.move(str(cwd_gmon), str(GMON_OUT))
        print(f"Moved gmon.out from {cwd_gmon} to {GMON_OUT}")
    
    if not GMON_OUT.exists():
        print("Error: gmon.out was not generated. Did the executable run successfully?")
        sys.exit(1)
        
    print(f"Profile data generated: {GMON_OUT}")

def parse_gprof_output(output):
    """
    Parses gprof output to categorize functions into "Overhead" vs "Computation".
    """
    lines = output.splitlines()
    
    # Categories patterns
    computation_patterns = [
        r"Processor",
        r"PreProcessor",
        r"Cartographer",
        r"Survey",
        r"Map",
        r"Pixel",
        r"Analysis",
        r"Grid",
        r"Kernel",
        r"Source",
        r"RCR",
        r"Composite"
    ]
    
    overhead_patterns = [
        r"Parser",
        r"OutputFile",
        r"main",
        r"_init",
        r"_fini",
        r"std::",
        r"__",
        r"frame_dummy",
        r"Thread",
        r"mutex",
        r"filesystem"
    ]

    category_times = {"Computation": 0.0, "Overhead": 0.0, "Uncategorized": 0.0}
    
    # Gprof flat profile format usually looks like:
    #  %   cumulative   self              self     total           
    # time   seconds   seconds    calls  ms/call  ms/call  name    
    # 25.00      0.01     0.01        1    10.00    10.00  foo()
    
    # Skip header
    start_processing = False
    
    for line in lines:
        if line.strip() == "":
            continue
        if "time" in line and "seconds" in line and "name" in line:
            start_processing = True
            continue
        if not start_processing:
            continue
        if line.startswith("\x0c"): # Form feed, end of flat profile
            break
            
        parts = line.split()
        if len(parts) < 7:
            continue
            
        try:
            # part[0] is % time, part[2] is self seconds
            self_time = float(parts[2])
            func_name = " ".join(parts[6:])
            
            matched = False
            
            # Check computation first
            for pat in computation_patterns:
                if re.search(pat, func_name):
                    category_times["Computation"] += self_time
                    matched = True
                    break
            
            if not matched:
                for pat in overhead_patterns:
                    if re.search(pat, func_name):
                        category_times["Overhead"] += self_time
                        matched = True
                        break
            
            if not matched:
                category_times["Uncategorized"] += self_time
                # print(f"Uncategorized: {func_name}") # Debug
                
        except ValueError:
            continue

    total_time = sum(category_times.values())
    
    summary_lines = []
    summary_lines.append("Profile Summary")
    summary_lines.append("===============\n")
    
    if total_time == 0:
        summary_lines.append("No significant time recorded in profile (Total time ~ 0s).")
    else:
        summary_lines.append(f"Total Recorded Time: {total_time:.4f} s\n")
        summary_lines.append("Category Breakdown:")
        for cat, time in category_times.items():
            percentage = (time / total_time) * 100
            summary_lines.append(f"  {cat}: {time:.4f} s ({percentage:.2f}%)")
            
    return "\n".join(summary_lines), total_time

def analyze_profile():
    """Runs gprof and analyzes the output."""
    print("Analyzing profile data...")
    
    cmd_gprof = ["gprof", str(EXECUTABLE_PATH), str(GMON_OUT)]
    
    try:
        result = subprocess.run(cmd_gprof, capture_output=True, text=True, check=True)
        gprof_output = result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running gprof: {e}")
        # Gprof might return non-zero if no data
        return

    # Save raw gprof output
    PROFILE_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    import datetime
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    raw_output_file = PROFILE_RESULTS_DIR / f"profile_raw_{timestamp}.txt"
    
    with open(raw_output_file, "w") as f:
        f.write(gprof_output)
    print(f"Raw gprof output saved to: {raw_output_file}")

    # Parse and categorize
    summary, total_time = parse_gprof_output(gprof_output)
    
    summary_file = PROFILE_RESULTS_DIR / f"profile_summary_{timestamp}.txt"
    with open(summary_file, "w") as f:
        f.write(summary)
    print(f"Profile summary saved to: {summary_file}")
    print("\n" + summary)
    
    # Append to a history log for comparison
    history_file = PROFILE_RESULTS_DIR / "history.log"
    with open(history_file, "a") as f:
        f.write(f"{timestamp}: Total={total_time:.4f}s\n")

def main():
    parser = argparse.ArgumentParser(description="Profile Radio Cartographer C++ code.")
    parser.add_argument("--rebuild", action="store_true", help="Clean and rebuild with profiling enabled.")
    args = parser.parse_args()

    # Check if executable exists, if not, force rebuild
    if not EXECUTABLE_PATH.exists():
        print("Executable not found. Forcing build.")
        args.rebuild = True

    if args.rebuild:
        print("Rebuilding project...")
        configure_and_build(rebuild=True)
        
    run_profiling_session()
    analyze_profile()

if __name__ == "__main__":
    main()
