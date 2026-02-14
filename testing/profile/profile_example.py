import sys
import os
import cProfile
import pstats
import contextlib
import tempfile
from pathlib import Path
from datetime import datetime

# Adjust path to find run_example
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
EXAMPLES_DIR = REPO_ROOT / "pyrc" / "examples"
sys.path.append(str(EXAMPLES_DIR))

import run_example

def main():
    # Determine executable path (try local project root first)
    executable = REPO_ROOT / "radio-cartographer"
    # if not executable.exists():
    #      executable = REPO_ROOT / "build" / "radio-cartographer"
    # if not executable.exists():
    #      executable = Path("/skynet/radio-cartographer/radio-cartographer")

    config_file = "run_config.toml" # Rel to run_example.py
    
    # Output file for profile results
    output_date_time_str = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_file = Path(__file__).parent / "results" / f"profile_{output_date_time_str}.txt"

    # Use a temporary directory to avoid cluttering the workspace with output files
    with tempfile.TemporaryDirectory() as temp_dir:
        original_cwd = os.getcwd()
        os.chdir(temp_dir)
        try:
            print(f"Profiling run_example with executable: {executable}")
            print(f"Running in temporary directory: {temp_dir}")
            
            with contextlib.redirect_stdout(open(os.devnull, 'w')):
                # We might also want to redirect stderr if needed
                # with contextlib.redirect_stderr(open(os.devnull, 'w')):
                profiler = cProfile.Profile()
                profiler.enable()
                try:
                    # run_example expects executable_path and config_file
                    # The config path resolution in run_example.py handles relative paths correctly
                    run_example.run_example(str(executable), config_file)
                except SystemExit:
                    pass
                profiler.disable()
        finally:
            os.chdir(original_cwd)

    print(f"Writing profile results to {results_file}")
    with open(results_file, 'w') as f:
        stats = pstats.Stats(profiler, stream=f)
        stats.strip_dirs()
        stats.sort_stats('cumulative')
        stats.print_stats()

if __name__ == "__main__":
    main()
