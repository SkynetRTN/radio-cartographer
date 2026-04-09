import os
import glob
import json
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gain_calibration_validation import Validation
from spectrum import Spectrum
import stats


REPO_ROOT = Path("/skynet/radio-cartographer")
SCRIPT_DIR = REPO_ROOT / "pyrc" / "pyrc"
TEST_FILES_DIR = SCRIPT_DIR / "spectral_testing"
EXECUTABLE_PATH = REPO_ROOT / "build" / "debug" / "spectral-cleaner"


def process_file(fits_file):
    fits_file = Path(fits_file)
    print(f"Processing {fits_file}...")

    # 1. Validation
    try:
        validator = Validation(str(fits_file))
        validated_filepath = Path(validator.validate())
        print(f"Validated file created: {validated_filepath}")
    except Exception as e:
        print(f"Failed to validate {fits_file}: {e}")
        return

    # 2. Spectral RFI Removal
    including_frequency_ranges = None
    excluding_frequency_ranges = None
    including_time_ranges = None
    excluding_time_ranges = None

    try:
        s = Spectrum(
            str(validated_filepath),
            0,
            1,
            including_frequency_ranges,
            excluding_frequency_ranges,
            including_time_ranges,
            excluding_time_ranges,
        )
        spectrum_data = s.spectrum()
    except Exception as e:
        print(f"Failed to generate spectrum for {validated_filepath}: {e}")
        if validated_filepath.exists():
            validated_filepath.unlink()
        return

    freqs = spectrum_data[0]
    intensities = spectrum_data[1]

    # C++ Background Algorithm assumes angular distance (frequency) is monotonically increasing
    sort_idx = np.argsort(freqs)
    freqs = freqs[sort_idx]
    intensities = intensities[sort_idx]

    scatter = np.std(np.diff(intensities)) / np.sqrt(2.0)
    baseline_scale = abs((freqs[-1] - freqs[0]) * 0.05)

    input_data = {
        "frequencies": freqs.tolist(),
        "intensities": intensities.tolist(),
        "scatter": scatter,
        "baseline": baseline_scale,
    }

    # Use per-file temp JSON names so multiple runs/files do not collide
    temp_stem = fits_file.stem
    in_json = SCRIPT_DIR / f"{temp_stem}_temp_in.json"
    out_json = SCRIPT_DIR / f"{temp_stem}_temp_out.json"

    with open(in_json, "w") as f:
        json.dump(input_data, f)

    try:
        subprocess.run(
            [
                str(EXECUTABLE_PATH),
                str(in_json),
                str(out_json),
            ],
            check=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Failed to run spectral cleaner: {e}")
        if validated_filepath.exists():
            print("removing validated file")
            try:
                validated_filepath.unlink()
            except OSError:
                pass
        if in_json.exists():
            in_json.unlink()
        return

    if not out_json.exists():
        print(f"Output json {out_json} not found. Skipping plot.")
        if validated_filepath.exists():
            validated_filepath.unlink()
        if in_json.exists():
            in_json.unlink()
        return

    with open(out_json, "r") as f:
        out_data = json.load(f)

    modeled_spectrum = np.array(out_data["modeled_spectrum"])
    diff = intensities - modeled_spectrum

    try:
        mask, mu, gamma = stats.chauvenet(diff, clip_lo = False)
    except Exception as e:
        print(f"Error performing RCR chauvenet: {e}")
        mask = np.zeros_like(diff, dtype=bool)

    # 3. Generating the plots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Top Plot - Overlay
    ax1.plot(freqs, intensities, color="black", label="Raw Spectrum", linewidth=0.5)
    ax1.plot(freqs, modeled_spectrum, color="red", label="Modeled Background")
    ax1.set_ylabel("Signal [Dimensionless Units]")
    ax1.set_title("Spectrum with Local Background Model")
    ax1.legend()

    # Bottom Plot - Difference
    ax2.plot(freqs, diff, color="black", label="Difference (Raw - Model)", linewidth=0.5)

    rejected_freqs = freqs[mask]
    rejected_diffs = diff[mask]
    if len(rejected_freqs) > 0:
        ax2.scatter(
            rejected_freqs,
            rejected_diffs,
            color="red",
            label=f"Rejected points ({len(rejected_freqs)})",
            zorder=5,
        )

    ax2.set_xlabel("Frequency (MHz)")
    ax2.set_ylabel("Difference")
    ax2.set_title("Difference with Rejected Points Highlighted")
    ax2.legend()

    plt.tight_layout()

    root_name = fits_file.stem
    plot_name = f"{root_name}_plotting.png"

    # Save plots to spectral_testing/plots
    plot_dir = fits_file.parent / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    plot_path = plot_dir / plot_name
    plt.savefig(plot_path)
    plt.close(fig)
    print(f"Saved plot: {plot_path}")

    # 4. Clean up
    if in_json.exists():
        in_json.unlink()
    if out_json.exists():
        out_json.unlink()

    if validated_filepath.exists():
        try:
            validated_filepath.unlink()
            print(f"Cleaned up validated file: {validated_filepath}")
        except OSError:
            pass


def main():
    if not TEST_FILES_DIR.is_dir():
        print(f"Directory {TEST_FILES_DIR} does not exist. Creating it.")
        try:
            TEST_FILES_DIR.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            print(f"Could not create directory {TEST_FILES_DIR}: {e}")

    fits_files = glob.glob(str(TEST_FILES_DIR / "*.fits"))
    # Filter out files that already have "_validated"
    fits_files = [f for f in fits_files if "_validated" not in f]

    print(f"Found {len(fits_files)} fits files in {TEST_FILES_DIR}")

    for fits_file in fits_files:
        process_file(fits_file)


if __name__ == "__main__":
    main()
