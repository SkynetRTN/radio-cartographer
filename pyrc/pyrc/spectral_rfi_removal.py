import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.stats import linregress
import matplotlib.pyplot as plt
import rcr
import utils
from spectrum import Spectrum
import subprocess
import json
import os
from pathlib import Path
import stats
import sys
from gain_calibration_validation import Validation

including_frequency_ranges = None
excluding_frequency_ranges = None
including_time_ranges = None
excluding_time_ranges = None

repo_root = Path("/skynet/radio-cartographer")
script_dir = repo_root / "pyrc" / "pyrc"

file_path = script_dir / "spectral_testing" / "0152191.fits"
#file_path = "testing/test_files/Skynet_56870_survey_55_ra2_9650_10467.cybmerge2_interpolated.fits"
print(f"Processing {file_path}...")

try:
    validator = Validation(str(file_path))
    validated_filepath = Path(validator.validate())
    print(f"Validated file created: {validated_filepath}")
except Exception as e:
    print(f"Failed to validate {file_path}: {e}")
    sys.exit(1)

s = Spectrum(str(validated_filepath), 0, 0,
             including_frequency_ranges, excluding_frequency_ranges,
             including_time_ranges, excluding_time_ranges)
spectrum = s.spectrum()

freqs = spectrum[0]
intensities = spectrum[1]

sort_idx = np.argsort(freqs)
freqs = freqs[sort_idx]
intensities = intensities[sort_idx]

scatter = np.std(np.diff(intensities)) / np.sqrt(2.0)
#baseline_scale = abs((freqs[-1] - freqs[0]) * 0.05) #lband scale
baseline_scale = abs((freqs[-1] - freqs[0]) * 0.02)
input_data = {
    "frequencies": freqs.tolist(),
    "intensities": intensities.tolist(),
    "scatter": scatter,
    "baseline": baseline_scale
}

in_json = script_dir / "temp_in.json"
out_json = script_dir / "temp_out.json"
executable_path = repo_root / "build" / "debug" / "spectral-cleaner"

with open(in_json, "w") as f:
    json.dump(input_data, f)

subprocess.run(
    [
        str(executable_path),
        str(in_json),
        str(out_json),
    ],
    check=True
)

with open(out_json, "r") as f:
    out_data = json.load(f)
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
modeled_spectrum = np.array(out_data["modeled_spectrum"])
#raw spectrum
plt.figure(figsize = (10, 6))
plt.plot(freqs, intensities, color="black", label="Raw Spectrum", linewidth=0.5)
plt.savefig("RawSpectrumxband")
plt.close()


plt.figure(figsize = (15, 11))
plt.plot(freqs, intensities, color="black", label="Raw Spectrum", linewidth = 0.7)
plt.plot(freqs, modeled_spectrum, color="red", label="Modeled Background")
plt.xlabel("Frequency (MHz)", fontsize = 25)
plt.ylabel("Signal [Dimensionless Units]", fontsize = 25)
plt.title("Spectrum with Background Model", fontsize = 30)
plt.legend(fontsize = 25)
plt.tight_layout()
plt.savefig("testspectrumxband.png")
plt.close()
print("saved spectrum")

diff = intensities - modeled_spectrum

try:
    mask, mu, gamma = stats.chauvenet(diff, clip_lo = False)
except Exception as e:
    print(f"Error performing RCR chauvenet: {e}")
    mask = np.zeros_like(diff, dtype=bool)

rejected_freqs = freqs[mask]
rejected_diffs = diff[mask]
plt.figure(figsize = (15, 9))
plt.plot(freqs, diff, color="black", label="Difference (Raw - Model)", )
if len(rejected_freqs) > 0:
    plt.scatter(
        rejected_freqs,
        rejected_diffs,
        color="red",
        label=f"Rejected points ({len(rejected_freqs)})",
        zorder=5,
        s = 55
    )
print(f"rejected freqs: {rejected_freqs}")
print(len(rejected_freqs))
plt.xlabel("Frequency (MHz)", fontsize = 25)
plt.ylabel("Difference", fontsize = 25)
plt.title("Difference with Rejected Points Highlighted", fontsize = 30)
plt.legend(fontsize = 25)
plt.tight_layout()
plt.savefig("xband.png")

plt.close()
if in_json.exists():
    in_json.unlink()
if out_json.exists():
    out_json.unlink()
if 'validated_filepath' in locals() and validated_filepath.exists():
    try:
        validated_filepath.unlink()
        print(f"Cleaned up validated file: {validated_filepath}")
    except OSError:
        pass