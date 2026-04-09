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

including_frequency_ranges = None
excluding_frequency_ranges = None
including_time_ranges = None
excluding_time_ranges = None

repo_root = Path("/skynet/radio-cartographer")
script_dir = repo_root / "pyrc" / "pyrc"

file_path = script_dir / "spectral_testing" / "0149927_validated.fits"
s = Spectrum(str(file_path), 0, 1,
             including_frequency_ranges, excluding_frequency_ranges,
             including_time_ranges, excluding_time_ranges)
spectrum = s.spectrum()

freqs = spectrum[0]
intensities = spectrum[1]

sort_idx = np.argsort(freqs)
freqs = freqs[sort_idx]
intensities = intensities[sort_idx]

scatter = np.std(np.diff(intensities)) / np.sqrt(2.0)
baseline_scale = abs((freqs[-1] - freqs[0]) * 0.05)

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

modeled_spectrum = np.array(out_data["modeled_spectrum"])

plt.figure(figsize=(10, 6))
plt.plot(freqs, intensities, color="black", label="Raw Spectrum", linewidth=0.5)
plt.plot(freqs, modeled_spectrum, color="red", label="Modeled Background")
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal [Dimensionless Units]")
plt.title("Spectrum with Local Background Model")
plt.legend()
plt.tight_layout()
plt.savefig("testspectrum.png")
print("saved spectrum")

if in_json.exists():
    in_json.unlink()
if out_json.exists():
    out_json.unlink()