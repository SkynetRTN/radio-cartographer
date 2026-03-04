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

including_frequency_ranges = None
excluding_frequency_ranges = None
including_time_ranges = None
excluding_time_ranges = None
file_path = "../../testing/test_files/0144845_validated.fits"
s = Spectrum(file_path, 0, 1, including_frequency_ranges, excluding_frequency_ranges, including_time_ranges, excluding_time_ranges)
spectrum = s.spectrum()

freqs = spectrum[0]
intensities = spectrum[1]

# C++ Background Algorithm assumes angular distance (frequency) is monotonically increasing
sort_idx = np.argsort(freqs)
freqs = freqs[sort_idx]
intensities = intensities[sort_idx]

# Estimate scatter as standard deviation of differences divided by sqrt(2)
# and pick a simple baseline scale (in frequency units)
scatter = np.std(np.diff(intensities)) / np.sqrt(2.0)
baseline_scale = (freqs[-1] - freqs[0]) * 0.05 # 5% of total bandwidth
if baseline_scale < 0:
    baseline_scale = -baseline_scale

input_data = {
    "frequencies": freqs.tolist(),
    "intensities": intensities.tolist(),
    "scatter": scatter,
    "baseline": baseline_scale
}

in_json = "temp_in.json"
out_json = "temp_out.json"

with open(in_json, "w") as f:
    json.dump(input_data, f)

executable_path = "/skynet/radio-cartographer/cmake-build-debug-docker/spectral-cleaner"

subprocess.run(["bash", "../../docker-shell.sh", executable_path, "pyrc/pyrc/" + in_json, "pyrc/pyrc/" + out_json], check=True)

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

# Clean up
if os.path.exists(in_json):
    os.remove(in_json)
if os.path.exists(out_json):
    os.remove(out_json)