"""rcpy — native Python port of Radio Cartographer.

Implements the single-dish radio mapping pipeline from Martin et al. 2018,
operating on the same time-ordered measurements as the C++ implementation
but expressed via numpy/scipy and the already-ported RCR2 library.

The public surface mirrors the C++ stage names so the port can be validated
stage-by-stage against the C++ binary used as a parity oracle.

Pipeline stages (see ../PIPELINE_STAGES.txt for the full spec):

    0.  rcpy.io.read_sdfits         — read SDFITS into a Survey
    1.  rcpy.coordinates.project    — sinusoidal cos-Dec projection
    2.  rcpy.calibration.gain       — gain calibration (§3.1)
    3.  rcpy.noise.measure_1d       — point-to-point along-scan noise (§3.2)
    4.  rcpy.background.subtract    — 1D background subtraction (§3.3)
    5.  rcpy.timeshift.correct      — time-delay correction (§3.4)
    6.  rcpy.noise.measure_2d       — across-scan noise (§3.5)
    7.  rcpy.rfi.subtract           — 2D RFI subtraction (§3.6)
    8.  rcpy.thetagap.compute       — bubble-blowing density (§3.7)
    9.  rcpy.surface.regrid         — local-polynomial regridding (§3.7)
    10. rcpy.photometry.aperture    — aperture photometry (§4)

    rcpy.pipeline.run               — orchestrates all of the above
"""

from rcpy.types import (
    Survey,
    Scan,
    Composite,
    PartitionSet,
    Channel,
    CalMethod,
    Coordinates,
    MapType,
    Telescope,
)

__version__ = "0.1.0.dev0"

__all__ = [
    "Survey",
    "Scan",
    "Composite",
    "PartitionSet",
    "Channel",
    "CalMethod",
    "Coordinates",
    "MapType",
    "Telescope",
]
