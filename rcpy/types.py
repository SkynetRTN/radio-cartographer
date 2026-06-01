"""Data model for the rcpy pipeline.

The C++ uses a heavyweight Scan class with ~200 getters and ~10 parallel
coordinate variants accumulated by in-place mutation. We collapse that into
a single Scan dataclass holding NumPy arrays, with stages adding fields
rather than overwriting them. A Survey is a list of Scans plus metadata; a
Composite bundles multiple Surveys for joint RFI subtraction and regridding.

Field naming convention: snake_case scalar-or-array attributes. Per-sample
arrays are 1D, length N_samples_in_scan. Optional fields are None until the
relevant stage populates them — making it explicit which fields each stage
consumes vs produces.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, List
import numpy as np


# ---------------------------------------------------------------------------
# Enums (mirror C++ enum names from include/Structures.h)
# ---------------------------------------------------------------------------

class Channel(Enum):
    LEFT = "left"
    RIGHT = "right"
    COMPOSITE = "composite"


class CalMethod(Enum):
    PRE = "pre"
    POST = "post"
    INTERPOLATED = "interpolated"
    CONSTANT = "constant"  # mean(start, end) — guards against spurious cal drift
    NONE = "none"


class Coordinates(Enum):
    EQUATORIAL = "equatorial"
    GALACTIC = "galactic"


class MapType(Enum):
    RASTER = "raster"
    NODDING = "nodding"
    DAISY = "daisy"


class Telescope(Enum):
    FOURTY_FOOT = "40ft"
    TWENTY_METER = "20m"
    GBT = "gbt"


class SurfaceModel(Enum):
    """Polynomial order fallback chain for surface modeling (§3.7)."""
    POLY_3RD = "poly_3rd"   # 10 coeffs, prefer this
    POLY_2ND = "poly_2nd"   # 6 coeffs, if insufficient data for 3rd
    POLY_1ST = "poly_1st"   # 3 coeffs, last resort
    EXCISED = "excised"     # not enough data even for a plane


# ---------------------------------------------------------------------------
# Geometry / partition
# ---------------------------------------------------------------------------

@dataclass
class PartitionSet:
    """Geometric envelope of a Survey (or Composite). Mirrors C++
    PartitionSet from include/Structures.h:54.

    All coordinates are in the projected/internal frame (post cos-Dec)
    except center_ra_deg / center_dec_deg which are the un-projected
    sky coords used to invert the projection for FITS output.
    """
    map_type: MapType

    # Bounds in projected frame
    min_ra: float = 0.0
    max_ra: float = 0.0
    min_dec: float = 0.0
    max_dec: float = 0.0

    # Sky-frame center (degrees) — used to invert the cos-Dec projection
    center_ra_deg: float = 0.0
    center_dec_deg: float = 0.0

    # Median in internal coords (from RCR on all sample positions)
    median_ra: float = 0.0
    median_dec: float = 0.0

    # Edge geometry. For rasters/noddings these are four linear edge
    # polynomials (each is a list of polynomial coefficients).
    # For daisies, edge_radius is the only edge.
    edge_one: Optional[np.ndarray] = None
    edge_two: Optional[np.ndarray] = None
    edge_three: Optional[np.ndarray] = None
    edge_four: Optional[np.ndarray] = None
    edge_radius: float = 0.0

    tracking: bool = False
    trim_size: float = 0.0


# ---------------------------------------------------------------------------
# Scan
# ---------------------------------------------------------------------------

@dataclass
class Scan:
    """One slew of the telescope between direction reversals (or one daisy
    petal, or one nodding leg).

    Mirrors C++ Scan class from include/Scan.h. Where C++ uses ~10 parallel
    coordinate variants accumulated by in-place mutation, we keep raw
    coordinates in `ra`, `dec` and add optional fields as the pipeline
    progresses (`ra_projected`, `ra_timeshifted`, etc.).

    All array fields are 1D of length N (samples in this scan), except
    where noted. Optional arrays are None until the relevant stage runs.
    """
    # ----- core (always present after ingest) -----
    scan_index: int            # position of this scan in the parent survey
    time: np.ndarray           # UT seconds, shape (N,)
    ra: np.ndarray             # right ascension, degrees, shape (N,)
    dec: np.ndarray            # declination, degrees, shape (N,)
    elevation: np.ndarray      # elevation angle, degrees, shape (N,)
    azimuth: np.ndarray        # azimuth, degrees, shape (N,)

    # Raw per-channel flux + integration weights
    flux_l: np.ndarray         # left polarization, shape (N,)
    flux_r: np.ndarray         # right polarization, shape (N,)
    dumps: np.ndarray          # integration count per sample, shape (N,)
    cal_state: np.ndarray      # 0=off, 1=on, etc.; calibration-diode flag

    # Composite (L + R, post-calibration). Populated by Stage 2.
    flux_composite: Optional[np.ndarray] = None

    # ----- projected coordinates (set by Stage 1) -----
    ra_proj: Optional[np.ndarray] = None
    dec_proj: Optional[np.ndarray] = None

    # ----- after time-delay correction (Stage 5) -----
    ra_ts: Optional[np.ndarray] = None
    dec_ts: Optional[np.ndarray] = None

    # ----- calibrated flux (the "working channel" in C++) -----
    flux: Optional[np.ndarray] = None

    # ----- noise + background -----
    noise_1d: Optional[np.ndarray] = None      # per-sample sigma (Stage 3)
    noise_2d: Optional[np.ndarray] = None      # per-sample sigma (Stage 6)
    background: Optional[np.ndarray] = None    # subtracted background (Stage 4)
    flux_bg: Optional[np.ndarray] = None       # flux - background

    # ----- RFI / theta-gap -----
    flux_rfi: Optional[np.ndarray] = None      # RFI-subtracted flux (Stage 7)
    rfi_keep_mask: Optional[np.ndarray] = None # bool, False = excised
    theta_gap: Optional[np.ndarray] = None     # per-sample density (Stage 8)
    theta_corr: Optional[np.ndarray] = None    # correlation length for photometry
    gm_weight: Optional[np.ndarray] = None     # GMWeight (sum LMWeight/rfiCount over anchors)
    edge_trim_mask: Optional[np.ndarray] = None # bool, False = trimmed turning-edge sample
    min_rcr_theta_gap: Optional[np.ndarray] = None  # per-sample minimum-distance filter
                                                     # for theta_gap neighbour collection
                                                     # (C++ Scan::minRCRThetaGap)

    # ----- along-scan angular distance (derived once after projection) -----
    ang_dist: Optional[np.ndarray] = None

    # ----- misc flags -----
    edge_flag: Optional[np.ndarray] = None     # near end-of-scan
    turning_point_flag: Optional[np.ndarray] = None
    rejected_flag: Optional[np.ndarray] = None # per-stage exclusion

    # ----- meta -----
    scan_in_ra: bool = True  # True if scan progresses primarily in RA

    @property
    def size(self) -> int:
        return int(self.time.size)

    def working_coords(self) -> tuple[np.ndarray, np.ndarray]:
        """Return the most-processed (ra, dec) pair available.

        Precedence: time-shifted > projected > raw. Most pipeline stages
        should call this rather than reaching for specific fields, so the
        port doesn't bake ordering assumptions into every site.
        """
        if self.ra_ts is not None and self.dec_ts is not None:
            return self.ra_ts, self.dec_ts
        if self.ra_proj is not None and self.dec_proj is not None:
            return self.ra_proj, self.dec_proj
        return self.ra, self.dec

    def working_flux(self) -> np.ndarray:
        """Return the most-processed flux array available.

        Precedence: rfi-subtracted > bg-subtracted > calibrated > L-channel.
        """
        for candidate in (self.flux_rfi, self.flux_bg, self.flux, self.flux_l):
            if candidate is not None:
                return candidate
        raise ValueError(f"Scan {self.scan_index} has no flux data")


# ---------------------------------------------------------------------------
# Survey
# ---------------------------------------------------------------------------

@dataclass
class Survey:
    """One observation. Owns scans + metadata + computed partition.

    Mirrors C++ Survey from include/Survey.h. A Survey is a single SDFITS
    file (or one logical observation). Multiple Surveys can be bundled
    into a Composite for joint processing.
    """
    scans: List[Scan]

    # Telescope / observation metadata
    telescope: Telescope
    frequency_ghz: float
    mjd: float
    map_type: MapType
    channel: Channel = Channel.COMPOSITE
    cal_method: CalMethod = CalMethod.INTERPOLATED
    p_coordinate: Coordinates = Coordinates.EQUATORIAL  # processing frame
    m_coordinate: Coordinates = Coordinates.EQUATORIAL  # output frame

    # Beam size (degrees). Set per (telescope, frequency) at ingest time.
    psf_fwhm: float = 0.0

    # Tracking (object-centered for daisies on moving objects)
    tracking: bool = False

    # Commanded source position from the FITS primary header (degrees).
    # For moving-object observations where the science samples don't
    # cover the source (telescope parks at the source then scans away),
    # this is the only way to anchor the projection.
    target_ra_deg: float = float("nan")
    target_dec_deg: float = float("nan")

    # Geometry, populated by Stage 1
    partition_sss: Optional[PartitionSet] = None
    partition_lss: Optional[PartitionSet] = None

    # Calibration samples (set by io.read_sdfits, consumed by calibration.gain).
    # These are the SWPVALID=0 samples — diode-on (CALSTATE=1) and turnaround
    # diode-off (CALSTATE=0). Skynet uses turnaround samples as the baseline.
    cal_time: Optional[np.ndarray] = None       # UTSECS of cal samples
    cal_state: Optional[np.ndarray] = None      # CALSTATE per cal sample
    cal_dumps: Optional[np.ndarray] = None      # NSAMPS per cal sample
    cal_flux_l: Optional[np.ndarray] = None     # band-integrated L flux
    cal_flux_r: Optional[np.ndarray] = None     # band-integrated R flux

    # Calibration deltas measured by Stage 2 — one per channel
    gain_delta_l: Optional[np.ndarray] = None  # shape (n_intervals,)
    gain_delta_r: Optional[np.ndarray] = None
    gain_delta_times: Optional[np.ndarray] = None  # midpoint times of each interval

    # Per-survey noise model summary (linear fit across scan number)
    noise_1d_slope: float = 0.0
    noise_1d_intercept: float = 0.0
    noise_2d_slope: float = 0.0
    noise_2d_intercept: float = 0.0

    # Time-delay correction value (seconds)
    time_shift: float = 0.0

    # Smallest along-scan inter-sample distance (after RCR rejection).
    # C++ Survey::setStandardThetaGap stores this as `minGapThreshold`.
    # Used as the per-sample lower bound on RCRMinThetaGap, which in
    # turn filters too-close neighbours out of the theta_gap calculation
    # (ProcessorThetaGap.cpp:485). Defaults to 0 (no filtering) until set.
    min_gap_threshold: float = 0.0

    # Path to source file (for provenance)
    source_path: str = ""

    survey_number: int = 0

    @property
    def n_scans(self) -> int:
        return len(self.scans)

    @property
    def n_samples(self) -> int:
        return sum(s.size for s in self.scans)


# ---------------------------------------------------------------------------
# Composite
# ---------------------------------------------------------------------------

@dataclass
class Composite:
    """One or more Surveys appended for joint downstream processing.

    After Stage 7 (per-survey edge calculation), surveys get bundled here.
    From this point on, RFI subtraction, theta-gap, and surface modeling
    operate on the combined set.
    """
    surveys: List[Survey]
    partition_sss: PartitionSet
    partition_lss: Optional[PartitionSet] = None

    @property
    def all_scans(self) -> List[Scan]:
        """Concatenate scans across surveys (preserving original order)."""
        out: list[Scan] = []
        for s in self.surveys:
            out.extend(s.scans)
        return out

    @property
    def psf_fwhm(self) -> float:
        # All surveys in a composite should share a beam size; if they
        # don't, take the first (and we should log a warning in pipeline).
        return self.surveys[0].psf_fwhm if self.surveys else 0.0

    @property
    def map_type(self) -> MapType:
        return self.surveys[0].map_type if self.surveys else MapType.RASTER


# ---------------------------------------------------------------------------
# Map (output container)
# ---------------------------------------------------------------------------

@dataclass
class Map:
    """Output 2D grids from Stage 9 surface modeling.

    Mirrors C++ Map from include/Map.h. Each field is a 2D ndarray of
    shape (n_dec_pixels, n_ra_pixels) — origin at (min_dec, min_ra).
    """
    # Coordinate metadata
    min_ra: float
    max_ra: float
    min_dec: float
    max_dec: float
    center_ra_deg: float
    center_dec_deg: float
    resolution: float           # pixel size in same units as coords

    # Final flux map (the main image)
    flux: np.ndarray            # shape (ny, nx)

    # Per-pixel diagnostics
    weight: np.ndarray          # weighted sample count per pixel
    weight_corr: np.ndarray     # correlation-corrected weight
    scale: np.ndarray           # theta_w used at this pixel
    correlation: np.ndarray     # combined theta_corr / theta_w
    path: np.ndarray            # scan-pattern visualization

    # Optional raw map (no cleaning applied)
    raw_flux: Optional[np.ndarray] = None

    @property
    def shape(self) -> tuple[int, int]:
        return self.flux.shape
