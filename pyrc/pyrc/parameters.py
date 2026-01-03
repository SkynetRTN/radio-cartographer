from dataclasses import dataclass, field
from enum import Enum
from typing import List, Union

class Channel(Enum):
    LEFT = "left"
    RIGHT = "right"
    COMPOSITE = "composite"

class CalibrationMethod(Enum):
    PRE = "pre"
    POST = "post"
    INTERPOLATED = "interpolated"
    NONE = "none"

class CoordinateSystem(Enum):
    EQUATORIAL = "equatorial"
    GALACTIC = "galactic"

class Receiver(Enum):
    HI = 1
    LO = 2

class CentroidType(Enum):
    CENTER = "center"
    BRIGHTEST = "brightest"
    COORDINATES = "coordinates"

class TimeShiftMode(Enum):
    AUTO = "auto"
    CUSTOM = "custom"
    OFF = "off"

@dataclass
class RadioCartographerConfig:
    """
    Configuration for a Radio Cartographer job.
    Mirrors the input requirements of the C++ executable.
    """
    # Survey Parameters
    channel: Channel = Channel.COMPOSITE
    calibration_method: CalibrationMethod = CalibrationMethod.INTERPOLATED
    coordinate_system: CoordinateSystem = CoordinateSystem.EQUATORIAL
    trim_size: float = 0.0

    # Spectral Parameters
    receiver: Receiver = Receiver.HI
    min_freq: float = 1350.0
    max_freq: float = 1435.0
    exclusion_bands: List[float] = field(default_factory=list)

    # Processor / Map Parameters
    rfi_scale: float = 0.8  # Corresponds to argv[12]
    weight_scale: float = 1.0  # Corresponds to argv[14]
    bg_scale: float = 0.7  # Corresponds to argv[11]
    
    # Flags
    photometry_enabled: bool = False  # argv[15]
    m10_plus_processing: bool = False  # argv[13]
    lss_mapping: bool = False  # argv[21]
    generate_raw_map: bool = False  # argv[8]

    # Time Shift
    time_shift_mode: TimeShiftMode = TimeShiftMode.AUTO
    time_shift_value: float = 0.0

    # Skip Flags
    skip_time_shift: bool = False
    skip_background_subtraction: bool = False
    skip_rfi_removal: bool = False
    skip_surface_modeling: bool = False

    # Photometry Parameters
    photo_inner_radius: float = 1.25
    photo_outer_radius: float = 5.0
    photo_centroid_type: CentroidType = CentroidType.BRIGHTEST

    def to_args(self) -> List[str]:
        """
        Generates the list of CLI arguments expected by the C++ executable.
        Note: The filename (argv[1]) is NOT included here and must be prepended.
        """
        args = ["" for _ in range(22)] # Allocate space up to argv[21] + 1 for safety? 
        # Actually Source.cpp accesses specific indices.
        # Let's map them explicitly to a list of the correct size/content.
        # Indices are 1-based in Source.cpp (argv[1] is filename).
        # We will return list where index 0 maps to argv[2] (since argv[1] is separate).
        # WAIT. It's safer to return the WHOLE list excluding argv[0] (executable) and argv[1] (filename).
        # But Source.cpp indices are fixed.
        # argv[2] is the first config arg.
        
        # We need to constructing a list of strings such that:
        # result[0] will be argv[2]
        # result[1] will be argv[3]
        # ...
        
        # Let's build a dict map first for clarity, then convert to list.
        # Only filling indices 2 through 21, plus extra for exclusion bands.
        
        arg_map = {}
        
        # Survey
        arg_map[2] = self.channel.value
        # argv[3] is Receiver. 2 -> LO, else HI.
        arg_map[3] = "2" if self.receiver == Receiver.LO else "1"
        arg_map[4] = self.calibration_method.value
        arg_map[5] = self.coordinate_system.value
        
        # Processor
        # argv[6]: Time Shift String
        if self.time_shift_mode == TimeShiftMode.AUTO:
            arg_map[6] = "auto"
        elif self.time_shift_mode == TimeShiftMode.OFF:
            arg_map[6] = "off"
        else:
            arg_map[6] = "custom" # Or maybe just the number?
            # Source.cpp: procParams.timeShift = argv[6];
            # if (procParams.timeShift == "auto") headerInfo["RCTS"] = "2"
            # else headerInfo["RCTS"] = "0"
            # It seems Source.cpp stores the string. 
            # `input.py` does: if is_number(self.time_delay) -> CUSTOM.value (numeric?) No `CUSTOM` is a value.
            
            # Re-checking input.py:
            # timeshift_option = TimeShift.CUSTOM.value if is_number(self.time_delay) else self.time_delay.value
            # TimeShift.CUSTOM = "custom"
            # So if it's a number, it passes "custom".
            
        arg_map[7] = str(self.time_shift_value)
        arg_map[8] = "1" if self.generate_raw_map else "0"
        
        # Spectral
        arg_map[9] = str(self.min_freq)
        arg_map[10] = str(self.max_freq)
        
        # Processor / Map
        arg_map[11] = str(self.bg_scale)
        arg_map[12] = str(self.rfi_scale)
        arg_map[13] = "1" if self.m10_plus_processing else "0"
        arg_map[14] = str(self.weight_scale)
        
        # Flags / Photo
        arg_map[15] = "1" if self.photometry_enabled else "0"
        arg_map[16] = str(self.photo_inner_radius)
        arg_map[17] = str(self.photo_outer_radius)
        arg_map[18] = self.photo_centroid_type.value
        arg_map[19] = "0" # Unused placeholder?
        
        # Survey / Map (cont)
        arg_map[20] = str(self.trim_size)
        arg_map[21] = "1" if self.lss_mapping else "0"

        # Skip Flags
        arg_map[22] = "1" if self.skip_time_shift else "0"
        arg_map[23] = "1" if self.skip_background_subtraction else "0"
        arg_map[24] = "1" if self.skip_rfi_removal else "0"
        arg_map[25] = "1" if self.skip_surface_modeling else "0"
        
        # Convert to list.
        # argv[2] is at index 0 of our result.
        # Check range 2 to 25.
        args_list = []
        for i in range(2, 26):
            val = arg_map.get(i, "0")
            args_list.append(str(val))
            
        # Append exclusion bands (argv[22]...)
        for band_edge in self.exclusion_bands:
            args_list.append(str(band_edge))
            
        return args_list
