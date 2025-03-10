import os
import pathlib
from enum import Enum


class RadioCartographerSetting():
    channel = None
    gain_calibration = None
    image_coordinate_system = None
    time_delay = None
    include_raw = None

    min_frequency = None
    max_frequency = None
    skip_frequencies = []

    two_d_background_subtraction_scale = None
    one_d_background_subtraction_scale = None
    two_d_surface_model = None
    surface_model_weighting_scale = None

    def set_basics(self, channel, gain_calibration, image_coordinate_system, time_delay, include_raw):
        self.channel = channel
        self.gain_calibration = gain_calibration
        self.image_coordinate_system = image_coordinate_system
        self.time_delay = time_delay
        self.include_raw = include_raw

    def set_frequency(self, min_frequency, max_frequency, skip_frequencies):
        self.min_frequency = min_frequency
        self.max_frequency = max_frequency
        self.skip_frequencies = skip_frequencies

    def set_rfi(self, one_d_background_subtraction_scale, preset):
        if preset == RfiPreset.BRIGHT:
            self.set_rfi_custom(one_d_background_subtraction_scale, 0.35, SurfaceModel.WITH_NOISE_PRIOR, 0.33)
        elif preset == RfiPreset.FAINT:
            self.set_rfi_custom(one_d_background_subtraction_scale, 0.7, SurfaceModel.WITHOUT_NOISE_PRIOR, 0.33)

    def set_rfi_custom(self, one_d_background_subtraction_scale, two_d_background_subtraction_scale,
                two_d_surface_model, surface_model_weighting_scale):
        self.two_d_background_subtraction_scale = two_d_background_subtraction_scale
        self.one_d_background_subtraction_scale = one_d_background_subtraction_scale
        self.two_d_surface_model = two_d_surface_model
        self.surface_model_weighting_scale = surface_model_weighting_scale

    def get_command(self):
        if (self.channel is None or self.gain_calibration is None or self.image_coordinate_system is None
                or self.time_delay is None or self.include_raw is None):
            raise ValueError("Basic Parameters must be set")
        if self.min_frequency is None or self.max_frequency is None:
            raise ValueError("Frequency Parameters must be set")
        for freq in self.skip_frequencies:
            if freq[0] >= freq[1] or freq[0] < self.min_frequency or freq[1] > self.max_frequency:
                raise ValueError("Invalid skip frequency")
        if (self.two_d_background_subtraction_scale is None or self.one_d_background_subtraction_scale is None
                or self.two_d_surface_model is None or self.surface_model_weighting_scale is None):
            raise ValueError("RFI Parameters must be set")

        timeshift_option = TimeShift.CUSTOM.value if is_number(self.time_delay) else self.time_delay.value
        timeshift_value = self.time_delay if is_number(self.time_delay) else None
        # composite 1 none equatorial auto None 0 1350 1435 20 0.7 1 0.6666 0 1.25 5.0 brightest 0 0.0 0
        return (f"{self.channel.value} 1 {self.gain_calibration.value} {self.image_coordinate_system.value}"
                f" {timeshift_option} {timeshift_value}"
                f" {self.include_raw.value} {self.min_frequency} {self.max_frequency}"
                f" {self.one_d_background_subtraction_scale} {self.two_d_background_subtraction_scale}"
                f" {self.two_d_surface_model.value} {self.surface_model_weighting_scale}"
                f" 0 1.25 5.0 brightest 0 0.0 0"  # photometry parameters
                f" {' '.join([f'{freq[0]} {freq[1]}' for freq in self.skip_frequencies])}")


class CoordinateSystem(Enum):
    EQUATORIAL = "equatorial"
    GALACTIC = "galactic"


class JobChannel(Enum):
    LEFT = "left"
    RIGHT = "right"
    COMPOSITE = "composite"


class Receiver(Enum):
    ONE = "1"
    TWO = "2"


class Calibration(Enum):
    PRE = "pre"
    POST = "post"
    INTERPOLATED = "interpolated"
    NONE = "none"


class TimeShift(Enum):
    AUTO = "auto"
    CUSTOM = "custom"
    OFF = "off"


class RfiPreset(Enum):
    BRIGHT = 1
    FAINT = 0
    NA = -1


class SurfaceModel(Enum):
    WITH_NOISE_PRIOR = 1
    WITHOUT_NOISE_PRIOR = 0


class IncludeRaw(Enum):
    TRUE = 1
    FALSE = 0


def is_number(s):
    return isinstance(s, int) or isinstance(s, float)


if __name__ == "__main__":
    filename = "/skynet/radio-cartographer/testing/Temporary/ra_zero.fits"

    settings = RadioCartographerSetting()

    settings.set_basics(JobChannel.COMPOSITE,
                        Calibration.INTERPOLATED,
                        CoordinateSystem.EQUATORIAL,
                        TimeShift.AUTO, IncludeRaw.TRUE)
    settings.set_frequency(1355, 1455, [])
    settings.set_rfi(0.7, RfiPreset.FAINT)

    full_command = f"{filename} {settings.get_command()}"
    print(full_command)
    f = open(os.path.join(pathlib.Path(__file__).parent.parent.absolute(), "Temporary", "command"), "w")
    f.write(full_command)
    f.close()
