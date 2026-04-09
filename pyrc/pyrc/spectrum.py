import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
import utils

class Spectrum:
    def __init__(self, file_path: str, ifnum, plnum, including_frequency_ranges, excluding_frequency_ranges, including_time_ranges, excluding_time_ranges):
        self.filepath = file_path

        with fits.open(self.filepath) as hdul:
            self.header = hdul[0].header
            self.data = Table(hdul[1].data)

            # Find total number of feeds and channels
            ifnums = np.unique(self.data['IFNUM'])
            plnums = np.unique(self.data['PLNUM'])

            # Find total number of channels
            self.channel_count = len(ifnums) * len(plnums)

            self.data = self.data[
                (self.data['IFNUM'] == ifnum) &
                (self.data['PLNUM'] == plnum) &
                (self.data['CALSTATE'] == 0) &
                (self.data['SWPVALID'] == 0)
                ]

            self.ifnum = ifnum
            self.plnum = plnum

            # Accept frequency ranges
            self.including_frequency_ranges = including_frequency_ranges
            self.excluding_frequency_ranges = excluding_frequency_ranges

            # Accept time ranges
            self.including_time_ranges = including_time_ranges
            self.excluding_time_ranges = excluding_time_ranges

    def spectrum(self):
        if self.including_time_ranges or self.excluding_time_ranges:
            self.data = utils.filter_time_ranges(self.header, self.data, self.including_time_ranges, self.excluding_time_ranges)
        if self.including_frequency_ranges or self.excluding_frequency_ranges:
            frequencies, self.data['DATA'] = utils.filter_frequency_ranges(self.header, self.data, self.ifnum, self.including_frequency_ranges, self.excluding_frequency_ranges)
        else:
            frequencies = utils.get_frequency_range(self.header, self.ifnum)
            frequencies = np.linspace(frequencies[1], frequencies[0], frequencies[2])

        data_start_index, post_cal_start_index, off_start_index = utils.find_calibrations(self.header, self.data, self.channel_count)
        self.off_start_index = off_start_index

        if self.off_start_index:
            on_spectrum = utils.integrate_data(self.header, self.data['DATA'][:self.off_start_index], "spectrum")
            off_spectrum = utils.integrate_data(self.header, self.data['DATA'][self.off_start_index:], "spectrum")

            spectrum = on_spectrum - off_spectrum
        else:
            spectrum = utils.integrate_data(self.header, self.data['DATA'], "spectrum")

        return [frequencies, spectrum]


if __name__ == "__main__":
    filepath = "C:/Users/starb/Downloads/0144767_validated_corrected.fits"

    s0 = Spectrum(filepath, 0, 0, None, None, None, None)
    spectrum0 = s0.spectrum()

    s1 = Spectrum(filepath, 0, 1, None, None, None, None)
    spectrum1 = s1.spectrum()