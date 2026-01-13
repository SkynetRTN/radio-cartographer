import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.stats import linregress
import rcr
from . import utils


class Gain_Calibration:
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
                (self.data['PLNUM'] == plnum)
                ]

            self.ifnum = ifnum
            self.plnum = plnum

            # Accept frequency ranges
            self.including_frequency_ranges = including_frequency_ranges
            self.excluding_frequency_ranges = excluding_frequency_ranges

            # Accept time ranges
            self.including_time_ranges = including_time_ranges
            self.excluding_time_ranges = excluding_time_ranges

    def parse_calibration_spike(self, data):
        print(data["CALSTATE"])
        on_mask = data[
            (data['CALSTATE'] == 1) &
            (data['SWPVALID'] == 0)
            ]

        off_mask = data[
            (data['CALSTATE'] == 0) &
            (data['SWPVALID'] == 0)
            ]

        return on_mask, off_mask

    def linear(self, x, params): # model function
        return params[0] + x * params[1]

    def d_linear_1(self, x, params): # first model parameter derivative
        return 1

    def d_linear_2(self, x, params): # second model parameter derivative
        return x

    def perform_rcr(self, array):
        x = array[0].copy()
        x -= np.average(x)

        y = array[1].copy()

        result = linregress(x, y)
        guess = [result.slope, result.intercept]

        model = rcr.FunctionalForm(self.linear,
                                   x,
                                   y,
                                   [self.d_linear_1, self.d_linear_2],
                                   guess
                                   )

        r = rcr.RCR(rcr.SS_MEDIAN_DL)
        r.setParametricModel(model)
        r.performBulkRejection(y)

        # Fetch indices
        indices = r.result.indices

        # Keep on valid indices
        x = np.array([x[i] for i in indices])
        y = np.array([y[i] for i in indices])
        best_fit_parameters = model.result.parameters

        sigma = (1 / (len(x) - 2)) * np.sum((y - best_fit_parameters[1] * x - best_fit_parameters[0]) ** 2)
        m_sd = np.sqrt(sigma / np.sum((x - np.mean(x)) ** 2))
        b_sd = np.sqrt(sigma * ((1 / len(x)) + ((np.mean(x) ** 2) / np.sum((x - np.mean(x)) ** 2))))
        uncertainties = (b_sd, m_sd)

        return best_fit_parameters, uncertainties

    def calculate_calibration_height(self, calibration):
        diode_on, diode_off = self.parse_calibration_spike(calibration)

        # Check that on and off sections are greater than 2 points to perform fitting
        if len(diode_on) >= 4 and len(diode_off) >= 4:
            diode_on_array = utils.integrate_data(self.header, diode_on, "continuum")
            diode_off_array = utils.integrate_data(self.header, diode_off, "continuum")

            diode_on_best_fit_parameters, diode_on_uncertainties = self.perform_rcr(diode_on_array)
            diode_off_best_fit_parameters, diode_off_uncertainties = self.perform_rcr(diode_off_array)

            evaluation_time = (np.average(diode_on_array[0]) + np.average(diode_off_array[0])) / 2
            diode_on_evaluation_time = evaluation_time - np.average(diode_on_array[0])
            diode_off_evaluation_time = evaluation_time - np.average(diode_off_array[0])

            diode_on_y = diode_on_evaluation_time * diode_on_best_fit_parameters[1] + diode_on_best_fit_parameters[0]
            diode_off_y = diode_off_evaluation_time * diode_off_best_fit_parameters[1] + diode_off_best_fit_parameters[0]

            calibration_delta = diode_on_y - diode_off_y
            calibration_uncertainty = np.sqrt(diode_on_uncertainties[0]**2 + diode_off_uncertainties[0]**2 + (diode_on_uncertainties[1] * diode_on_evaluation_time)**2 + (diode_off_uncertainties[1] * diode_off_evaluation_time)**2)

            return calibration_delta, calibration_uncertainty
        else:
            return None, None

    def Gain_calibration(self):
        if self.including_time_ranges or self.excluding_time_ranges:
            self.data = utils.filter_time_ranges(self.header, self.data, self.including_time_ranges, self.excluding_time_ranges)
        if self.including_frequency_ranges or self.excluding_frequency_ranges:
            frequencies, self.data['DATA'] = utils.filter_frequency_ranges(self.header, self.data, self.ifnum, self.including_frequency_ranges, self.excluding_frequency_ranges)
        else:
            frequencies = utils.get_frequency_range(self.header, self.ifnum)
            frequencies = np.linspace(frequencies[1], frequencies[0], frequencies[2])

        data_start_index, post_cal_start_index, off_start_index = utils.find_calibrations(self.header, self.data, self.channel_count)
        self.data_start_index = data_start_index
        self.post_cal_start_index = post_cal_start_index
        self.off_start_index = off_start_index

        pre_calibration = self.data[:self.data_start_index]
        post_calibration = self.data[self.post_cal_start_index:]

        pre_calibration_intensity = None
        post_calibration_intensity = None

        pre_calibration_intensity, pre_calibration_uncertainty = self.calculate_calibration_height(pre_calibration)
        post_calibration_intensity, post_calibration_uncertainty = self.calculate_calibration_height(post_calibration)

        continuum = utils.integrate_data(self.header, self.data[self.data_start_index:self.post_cal_start_index], "continuum")
        deltapre = "None"
        deltapost = "None"
        cal_method = "None"
        if pre_calibration_intensity and post_calibration_intensity:
            z_score = abs(pre_calibration_intensity - post_calibration_intensity) / np.sqrt(pre_calibration_uncertainty ** 2 + post_calibration_uncertainty ** 2)
            deltapre = pre_calibration_intensity
            deltapost = post_calibration_intensity
            if z_score >= 1.96:
                cal_method = "interpolated"
            else:
                cal_method = "average"
        elif pre_calibration_intensity:
            continuum[1] /= pre_calibration_intensity
            deltapre = pre_calibration_intensity
            deltapost = "None"
            method = "pre"
        elif post_calibration_intensity:
            continuum[1] /= post_calibration_intensity
            deltapost = post_calibration_intensity
            deltapre = "None"
            method = "post"


        return deltapre, deltapost, method
