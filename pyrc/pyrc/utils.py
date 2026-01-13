import os
from contextlib import contextmanager
import re
from astropy.time import Time
import astropy.units as u
import numpy as np
from .gain_calibration_validation import Validation
@contextmanager
def validated_temp_path(filepath: str):
    v = Validation(filepath)
    validated_path = v.validate()  # temp path

    yield validated_path
    '''finally:
        try:
            os.remove(validated_path)
        except FileNotFoundError:
            pass'''

def parse_history(header):
    # Find all instances of the HISTORY card in the FITS header
    entries = header.get('HISTORY', [])
    if isinstance(entries, str):
        entries = [entries]

    parsed = {}
    extra_lines = []

    # Parse each HISTORY card
    for entry in entries:
        # Remove inline comments
        clean_entry = entry.split('/')[0].strip()

        # Definte the parsing strategy
        match = re.match(r'^\s*([A-Za-z0-9_,]+(?: [A-Za-z0-9_,]+)*)\s+(.*)', clean_entry)
        if match:
            key = match.group(1).strip()
            val_str = match.group(2).strip()

            # Handle underscore-separated numeric values like "1355_1435"
            if re.fullmatch(r'\d+_\d+', val_str):
                a, b = val_str.split('_')
                parsed[key] = (float(a), float(b))
                continue

            # Handle comma/space-separated numeric values
            parts = val_str.replace(',', ' ').split()
            try:
                if all(re.fullmatch(r'-?\d+(\.\d+)?', p) for p in parts):
                    vals = [float(p) for p in parts]
                    parsed[key] = vals if len(vals) > 1 else vals[0]
                else:
                    parsed[key] = val_str
            except ValueError:
                parsed[key] = val_str
        elif clean_entry:
            extra_lines.append(entry.strip())

    if extra_lines:
        parsed["_extra"] = extra_lines

    return parsed

def get_frequency_range(header, ifnum):
    # Use the parse history function to get dictionary of sub FITS header
    history = parse_history(header)

    # Determine whether the file is LOW or HIGH resolution
    datamode = history.get('DATAMODE')

    # Find the valid channel range within the file
    low_channel = int(history.get('START,STOP channels')[0])
    high_channel = int(history.get('START,STOP channels')[1])

    # Find the number of channels used across the file
    channel_count = high_channel - low_channel + 1

    if datamode == 'HIRES':
        # For HIGH resolution files, retrieve the proper band center from the history field
        band_center = history.get('HIRES bands')[ifnum]
        # Process to find the band width
        band_width = header['OBSBW']

        # Find the lowest and highest frequencies of the file
        low_frequency = band_center - (band_width / 2)
        high_frequency = band_center + (band_width / 2)

        return low_frequency, high_frequency, channel_count

    elif datamode == 'LOWRES':
        # For LOW resolution files, retrieve the proper band center from the FITS header
        band_center = header['OBSFREQ']
        # Process to find the band width
        band_width = header['OBSBW']

        # Pull the low and high frequencies from the provided header
        low_frequency = history.get('RFFILTER')[0] # band_center - (band_width / 2)
        high_frequency = history.get('RFFILTER')[1] # band_center + (band_width / 2)

        return low_frequency, high_frequency, channel_count

    else:
        # If the file is not LOW or HIGH resolution then raise an exception
        # TODO add graceful error handling
        raise ValueError(f"Unknown datamode: {datamode}")

def integrate_data(header, data, mode, frequencies=None):
    if mode == "continuum":
        intensities = np.array(data['DATA'])
        intensities = np.sum(intensities, axis=1)

        times = Time(data["DATE-OBS"], format='isot')
        t0 = Time(header["DATE"], format="isot")
        time_rel = (times - t0)

        return [time_rel.sec, intensities]

    elif mode == "spectrum":
        intensities = np.array(data)
        intensities = np.sum(intensities, axis=0)

        return intensities
    return None


def find_calibrations(header, data, channel_count):
    '''
    Calibration spikes must be systematically located using
    CALSTATE and SWPVALID to determine location and OBSMODE
    to determine if an ON/OFF transition occurs.
    '''

    # Initialize necessary indices
    data_start_ind = None
    post_cal_start_ind = None
    off_start_index = None

    # Create a counter for valid data
    counter = 0

    # Initialize confirmation Booleans
    cal_started = False
    pre_cal_complete = False

    for ind, i in enumerate(data):
        # If the CALSTATE is equal to 1 then state that calibration has started
        if i['CALSTATE'] == 1:
            cal_started = True

        # If data begins being collected (CALSTATE = 0 and SWPVALID = 1) then state that pre calibration has started and the first data index can be recorded
        if cal_started and i["CALSTATE"] == 0 and i["SWPVALID"] == 1 and not pre_cal_complete:
            data_start_ind = ind
            pre_cal_complete = True

        # If the pre calibration is complete and the sweep is no longer valid then keep track of this as the post calibration start
        if ind > 0 and pre_cal_complete and i["SWPVALID"] == 0 and data[ind - 1]["SWPVALID"] == 0:
            if post_cal_start_ind is None:
                post_cal_start_ind = ind - 1

        # Reset the post cal index to None if the above condition is False.
        # This allows for sweeps to be invalid within the observation such as during on/off transition or invalid data blips.
        else:
            post_cal_start_ind = None

        # Keep track of contiguous data points
        if pre_cal_complete and i['CALSTATE'] == 0 and i['SWPVALID'] == 1:
            counter += 1

        # If 3 or less valid data points across all channels have been collected and the sweep becomes
        # invalid then treat this section of data as invalid and reset necessary params.
        if counter <= 3 * channel_count and i['SWPVALID'] == 0 and data_start_ind:
            data_start_ind = None
            pre_cal_complete = False

        # When pre calibration is complete and a new cal spike begins,
        # break the loop (post cal index will already be recorded if available)
        if pre_cal_complete and i['SWPVALID'] == 0 and i['CALSTATE'] == 1:
            break

    if not pre_cal_complete:
        pre_cal_complete = True
        data_start_ind = 0

        for ind, i in enumerate(data):
            # If the pre calibration is complete and the sweep is no longer valid then keep track of this as the post calibration start
            if ind > 0 and pre_cal_complete and i["SWPVALID"] == 0 and data[ind - 1]["SWPVALID"] == 0:
                if post_cal_start_ind is None:
                    post_cal_start_ind = ind - 1

            # Reset the post cal index to None if the above condition is False.
            # This allows for sweeps to be invalid within the observation such as during on/off transition or invalid data blips.
            else:
                post_cal_start_ind = None

            # When pre calibration is complete and a new cal spike begins,
            # break the loop (post cal index will already be recorded if available)
            if i['SWPVALID'] == 0 and i['CALSTATE'] == 1:
                break

    if not post_cal_start_ind:
        post_cal_start_ind = len(data) - 1

    # If the file is an on/off file then find when the transition occurs and store an additional index
    if header['OBSMODE'] == 'onoff':
        for ind, i in enumerate(data):
            target = 'onoff:off'

            if target in i['OBSMODE']:
                offstart = ind
                off_start_index = offstart

                break

    return data_start_ind, post_cal_start_ind, off_start_index


def filter_time_ranges(header, data, including_time_ranges, excluding_time_ranges):
    # Create time array to be filtered
    t0 = Time(header["DATE"], format="isot")
    dt = Time(data["DATE-OBS"], format="isot") - t0
    times = dt.to_value(u.s)

    # If there are ranges to include then filter all times not within those ranges
    if including_time_ranges:
        include_mask = np.zeros(len(times), dtype=bool)

        # Use a mask to get all valid indices
        for start_time, end_time in including_time_ranges:
            include_mask |= (start_time < times) & (times < end_time)

        # Apply the mask
        data = data[include_mask]

    # If there are ranges to exclude then filter all times within those ranges
    if excluding_time_ranges:
        exclude_mask = np.ones(len(times), dtype=bool)

        # Use a mask to get all valid indices
        for start_time, end_time in excluding_time_ranges:
            exclude_mask &= ~((start_time < times) & (times < end_time))

        # Apply the mask
        data = data[exclude_mask]

    return data

def filter_frequency_ranges(header, data, ifnum, including_frequency_ranges, excluding_frequency_ranges):
    # Get the necessary frequency data
    low_frequency, high_frequency, n_channels = get_frequency_range(header, ifnum)

    # Create an array of frequencies from the highest frequency to the lowest frequency of length of total channels
    frequencies = np.linspace(high_frequency, low_frequency, n_channels)

    # If there are ranges to include then filter all frequencies not within those ranges
    if including_frequency_ranges:
        include_freq_mask = np.zeros(len(frequencies), dtype=bool)

        # Use a mask to get all valid indices
        for fmin, fmax in including_frequency_ranges:
            low, high = sorted((fmin, fmax))
            include_freq_mask |= (frequencies > low) & (frequencies < high)

        # Update the frequencies and data to reflect the mask
        frequencies = frequencies[include_freq_mask]
        data['DATA'] = [row[include_freq_mask] for row in data['DATA']]

    if excluding_frequency_ranges:
        exclude_freq_mask = np.ones(len(frequencies), dtype=bool)

        # Use a mask to get all valid indices
        for fmin, fmax in excluding_frequency_ranges:
            low, high = sorted((fmin, fmax))
            exclude_freq_mask &= ~((frequencies > low) & (frequencies < high))

        # Update the frequencies and data to reflect the mask
        frequencies = frequencies[exclude_freq_mask]
        data['DATA'] = [row[exclude_freq_mask] for row in data['DATA']]

    return frequencies, data['DATA']
