import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

class Validation:
    def __init__(self, file_path: str):
        '''
        Initialization function for provided file. Responsible for
        opening the SDFITS file's header and data and initializing
        necessary params.
        '''

        self.filepath = file_path

        with fits.open(self.filepath) as hdul:
            # Use astropy's built in verification methods
            hdul.verify('exception')

            self.header = hdul[0].header
            self.data = Table(hdul[1].data)

    def mask_nan_values(self): # TODO revisit to see if necessary
        '''
        Masks values of data that are nonphysical to prevent errors.
        '''

        # Initialize data cube in SDFITS file
        data_array = np.ma.array(self.data['DATA'])

        # Ensure no NaN values in data cube
        nan_mask = np.isnan(data_array)
        if np.any(nan_mask):
            # Mask values if NaN values exist
            self.data['DATA'].mask = nan_mask

    def validate_time(self):
        '''
        Ensure necessary time keeping parameters function properly.
        '''

        try:
            t0 = Time(self.header["DATE"], format="isot")
            dt = Time(self.data["DATE-OBS"], format="isot") - t0

        except Exception as e:
            print("Could not parse observation times!")

    def validate_physical_values(self):
        '''
        Ensure physical values are indeed physical.
        '''

        columns = ["DURATION", "EXPOSURE", "TSYS", "TCAL", "LST", "ELEVATIO", "TAMBIENT", "PRESSURE",
                   "HUMIDITY", "RESTFREQ", "FREQRES", "TRGTLONG", "MJD", "UTSECS"]

        for column in columns:
            try:
                # Initialize the array for the given column
                column_data_array = self.data[column]

                # Ensure that columns that should not have negative values do not have negative values
                if np.any(column_data_array < 0):
                    mask = column_data_array >= 0
                    self.data = self.data[mask]

            except Exception as e:
                print(f"Nonphysical value detected for {column}!")

    def get_channels(self):
        '''
        SDFITS files are marked with channels not to be used due
        to poor signal or band-pass roll off. These should be
        removed for accurate measurements.
        '''

        for key, value in self.header.items():
            if key == ("HISTORY"):
                if value.startswith("START,STOP"):
                    # Extract all integers from the string
                    value = value.replace(",", " ").strip()
                    value = value.split(" ")

                    channels = []
                    for k in value:
                        k = str(k).strip()
                        try:
                            # If it's an integer add it to the channels list
                            k = int(k)
                            channels.append(k)

                        except ValueError:
                            # If it can't become a integer, skip it
                            continue

                    # Remove all string type characters from the channels list
                    start_channel = int(channels[0])
                    stop_channel = int(channels[1])

        self.data['DATA'] = [row[start_channel:stop_channel + 1] for row in self.data['DATA']]

    def validate(self):
        '''
        Validates the data in a file. Ensures all date cards
        comply to the datetime library standard and that
        recorded measurements are physical.

        Finds the pre- and post- calibration spikes and
        removes invalid channels. Saves the polished file.
        '''

        # Mask nan values
        self.mask_nan_values()

        # Validate necessary time elements
        self.validate_time()

        # Validate physical values
        self.validate_physical_values()

        # Remove poor data channels
        self.get_channels()