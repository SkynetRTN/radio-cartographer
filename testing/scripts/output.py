import os

from astropy.io import fits
import numpy as np


def get_data(filename):
    image_data = np.genfromtxt(filename + ".txt", delimiter='\t')
    os.remove(filename + ".txt")
    return image_data


def create_file():
    files = ['main'] #, 'path', 'scale', 'weight', 'correlation', 'raw']
    hdu_list = fits.HDUList()
    for file in files:
        if not os.path.exists(f'{file}.txt'):
            continue
        with open(f'{file}.txt', 'r') as f:
            data = get_data(file)
            if file == 'main':
                hdu_list.append(fits.PrimaryHDU(data))
            else:
                hdu_list.append(fits.ImageHDU(data))

    hdu_list.writeto('output.fits', overwrite=True)


if __name__ == "__main__":
    create_file()
