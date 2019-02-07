from astropy.io import fits
import numpy as np
import pandas as pd


hdul = fits.open("F0056+58_1132.bb/Stack1.fits")
hdr = hdul[0].header

data = hdul[0].data
for key, val in hdr.items():
    print(key, val)
