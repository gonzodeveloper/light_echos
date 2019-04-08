from astropy.io import fits
from astropy.wcs import WCS
from regions import write_ds9
from regions import PixCoord, LinePixelRegion
from multiprocessing import Pool
import matplotlib.pyplot as plt
from glob import glob
import numpy as np
import pandas as pd
import os
from pathlib import Path


def block_array(arr, nrows, ncols):
    """
    Thank you unutbu @ stack overflow
    """
    h, w = arr.shape
    return (arr.reshape(h//nrows, nrows, -1, ncols)
               .swapaxes(1,2)
               .reshape(-1, nrows, ncols))


def plot_cut(cut, field_name, cut_idx):
    title = "{}_{}".format(field_name, cut_idx)

    plt.figure(figsize=(10,10))
    plt.title(title)
    plt.imshow(np.flipud(cut), cmap='gray', vmin=0, vmax=3.5)
    ax = plt.gca()
    ax.set_xticks(np.arange(0, 128, 32))
    ax.set_yticks(np.arange(0, 128, 32))
    plt.grid(color='r', linewidth='4')
    plt.show()


def load_print_cut(filename):
    echo = np.memmap(filename, dtype='float32', mode='r', shape=(16, 32, 32))

    for idx, cut in enumerate(echo):
        plt.subplot(4,4,idx+1)
        plt.imshow(cut, cmap='gray', vmin=0, vmax=3.5)

    plt.show()


def cut_echo(data, midpoint):

    # 128 x 128 box around midpoint of region
    size = 128
    # Offset cuts, exact centering will split echo weirdly
    offset = 16

    # Make a bounding box
    x1 = int(midpoint[0] - (size / 2))
    x2 = int(midpoint[0] + (size / 2))
    y1 = int(midpoint[1] - (size / 2))
    y2 = int(midpoint[1] + (size / 2))

    # Isolate echo
    echo = data[x1 + offset:x2 + offset, y1 + offset:y2 + offset]

    # Cut echo into 32x32 sub-arrays for training
    cuts = block_array(echo, 32, 32)

    boundary = (x1 + offset, x2 + offset, y1 + offset, y2 + offset)

    return cuts, boundary


def write_cut_region(boundaries, write_path, wcs):
    Y1, Y2, X1, X2 = boundaries
    region_lines = []
    # Draw x lines
    for x_coord in np.arange(X1, X2+1, 32):
        line = LinePixelRegion(start=PixCoord(x=x_coord, y=Y1), end=PixCoord(x=x_coord, y=Y2))
        region_lines.append(line.to_sky(wcs))

    # Draw y lines
    for y_coord in np.arange(Y1, Y2+1, 32):
        line = LinePixelRegion(start=PixCoord(x=X1, y=y_coord), end=PixCoord(x=X2, y=y_coord))
        region_lines.append(line.to_sky(wcs))

    region_file = "{}.reg".format(write_path)
    write_ds9(region_lines, region_file)


def process_file(region_file, path, results_path):

    regions = pd.read_csv(region_file, sep=' ', header=None, index_col=False,
                          names=['RA1', 'Dec1', 'RA2', 'Dec2', 'status'], skiprows=[0])

    # Only evaluate files with echoes marked as definite hits
    mask = regions['status'] == "definitely"
    definite = regions.loc[mask]

    if len(definite) == 0:
        return 0, 0

    # Get the field name and make a results sub-folder for that field
    field = region_file.split("/")[3].split(".")[0]\

    if not os.path.exists("{}/{}".format(results_path, field)):
        os.mkdir("{}/{}".format(results_path, field))

    # Get all subtractions for the field
    subtractions = glob("{}/{}/sub*".format(path, field))

    # Iterate through subtractions for individual FITS file
    for idx_s, sub in enumerate(subtractions):
        # Get data from fits file and a world coordinate system
        hdul = fits.open(sub)
        data = hdul[0].data
        wcs = WCS(sub)

        # Make a symbolic link to subtraction in results sub-dir for field
        sub_name = "sub_{}.fits".format(idx_s)
        os.link(sub, "{}/{}/{}".format(results_path, field, sub_name))

        # Iterate through 'definite' hits as defined by REG file (each is marked by a line)
        for idx_e, echo in definite.iterrows():
            # Get the xy midpoint for the echo
            x1, y1 = wcs.wcs_world2pix(echo['RA1'], echo['Dec1'], 0)
            x2, y2 = wcs.wcs_world2pix(echo['RA2'], echo['Dec2'], 0)
            midpoint = ((x1 + x2) / 2,
                        (y1 + y2) / 2)

            # Cuts = x32 numpy arrays of cut up image
            # Boundary = X1 X2 Y1 Y2 boundaries for image
            cuts, boundary = cut_echo(data, midpoint)

            # Write out boundary region
            cut_path = "{}/{}/sub_{}_cut_{}".format(results_path, field, idx_s, idx_e)
            write_cut_region(boundary, cut_path, wcs)
            # Write out cut stack
            cut_fname = "{}.dat".format(cut_path)
            fp = np.memmap(cut_fname, dtype='float32', mode='w+', shape=cuts.shape)
            fp[:] = cuts[:]
            del fp
            # Touch a file for labels
            Path("{}.lab".format(cut_path))
    return len(definite), len(subtractions)


if __name__ == "__main__":

    path = "./echoes"
    results_path = "{}/finds".format(path)
    files = glob("{}/echo_regions/*.deg.reg".format(path))

    echo_sum = 0
    sub_sum = 0

    total = len(files)
    with Pool(processes=2) as pool:
        results = [pool.apply_async(process_file, args=(region, path, results_path,)) for region in files]
        for idx, r in enumerate(results):
            sums = r.get()
            n_echos, n_subs = sums
            echo_sum += n_echos
            sub_sum += n_subs
            print('Completed processing %s of %s files' % (idx+1, total), end='\n')

        print("Found {} echoes in {} subtractions!".format(echo_sum, sub_sum))


