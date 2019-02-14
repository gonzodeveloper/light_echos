from astropy.io import fits
from multiprocessing import Pool
import matplotlib.pyplot as plt
from glob import glob
import numpy as np
import pandas as pd
import pickle

def block_array(arr, nrows, ncols):
    """
    Thank you unutbu @ stack overflow
    """
    h, w = arr.shape
    return (arr.reshape(h//nrows, nrows, -1, ncols)
               .swapaxes(1,2)
               .reshape(-1, nrows, ncols))

def plot_cut(cut, field_name, cut_idx, path):
    title = "{}_{}".format(field_name, cut_idx)

    plt.figure(figsize=(10,10))
    plt.title(title)
    plt.imshow(cut, cmap='gray', vmin=0, vmax=3.5)
    ax = plt.gca()
    ax.set_xticks(np.arange(0, 128, 32))
    ax.set_yticks(np.arange(0, 128, 32))
    plt.grid(color='r', linewidth='4')
    plt.savefig("{}.png".format(title), dpi=100)
    plt.close()


def cut_echo(data, coords):

    midpoint = ((coords['x1'] + coords['x2']) / 2,
                (coords['y1'] + coords['y2']) / 2)

    # 128 x 128 box around midpoint of region
    size = 128
    # Offset cuts, exact centering will split echo weirdly
    offset = 16

    # Make a bounding box
    x1 = midpoint[0] - (size / 2)
    x2 = midpoint[0] + (size / 2)
    y1 = midpoint[1] - (size / 2)
    y2 = midpoint[1] + (size / 2)

    cuts = []
    cuts[0] = data[x1 + offset:x2 + offset, y1 + offset:y2 + offset]
    cuts[1] = data[x1 + offset:x2 + offset, y1 - offset:y2 - offset]
    cuts[2] = data[x1 - offset:x2 - offset, y1 + offset:y2 + offset]
    cuts[3] = data[x1 - offset:x2 - offset, y1 - offset:y2 - offset]

    return cuts


def process_file(region_file, echo_path):

    regions = pd.read_csv(region_file, sep=' ', header=None, names=['x1', 'y1', 'x2', 'y2', 'status'], skiprows=[0])
    mask = regions['status'] == "definitely"
    definite = regions.loc[mask]

    if len(definite) == 0:
        return -1

    field = region_file.split(".")[0]
    subtractions = glob("{}/{}/sub*".format(echo_path, field))

    for sub in subtractions:
        hdul = fits.open("../{}".format(field))
        data = hdul[0].data

        for idx, echo in definite.iterrows():

            cuts = cut_echo(data, echo)

            # Make a stack of all 32x32 cuts of echo including 4x shift
            cut_stack = [blocks for arr in cuts for blocks in block_array(arr, 32, 32)]




if __name__ == "__main__":

    echo_path = "./echos"

    files = glob("regions/*.xy.reg")

    total = len(files)
    with Pool(processes=16) as pool:
        results = [pool.apply_async(process_file, args=(region, echo_path,)) for region in files]
        for idx, r in enumerate(results):
            print('Completed processing % of % files' % (idx+1, total), end='\r')

        print()


