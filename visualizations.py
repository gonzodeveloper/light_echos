from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import pickle


def main_plot(data, echo):
    plt.suptitle("Tycho Light Echo")
    plt.subplot(121)
    plt.imshow(data, cmap='gray', vmin=0, vmax=3.5)
    plt.subplot(122)
    plt.imshow(echo, cmap='gray', vmin=0, vmax=3.5)
    ax = plt.gca()
    ax.set_xticks(np.arange(0,128,16))
    ax.set_yticks(np.arange(0,128,16))
    plt.grid(color='r')
    plt.show()


def cut_echo(echo, size=16):
    cuts = []
    for x in np.arange(0, len(echo), size):
        for y in np.arange(0, len(echo), size):
            cut = echo[x:x + size, y:y + size]
            cuts.append(cut)
    return cuts


def plot_cuts(cuts, label):
    dim = np.sqrt(len(cuts))
    for i in range(len(cuts)):
        plt.subplot(dim, dim, i+1)
        plt.imshow(cuts[i], cmap='gray', vmin=0, vmax=3.6)
    plt.suptitle(label)
    plt.show()


if __name__ == "__main__":
    hdul = fits.open("F0014+62_1130.bd/sub_0-1.fits")
    hdr = hdul[0].header

    data = hdul[0].data
    x0 = 620
    x1 = 748
    y0 = 620
    y1 = 748

    echo_0 = data[x0:x1, y0:y1]
    echo_1 = data[x0+4:x1+4, y0:y1]
    echo_2 = data[x0:x1, y0+4:y1+4]
    echo_3 = data[x0+4:x1+4, y0+4:y1+4]

    echos = [echo_0, echo_1, echo_2, echo_3]
    cuts = [cut_echo(x) for x in echos]
    #plot_cuts(cuts[0], "echo_0")
    main_plot(data, echo_0)
    exit()
    plot_cuts(cuts[1], "echo_1")
    plot_cuts(cuts[2], "echo_2")
    plot_cuts(cuts[3], "echo_3")
    cut_stack = [y for x in cuts for y in x]
    print(len(cut_stack))
    print(cut_stack[1].shape)
    with open("data/F0014+62_cuts_01.pkl", "wb+") as file:
        pickle.dump(cuts, file)






