#!/usr/bin/env python3

import os
import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib.patches import FancyArrowPatch


KEYS = {'x': 0,
        'y': 1,
        'fluid': 2,
        'level': 3,
        'vorticity': 4,
        'speed': 5,
        'u': 6,
        'v': 7,
        'pressure': 8,
        'h': 9}


def get_frame(i, Ly=2.0):
    """
    Get the 1D and 2D data for the ith output frame
    """
    fname1 = f"out/data-1-{i:010d}.dat"
    fname2 = f"out/data-2-{i:010d}.dat"

    # read in data from files
    data1 = np.loadtxt(fname1)
    data2 = np.loadtxt(fname2)

    # reshape 2D data
    ny = 0
    while data2[ny, 0] == 0.0:
        ny = ny + 1
    data2 = np.reshape(data2, (ny, -1, 10), order='F')

    # clip y range to [0,Ly]
    iy = 0
    while data2[iy, 0, 1] <= Ly:
        iy = iy + 1
        if iy == data2.shape[0]-1:
            break
    data2 = data2[0:iy, :, :]

    return data1, data2


def get_interface(i):
    return np.loadtxt(f"out/data-1-{i:010d}.dat")


def plot_single(i, save=True):
    """
    Plots the ith frame and saves or shows it
    """
    data1, data2 = get_frame(i)

    # plot 1D interface
    plt.plot(data1[:, 0], data1[:, 1])
    plt.xlim([0, 64])
    plt.ylim([0, 2])
    plt.gca().set_aspect(8.0)

    if save:
        plt.savefig(f"interface-{i:010d}.png", dpi=300)
    else:
        plt.show()
    plt.clf()

    # plot 2D fluid
    plt.imshow(data2[:, :, 2], interpolation='bilinear',
               origin='lower', extent=(0, 64, 0, 2), aspect=8.0)
    plt.plot(data1[:, 0], data1[:, 1])

    if save:
        plt.savefig(f"fluid-{i:010d}.png", dpi=300)
    else:
        plt.show()
    plt.clf()

    # plot level
    plt.imshow(data2[:, :, 3], interpolation='bilinear',
               origin='lower', extent=(0, 64, 0, 2), aspect=8.0)
    plt.plot(data1[:, 0], data1[:, 1])

    if save:
        plt.savefig(f"level-{i:010d}.png", dpi=300)
    else:
        plt.show()
    plt.clf()


def plot_series(data='fluid', Ly=2.0, n=0, track=False, rate=1):
    """
    Plots an animated series
    """
    try:
        i = KEYS[data]
    except KeyError:
        if data == 'interface':
            i = -1
        elif data == 'quiver':
            i = -2
        elif data == 'stream':
            i = -3

    ims = []
    fig, ax = plt.subplots()

    if n == 0:
        n = int(len([name for name in os.listdir('./out')
                    if os.path.isfile(os.path.join('./out', name))])/2)
    for k in range(0, n):
        if k % rate != 0:
            continue

        ax.set_xlim([0, 64])
        ax.set_ylim([0, Ly])
        ax.set_aspect(8.0)
        if i >= 0:
            data1, data2 = get_frame(k, Ly)
            N = data1.shape[0]

            # track the peak and keep it in the centre of the frame
            if track:
                # roll so that the peak is at the start
                imax = np.argmax(data1[:, 1])
                data1[:, 1:2] = np.roll(data1[:, 1:2], int(N*0.75)-imax, axis=0)
                data2[:, :, 2:] = np.roll(data2[:, :, 2:], int(N*0.75)-imax, axis=1)

            im1 = ax.imshow(data2[:, :, i], interpolation='bilinear',
                            origin='lower', extent=(0, 64, 0, Ly), aspect=8)
            im2, = ax.plot(data1[:, 0], data1[:, 1], color='black')
            title = plt.title(data, weight='bold', fontsize=14)
            text = plt.text(0.5, -1, s=("t = " + str(f"{k*0.01:.3f}")),
                            weight='bold', fontsize=12, ha='center',
                            transform=ax.transAxes)
            ims.append([im1, im2, title, text])
        elif i == -1:
            data1 = get_interface(k)
            N = data1.shape[0]

            # track the peak and keep it in the centre of the frame
            if track:
                # roll so that the peak is at the start
                imax = np.argmax(data1[:, 1])
                data1[:, 1:2] = np.roll(data1[:, 1:2], int(N*0.75)-imax, axis=0)

            im1, = ax.plot(data1[:, 0], data1[:, 1], color='black')
            im2, = ax.plot(data1[:, 0], data1[:, 2], color='red')
            title = plt.title(data, weight='bold', fontsize=14)
            text = plt.text(0.5, -1, s=("t = " + str(f"{k*0.01:.3f}")),
                            weight='bold', fontsize=12, ha='center',
                            transform=ax.transAxes)
            ims.append([im1, im2, title, text])
        elif i == -2:
            data1, data2 = get_frame(k, Ly)
            N = data1.shape[0]

            # track the peak and keep it in the centre of the frame
            if track:
                # roll so that the peak is at the start
                imax = np.argmax(data1[:, 1])
                data1[:, 1:2] = np.roll(data1[:, 1:2], int(N*0.75)-imax, axis=0)
                data2[:, :, 2:] = np.roll(data2[:, :, 2:], int(N*0.75)-imax, axis=1)

            im1 = ax.quiver(data2[:, ::8, 0], data2[:, ::8, 1], data2[:, ::8, 6], data2[:, ::8, 7])
            im2, = ax.plot(data1[:, 0], data1[:, 1], color='black')
            title = plt.title(data, weight='bold', fontsize=14)
            text = plt.text(0.5, -1, s=("t = " + str(f"{k*0.01:.3f}")),
                            weight='bold', fontsize=12, ha='center',
                            transform=ax.transAxes)
            ims.append([im1, im2, title, text])
        elif i == -3:
            data1, data2 = get_frame(k, Ly)
            N = data1.shape[0]

            # track the peak and keep it in the centre of the frame
            if track:
                # roll so that the peak is at the start
                imax = np.argmax(data1[:, 1])
                data1[:, 1:2] = np.roll(data1[:, 1:2], int(N*0.75)-imax, axis=0)
                data2[:, :, 2:] = np.roll(data2[:, :, 2:], int(N*0.75)-imax, axis=1)

            lines, arrows = ax.streamplot(0.125*np.round(data2[0, :, 0]/0.125),
                                          0.125*np.round(data2[:, 0, 1]/0.125),
                                          data2[:, :, 6], data2[:, :, 7])
            im2, = ax.plot(data1[:, 0], data1[:, 1], color='black')
            title = plt.title(data, weight='bold', fontsize=14)
            text = plt.text(0.5, -1, s=("t = " + str(f"{k*0.01:.3f}")),
                            weight='bold', fontsize=12, ha='center',
                            transform=ax.transAxes)
            ims.append([lines, im2, title, text])
    if i >= 0:
        fig.colorbar(im1, location='bottom')
    fig.add_axes(ax)
    final_an = anim.ArtistAnimation(fig, ims, interval=1,
                                    blit=True, repeat=False)
    final_an.save(f'{data}.gif', writer=anim.PillowWriter(fps=100))
    plt.clf()


def plot_stream(Ly=2.0, n=0, track=False, rate=1):
    data1, data2 = get_frame(0, Ly)
    N = data1.shape[0]
    x = 0.125*np.round(data2[0, :, 0]/0.125)
    y = 0.125*np.round(data2[:, 0, 1]/0.125)
    s = np.sqrt(data2[:, :, 6]*data2[:, :, 6] + data2[:, :, 7]*data2[:, :, 7])

    fig, ax = plt.subplots()
    stream = ax.streamplot(x, y, data2[:, :, 6], data2[:, :, 7], density=4, linewidth=s/s.max())
    # line = ax.plot(data1[:, 0], data1[:, 1], color='black')

    def animate(i):
        ax.collections = [] # clear lines streamplot

        # Clear arrowheads streamplot.
        for artist in ax.get_children():
            if isinstance(artist, FancyArrowPatch):
                artist.remove()

        data1, data2 = get_frame(i, Ly)
        s = np.sqrt(data2[:, :, 6]*data2[:, :, 6] + data2[:, :, 7]*data2[:, :, 7])
        stream = ax.streamplot(x, y, data2[:, :, 6], data2[:, :, 7], density=4, linewidth=s/s.max())
        return stream

    an = anim.FuncAnimation(fig, animate, frames=100, interval=1, blit=True, repeat=False)
    an.save('stream.gif', writer=anim.PillowWriter(fps=100))


def plot_convergence(n=0):
    """
    plots the convergence of the interface to a travelling wave
    """
    if n == 0:
        n = int(len([name for name in os.listdir('./out')
                    if os.path.isfile(os.path.join('./out', name))])/2)

    d = np.zeros([n-1])
    h = np.zeros([n-1])

    data0 = get_interface(0)
    imax = np.argmax(data0[:, 1])
    hmax0 = data0[imax, 1]
    data0[:, 1:2] = np.roll(data0[:, 1:2], -imax, axis=0)
    for k in range(1, n):
        data1 = get_interface(k)

        imax = np.argmax(data1[:, 1])
        hmax1 = data1[imax, 1]
        data1[:, 1:2] = np.roll(data1[:, 1:2], -imax, axis=0)

        # find the difference
        diff = 64/512 * (data1[:, 1]-data0[:, 1])
        d[k-1] = np.sqrt(np.sum(diff*diff))
        h[k-1] = hmax1 - 1.0

        data0 = data1

    plt.plot(h)
    plt.savefig("hmax.png", dpi=300)
    plt.clf()


plot_series(data='interface', n=4, track=True, rate=1)
plot_series(data='vorticity', Ly=2.0, n=4, track=True, rate=1)
plot_series(data='quiver', Ly=2.0, n=4, track=True, rate=1)
plot_series(data='u', Ly=2.0, n=4, track=True, rate=1)
plot_series(data='v', Ly=2.0, n=4, track=True, rate=1)
plot_series(data='level', Ly=2.0, n=4, track=True, rate=1)
# plot_stream(Ly=2.0, n=3, track=True, rate=50)
# plot_convergence(n=4140)
# plot_series(data='level', Ly=4.0)
# plot_series(data='fluid', Ly=2.0)
# plot_series(data='u', Ly=2.0)
# plot_series(data='v', Ly=2.0)
# plot_series(data='speed', Ly=2.0)
# plot_series(data='vorticity', Ly=2.0)
