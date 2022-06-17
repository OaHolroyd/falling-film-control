#!/usr/bin/env python3
import os
import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

INDICES = {'x': 0,
           'y': 1,
           'fluid': 2,
           'level': 3,
           'vorticity': 4,
           'speed': 5,
           'u': 6,
           'v': 7,
           'pressure': 8,
           'h': 9}

# COLORMAPS
CMAP_LEVEL = ListedColormap([[0.26, 0.51, 0.92],
                            [0.00, 1.00, 1.00],
                            [0.29, 1.00, 0.45],
                            [0.98, 1.00, 0.00],
                            [1.00, 0.49, 0.00],
                            [1.00, 0.00, 0.00],
                            [0.50, 0.00, 0.00],
                            [0.70, 0.00, 0.29],
                            [1.00, 0.0, 0.93]])
CMAP_SINGLE = LinearSegmentedColormap("", {'red':   [(0.0,  1.00, 1.00),
                                                     (1.0,  0.70, 0.70)],
                                           'green': [(0.0,  1.00, 1.00),
                                                     (1.0,  0.20, 0.20)],
                                           'blue':  [(0.0,  1.00, 1.00),
                                                     (1.0,  0.30, 0.30)]})
CMAP_DOUBLE = LinearSegmentedColormap("", {'red':   [(0.0,  0.20, 0.20),
                                                     (0.5,  1.00, 1.00),
                                                     (1.0,  0.70, 0.70)],
                                           'green': [(0.0,  0.40, 0.40),
                                                     (0.5,  1.00, 1.00),
                                                     (1.0,  0.20, 0.20)],
                                           'blue':  [(0.0,  0.70, 0.70),
                                                     (0.5,  1.00, 1.00),
                                                     (1.0,  0.30, 0.30)]})
CMAP_PHASE = LinearSegmentedColormap("", {'red':   [(0.0,  1.00, 1.00),
                                                    (1.0,  0.20, 0.20)],
                                          'green': [(0.0,  1.00, 1.00),
                                                    (1.0,  0.40, 0.40)],
                                          'blue':  [(0.0,  1.00, 1.00),
                                                    (1.0,  0.70, 0.70)]})


def get_n():
    """
    Returns the largest possible n to use
    """
    return int(len([name for name in os.listdir('./out')
                    if os.path.isfile(os.path.join('./out', name))])/2)


def get_1D_data(i):
    """
    Extracts the 1D data from the ith output. Returns the time and the data
    """
    fname = f"out/data-1-{i:010d}.dat"

    with open(fname, encoding="utf-8") as f:
        first_line = f.readline()
    t = float(first_line[5:])

    data1 = np.loadtxt(fname)
    return t, data1


def get_2D_data(i, Ly=None):
    """
    Extracts the 2D data from the ith output. Returns the time and the data.
    Optionally clips the data to the desired y value
    """
    fname = f"out/data-2-{i:010d}.dat"

    with open(fname, encoding="utf-8") as f:
        first_line = f.readline()
    t = float(first_line[5:])

    data2 = np.loadtxt(fname)

    # reshape data into 3D array
    ny = 0
    x0 = data2[0, 0]
    while data2[ny, 0] == x0:
        ny = ny + 1
    data2 = np.reshape(data2, (ny, -1, 10), order='F')

    if Ly is None:
        return t, data2

    # optionally clip y range to [0,clip]
    # clip y range to [0,Ly]
    iy = 0
    while data2[iy, 0, 1] <= Ly:
        iy = iy + 1
        if iy == data2.shape[0]-1:
            break
    data2 = data2[0:iy, :, :]

    return t, data2


def track_peak(data1, data2=None, p=0.75):
    """
    Shifts the data so that the peak is at p*Lx
    """
    N = data1.shape[0]
    imax = np.argmax(data1[:, 1])
    data1[:, 1:2] = np.roll(data1[:, 1:2], int(N*p)-imax, axis=0)
    if data2 is not None:
        data2[:, :, 2:] = np.roll(data2[:, :, 2:], int(N*p)-imax, axis=1)
        return data1, data2
    return data1


def get_extent(i, data1, data2, clip=False):
    """
    Returns the extent (ie max(abs(X))) for the data, optionally clipped
    """
    v = np.abs(data2[:, :, i]).max()
    if clip:
        v = 0.0
        for k in range(0, data1[:, 0].shape[0]):
            for j in range(0, data2[:, 0, i].shape[0]):
                if data2[j, k, 1] > data1[k, 1]:
                    break
                v = max(v, data2[j, k, i])
    return v


def get_extent_series(i, series, Ly=2.0, clip=False):
    """
    Returns the extent (ie max(abs(X))) for the data, optionally clipped over
    the entire series
    """
    v = 0.0
    for j in series:
        _, data1 = get_1D_data(j)
        _, data2 = get_2D_data(j, Ly)
        v = max(v, get_extent(i, data1, data2, clip=False))
    return v


def get_im(ax, data2, key="fluid", Lx=64.0, Ly=2.0, v=1.0):
    """
    plots the ith data on the given axis
    """
    if key == "fluid":
        # plot fluid and interface
        im = ax.imshow(data2[:, :, 2], interpolation='bilinear',
                       origin='lower', aspect=8.0, cmap=CMAP_PHASE,
                       extent=(0, Lx, 0, Ly))
    elif key == "level":
        # plot level and interface
        im = ax.imshow(data2[:, :, 3], interpolation='nearest',
                       origin='lower', aspect=8.0, cmap=CMAP_LEVEL,
                       extent=(0, Lx, 0, Ly), vmin=3.5, vmax=12.5)
    elif key == "vorticity":
        # plot vorticity and interface
        im = ax.imshow(data2[:, :, 4], interpolation='bilinear',
                       origin='lower', aspect=8.0, cmap=CMAP_DOUBLE,
                       extent=(0, Lx, 0, Ly), vmin=-v, vmax=v)
    elif key == "speed":
        # plot speed and interface
        im = ax.imshow(data2[:, :, 5], interpolation='bilinear',
                       origin='lower', aspect=8.0, cmap=CMAP_SINGLE,
                       extent=(0, Lx, 0, Ly), vmax=v)
    elif key == "u":
        # plot u and interface
        im = ax.imshow(data2[:, :, 6], interpolation='bilinear',
                       origin='lower', aspect=8.0, cmap=CMAP_DOUBLE,
                       extent=(0, Lx, 0, Ly), vmin=-v, vmax=v)
    elif key == "v":
        # plot v and interface
        im = ax.imshow(data2[:, :, 7], interpolation='bilinear',
                       origin='lower', aspect=8.0, cmap=CMAP_DOUBLE,
                       extent=(0, Lx, 0, Ly), vmin=-v, vmax=v)
    elif key == "pressure":
        # plot pressure and interface
        im = ax.imshow(data2[:, :, 8], interpolation='bilinear',
                       origin='lower', aspect=8.0, cmap=CMAP_DOUBLE,
                       extent=(0, Lx, 0, Ly), vmin=-v, vmax=v)
    return im


def plot_frame(i, key="interface", save=True, Ly=2.0, track=False, clip=False):
    """
    Plots the ith frame with the following possible keys:
      interface [default], control, fluid, level, vorticity, speed, u, v,
      or pressure
    and either saves [default] or shows the plot
    """
    t, data1 = get_1D_data(i)

    Lx = data1[:, 0].max()

    # get 2D data if required
    if key not in ("interface", "control"):
        _, data2 = get_2D_data(i, Ly)

    # track if required
    if track:
        if key not in ("interface", "control"):
            data1, data2 = track_peak(data1, data2)
        else:
            data1 = track_peak(data1)

    # plot
    fig, ax = plt.subplots()
    if key == "interface":
        ax.plot(data1[:, 0], data1[:, 1], color="blue")

        # keys and axes
        ax.set_xlim([0, Lx])
        ax.set_xlim([0, Ly])
        ax.gca().set_aspect(8.0)
    elif key == "control":
        ax.plot(data1[:, 0], data1[:, 1], color="blue")
        ax.plot(data1[:, 0], data1[:, 2], color="red")

        # keys and axes
        ax.set_xlim([0, Lx])
        ax.set_ylim([0, Ly])
        ax.gca().set_aspect(8.0)
    else:  # 2D plots
        if key in ("fluid", "level"):  # with a fixed upper/lower bound
            im = get_im(ax, data2, key=key, Lx=Lx, Ly=Ly)
        else:  # with an extent
            v = get_extent(INDICES[key], data1, data2, clip=clip)
            im = get_im(ax, data2, key=key, Lx=Lx, Ly=Ly, v=v)
        if clip:
            ax.fill_between(data1[:, 0], data1[:, 1], Ly, color="white")
        ax.plot(data1[:, 0], data1[:, 1], color="black")  # interface

        # keys and axes
        fig.colorbar(im, location='bottom')
        ax.set_xlim([0, Lx])
        ax.set_ylim([0, Ly])

    ax.set_title(f"{key} (t = {t})", weight='bold', fontsize=14)
    fig.add_axes(ax)

    # save or display
    if save:
        fig.savefig(f"plots/{key}-{i:010d}.png", dpi=300)
    else:
        plt.show()


def plot_series(n=None, key="interface", save=True, Ly=2.0, track=False,
                clip=False, rate=1):
    """
    Plots an animation with the following possible keys:
      interface [default], control, fluid, level, vorticity, speed, u, v,
      or pressure
    and either saves [default] or shows the animation
    """
    # get max value of n if it is not supplied
    if n is None:
        n = get_n()

    t, data1 = get_1D_data(0)

    Lx = data1[:, 0].max()
    N = data1[:, 0].shape[0]

    # get 2D data if required
    if key not in ("interface", "control"):
        _, data2 = get_2D_data(0, Ly)

    # track if required
    if track:
        if key not in ("interface", "control"):
            data1, data2 = track_peak(data1, data2)
        else:
            data1 = track_peak(data1)

    # plot
    fig, ax = plt.subplots()
    if key == "interface":
        l1, = ax.plot(data1[:, 0], data1[:, 1], color="blue")

        # keys and axes
        ax.set_xlim([0, Lx])
        ax.set_ylim([0, Ly])
        ax.set_aspect(8.0)
    elif key == "control":
        l1, = ax.plot(data1[:, 0], data1[:, 1], color="blue")
        l2, = ax.plot(data1[:, 0], data1[:, 2], color="red")

        # keys and axes
        ax.set_xlim([0, Lx])
        ax.set_ylim([0, Ly])
        ax.set_aspect(8.0)
    else:  # 2D plots
        if key in ("fluid", "level"):  # with a fixed upper/lower bound
            im = get_im(ax, data2, key=key, Lx=Lx, Ly=Ly)
        else:  # with an extent
            v = get_extent_series(INDICES[key], range(0, n, rate),
                                  Ly=2.0, clip=False)
            im = get_im(ax, data2, key=key, Lx=Lx, Ly=Ly, v=v)
        if clip:
            f1 = ax.fill_between(data1[:, 0], data1[:, 1], Ly, color="white")
        l1, = ax.plot(data1[:, 0], data1[:, 1], color="black")  # interface

        # keys and axes
        fig.colorbar(im, location='bottom')
        ax.set_xlim([0, Lx])
        ax.set_ylim([0, Ly])

    ax.set_title(f"{key} (t = {t})", weight='bold', fontsize=14)
    fig.add_axes(ax)

    def step(i):
        """
        Animation stepper function
        """
        t, data1 = get_1D_data(i)
        ax.set_title(f"{key} (t = {t})", weight='bold', fontsize=14)

        # get 2D data if required
        if key not in ("interface", "control"):
            _, data2 = get_2D_data(i, Ly)

        # track if required
        if track:
            if key not in ("interface", "control"):
                data1, data2 = track_peak(data1, data2)
            else:
                data1 = track_peak(data1)

        # plot
        if key == "interface":
            l1.set_data(data1[:, 0], data1[:, 1])
            return (l1)
        if key == "control":
            l1.set_data(data1[:, 0], data1[:, 1])
            l2.set_data(data1[:, 0], data1[:, 2])
            return (l1, l2)
        # 2D plots
        im.set_data(data2[:, :, INDICES[key]])
        if clip:
            path = f1.get_paths()[0]
            path.vertices[1:N+1, 1] = data1[:, 1]
        l1.set_data(data1[:, 0], data1[:, 1])  # interface
        return (im, l1)
    ani = anim.FuncAnimation(fig, step, frames=range(0, n, rate))

    if save:
        ani.save(f'plots/{key}.gif', writer=anim.PillowWriter(fps=100))
    else:
        plt.show()


def plot_hstats(n=None, rate=1, save=True):
    """
    Plots statistics about the interface (hmin, hmax, mean deviation)
    """
    if n is None:
        n = get_n()

    _, data1 = get_1D_data(0)
    dx = data1[1, 0] - data1[0, 0]
    L = data1[:, 0].max()

    # get data
    t = np.zeros(n)
    hmin = np.zeros(n)
    hmax = np.zeros(n)
    dh = np.zeros(n)
    for i in range(0, n, rate):
        t[i], data1 = get_1D_data(i)
        hmin[i] = 1.0 - data1[:, 1].min()
        hmax[i] = data1[:, 1].max() - 1.0
        dh[i] = dx * sum((data1[:, 1] - 1.0)**2) / L

    # plot data
    plt.plot(t, hmin, label="hmin")
    plt.plot(t, hmax, label="hmax")
    plt.plot(t, dh, label="dh")
    plt.gca().set_ylim(bottom=0)
    plt.xlabel("t")
    plt.ylabel("deviation")
    plt.legend()

    # save or display
    if save:
        plt.savefig("plots/hstats.png", dpi=300, bbox_inches='tight')
    else:
        plt.show()
    plt.clf()


def plot_umax(n=None, rate=1, save=True):
    """
    Plots umax vs t
    """
    if n is None:
        n = get_n()

    # get data
    t = np.zeros(n)
    umax = np.zeros(n)
    for i in range(0, n, rate):
        t[i], data1 = get_1D_data(i)
        _, data2 = get_2D_data(i, Ly=5.0)
        umax[i] = get_extent(INDICES["u"], data1, data2, clip=True)

    # plot data
    plt.plot(t, umax)
    plt.gca().set_ylim(bottom=0)
    plt.xlabel("t")
    plt.ylabel("umax")

    # save or display
    if save:
        plt.savefig("plots/umax.png", dpi=300, bbox_inches='tight')
    else:
        plt.show()
    plt.clf()
