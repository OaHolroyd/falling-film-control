import numpy as np
import matplotlib.pyplot as plt


def main():
    # Get 1D data
    data = np.loadtxt(f'out/ns-0.dat')
    t = data[:, 0]
    dh = data[:, 1]
    de = data[:, 2]
    dc = data[:, 3]
    c = data[:, 4]


    # Plot 1D data
    fig, ax = plt.subplots()
    ax.semilogy(t, dh)
    ax.semilogy(t, dc)
    ax.semilogy(t, de)
    fig.savefig(f"plots/lines.png")
    plt.close(fig)


    # Plot 2D frames
    # plot dummy data
    fig, ax = plt.subplots()
    data = np.loadtxt(f'out/ns-1-{0:010d}.dat')
    x = data[:, 0]

    hplot, = plt.plot(x, x)
    fplot, = plt.plot(x, x)
    # zplot, = plt.plot(x, x)

    plt.axis([0, 30, -2, 2])


    for i in range(len(t)):
        # Get 2D data for the ith step
        data = np.loadtxt(f'out/ns-1-{i:010d}.dat')
        x = data[:, 0]
        h = data[:, 1]
        f = data[:, 2]
        z = data[:, 3]
        q = data[:, 4]

        # TODO: could be faster if we just change the ydata
        hplot.set_ydata(h-1)
        fplot.set_ydata(f)
        # zplot.set_ydata(z)
        plt.title(f'time {t[i]} [step {i}]')

        fig.savefig(f"plots/{i}.png")

    # Turn plots into a gif
    # TODO


if __name__ == '__main__':
    main()
