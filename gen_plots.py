import numpy as np
import matplotlib.pyplot as plt


def main():
    pde = 'ns'

    # Get 1D data
    data = np.loadtxt(f'out/{pde}-0.dat')
    t = data[:, 0]
    dh = data[:, 1]
    de = data[:, 2]
    dc = data[:, 3]
    c = data[:, 4]


    # Plot 1D data
    fig, ax = plt.subplots()
    ax.semilogy(t, dh)
    ax.semilogy(t, de)
    fig.savefig(f"plots/lines.png")
    plt.close(fig)


    # Plot 2D frames
    # plot dummy data
    fig, ax = plt.subplots()
    data = np.loadtxt(f'out/{pde}-1-{0:010d}.dat')
    x = data[:, 0]
    dx = x[1] - x[0]

    hplot, = plt.plot(x, x)
    # fplot, = plt.plot(x, x)
    zplot, = plt.plot(x, x)

    plt.axis([0, 30, 0.75, 1.25])


    for i in range(len(t)):
        # Get 2D data for the ith step
        data = np.loadtxt(f'out/{pde}-1-{i:010d}.dat')
        x = data[:, 0]
        h = data[:, 1]
        # f = data[:, 2]
        z = data[:, 3]
        # q = data[:, 4]

        de[i] = np.sqrt(np.sum((h-z)*(h-z))*dx)

        # TODO: could be faster if we just change the ydata
        hplot.set_ydata(h)
        # fplot.set_ydata(f)
        zplot.set_ydata(z)
        plt.title(f'time {t[i]} [step {i}]')

        fig.savefig(f"plots/{i}.png")

    # Plot 1D data
    fig, ax = plt.subplots()
    ax.semilogy(t, dh)
    ax.semilogy(t, de)
    fig.savefig(f"plots/lines2.png")
    plt.close(fig)

    # Turn plots into a gif
    # TODO


if __name__ == '__main__':
    main()
