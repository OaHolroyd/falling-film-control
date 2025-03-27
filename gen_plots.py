import numpy as np
import matplotlib.pyplot as plt


def main():
    try:
        # plot columns of L
        L = np.loadtxt(f'out/L.dat')
        n = len(L[:, 0]) // 2
        x = np.arange(2*n)
        fig, ax = plt.subplots()
        plt.plot(x, L)
        fig.savefig("plots/L.png")
        plt.close(fig)
    except:
        pass

    try:
        # plot rows of K
        K = np.loadtxt(f'out/K.dat')
        x = np.arange(max(*K.shape))
        fig, ax = plt.subplots()
        plt.plot(x, K.T)
        fig.savefig("plots/K.png")
        plt.close(fig)
    except:
        pass

    try:
        # plot columns of B
        B = np.loadtxt(f'out/Bcf.dat')
        x = np.arange(2*n)
        fig, ax = plt.subplots()
        plt.plot(x, B)
        fig.savefig("plots/B.png")
        plt.close(fig)
    except:
        pass

    try:
        # sparsity patterns of A
        A = np.loadtxt(f'out/Acf.dat')
        Awr = np.loadtxt(f'out/A_wr.dat')
        fig, ax = plt.subplots(1, 2)
        ax[0].spy(A)
        ax[0].plot([-0.5, 2*n-0.5], [n-0.5, n-0.5])
        ax[0].plot([n-0.5, n-0.5], [-0.5, 2*n-0.5])
        ax[1].spy(Awr)
        ax[1].plot([-0.5, 2*n-0.5], [n-0.5, n-0.5])
        ax[1].plot([n-0.5, n-0.5], [-0.5, 2*n-0.5])
        fig.savefig("plots/A.png")
        plt.close(fig)
    except:
        pass


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
    fig.savefig("plots/lines.png")
    plt.close(fig)


    # Plot 2D frames
    # plot dummy data
    fig, ax = plt.subplots()
    data = np.loadtxt(f'out/{pde}-1-{0:010d}.dat')
    x = data[:, 0]
    dx = x[1] - x[0]

    hplot, = plt.plot(x, x)
    fplot, = plt.plot(x, x)
    zplot, = plt.plot(x, x)

    plt.axis([0, 30, 0.95, 1.05])


    for i in range(len(t)):
        # Get 2D data for the ith step
        data = np.loadtxt(f'out/{pde}-1-{i:010d}.dat')
        x = data[:, 0]
        h = data[:, 1]
        f = data[:, 2]
        z = data[:, 3]
        # q = data[:, 4]

        # TODO: could be faster if we just change the ydata
        hplot.set_ydata(h)
        fplot.set_ydata(1+f)
        zplot.set_ydata(z)
        plt.title(f'time {t[i]} [step {i}]')

        fig.savefig(f"plots/{i}.png")

    # Turn plots into a gif
    # TODO


if __name__ == '__main__':
    main()
