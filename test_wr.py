import numpy as np
import matplotlib.pyplot as plt


def der_mat(o, n, dx):
    """Compute derivative matrices"""
    D = np.zeros((n, n))

    def wrap(i):
        # i is the row index (z), j is the col index (x)
        return (i + n) % n

    # x
    if o == 0:
        for i in range(n):
            D[i, i] = 0.0

    # x
    elif o == 1:
        for i in range(n):
            D[i, wrap(i-1)] = -0.5
            D[i, wrap(i+1)] = 0.5

    # xxx
    elif o == 3:
        for i in range(n):
            D[i, wrap(i-2)] = -0.5
            D[i, wrap(i-1)] = 1.0
            D[i, wrap(i+1)] = -1.0
            D[i, wrap(i+2)] = 0.5

    return D / (dx ** o)


def main():
    L = np.loadtxt('out/L.dat')
    L0 = L[:, 0]
    L1 = L[:, 1]
    x = np.arange(len(L0))

    # Plot
    fig, ax = plt.subplots()

    ax.plot(x, L0)
    ax.plot(x, L1)

    plt.show()
    plt.close()


def main1():
    # grid parameters
    n = 128
    l = 30
    dx = l / n

    def wrap(i):
        # i is the row index (z), j is the col index (x)
        return (i + n) % n

    # PDE parameters
    theta = np.pi / 4
    Re = 15.0
    Ca = 0.05
    W = 0.1
    beta = 1.0 / np.tan(theta)

    # create the grid
    x = dx * np.arange(n).reshape((n, 1)) + 0.5 * dx

    # derivative matrices
    I = der_mat(0, n, dx)
    D1 = der_mat(1, n, dx)
    D3 = der_mat(3, n, dx)

    def compute_residual(dt, h, h0, q, q0):
        hx = D1 @ h
        h0x = D1 @ h0
        hxxx = D3 @ h
        h0xxx = D3 @ h0
        qx = D1 @ q
        q0x = D1 @ q0

        res_h = 2.0 * h + dt * qx - 2.0 * h0 + dt * q0x
        res_q = 2.0 * q - 9.0/7.0*dt*q*q/h/h*hx + 5.0*dt*beta/3.0/Re*h*hx + 17.0*dt/7.0*q/h*qx - 5.0*dt/3.0/Re*h - 5.0*dt/6.0/Ca/Re*h*hxxx + 5.0*dt/2.0/Re*q/h/h - 2.0 * q0 - 9.0/7.0*dt*q0*q0/h0/h0*h0x + 5.0*dt*beta/3.0/Re*h0*h0x + 17.0*dt/7.0*q0/h0*q0x - 5.0*dt/3.0/Re*h0 - 5.0*dt/6.0/Ca/Re*h0*h0xxx + 5.0*dt/2.0/Re*q0/h0/h0

        res = np.concatenate((res_h, res_q))
        res_norm_2 = np.sum(res*res)

        return res_norm_2, res

    def compute_jacobian(dt, h, q):
        hx = D1 @ h
        hxxx = D3 @ h
        qx = D1 @ q

        J = np.zeros((2*n, 2*n))

        # top left (dFh/dH)
        for i in range(n):
            J[i][i] = 2.0

        # top right (dFh/dq)
        for i in range(n):
            J[i][n+wrap(i-1)] = dt * (-0.5 / dx)
            J[i][n+wrap(i+1)] = dt * (0.5 / dx)

        # bottom left (dFq/dh)
        for i in range(n):
            c1 = dt * (9.0/7.0*q[i]*q[i]/h[i]/h[i] + 5.0*beta/3.0/Re*h[i]).item()
            c3 = (-dt * 5.0/6.0/Ca/Re*h[i]).item()
            J[n+i][wrap(i-2)] = (-0.5 / dx / dx / dx) * c3
            J[n+i][wrap(i-1)] = (1.0 / dx / dx / dx) * c3 + (-0.5 / dx) * c1
            J[n+i][wrap(i+0)] = dt * (18.0/7.0*q[i]*q[i]/h[i]/h[i]/h[i]*hx[i] + 5.0*beta/3.0/Re*hx[i] - 17.0/7.0*q[i]/h[i]/h[i]*qx[i] - 5.0/3.0/Re - 5.0/6.0/Ca/Re*hxxx[i] - 5.0/Re*q[i]/h[i]/h[i]/h[i]).item()
            J[n+i][wrap(i+1)] = (-1.0 / dx / dx / dx) * c3 + (0.5 / dx) * c1
            J[n+i][wrap(i+2)] = (0.5 / dx / dx / dx) * c3

        # bottom right (dFq/dq)
        for i in range(n):
            c1 = (dt * 17.0/7.0*q[i]/h[i]).item()
            J[n+i][n+wrap(i-1)] = (-0.5 / dx) * c1
            J[n+i][n+wrap(i+0)] = 2.0 + dt * (-18.0/7.0*q[i]/h[i]/h[i]*hx[i] + 17.0/7.0/h[i]*qx[i] + 5.0/2.0/Re/h[i]/h[i]).item()
            J[n+i][n+wrap(i+1)] = (0.5 / dx) * c1

        return J

    # initial condition
    h = 1.0 + 0.01 * np.sin(2.0 * np.pi * x / l)
    h0 = h
    q = 2.0 / 3.0 + 0 * x
    q0 = q

    # Plot 2D frames
    fig, ax = plt.subplots()
    hplot, = plt.plot(x, h)
    plt.axis([0, l, 0, 2])

    # time loop
    iter_max = 100
    res_tol_2 = 1.0e-10
    dt_base = 1.0 / 30.0
    t = 0
    t_end = 2.0
    dt_out = 0.5
    t_out = 0.0
    out_step = 0
    tmax = []
    hmax = []
    while t < t_end:
        # decide dt
        dt = dt_base
        output = False
        if (t + dt) > t_out:
            dt = t_out - t
            output = True
            t_out += dt_out

        # break if h is too small or nan
        if np.any(h < 0.0) or np.any(np.isnan(h)):
            print(f"blow up at t = {t}")
            break

        for k in range(iter_max):
            # compute the residual
            res_norm_2, res = compute_residual(dt, h, h0, q, q0)

            # finish early if converged
            if res_norm_2 < res_tol_2:
                break

            # compute Jacobian and solve linear system
            J = compute_jacobian(dt, h, q)
            dhq = np.linalg.solve(J, res)

            # update variables
            h = h - dhq[0:n]
            q = q - dhq[n:2*n]

        tmax.append(t)
        hmax.append(np.max(h0))
        h0[:] = h[:]
        q0[:] = q[:]
        t += dt

        if output:
            print(f"t = {t}")
            hplot.set_ydata(h)
            plt.title(f'time {t}')
            fig.savefig(f"plots/{out_step}.png")
            out_step += 1

    def growth_rate(k):
        # TODO: rescale k using l
        k = (2*np.pi / l) * k
        lam = -17/12 * 1j - 5/4/Re - 1j * np.sqrt(592*Ca*Ca*Re*Re*k*k + 11760*Ca*Ca*Re*k*k/np.tan(theta) + 5880*Ca*Re*k*k*k*k + 21000 * 1j *Ca*Ca*Re*k - 11025*Ca*Ca)/84/Ca/Re
        return lam
    k0 = np.sqrt(Ca * (8/5*Re - 2/np.tan(theta)))

    # plot growth rate
    fig, ax = plt.subplots()
    plt.semilogy(tmax, hmax)
    plt.xlabel('t')
    plt.ylabel('hmax')
    plt.title('growth rate')
    fig.savefig("plots/growth.png")

    fig, ax = plt.subplots()
    ks = [0.1 * k for k in range(100)]
    lams = [np.real(growth_rate(k)) for k in ks]
    plt.plot(ks, lams)
    plt.scatter(k0, 0.0)
    plt.xlabel('k')
    plt.ylabel('lam')
    plt.title('growth rate')
    fig.savefig("plots/lam.png")


if __name__ == '__main__':
    main()
