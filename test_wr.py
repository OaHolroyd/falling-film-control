import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_array


def _pentadiagonal_lu_factorise_nonperiodic(
        a: np.ndarray,
        b: np.ndarray,
        c: np.ndarray,
        d: np.ndarray,
        e: np.ndarray,
        n: int
) -> None:
    """
    Factorise a pentadiagonal matrix A = LU, where L is lower triangular with three bands and U is upper unit triangular
    with three bands.

    See section D, 13a and 13b in:
        https://arxiv.org/pdf/1807.07382
    for more details.

    Args:
        a: second lower diagonal band (OVERWRITTEN with the main diagonal band of L)
        b: first lower diagonal band (OVERWRITTEN with the first lower diagonal band of L)
        c: main diagonal band (OVERWRITTEN with the second lower diagonal band of L)
        d: first upper diagonal band (OVERWRITTEN with first upper diagonal band of U)
        e: second upper diagonal band (OVERWRITTEN with second upper diagonal band of U)
        n: size of the matrix (square)
    """
    # Compute the entries of L and U
    d[0] /= c[0]
    e[0] /= c[0]

    c[1] -= b[1] * d[0]
    d[1] = (d[1] - b[1] * e[0]) / c[1]
    e[1] /= c[1]

    for i in range(2, n - 2):
        b[i] -= a[i] * d[i - 2]
        c[i] -= a[i] * e[i - 2] + b[i] * d[i - 1]
        d[i] = (d[i] - b[i] * e[i - 1]) / c[i]
        e[i] /= c[i]

    b[n - 2] -= a[n - 2] * d[n - 4]
    c[n - 2] -= a[n - 2] * e[n - 4] + b[n - 2] * d[n - 3]
    d[n - 2] = (d[n - 2] - b[n - 2] * e[n - 3]) / c[n - 2]

    b[n - 1] -= a[n - 1] * d[n - 3]
    c[n - 1] -= a[n - 1] * e[n - 3] + b[n - 1] * d[n - 2]


def _pentadiagonal_lu_solve_nonperiodic(
        al: np.ndarray,
        be: np.ndarray,
        ep: np.ndarray,
        ga: np.ndarray,
        de: np.ndarray,
        f: np.ndarray,
        n: int
) -> None:
    """
    Solve a pentadiagonal system of equations that have been factorised as A = LU.

    Args:
        al: main diagonal band of L
        be: first lower diagonal band of L
        ep: second lower diagonal band of L
        ga: first upper diagonal band of U
        de: second upper diagonal band of U
        f: right-hand side (OVERWRITTEN with the solution on output)
        n: size of the matrix
    """
    # solve Ly = f via forward substitution
    y = f  # reuse the f array
    y[0] /= al[0]
    y[1] = (y[1] - be[1] * y[0]) / al[1]
    for i in range(2, n):
        y[i] = (y[i] - be[i] * y[i - 1] - ep[i] * y[i - 2]) / al[i]

    # solve Ux = y via backward substitution
    x = y  # reuse the y array (which is the f array)
    x[n - 2] -= ga[n - 2] * x[n - 1]
    for i in range(n - 3, -1, -1):
        x[i] -= ga[i] * x[i + 1] + de[i] * x[i + 2]


def pentadiagonal_solve_nonperiodic(
        a: np.ndarray,
        b: np.ndarray,
        c: np.ndarray,
        d: np.ndarray,
        e: np.ndarray,
        f: np.ndarray,
        n: int
) -> None:
    """
    Solve a pentadiagonal system of equations.

    This function solves Ax = f for A = pentadiagonal(a, b, c, d, e, periodic=False).

    There are three steps:
      1. Factor A = LU, where L is lower triangular with three bands and U is upper unit triangular with three bands.
      2. Solve Ly = f.
      3. Solve Ux = y.

    This method can be found in:
      Gisela Engeln-Müllges and Frank Uhlig. Numerical Algorithms with C. Springer-Verlag, Berlin, Heidelberg, 1996.
    and also in section D of:
      https://arxiv.org/pdf/1807.07382

    Args:
        a: second lower diagonal band
        b: first lower diagonal band
        c: main diagonal band
        d: first upper diagonal band
        e: second upper diagonal band
        f: right-hand side (OVERWRITTEN with the solution to the system of equations)
        n: size of the matrix (square)
    """
    # compute the LU factorisation of A
    # NOTE: in a real implementation, we overwrite a, b, c, d, e with the LU factors to save memory, and so they can be
    #       reused in further solves.
    _pentadiagonal_lu_factorise_nonperiodic(a, b, c, d, e, n)

    # solve the system of equations
    _pentadiagonal_lu_solve_nonperiodic(c, b, a, d, e, f, n)


def _pentadiagonal_lu_factorise_periodic(
        a: np.ndarray,
        b: np.ndarray,
        c: np.ndarray,
        d: np.ndarray,
        e: np.ndarray,
        k0: np.ndarray,
        k1: np.ndarray,
        n: int,
) -> None:
    """
    Compute the subsystem LU factorisation and prepares a periodic pentadiagonal matrix for solving.

    See `pentadiagonal_solve_periodic` for notation.

    Args:
        a: second lower diagonal band (OVERWRITTEN with the main diagonal band of E)
        b: first lower diagonal band (OVERWRITTEN with the firs lower diagonal band of E)
        c: main diagonal band (OVERWRITTEN with the second lower diagonal band of E)
        d: first upper diagonal band (OVERWRITTEN with the first upper diagonal band of E)
        e: second upper diagonal band (OVERWRITTEN with the second upper diagonal band of E)
        k0: OVERWRITTEN with the first column of E^-1 K
        k1: OVERWRITTEN with the second column of E^-1 K
        n: size of the matrix (square)
    """
    # set K = [k0 | k1]
    k0[0] = a[0]
    for i in range(1, n - 4):
        k0[i] = 0.0
    k0[n - 4] = e[n - 4]
    k0[n - 3] = d[n - 3]

    k1[0] = b[0]
    k1[1] = a[1]
    for i in range(2, n - 3):
        k1[i] = 0.0
    k1[n - 3] = e[n - 3]

    # compute the LU factorisation of E
    _pentadiagonal_lu_factorise_nonperiodic(a, b, c, d, e, n - 2)

    # solve E \ K
    _pentadiagonal_lu_solve_nonperiodic(c, b, a, d, e, k0, n - 2)
    _pentadiagonal_lu_solve_nonperiodic(c, b, a, d, e, k1, n - 2)

    # compute the 2x2 matrix C - H E^-1 K for eqn 12
    c[n - 2] -= e[n - 2] * k0[0] + a[n - 2] * k0[n - 4] + b[n - 2] * k0[n - 3]
    d[n - 2] -= e[n - 2] * k1[0] + a[n - 2] * k1[n - 4] + b[n - 2] * k1[n - 3]
    b[n - 1] -= d[n - 1] * k0[0] + e[n - 1] * k0[1] + a[n - 1] * k0[n - 3]
    c[n - 1] -= d[n - 1] * k1[0] + e[n - 1] * k1[1] + a[n - 1] * k1[n - 3]


def _pentadiagonal_lu_solve_periodic(
        a: np.ndarray,
        b: np.ndarray,
        c: np.ndarray,
        d: np.ndarray,
        e: np.ndarray,
        k0: np.ndarray,
        k1: np.ndarray,
        f: np.ndarray,
        n: int
) -> None:
    """
    Starting with a system proccesed by `_pentadiagonal_lu_factorise_periodic`, solve Ax=f.

    See `pentadiagonal_solve_periodic` for notation.

    Args:
        a: the main diagonal band of E
        b: the first lower diagonal band of E
        c: the second lower diagonal band of E
        d: the first upper diagonal band of E
        e: the second upper diagonal band of E
        k0: the first column of E^-1 K
        k1: the second column of E^-1 K
        n: size of the matrix (square)
    """
    # solve E \ f[:-2]
    _pentadiagonal_lu_solve_nonperiodic(c, b, a, d, e, f, n - 2)

    # compute rhs vector of eqn 12
    f[n - 2] -= e[n - 2] * f[0] + a[n - 2] * f[n - 4] + b[n - 2] * f[n - 3]
    f[n - 1] -= d[n - 1] * f[0] + e[n - 1] * f[1] + a[n - 1] * f[n - 3]

    # solve for the final two elements of the solution vector (eq 12)
    det = c[n - 2] * c[n - 1] - d[n - 2] * b[n - 1]
    tmp = (c[n - 1] * f[n - 2] - d[n - 2] * f[n - 1]) / det
    f[n - 1] = (c[n - 2] * f[n - 1] - b[n - 1] * f[n - 2]) / det
    f[n - 2] = tmp

    # update the rhs vector for the final solve (eq 11)
    for i in range(n - 2):
        f[i] -= k0[i] * f[n - 2] + k1[i] * f[n - 1]


def pentadiagonal_solve_periodic(
        a: np.ndarray,
        b: np.ndarray,
        c: np.ndarray,
        d: np.ndarray,
        e: np.ndarray,
        f: np.ndarray,
        n: int
) -> None:
    """
    Solve a periodic pentadiagonal system of equations.

    This function solves Ax = f for A = pentadiagonal(a, b, c, d, e, periodic=True).

    The method is adapted from
        https://people.sc.fsu.edu/~inavon/pubs1/navon1987.pdf
    and
        https://arxiv.org/pdf/1807.07382
    to handle non-Toeplitz matrices. Equation numbers mentioned in the code refer to the second paper.

    The input diagonals wrap row-wise, so that for n = 8 we have
        [ c0 d0 e0  0  0  0 a0 b0 ]
        [ b1 c1 d1 e1  0  0  0 a1 ]
        [ a2 b2 c2 d2 e2  0  0  0 ]
        [  0 a3 b3 c3 d3 e3  0  0 ]
        [  0  0 a4 b4 c4 d4 e4  0 ]
        [  0  0  0 a5 b5 c5 d5 e5 ]
        [ e6  0  0  0 a6 b6 c6 d6 ]
        [ d7 e7  0  0  0 a7 b7 c7 ]

    Within the algorithm the matrix is split into blocks:
        [ E  K ]
        [ H  C ]
    where
      - E is the (n-2) x (n-2) nonperiodic pentadiagonal submatrix in the upper left corner of A
      - H is the 2 x (n-2) matrix formed by the start of the last two rows of A
      - K is the (n-2) x 2 matrix formed by the start of the last two columns of A
      - C is the 2 x 2 matrix formed by the bottom right corner of A

    Args:
        a: second lower diagonal band
        b: first lower diagonal band
        c: main diagonal band
        d: first upper diagonal band
        e: second upper diagonal band
        f: right-hand side (OVERWRITTEN with the solution to the system of equations)
        n: size of the matrix (square)
    """
    # factorise and prepare K
    k0 = np.zeros((n - 2,))
    k1 = np.zeros((n - 2,))
    _pentadiagonal_lu_factorise_periodic(a, b, c, d, e, k0, k1, n)

    # solve the factorised system
    _pentadiagonal_lu_solve_periodic(a, b, c, d, e, k0, k1, f, n)


# grid parameters
N = 128
L = 30.0
DX = L / N

# PDE parameters
THETA = np.pi / 4.0
RE = 15.0
CA = 0.05
W = 0.1
BETA = 1.0 / np.tan(THETA)


def WRAP(i):
    """Periodically wrap an index"""
    return (i + N) % N


def der_mat_centre(order):
    """
    Compute a derivative matrix of the desired order.

    This evaluates the derivative at a point from a vector of relatively
    cell-centred values.
    """
    D = np.zeros((N, N))

    # I
    if order == 0:
        for i in range(N):
            D[i, i] = 1.0

    # Dx
    elif order == 1:
        for i in range(N):
            D[i, WRAP(i-1)] = -0.5
            D[i, WRAP(i+1)] = 0.5

    # Dxxx
    elif order == 3:
        for i in range(N):
            D[i, WRAP(i-2)] = -0.5
            D[i, WRAP(i-1)] = 1.0
            D[i, WRAP(i+1)] = -1.0
            D[i, WRAP(i+2)] = 0.5

    return D / (DX ** order)


def der_mat_left(order):
    """
    Compute a derivative matrix of the desired order.

    This evaluates the derivative at a point from a vector of relatively
    face-centred values, where the ith face is the left face of the ith cell.
    """
    D = np.zeros((N, N))

    # I
    if order == 0:
        for i in range(N):
            D[i, i] = 0.5
            D[i, WRAP(i+1)] = 0.5

    # Dx
    elif order == 1:
        for i in range(N):
            D[i, i] = -1.0
            D[i, WRAP(i+1)] = 1.0

    # Dxxx
    elif order == 3:
        for i in range(N):
            D[i, WRAP(i-1)] = -1.0
            D[i, i] = 3.0
            D[i, WRAP(i+1)] = -3.0
            D[i, WRAP(i+2)] = 1.0

    return D / (DX ** order)


def der_mat_right(order):
    """
    Compute a derivative matrix of the desired order.

    This evaluates the derivative at a point from a vector of relatively
    face-centred values, where the ith face is the right face of the ith cell.
    """
    D = np.zeros((N, N))

    # I
    if order == 0:
        for i in range(N):
            D[i, WRAP(i-1)] = 0.5
            D[i, i] = 0.5

    # Dx
    elif order == 1:
        for i in range(N):
            D[i, WRAP(i-1)] = -1.0
            D[i, i] = 1.0

    # Dxxx
    elif order == 3:
        for i in range(N):
            D[i, WRAP(i-2)] = -1.0
            D[i, WRAP(i-1)] = 3.0
            D[i, i] = -3.0
            D[i, WRAP(i+1)] = 1.0

    return D / (DX ** order)


def main():
    # create the grid
    xf = DX * np.arange(N).reshape((N,))
    xc = xf + 0.5 * DX

    # derivative matrices
    D1l = coo_array(der_mat_left(1))
    D1c = coo_array(der_mat_centre(1))
    D0r = coo_array(der_mat_right(0))
    D1r = coo_array(der_mat_right(1))
    D3r = coo_array(der_mat_right(3))

    def compute_residual_blocks(dt, hc, h0c, qf, q0f):
        qxc = D1l @ qf
        q0xc = D1l @ q0f

        hf = D0r @ hc
        h0f = D0r @ hc
        hxf = D1r @ hc
        h0xf = D1r @ hc
        hxxxf = D3r @ hc
        h0xxxf = D3r @ hc
        qxf = D1c @ qf
        q0xf = D1c @ q0f

        res_h = hc + 0.5 * dt * qxc - h0c + 0.5 * dt * q0xc
        res_q = qf \
              - 9.0/14.0*dt*qf*qf/hf/hf*hxf \
              + 5.0*dt*BETA/6.0/RE*hf*hxf \
              + 17.0*dt/14.0*qf/hf*qxf \
              - 5.0*dt/6.0/RE*hf \
              - 5.0*dt/12.0/CA/RE*hf*hxxxf \
              + 5.0*dt/4.0/RE*qf/hf/hf \
              - q0f \
              - 9.0/14.0*dt*q0f*q0f/h0f/h0f*h0xf \
              + 5.0*dt*BETA/6.0/RE*h0f*h0xf \
              + 17.0*dt/14.0*q0f/h0f*q0xf \
              - 5.0*dt/6.0/RE*h0f \
              - 5.0*dt/12.0/CA/RE*h0f*h0xxxf \
              + 5.0*dt/4.0/RE*q0f/h0f/h0f

        res_norm_2 = np.sum(res_h*res_h) + np.sum(res_q*res_q)

        return res_norm_2, res_h, res_q

    def fast_solve(dt, hc, qf, a, b):
        hf = D0r @ hc
        hxf = D1r @ hc
        hxxxf = D3r @ hc
        qxf = D1c @ qf

        # C interim coefficients
        c0 = 9.0/7.0*dt*qf*qf/hf/hf/hf*hxf + 5.0*dt*BETA/6.0/RE*hxf - 17.0*dt/14.0*qf/hf/hf*qxf - 5.0*dt/6.0/RE - 5.0*dt/12.0/CA/RE*hxxxf - 2.5*dt/RE*qf/hf/hf/hf
        c1 = -9.0/14.0*dt*qf*qf/hf/hf + 5.0*dt*BETA/6.0/RE*hf
        c3 = -5.0*dt/12.0/CA/RE*hf

        # construct the diagonals of C
        c_l2 = np.zeros((N,))
        c_l1 = np.zeros((N,))
        c_d0 = np.zeros((N,))
        c_u1 = np.zeros((N,))
        for i in range(N):
            c_l2[i] = (-1.0 / DX/DX/DX) * c3[i]
            c_l1[i] = (3.0 / DX/DX/DX) * c3[i] + (-1.0 / DX) * c1[i] + (0.5) * c0[i]
            c_d0[i] = (-3.0 / DX/DX/DX) * c3[i] + (1.0 / DX) * c1[i] + (0.5) * c0[i]
            c_u1[i] = (1.0 / DX/DX/DX) * c3[i]

        # Schur complement, D - C A \ B = D - C B
        d0 = 1.0 - 9.0/7.0*dt*qf/hf/hf*hxf + 17.0*dt/14.0/hf*qxf + 5.0*dt/RE/hf/hf
        d1 = 17.0*dt/14.0*qf/hf

        s_l2 = np.zeros((N,))
        s_l1 = np.zeros((N,))
        s_d0 = np.zeros((N,))
        s_u1 = np.zeros((N,))
        s_u2 = np.zeros((N,))
        for i in range(N):
            s_l1[i] = (-0.5 / DX) * d1[i]
            s_d0[i] = d0[i]
            s_u1[i] = (0.5 / DX) * d1[i]

            s_l2[i] += c_l2[i] * 0.5 * dt / DX
            s_l1[i] += (c_l1[i] - c_l2[i]) * 0.5 * dt / DX
            s_d0[i] += (c_d0[i] - c_l1[i]) * 0.5 * dt / DX
            s_u1[i] += (c_u1[i] - c_d0[i]) * 0.5 * dt / DX
            s_u2[i] += -c_u1[i] * 0.5 * dt / DX

        # z = C a
        z = np.zeros(a.shape)
        for i in range(N):
            z[i] += a[WRAP(i-2)] * c_l2[i]
            z[i] += a[WRAP(i-1)] * c_l1[i]
            z[i] += a[WRAP(i+0)] * c_d0[i]
            z[i] += a[WRAP(i+1)] * c_u1[i]

        # [b, z] = S \ [b, z]
        k0 = np.zeros((N,))
        k1 = np.zeros((N,))
        _pentadiagonal_lu_factorise_periodic(s_l2, s_l1, s_d0, s_u1, s_u2, k0, k1, N)
        _pentadiagonal_lu_solve_periodic(s_l2, s_l1, s_d0, s_u1, s_u2, k0, k1, b, N)
        _pentadiagonal_lu_solve_periodic(s_l2, s_l1, s_d0, s_u1, s_u2, k0, k1, z, N)

        # b = b - z
        b = b - z

        # a = a - B b
        for i in range(N):
            a[i] += (b[i] - b[WRAP(i+1)]) * 0.5 * dt / DX

        return a, b

    # initial condition
    h = 1.0 + 0.01 * np.sin(2.0 * np.pi * xc / L)
    h0 = h
    q = 2.0 / 3.0 + 0 * xf
    q0 = q

    # Plot 2D frames
    fig, _ = plt.subplots()
    hplot, = plt.plot(xc, h)
    qplot, = plt.plot(xf, q)
    plt.axis([0, L, 0, 2])

    # time loop
    iter_max = 100
    res_tol_2 = 1.0e-10
    dt_base = 1.0 / 30.0
    t = 0
    t_end = 100.0
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
            # res_norm_2, res = compute_residual(dt, h, h0, q, q0)
            res_norm_2, res_h, res_q = compute_residual_blocks(dt, h, h0, q, q0)


            # finish early if converged
            if res_norm_2 < res_tol_2:
                break

            # update h and q
            # dh, dq = block_solve(dt, h, q, res_h, res_q)
            dh, dq = fast_solve(dt, h, q, res_h, res_q)
            h = h - dh
            q = q - dq

        tmax.append(t)
        hmax.append(np.max(h0))
        h0[:] = h[:]
        q0[:] = q[:]
        t += dt

        if output:
            print(f"t = {t} [{k}]")
            hplot.set_ydata(h)
            qplot.set_ydata(q)
            plt.title(f'time {t}')
            fig.savefig(f"plots/{out_step}.png")
            out_step += 1


if __name__ == '__main__':
    main()
