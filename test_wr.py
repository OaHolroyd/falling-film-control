import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lu_factor, lu_solve


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


def block_schur_solve(B, C, S, a, b):
    """
    This is a code for the block-solution to

        [      ][   ] [   ]
        [ I  B ][ x ] [ a ]
        [      ][   ]=[   ]
        [ C  D ][ y ] [ b ]
        [      ][   ] [   ]

    With S = D - C B
    """
    # z = C a
    z = C @ a

    # [b, z] = S \ [b, z]
    # S = lu_factor(S)
    # b = lu_solve(S, b)
    # z = lu_solve(S, z)
    b = np.linalg.solve(S, b)
    z = np.linalg.solve(S, z)

    # b = b - z
    b = b - z

    # z = B b
    z = B @ b

    # a = a - z
    a = a - z

    return a, b


def main():
    # create the grid
    xf = DX * np.arange(N).reshape((N, 1))
    xc = xf + 0.5 * DX

    # derivative matrices
    D0c = der_mat_centre(0)
    D1c = der_mat_centre(1)
    D3c = der_mat_centre(3)
    D0l = der_mat_left(0)
    D1l = der_mat_left(1)
    D3l = der_mat_left(3)
    D0r = der_mat_right(0)
    D1r = der_mat_right(1)
    D3r = der_mat_right(3)

    def compute_residual(dt, hc, h0c, qf, q0f):
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

        res = np.concatenate((res_h, res_q))
        res_norm_2 = np.sum(res*res)

        return res_norm_2, res

    def compute_jacobian(dt, hc, qf):
        qxc = D1l @ qf

        hf = D0r @ hc
        hxf = D1r @ hc
        hxxxf = D3r @ hc
        qxf = D1c @ qf

        J = np.zeros((2*N, 2*N))

        # top left (dFh/dH)
        for i in range(N):
            J[i][i] = 1.0

        # top right (dFh/dq)
        for i in range(N):
            c0 = 0.5 * dt
            J[i][N+WRAP(i+0)] = (-1.0 / DX) * c0
            J[i][N+WRAP(i+1)] = (1.0 / DX) * c0

        # bottom left (dFq/dh)
        c0 = 9.0/7.0*dt*qf*qf/hf/hf/hf*hxf + 5.0*dt*BETA/6.0/RE*hxf - 17.0*dt/14.0*qf/hf/hf*qxf - 5.0*dt/6.0/RE - 5.0*dt/12.0/CA/RE*hxxxf - 2.5*dt/RE*qf/hf/hf/hf
        c1 = -9.0/14.0*dt*qf*qf/hf/hf + 5.0*dt*BETA/6.0/RE*hf
        c3 = -5.0*dt/12.0/CA/RE*hf
        for i in range(N):
            J[N+i, WRAP(i-1)] += (0.5) * c0[i].item()
            J[N+i, WRAP(i+0)] += (0.5) * c0[i].item()

            J[N+i, WRAP(i-1)] += (-1.0 / DX) * c1[i].item()
            J[N+i, WRAP(i+0)] += (1.0 / DX) * c1[i].item()

            J[N+i, WRAP(i-2)] = (-1.0 / DX/DX/DX) * c3[i].item()
            J[N+i, WRAP(i-1)] = (3.0 / DX/DX/DX) * c3[i].item()
            J[N+i, WRAP(i+0)] = (-3.0 / DX/DX/DX) * c3[i].item()
            J[N+i, WRAP(i+1)] = (1.0 / DX/DX/DX) * c3[i].item()

        # bottom right (dFq/dq)
        c0 = 1.0 - 9.0/7.0*dt*qf/hf/hf*hxf + 17.0*dt/14.0/hf*qxf + 5.0*dt/RE/hf/hf
        c1 = 17.0*dt/14.0*qf/hf
        for i in range(N):
            J[N+i][N+WRAP(i-1)] = (-0.5 / DX) * c1[i].item()
            J[N+i][N+WRAP(i+0)] = c0[i].item()
            J[N+i][N+WRAP(i+1)] = (0.5 / DX) * c1[i].item()

        return J

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

    def compute_jacobian_blocks(dt, hc, qf):
        qxc = D1l @ qf

        hf = D0r @ hc
        hxf = D1r @ hc
        hxxxf = D3r @ hc
        qxf = D1c @ qf

        A = np.zeros((N, N))
        B = np.zeros((N, N))
        C = np.zeros((N, N))
        D = np.zeros((N, N))

        # top left (dFh/dH)
        for i in range(N):
            A[i][i] = 1.0

        # top right (dFh/dq)
        for i in range(N):
            c0 = 0.5 * dt
            B[i][WRAP(i+0)] = (-1.0 / DX) * c0
            B[i][WRAP(i+1)] = (1.0 / DX) * c0

        # bottom left (dFq/dh)
        c0 = 9.0/7.0*dt*qf*qf/hf/hf/hf*hxf + 5.0*dt*BETA/6.0/RE*hxf - 17.0*dt/14.0*qf/hf/hf*qxf - 5.0*dt/6.0/RE - 5.0*dt/12.0/CA/RE*hxxxf - 2.5*dt/RE*qf/hf/hf/hf
        c1 = -9.0/14.0*dt*qf*qf/hf/hf + 5.0*dt*BETA/6.0/RE*hf
        c3 = -5.0*dt/12.0/CA/RE*hf
        for i in range(N):
            C[i, WRAP(i-1)] += (0.5) * c0[i].item()
            C[i, WRAP(i+0)] += (0.5) * c0[i].item()
            C[i, WRAP(i-1)] += (-1.0 / DX) * c1[i].item()
            C[i, WRAP(i+0)] += (1.0 / DX) * c1[i].item()
            C[i, WRAP(i-2)] = (-1.0 / DX/DX/DX) * c3[i].item()
            C[i, WRAP(i-1)] = (3.0 / DX/DX/DX) * c3[i].item()
            C[i, WRAP(i+0)] = (-3.0 / DX/DX/DX) * c3[i].item()
            C[i, WRAP(i+1)] = (1.0 / DX/DX/DX) * c3[i].item()

        # bottom right (dFq/dq)
        c0 = 1.0 - 9.0/7.0*dt*qf/hf/hf*hxf + 17.0*dt/14.0/hf*qxf + 5.0*dt/RE/hf/hf
        c1 = 17.0*dt/14.0*qf/hf
        for i in range(N):
            D[i][WRAP(i-1)] = (-0.5 / DX) * c1[i].item()
            D[i][WRAP(i+0)] = c0[i].item()
            D[i][WRAP(i+1)] = (0.5 / DX) * c1[i].item()

        return A, B, C, D

    def compute_jacobian_blocks_schur(dt, hc, qf):
        qxc = D1l @ qf

        hf = D0r @ hc
        hxf = D1r @ hc
        hxxxf = D3r @ hc
        qxf = D1c @ qf

        B = np.zeros((N, N))
        C = np.zeros((N, N))
        D = np.zeros((N, N))

        # top right (dFh/dq)
        for i in range(N):
            c0 = 0.5 * dt
            B[i][WRAP(i+0)] = (-1.0 / DX) * c0
            B[i][WRAP(i+1)] = (1.0 / DX) * c0

        # bottom left (dFq/dh)
        c0 = 9.0/7.0*dt*qf*qf/hf/hf/hf*hxf + 5.0*dt*BETA/6.0/RE*hxf - 17.0*dt/14.0*qf/hf/hf*qxf - 5.0*dt/6.0/RE - 5.0*dt/12.0/CA/RE*hxxxf - 2.5*dt/RE*qf/hf/hf/hf
        c1 = -9.0/14.0*dt*qf*qf/hf/hf + 5.0*dt*BETA/6.0/RE*hf
        c3 = -5.0*dt/12.0/CA/RE*hf
        for i in range(N):
            C[i, WRAP(i-1)] += (0.5) * c0[i].item()
            C[i, WRAP(i+0)] += (0.5) * c0[i].item()

            C[i, WRAP(i-1)] += (-1.0 / DX) * c1[i].item()
            C[i, WRAP(i+0)] += (1.0 / DX) * c1[i].item()

            C[i, WRAP(i-2)] += (-1.0 / DX/DX/DX) * c3[i].item()
            C[i, WRAP(i-1)] += (3.0 / DX/DX/DX) * c3[i].item()
            C[i, WRAP(i+0)] += (-3.0 / DX/DX/DX) * c3[i].item()
            C[i, WRAP(i+1)] += (1.0 / DX/DX/DX) * c3[i].item()

        # bottom right (dFq/dq)
        c0 = 1.0 - 9.0/7.0*dt*qf/hf/hf*hxf + 17.0*dt/14.0/hf*qxf + 5.0*dt/RE/hf/hf
        c1 = 17.0*dt/14.0*qf/hf
        for i in range(N):
            D[i][WRAP(i-1)] = (-0.5 / DX) * c1[i].item()
            D[i][WRAP(i+0)] = c0[i].item()
            D[i][WRAP(i+1)] = (0.5 / DX) * c1[i].item()

        S = np.zeros((N, N))
        for i in range(N):
            S[i]

        S = D - C @ B

        return B, C, S

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

            # compute Jacobian and solve linear system
            B, C, S = compute_jacobian_blocks_schur(dt, h, q)

            dh, dq = block_schur_solve(B, C, S, res_h, res_q)
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
