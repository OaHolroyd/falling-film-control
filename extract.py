import json
import numpy as np
import matplotlib.pyplot as plt


num = 280
group_data = np.zeros((num, 4))

for i in range(num):
    try:
        # Get paths and file names
        run_dir = f"py-work-0001/run-{i:04d}"
        out_dir = f"{run_dir}/out"
        params_file = f"{run_dir}/params.json"
        data_file = f"{out_dir}/ns-0.dat"

        # Read parameters from file
        with open(params_file, 'r') as fp:
            params = json.load(fp)

        weight = params["CONTROL"]["del"]
        m = params["CONTROL"]["M"]

        # Read data
        data = np.loadtxt(data_file)
        t = data[:, 0]
        dh = data[:, 1]
        de = data[:, 2]
        dc = data[:, 3]
        c = data[:, 4]

        # Compute sum of cost
        dt = t[1]-t[0]
        total_cost = np.sum(dc) * dt

        tmax = t[-1]


        # # Plot 1D data
        # fig, ax = plt.subplots()

        # ax.semilogy(t, dh, label=f"dh {weight} [{total_cost}]")
        # ax.semilogy(t, dc, label=f"dc {weight} [{total_cost}]")

        # plt.legend()
        # plt.title(f"[{i}] m {m}, w {weight}")
        # fig.savefig(f"plots/lines-{i}.png")

        # plt.close(fig)

        if tmax < 50.0:
            total_cost = np.nan


        # Store overall metrics
        group_data[i, :] = [m, weight, total_cost, tmax]

        print(f"{i}: {group_data[i, :]}")

    except:
        print(f"something went wrong for i={i}")


# Plot group data
fig, ax = plt.subplots()
for m in list(set(list(group_data[:, 0]))):
    if m == 0:
        continue
    if m == 9:
        continue

    mask = np.logical_and(group_data[:, 0] == m, group_data[:, 1] >= 0.0)
    m_data = group_data[mask][:, 1:]

    m_data = m_data[m_data[:, 0].argsort()]

    ax.plot(m_data[:, 0], m_data[:, 1], label=f"{m}")


    figm, axm = plt.subplots()
    axm.plot(m_data[:, 0], m_data[:, 1], label=f"{m}")
    plt.legend()
    figm.savefig(f"plots/group_lines_{m}.png")
    plt.close(figm)

# ax.set_ylim(top=5.0)
plt.legend()
fig.savefig(f"plots/group_lines.png")

plt.close(fig)
