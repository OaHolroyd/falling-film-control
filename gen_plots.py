#!/usr/bin/env python3
import src.plot as plt

# static plots
plt.plot_hstats()
plt.plot_umax()

# animations
plt.plot_series(key="interface", save=True, Ly=2.0, track=True, clip=False, rate=5)
plt.plot_series(key="control", save=True, Ly=2.0, track=True, clip=False, rate=5)
plt.plot_series(key="fluid", save=True, Ly=2.0, track=True, clip=False, rate=5)
plt.plot_series(key="level", save=True, Ly=8.0, track=True, clip=False, rate=5)
plt.plot_series(key="vorticity", save=True, Ly=4.0, track=True, clip=False, rate=5)
plt.plot_series(key="speed", save=True, Ly=4.0, track=True, clip=False, rate=5)
plt.plot_series(key="u", save=True, Ly=4.0, track=True, clip=False, rate=5)
plt.plot_series(key="v", save=True, Ly=4.0, track=True, clip=False, rate=5)
plt.plot_series(key="pressure", save=True, Ly=4.0, track=True, clip=False, rate=5)
