import shutil
import subprocess
from contextlib import ExitStack
from pathlib import Path

import numpy as np
from PIL import Image
from matplotlib import pyplot as plt

from scripts.config import Config

EXE = "./film-ns"


# define a custom error for when a simulation cannot be run
class SimulationError(Exception):
    pass


def rm_tree(path: Path):
    """
    Recursively remove a directory and all its contents.
    """
    if path.is_dir():
        for child in path.iterdir():
            rm_tree(child)
        path.rmdir()
    elif path.is_file():
        path.unlink()


class Simulation:
    def __init__(self, config: Config | Path | str, base_dir: Path | str = ".", exe: Path | str | None = None):
        """
        Initialize the simulation with a configuration object.

        The simulation is housed in a directory:
            base_dir/
              out/
                <various output files>
              plots/
                <various plots>
              film-ns
              params.json
              output.txt

        Args:
            config: Configuration object containing simulation parameters or a path to a JSON file.
            base_dir: Base directory for the simulation files (will be created if it doesn't exist).
            exe: Path to the executable. If None, defaults to 'film-ns' in the base_dir.
        """
        if not isinstance(config, Config):
            config = Config.from_json(config)
        self.config = config
        self.base_dir = Path(base_dir)

        self.base_dir.mkdir(parents=True, exist_ok=True)

        self.exe = self.base_dir / EXE
        if exe is not None:
            # need to copy the executable to the base_dir
            shutil.copy(exe, self.exe)

        self.output_dir = self.base_dir / "out"
        self.plots_dir = self.base_dir / "plots"
        self.dump_dir = self.base_dir / "dump"
        self.params = self.base_dir / "params.json"
        self.output = self.base_dir / "output.txt"

        # check that exe is a valid path
        if not self.exe.exists():
            raise ValueError(f"Executable '{self.exe}' does not exist.")

        # create the output and plots directories if they do not exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.dump_dir.mkdir(parents=True, exist_ok=True)

    @property
    def has_completed(self):
        """
        Check if the simulation has completed successfully.
        """
        # If it has completed successfully, the final line of the log file should start with a #
        try:
            with open(self.output, "r") as fp:
                lines = fp.readlines()
        except FileNotFoundError:
            # No output file found, cannot have started
            return False

        return lines[-1][0] == '#'

    @property
    def has_converged(self):
        """
        Check if the simulation has converged successfully.
        """
        # Cannot have converged if it has not completed
        if not self.has_completed:
            return False

        # Check if the final film thickness is within the tolerance
        data = np.loadtxt(self.output_dir / 'ns-0.dat')
        dh = data[-1, 1]  # final interfacial deviation, ||h - 1||2
        tol  = 1e-3  # tolerance for convergence

        return dh < tol

    def dump_config(self):
        """
        Dump the configuration to a JSON file.
        """
        self.config.to_json(self.params)

    def run(self, force: bool = False, timeout: int = 21600):
        """
        Run the simulation.

        Args:
            force: If True, overwrite existing output and plots directories.
            timeout: Timeout for the simulation in seconds (default is 6 hours).
        """
        # if force run, clear everything before we start
        if force:
            rm_tree(self.output_dir)
            rm_tree(self.plots_dir)

            rm_tree(self.params)
            rm_tree(self.output)

            self.output_dir.mkdir(parents=True, exist_ok=True)
            self.plots_dir.mkdir(parents=True, exist_ok=True)

        # check if the output directory is empty
        if len(list(self.output_dir.iterdir())) > 0:
            raise SimulationError("Output directory is not empty. Use force=True to overwrite.")

        # check if the plots directory is empty
        if len(list(self.plots_dir.iterdir())) > 0:
            raise SimulationError("Plots directory is not empty. Use force=True to overwrite.")

        # check if the params file exists
        if self.params.exists():
            raise SimulationError("Params file already exists. Use force=True to overwrite.")

        # check if the output file exists
        if self.output.exists():
            raise SimulationError("Output file already exists. Use force=True to overwrite.")

        # dump the configuration to a JSON file to be read by the executable
        self.dump_config()

        # run the simulation
        try:
            with open(self.output, "w") as fp:
                subprocess.run([EXE], stderr=fp, stdout=fp, cwd=self.base_dir, timeout=timeout)
        except Exception as e:
            with open(self.output, "w") as fp:
                fp.write(str(e))

    def plot_0d(self):
        """
        Plot the results of the simulation in 0D (ie convergence plots etc).
        """
        n = self.config.n
        x = self.config.x

        title = f"Re = {self.config.re}, Ca = {self.config.ca}, m = {self.config.m}"
        if self.config.uses_estimator:
            title += f", p = {self.config.p}"

        # Load the data
        data = np.loadtxt(self.output_dir / 'ns-0.dat')
        t = data[:, 0]  # time
        dh = data[:, 1]  # interfacial deviation, ||h - 1||2
        de = data[:, 2]  # estimator error, ||h - z||2

        # Plot the data
        fig, ax = plt.subplots()
        ax.semilogy(t, dh, label="dh")

        if self.config.uses_estimator:
            ax.semilogy(t, de, label="de")

        ax.set_xlabel("t")
        ax.legend()
        ax.set_title(title)
        fig.savefig(self.plots_dir / "lines.png")

        if (self.config.strategy == "estimator") and (self.config.p > 0):
            # Load the L matrix
            L = np.atleast_2d(np.loadtxt(self.output_dir / 'L.dat'))
            lh = L[:n, 0]
            lq = L[n:, 0]

            # Plot the first column of L
            fig, ax = plt.subplots(2, 1)
            ax[0].plot(x, lh, label="L0 (interface)")
            ax[0].set_xlabel("x")
            ax[0].set_ylabel("L0")
            ax[0].set_title(title)
            ax[0].legend()

            ax[1].plot(x, lq, label="L0 (flux)")
            ax[1].set_xlabel("x")
            ax[1].set_ylabel("L0")
            ax[1].legend()

            fig.savefig(self.plots_dir / "L.png")
            plt.close(fig)

        if (self.config.strategy in ["lqr", "estimator"]) and (self.config.m > 0):
            # Load the K matrix
            K = np.atleast_2d(np.loadtxt(self.output_dir / 'K.dat'))
            kh = K[0, :n]
            kq = K[0, n:]

            # Plot the first column of L
            fig, ax = plt.subplots(2, 1)
            ax[0].plot(x, kh, label="K0 (interface)")
            ax[0].set_xlabel("x")
            ax[0].set_ylabel("K0")
            ax[0].set_title(title)
            ax[0].legend()

            ax[1].plot(x, kq, label="K0 (flux)")
            ax[1].set_xlabel("x")
            ax[1].set_ylabel("K0")
            ax[1].legend()

            fig.savefig(self.plots_dir / "K.png")
            plt.close(fig)

    def plot_1d(self):
        """
        Plot the results of the simulation in 1D (ie interface plots).
        """
        x = self.config.x

        # Save the interface plots to a subdirectory
        interface_dir = self.plots_dir / "interfaces"
        interface_dir.mkdir(parents=True, exist_ok=True)

        # Load 1D data to get the time
        data = np.loadtxt(self.output_dir / 'ns-0.dat')
        t = data[:, 0]
        nsteps = len(t)

        # Go through the output files and find the max/min interfacial (and estimator) heights
        hmax = 0
        hmin = 0
        for i in range(nsteps):
            # Extract the data
            filename = self.output_dir / f'ns-1-{i:010d}.dat'
            data = np.loadtxt(filename)
            h = data[:, 1] - 1.0
            z = data[:, 3] - 1.0

            hmax = max(hmax, np.max(h))
            hmin = min(hmin, np.min(h))
            hmax = max(hmax, np.max(z))
            hmin = min(hmin, np.min(z))
        hmax = np.ceil(hmax / 0.5) * 0.5
        hmin = np.floor(hmin / 0.5) * 0.5

        # plot some dummy data
        fig, ax = plt.subplots()

        ax.set_ylim(hmin, hmax)
        ax.set_xlim(0, self.config.lx)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title("TITLE")

        hplot, = ax.plot(x, x, label="h - 1")
        fplot, = ax.plot(x, x, label="f")
        zplot = None
        if self.config.uses_estimator:
            zplot, = ax.plot(x, x, label="z - 1")
        ax.legend(loc="lower right")

        # Plot the real data
        for i in range(nsteps):
            # Extract the data
            filename = self.output_dir / f'ns-1-{i:010d}.dat'
            data = np.loadtxt(filename)
            h = data[:, 1] - 1.0
            f = data[:, 2]
            z = data[:, 3] - 1.0

            # Update the plots
            ax.set_title(f'time {t[i]} [step {i}]')
            hplot.set_ydata(h)
            fplot.set_ydata(f)
            if zplot is not None:
                zplot.set_ydata(z)

            fig.savefig(interface_dir / f'interface-{i:010d}.png')

        # Create a GIF from the interface plots
        fps = 20
        gif_path = self.plots_dir / "interface.gif"
        with ExitStack() as stack:
            # lazily load images
            images = (
                stack.enter_context(Image.open(f))
                for f in sorted(interface_dir.glob('*.png'))
            )

            # extract first image from iterator
            image = next(images)

            image.save(
                fp=gif_path,
                format='GIF',
                append_images=images,
                save_all=True,
                duration=1000 / fps,
                loop=0,
                quality=99
            )

    def plot(self):
        """
        Plot the results of the simulation.
        """
        # check if the output directory is empty
        if len(list(self.output_dir.iterdir())) == 0:
            print("Output directory is empty. Run the simulation first.")
            return

        if self.config.output_dim >= 0:
            self.plot_0d()
        if self.config.output_dim >= 1:
            self.plot_1d()
