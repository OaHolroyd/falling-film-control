import itertools
from datetime import datetime
from pathlib import Path

import numpy as np
from mpi4py import MPI

from scripts.config import Config
from scripts.simulation import Simulation, SimulationError


class MultiSim:
    def __init__(self, base_dir: Path | str, config: Config, exe: Path | str, shift: int = 0):
        """
        Set up a parallel batch of simulation runs

        Args:
            base_dir: base directory for individual simulation subdirectories
            config: base configuration object
            exe: path to the executable
            shift: shift to apply to the folder naming (to avoid collisions if running more than one MultiSim instance
                with the same base_dir)
        """
        self.base_dir = Path(base_dir)
        self.base_config = config
        self.exe = Path(exe)
        self.shift = shift

        self.count = 0
        self.configs: list[(int, Config)] = []

        comm = MPI.COMM_WORLD
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()

        if self.rank == 0:
            self.base_dir.mkdir(parents=True, exist_ok=True)

    def reset(self, config: Config | None = None):
        """
        Reset the simulation configurations to a new base configuration.

        Args:
            config: New base configuration object.
        """
        if config is None:
            config = self.base_config
        self.base_config = config

        self.configs = []

    def add_config(self, config: Config):
        """
        Add a new configuration to the list of configurations.

        Args:
            config: Configuration object to add.
        """
        self.configs.append((self.count + self.shift, config))
        self.count += 1

    def add_variant(self, **kwargs):
        """
        Add a new configuration variant based on the base configuration.

        Args:
            **kwargs: Configuration parameters to modify.
        """
        config = self.base_config.copy()
        config.update(**kwargs)
        self.add_config(config)

    def add_variants(self, **kwargs):
        """
        Add multiple configuration variants based on the base configuration.

        Args:
            **kwargs: Configuration parameters to modify. Each entry should be a list of the same length.
        """
        # Get an arbitrary value from kwargs
        val = next(iter(kwargs.values()))

        # Check if the value is a list
        if not isinstance(val, list):
            raise ValueError("Expected list for variant parameters")
        n = len(val)

        # the values in kwargs should be lists with the same length
        for k, v in kwargs.items():
            if not isinstance(v, list):
                raise ValueError(f"Expected list for {k}, got {type(v)}")
            if len(v) != n:
                raise ValueError("All lists must have the same length")

        # for each value in the lists, create a new configuration
        for i in range(n):
            new_kwargs = {}
            for k, v in kwargs.items():
                new_kwargs[k] = v[i]
            self.add_variant(**new_kwargs)

    def add_variants_product(self, **kwargs):
        """
        Add multiple configuration variants based on the base configuration using a Cartesian product.

        Args:
            **kwargs: Configuration parameters to modify. Each entry should be a list.
        """
        # the values in kwargs should be lists
        for k, v in kwargs.items():
            if not isinstance(v, (list, tuple, np.ndarray)):
                raise ValueError(f"Expected list for {k}, got {type(v)}")

        keys = kwargs.keys()
        for instance in itertools.product(*kwargs.values()):
            new_kwargs = dict(zip(keys, instance))
            self.add_variant(**new_kwargs)

    def run(self, index: int, timeout: int = 21600, plot: bool = True):
        """
        Run a simulation with the given index.

        Args:
            index: Index of the simulation to run.
            timeout: Timeout for the simulation in seconds (default is 6 hours).
            plot: If True, generate plots after the simulation.
        """
        if index < 0 or index >= len(self.configs):
            raise ValueError("Index out of range")

        i, config = self.configs[index]
        sim = Simulation(config=config, base_dir=self.base_dir / f"run-{i}", exe=self.exe)

        try:
            sim.run(force=False, timeout=timeout)
        except SimulationError:
            # if the simulation has already been run, check if it ran to completion
            if not sim.has_completed:
                print(f"  rerunning {index}")
                sim.run(force=True, timeout=timeout)
            else:
                print(f"  simulation {index} already completed, skipping")
        if plot:
            sim.plot()

    def run_all(self, timeout: int = 21600, plot: bool = True):
        # decide which configurations to run
        n_configs = len(self.configs)
        indices = []
        for i in range(self.rank, n_configs, self.size):
            indices.append(i)

        for i, ind in enumerate(indices):
            tstart = datetime.now()
            print(f"[{self.rank:2d}] Starting simulation {ind} ({i + 1}/{len(indices)}) at {tstart}")
            self.run(ind, timeout=timeout, plot=plot)
            tend = datetime.now()
            print(f"[{self.rank:2d}] Finished simulation {ind} ({i + 1}/{len(indices)}) after {tend - tstart}")
