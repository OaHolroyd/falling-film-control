import itertools
import shutil
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

    def run(self, index: int, run_dir: Path | str, timeout: int = 21600, plot: bool = True):
        """
        Run a simulation with the given index.

        Args:
            index: Index of the simulation to run.
            run_dir: Base directory for the simulation.
            timeout: Timeout for the simulation in seconds (default is 6 hours).
            plot: If True, generate plots after the simulation.
        """
        if index < 0 or index >= len(self.configs):
            raise ValueError("Index out of range")

        _, config = self.configs[index]
        sim = Simulation(config=config, base_dir=run_dir, exe=self.exe)

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

    def divide_problems(self):
        """
        Divide the problems among the ranks.

        Returns:
            my_indices: Indices of the configs assigned to this rank.
            my_total: Total expected runtime for this rank.
            my_prop: Proportion of the expected runtime that this rank will do.
        """
        # try and balance each proc's expected total runtime
        indices = list(range(len(self.configs)))
        indices.sort(key=lambda i: self.configs[i][1].expected_runtime)
        totals = [0 for _ in range(self.size)]
        split_indices = [[] for _ in range(self.size)]
        while len(indices) > 0:
            index = indices.pop()

            # add the config to the list with the smallest total
            i = 0
            total = totals[0]
            for j in range(1, self.size):
                if totals[j] < total:
                    total = totals[j]
                    i = j
            split_indices[i].append(index)
            totals[i] += self.configs[index][1].expected_runtime

        # Compute the proportion of the expected runtime that this rank will do
        my_total = totals[self.rank]
        my_indices = split_indices[self.rank]
        total = 0
        for t in totals:
            total += t
        my_prop = my_total / total

        return my_indices, my_total, my_prop

    def run_all(self, timeout: int = 21600, plot: bool = True):
        """
        Run all simulations for this rank.

        Args:
            timeout: Timeout for the simulations in seconds (default is 6 hours).
            plot: If True, generate plots after each simulation.
        """
        # try and balance each proc's expected total runtime
        my_indices, my_total, my_prop = self.divide_problems()
        print(f"[{self.rank:2d}] Performing {100 * my_prop:.2f}% of the total")

        # Run the simulations for this rank
        tbegin = datetime.now()
        completed = 0
        for i, ind in enumerate(my_indices):
            tstart = datetime.now()
            print(f"[{self.rank:2d}] Starting simulation {ind} ({i + 1}/{len(my_indices)}) at {tstart}")
            self.run(ind, run_dir=self.base_dir / f"run-{i}", timeout=timeout, plot=plot)
            tend = datetime.now()

            completed += self.configs[ind][1].expected_runtime / my_total
            rem_time = (tend - tbegin) * (1.0 / completed - 1.0)

            print(
                f"[{self.rank:2d}] Finished simulation {ind} ({i + 1}/{len(my_indices)}) after {tend - tstart}\n     estimated time remaining {rem_time}"
            )

    def run_binary_search(self, mmax, index: int, base_dir: Path, timeout: int = 21600, plot: bool = True):
        """
        Run a binary search to find the minimum value of m that works.

        Args:
            mmax: Maximum value of m to search.
            index: Index of the simulation to run.
            base_dir: Base directory for the simulation.
            timeout: Timeout for the simulation in seconds (default is 6 hours).
            plot: If True, generate plots after the simulation.
        """
        base_config = self.configs[index][1]

        # define endpoints for binary search
        m0 = 0
        m1 = mmax

        # define function to check if the simulation works
        def works(config: Config, run_dir: Path | str):
            # Run the simulation
            sim = Simulation(config=config, base_dir=run_dir, exe=self.exe)
            try:
                sim.run(force=False, timeout=timeout)
            except SimulationError:
                # if the simulation has already been run, check if it ran to completion
                if not sim.has_completed:
                    sim.run(force=True, timeout=timeout)

            if plot:
                sim.plot()

            return sim.has_converged

        # check lower endpoint (ie no control)
        config = base_config.copy()
        config.m = m0
        if works(config, run_dir=base_dir / "run-0"):
            # don't need to do anything because m=0 works
            return

        # move dump files to the base directory so they can be reused
        if not config.uses_estimator:
            if not (base_dir / "dump").exists():
                shutil.move(base_dir / "run-0" / "dump", base_dir / "dump")

        # check upper endpoint (max controls)
        config = base_config.copy()
        config.m = m1
        if not works(config, run_dir=base_dir / "run-1"):
            # mmax doesn't work, so no need to search
            return

        # binary search for the minimum m that works
        def binary_search(a0, a1, iter):
            # we are done if a0 and a1 are adjacent
            if a0 == a1 - 1:
                return

            # check the midpoint
            a = round((a0 + a1) / 2)
            config = base_config.copy()
            config.m = a

            # check if the simulation works
            if works(config, run_dir=base_dir / f"run-{iter}"):
                # if it works, search the lower half
                return binary_search(a0, a, iter + 1)
            else:
                # if it doesn't work, search the upper half
                return binary_search(a, a1, iter + 1)

        binary_search(m0, m1, 2)

    def run_all_binary_search(self, mmax, timeout: int = 21600, plot: bool = True):
        """
        For each config assigned to this rank, run a binary search to find the minimum value of m that works.

        Should not be run across configs that differ only in m, as this will duplicate work.

        Args:
            mmax: Maximum value of m to search.
            timeout: Timeout for the simulations in seconds (default is 6 hours).
            plot: If True, generate plots after each simulation.
        """
        # try and balance each proc's expected total runtime
        my_indices, my_total, my_prop = self.divide_problems()
        print(f"[{self.rank:2d}] Performing {100 * my_prop:.2f}% of the total")

        # Run the simulations for this rank
        for i, ind in enumerate(my_indices):
            tstart = datetime.now()
            print(f"[{self.rank:2d}] Starting search {ind} ({i + 1}/{len(my_indices)}) at {tstart}")
            self.run_binary_search(mmax=mmax, index=ind, base_dir=self.base_dir / f"search-{self.rank}-{i}",
                                   timeout=timeout, plot=plot)
            tend = datetime.now()

            print(f"[{self.rank:2d}] Finished search {ind} ({i + 1}/{len(my_indices)}) after {tend - tstart}")
