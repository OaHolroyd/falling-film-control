import json
from pathlib import Path

import numpy as np


class Config:
    def __init__(
            self,
            re: float = 15.0,
            ca: float = 0.05,
            lx: float = 30.0,
            level: int = 8,
            tmax: float = 300.0,
            cstart: float = 200.0,
            m: int = 5,
            p: int = 10,
            dtout: float = 0.5,
            shift: float = 0.0,
            strategy: str = "lqr",
            model: str = "wr",
            exact_flux: bool = True,
            cost_weight: float = 1.0,
            t0: float = 0.0,
            theta: float = 1.047197551,
            rho_ratio: float = 1000.0,
            mu_ratio: float = 1.0,
            output_dim: int = 1,
            actuator_width: float = 0.1,
            strength: float = 1.0,
    ):
        """
        Create a controlled NS simulation configuration.

        Args:
            re: Reynolds number.
            ca: Capillary number.
            lx: Domain length.
            level: Grid level (so there are 2^level grid points).
            tmax: Maximum time.
            cstart: Time at which the control starts.
            m: Number of actuators.
            p: Number of observers.
            dtout: Output time step.
            shift: How much to shift the observers (normalised so that a shift of 1 shifts all observers upstream by one spacing).
            strategy: Control strategy to use. Options are "lqr", "static", "dynamic", "pair", or "estimator".
            model: Which thin film model to use when deriving the control. Options are "wr" or "benney".
            exact_flux: Whether to use the exact flux or an approximation.
            cost_weight: Weighting of the interfacial deviation (vs the actuation) in the cost function.
            t0: If this is > 0, run a (fast) WR simulation to t=t0 and use this as an initial condition for the full simulation.
            theta: Angle of the plate.
            rho_ratio: Ratio of the density of the lower fluid to the upper fluid.
            mu_ratio: Ratio of the dynamic viscosity of the lower fluid to the upper fluid.
            output_dim: Maximum dimensionality of the output data. (0 for integral statistics, 1 for interfaces, and 2 for full fields).
            actuator_width: Width parameter for the actuators (smaller means a narrower actuator).
            strength: Scaling to apply to the actuators.
        """
        self.re = re
        self.ca = ca
        self.lx = lx
        self.level = level
        self.tmax = tmax
        self.cstart = cstart
        self.m = m
        self.p = p
        self.dtout = dtout
        self.shift = shift
        self.strategy = strategy
        self.model = model
        self.exact_flux = exact_flux
        self.cost_weight = cost_weight
        self.t0 = t0
        self.theta = theta
        self.rho_ratio = rho_ratio
        self.mu_ratio = mu_ratio
        self.output_dim = output_dim
        self.actuator_width = actuator_width
        self.strength = strength

    @property
    def n(self):
        """
        Number of grid points in the domain.
        """
        return 2 ** self.level

    @property
    def dx(self):
        """
        Grid spacing in the domain.
        """
        return self.lx / self.n

    @property
    def x(self):
        """
        Grid points in the domain.
        """
        return self.dx * (0.5 + np.arange(self.n))

    @property
    def uses_estimator(self):
        """
        Whether the control strategy uses an estimator.
        """
        return self.strategy in ["estimator", "dynamic"]

    def update(self, **kwargs):
        """
        Update the configuration with new parameters.

        Args:
            kwargs: Keyword arguments for the parameters to update.
        """
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise ValueError(f"Invalid parameter '{key}'")

    def to_dict(self, flat: bool = False) -> dict:
        """
        Return the configuration as a dictionary.

        Args:
            flat: If True, return a flat dictionary. Otherwise, return a nested dictionary (used for the simulations).
        """
        if flat:
            return {
                "re": self.re,
                "ca": self.ca,
                "lx": self.lx,
                "level": self.level,
                "tmax": self.tmax,
                "cstart": self.cstart,
                "m": self.m,
                "p": self.p,
                "dtout": self.dtout,
                "shift": self.shift,
                "strategy": self.strategy,
                "model": self.model,
                "exact_flux": self.exact_flux,
                "cost_weight": self.cost_weight,
                "t0": self.t0,
                "theta": self.theta,
                "rho_ratio": self.rho_ratio,
                "mu_ratio": self.mu_ratio,
                "output_dim": self.output_dim,
                "actuator_width": self.actuator_width,
                "strength": self.strength
            }

        else:
            return {
                "DOMAIN": {
                    "Lx": self.lx,
                    "Ly": self.lx * 0.25,
                    "tmax": self.tmax,
                    "t0": self.t0
                },
                "PHYSICAL": {
                    "Re": self.re,
                    "Ca": self.ca,
                    "theta": self.theta,
                    "rho_ratio": self.rho_ratio,
                    "mu_ratio": self.mu_ratio,
                },
                "SOLVER": {
                    "level": self.level,
                    "dtout": self.dtout,
                    "output": self.output_dim,
                },
                "CONTROL": {
                    "M": self.m,
                    "P": self.p,
                    "start": self.cstart,
                    "width": self.actuator_width,
                    "alpha": self.strength,
                    "del": self.shift,
                    "mu": self.cost_weight,
                    "rom": self.model,
                    "strategy": self.strategy,
                    "exact_flux": self.exact_flux,
                }
            }

    @classmethod
    def from_dict(cls, config: dict):
        """
        Load the configuration from a (flat or nested) dictionary.

        Args:
            config: Dictionary to load from.
        """
        # check if the dict is flat or not
        if "DOMAIN" in config:
            return cls(
                re=config["PHYSICAL"]["Re"],
                ca=config["PHYSICAL"]["Ca"],
                lx=config["DOMAIN"]["Lx"],
                level=config["SOLVER"]["level"],
                tmax=config["DOMAIN"]["tmax"],
                cstart=config["CONTROL"]["start"],
                m=config["CONTROL"]["M"],
                p=config["CONTROL"]["P"],
                dtout=config["SOLVER"]["dtout"],
                shift=config["CONTROL"]["del"],
                strategy=config["CONTROL"]["strategy"],
                model=config["CONTROL"]["rom"],
                exact_flux=config["CONTROL"]["exact_flux"],
                cost_weight=config["CONTROL"]["mu"],
                t0=config["DOMAIN"]["t0"],
                theta=config["PHYSICAL"]["theta"],
                rho_ratio=config["PHYSICAL"]["rho_ratio"],
                mu_ratio=config["PHYSICAL"]["mu_ratio"],
                output_dim=config["SOLVER"]["output"],
                actuator_width=config["CONTROL"]["width"],
                strength=config["CONTROL"]["alpha"],
            )
        else:
            return cls(**config)

    def to_json(self, filename: Path | str, flat: bool = False):
        """
        Write the configuration to a JSON file.

        Args:
            filename: Filename to write to.
            flat: If True, write a flat JSON file. Otherwise, write a nested JSON file (used for the simulations).
        """
        with open(filename, "w") as fp:
            json.dump(self.to_dict(flat=flat), fp, indent=4)

    @classmethod
    def from_json(cls, filename: Path | str):
        """
        Load the configuration from a JSON file.

        Args:
            filename: Filename to load from.
        """
        with open(filename, "r") as fp:
            config = json.load(fp)
        return cls.from_dict(config)
