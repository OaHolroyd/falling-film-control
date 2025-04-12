from pathlib import Path

from scripts.config import Config
from scripts.multisim import MultiSim
from scripts.simulation import Simulation


def main():
    config = Config.from_json("params.json")

    multisim = MultiSim(base_dir="multirun", config=config, exe="film-ns")
    multisim.add_variants_product(re=[1.0, 2.0], ca=[0.05, 0.01])
    multisim.run_all()


if __name__ == "__main__":
    main()
