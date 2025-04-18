import numpy as np

from scripts.config import Config
from scripts.multisim import MultiSim


def main():
    config = Config.from_json("params.json")

    re = np.logspace(0, 2, 20)
    ca = np.logspace(-3, -1, 20)

    re = [5.0, 15.0]
    ca = [0.05]

    multisim = MultiSim(base_dir="multirun", config=config, exe="film-ns")
    multisim.add_variants_product(re=re, ca=ca)
    # multisim.run_all()

    multisim.run_all_binary_search(11)


if __name__ == "__main__":
    main()
