import numpy as np

from scripts.config import Config
from scripts.multisim import MultiSim


def main():
    config = Config.from_json("params.json")

    re = np.logspace(0, 2, 20)
    ca = np.logspace(-3, -1, 20)
    p = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]

    multisim = MultiSim(base_dir="multirun", config=config, exe="film-ns")
    multisim.add_variants_product(re=re, ca=ca, p=p)
    multisim.run_all()


if __name__ == "__main__":
    main()
