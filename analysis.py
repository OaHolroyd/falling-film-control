from argparse import ArgumentParser

import numpy as np

from scripts.config import Config
from scripts.multisim import MultiSim


def main():
    # # set up the argument parser
    # parser = ArgumentParser()
    # parser.add_argument("--re", nargs='+', type=float, default=None)
    # parser.add_argument("--ca", nargs='+', type=float, default=None)
    # parser.add_argument("--m", type=int, default=None)
    # parser.add_argument("--p", type=int, default=None)
    # parser.add_argument("--params", type=str, default="params.json")
    #
    # # unpack the arguments
    # args = parser.parse_args()._get_kwargs()
    # args = {k: v for k, v in args if v is not None}
    #
    # # remove non-parameter arguments
    # params_file = args.pop("params")
    #
    # # set up and run the sims
    # config = Config.from_json(params_file)
    # multisim = MultiSim(base_dir="multirun", config=config, exe="film-ns")
    # multisim.add_variants_product(**args)
    # multisim.run_all()

    config = Config.from_json("params.json")

    re = np.logspace(0, 2, 20)
    ca = np.logspace(-3, -1, 20)

    multisim = MultiSim(base_dir="multirun", config=config, exe="film-ns")
    multisim.add_variants_product(re=re, ca=ca)
    multisim.run_all()
    # multisim.run_all_binary_search(19)


if __name__ == "__main__":
    main()
