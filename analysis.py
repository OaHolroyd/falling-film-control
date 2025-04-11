from pathlib import Path

from scripts.config import Config
from scripts.simulation import Simulation


def main():
    config = Config.from_json("params.json")
    config.update(re=8.0)
    config.update(tmax=150.0)
    config.update(cstart=100.0)

    test_dir = Path("test")
    test_dir.mkdir(parents=True, exist_ok=True)

    sim = Simulation(config=config, base_dir=test_dir, exe="film-ns")
    # sim.run()
    sim.plot()


if __name__ == "__main__":
    main()
