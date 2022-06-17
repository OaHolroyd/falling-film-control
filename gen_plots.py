#!/usr/bin/env python3
import argparse
import src.plot as plt


def main(n, track, clip, rate, start):
    """
    Main function for CLI
    """
    # static plots
    plt.plot_hstats(n)
    print("hstats done")

    plt.plot_umax(n)
    print("umax done")

    # animations
    plt.plot_series(n=n, i0=start, key="interface", Ly=2.0,
                    track=track, clip=clip, rate=rate)
    print("interface done")
    plt.plot_series(n=n, i0=start, key="control", Ly=2.0,
                    track=track, clip=clip, rate=rate)
    print("control done")
    plt.plot_series(n=n, i0=start, key="fluid", Ly=2.0,
                    track=track, clip=clip, rate=rate)
    print("fluid done")
    plt.plot_series(n=n, i0=start, key="level", Ly=8.0,
                    track=track, clip=clip, rate=rate)
    print("level done")
    plt.plot_series(n=n, i0=start, key="vorticity", Ly=4.0,
                    track=track, clip=clip, rate=rate)
    print("vorticity done")
    plt.plot_series(n=n, i0=start, key="speed", Ly=4.0,
                    track=track, clip=clip, rate=rate)
    print("speed done")
    plt.plot_series(n=n, i0=start, key="u", Ly=4.0,
                    track=track, clip=clip, rate=rate)
    print("u done")
    plt.plot_series(n=n, i0=start, key="v", Ly=4.0,
                    track=track, clip=clip, rate=rate)
    print("v done")
    plt.plot_series(n=n, i0=start, key="pressure", Ly=4.0,
                    track=track, clip=clip, rate=rate)
    print("pressure done")


if __name__ == "__main__":
    # process args
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', action='store', type=int,
                        help='number of frames')
    parser.add_argument('-s', '--start', action='store', type=int, default=0,
                        help='start frame')
    parser.add_argument('-t', '--track', action='store_true',
                        help='track peak')
    parser.add_argument('-c', '--clip', action='store_true',
                        help='clip to fluid')
    parser.add_argument('-r', '--rate', action='store', type=int, default=5,
                        help='frame skips')
    args = parser.parse_args()

    main(args.n, args.track, args.clip, args.rate, args.start)
