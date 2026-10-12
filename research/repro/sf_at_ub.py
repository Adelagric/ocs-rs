"""Dense sexed prototype at an explicit cap. Prints only verified optima.

Usage: python3 research/repro/sf_at_ub.py wheat <ub> [<ub> ...]
Historical n-keyed /tmp exports require an explicit --legacy flag.
"""
import argparse
import time
from dense_support import load_panel, support_first


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("panel")
    parser.add_argument("ub", type=float, nargs="+")
    parser.add_argument("--legacy", action="store_true")
    args = parser.parse_args()
    K, bv, male, _ = load_panel(args.panel, args.legacy)
    for ub in args.ub:
        t0 = time.perf_counter()
        c, S, iters = support_first(K, bv, male, ub)
        dt = time.perf_counter() - t0
        print(f"support-first @ ub={ub:.6f}: gain={bv @ c:+.5f} c'Kc={c @ K @ c:.6f} "
              f"support={len(S)} iters={iters} verified=True temps={dt:.3f}s")


if __name__ == "__main__":
    main()
