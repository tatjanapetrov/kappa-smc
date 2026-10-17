#!/usr/bin/env python3
"""Reproduce Figure 4a using models/polimerisation.ka and src.KappaSMC.

Run from the repository root:
    python reproduce_polymerization.py --kappa-bin /path/to/kappa/bin

Requires numpy, matplotlib, kappy and the KaSim executables.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib.pyplot as plt
import numpy as np

from kappa_smc.bltl import AP, BLTLMonitor, F_leq
from kappa_smc.estimation import KappaSMC, hoeffding_sample_size
from kappa_smc.simulation import KappaSimulator, SimulationConfig


ROOT = Path(__file__).resolve().parent.parent


def chain_pattern(i: int) -> str:
    """Pattern |⊢ A_i −|: >=i+1 consecutive agents, free x endpoint.

    We intentionally leave the other site of the initial agent unspecified:
    the matching sequence may be part of a longer polymer.
    """
    if not 1 <= i <= 5:
        raise ValueError("i must be between 1 and 5")
    if i == 1:
        return "A(y[_], x[.])"  
    # agents = ["A(x[1])"]
    agents = ["A(y[_], x[1])"]
    for bond in range(1, i - 1):
        agents.append(f"A(y[{bond}], x[{bond + 1}])")
    agents.append(f"A(y[{i - 1}], x[.])")
    return ", ".join(agents)


def expanded_model(source: Path) -> str:
    """Keep the official model unchanged, except for five new observables."""
    model = source.read_text(encoding="utf-8")
    extra = "\n".join(f"%obs: 'C_{i}' |{chain_pattern(i)}|" for i in range(1, 6))
    return model.rstrip() + "\n\n// Figure 4a additional pattern counts\n" + extra + "\n"


def reproduce(model_path: Path, kappa_bin: str, runs: int, plot_period: float, out: Path) -> None:
    if not model_path.is_file():
        raise FileNotFoundError(f"Missing official model: {model_path}")
    out.mkdir(parents=True, exist_ok=True)

    # Use one combined formula just to make KappaSMC instrument ALL five APs
    # simultaneously; then check each F formula separately on the same trace.
    formula = 'F_leq(10, AP("C_1>4"))'
    for i in range(2, 6):
        formula = f'Or({formula}, F_leq(10, AP("C_{i}>4")))'

    simulator = KappaSimulator(
        kappa_bin_path=kappa_bin,
        config=SimulationConfig(max_time=10.0, plot_period=plot_period),
    )
    smc = KappaSMC(simulator)

    with TemporaryDirectory(prefix="kappa_fig4a_") as tmp:
        augmented = Path(tmp) / "polimerisation_fig4a.ka"
        augmented.write_text(expanded_model(model_path), encoding="utf-8")
        prepared = smc.prepare(model_path=str(augmented), formula_source=formula)

        # The fully instrumented model is saved so the input is reproducible.
        (out / "polimerisation_fig4a_instrumented.ka").write_text(
            prepared.instrumented_model, encoding="utf-8"
        )

        properties = [F_leq(10, AP(f"C_{i}>4", idx=i - 1)) for i in range(1, 6)]
        successes = np.zeros(5, dtype=int)

        for run in range(1, runs + 1):
            trace = smc.run_trace(prepared)
            monitor = BLTLMonitor(trace.times, trace.ap_values)
            for i, property_ in enumerate(properties):
                successes[i] += int(monitor.holds(property_))
            if run == 1 or run % max(1, runs // 10) == 0 or run == runs:
                print(f"{run}/{runs}  P(F): {np.round(successes / run, 3)}", flush=True)

    p_f = successes / runs
    p_g = 1 - p_f  # G_[0,10](C_i<5) = not F_[0,10](C_i>4)
    x = np.arange(1, 6)
    epsilon = 0.05  # the paper's displayed error tolerance, not a Wilson CI

    fig, ax = plt.subplots(figsize=(7, 4.8))
    for values, color, legend in (
        (p_f, "forestgreen", r"$F_{[0,10]}(C_i>4)$"),
        (p_g, "black", r"$G_{[0,10]}(C_i<5)$"),
    ):
        lower = np.minimum(epsilon, values)
        upper = np.minimum(epsilon, 1 - values)
        ax.errorbar(x, values, yerr=np.vstack([lower, upper]),
                    color=color, ecolor="indianred", lw=1.3,
                    marker="+", capsize=2, label=legend)

    ax.set_xticks(x, [rf"$|\vdash A_{i}-|$" for i in x], rotation=25)
    ax.set_xlim(0.8, 5.2)
    ax.set_ylim(-0.025, 1.025)
    ax.set_xlabel("Polymer chain length")
    ax.set_ylabel("Satisfaction probability")
    ax.legend()
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out / f"figure4a_polymerisation.{ext}", dpi=220)
    plt.close(fig)

    with (out / "figure4a_polymerisation.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["i", "pattern", "successes", "runs", "P_F", "P_G"])
        for j in range(5):
            writer.writerow([j + 1, chain_pattern(j + 1), successes[j], runs, p_f[j], p_g[j]])
    print(f"Results saved to {out.resolve()}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=ROOT / "models" / "polimerisation.ka")
    parser.add_argument("--kappa-bin", required=True, help="Directory containing KaSim binaries")
    parser.add_argument("--runs", type=int, default=1000)
    parser.add_argument("--plot-period", type=float, default=0.01)
    parser.add_argument("--out", type=Path, default=ROOT / "results" / "polymerisation")
    parser.add_argument("--epsilon", type=float, default=None)
    parser.add_argument("--confidence-level", type=float, default=0.99)
    args = parser.parse_args()
    if args.epsilon is not None:
        args.runs = hoeffding_sample_size(
            epsilon=args.epsilon,
            confidence_level=args.confidence_level,
        )
    if args.runs < 1 or args.plot_period <= 0:
        parser.error("--runs and --plot-period must be positive")
    reproduce(args.model, args.kappa_bin, args.runs, args.plot_period, args.out)


if __name__ == "__main__":
    main()