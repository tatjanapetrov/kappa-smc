import time
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from statistics import NormalDist

from .bltl import (BLTLMonitor, Formula, bind_ap_indices,
                   extract_atomic_propositions, formula_horizon,
                   parse_ap_string, parse_formula)
from .instrumentation import instrument_kappa_model
from .simulation import KappaSimulator, process_simulation_results

import math

def hoeffding_sample_size(
    epsilon: float,
    confidence_level: float,
) -> int:
    """Minimum sample size satisfying Hoeffding's two-sided bound."""
    if not 0 < epsilon < 1:
        raise ValueError("epsilon must be between 0 and 1.")
    if not 0 < confidence_level < 1:
        raise ValueError("confidence_level must be between 0 and 1.")

    alpha = 1 - confidence_level

    return math.ceil(
        math.log(2 / alpha) / (2 * epsilon**2)
    )

@dataclass(frozen=True)
class SMCConfig:
    runs: int = 100
    progress_every: int = 10
    confidence_level: float = 0.95

    @classmethod
    def from_hoeffding(
        cls,
        epsilon: float,
        confidence_level: float = 0.95,
        progress_every: int = 10,
    ) -> "SMCConfig":
        return cls(
            runs=hoeffding_sample_size(epsilon, confidence_level),
            progress_every=progress_every,
            confidence_level=confidence_level,
        )
    
    def __post_init__(self) -> None:
        if self.runs <= 0:
            raise ValueError("runs must be positive.")

        if self.progress_every <= 0:
            raise ValueError("progress_every must be positive.")

        if not 0.0 < self.confidence_level < 1.0:
            raise ValueError("confidence_level must be between 0 and 1.")

@dataclass(frozen=True)
class PreparedProperty:
    """
    Formula and model after all preprocessing has been completed.
    """
    formula: Formula
    atomic_propositions: tuple[str, ...]
    instrumented_model: str


@dataclass
class TraceResult:
    """
    Result of one stochastic simulation and one BLTL evaluation.
    """
    holds: bool
    times: np.ndarray
    ap_values: np.ndarray
    atomic_propositions: tuple[str, ...]
    legend: tuple[str, ...]
    series: np.ndarray

    def print_results_table(self, max_rows: int | None = None) -> None:
        """Print a simulation table in chronological order."""
        header = " | ".join(f"{column:>12}" for column in self.legend)

        print(header)
        print("-" * len(header))

        rows = self.series if max_rows is None else self.series[:max_rows]
        
        for row in rows:
            print(" | ".join(f"{value:12.3f}" for value in row))

        if (max_rows is not None and len(self.series) > max_rows):
            print(f"... {len(self.series) - max_rows} additional rows")
        
    def plot(
        self,
        variables: list[str] | None = None,
        *,
        title: str = "Kappa simulation trace",
        ylabel: str = "Value",
        step: bool = False,
    ):
        """
        Plot selected variables or observables over simulation time.

        Parameters
        ----------
        variables: names of columns in `legend`. If omitted, all columns except `[T]` are plotted.

        title: plot title.

        ylabel: label of the vertical axis.

        step: use a step plot instead of a line plot, useful for Boolean observables.

        Returns
        -------
        fig, ax: Matplotlib figure and axes.
        """
        legend = list(self.legend)

        if "[T]" not in legend:
            raise ValueError("The trajectory does not contain the '[T]' column.")

        if variables is None:
            variables = [name for name in legend if name != "[T]"]

        missing = [variable for variable in variables if variable not in legend]

        if missing:
            raise ValueError(f"Variables not found: {missing}\n Available columns: {legend}")

        time_index = legend.index("[T]")
        times = self.series[:, time_index]

        fig, ax = plt.subplots()

        for variable in variables:
            variable_index = legend.index(variable)
            values = self.series[:, variable_index]

            if step:
                ax.step(times,values,where="post",label=variable)
            else:
                ax.plot(times,values,label=variable)

        ax.set_xlabel("Simulation time")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3)

        fig.tight_layout()
        plt.show()

        return fig, ax


@dataclass(frozen=True)
class SMCResult:
    """
    Result of a fixed-sample statistical model checking experiment.
    """
    successes: int
    runs: int
    elapsed_seconds: float
    confidence_level: float = 0.95

    def __post_init__(self) -> None:
        if self.runs <= 0:
            raise ValueError("runs must be positive.")

        if not 0 <= self.successes <= self.runs:
            raise ValueError("successes must be between 0 and runs.")

        if not 0.0 < self.confidence_level < 1.0:
            raise ValueError("confidence_level must be between 0 and 1.")
        
    @property
    def estimated_probability(self) -> float:
        if self.runs == 0:
            return float("nan")

        return self.successes / self.runs

    @property
    def confidence_interval(self) -> tuple[float, float]:
        """
        Return the two-sided Wilson confidence interval.
        """
        p_hat = self.estimated_probability
        n = self.runs

        z = NormalDist().inv_cdf(0.5 + self.confidence_level / 2)
        z_squared = z**2

        denominator = 1.0 + z_squared / n

        centre = (p_hat+ z_squared / (2.0 * n)) / denominator

        half_width = (z/denominator * (p_hat * (1.0 - p_hat) / n + z_squared / (4.0 * n**2)) ** 0.5)

        lower = max(0.0, centre - half_width)
        upper = min(1.0, centre + half_width)

        return lower, upper

    @property
    def confidence_lower(self) -> float:
        return self.confidence_interval[0]

    @property
    def confidence_upper(self) -> float:
        return self.confidence_interval[1]

    def plot(
        self,
        *,
        label: str = "BLTL property",
        title: str = "Statistical model checking result",
    ):
        """
        Plot the estimated probability and its confidence interval.

        Returns
        -------
        fig, ax: Matplotlib figure and axes.
        """
        probability = self.estimated_probability
        lower, upper = self.confidence_interval

        lower_error = probability - lower
        upper_error = upper - probability

        fig, ax = plt.subplots(figsize=(5.5, 4.5))

        ax.errorbar(
            [0],
            [probability],
            yerr=[
                [lower_error],
                [upper_error],
            ],
            fmt="o",
            markersize=7,
            capsize=6,
        )

        ax.set_xticks([0])
        ax.set_xticklabels([label])
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(0.0, 1.0)

        ax.set_ylabel("Estimated satisfaction probability")
        ax.set_title(title)
        ax.grid(axis="y",alpha=0.3)

        confidence_percentage = (100 * self.confidence_level)

        ax.annotate(
            (
                rf"$\hat{{p}}={probability:.3f}$"
                "\n"
                f"{confidence_percentage:.0f}% CI: "
                f"[{lower:.3f}, {upper:.3f}]"
            ),
            xy=(0, probability),
            xytext=(15, 0),
            textcoords="offset points",
            va="center",
        )

        fig.tight_layout()

        return fig, ax

    
# Complete Kappa SMC workflow
class KappaSMC:
    """
    Coordinate formula parsing, model instrumentation,
    simulation and BLTL monitoring.
    """

    def __init__(self,simulator: KappaSimulator) -> None:
        self.simulator = simulator

    def prepare(self, *, model_path: str, formula_source: str) -> PreparedProperty:
        """
        Parse the property and instrument the Kappa model.

        This is performed once before running one or more
        simulations.
        """
        formula = parse_formula(formula_source)
        required_horizon = formula_horizon(formula)

        if self.simulator.config.max_time < required_horizon:
            raise ValueError(
                "The requested simulation horizon is too short.\n"
                f"Formula horizon: {required_horizon:g}\n"
                f"Simulation max_time: {self.simulator.config.max_time:g}"
            )

        ap_names = extract_atomic_propositions(formula)

        if not ap_names:
            raise ValueError(
                "The formula does not contain any atomic propositions."
            )

        ap_triples = [parse_ap_string(name) for name in ap_names]

        return PreparedProperty(
                    formula=bind_ap_indices(formula,ap_names),
                    atomic_propositions=tuple(ap_names),
                    instrumented_model=instrument_kappa_model(model_path,ap_triples),
                )

    def run_trace(
        self,
        prepared: PreparedProperty,
        *,
        verbose_monitor: bool = False,
    ) -> TraceResult:
        """
        Execute and monitor one stochastic trajectory.
        """
        raw_results = self.simulator.simulate(prepared.instrumented_model)
        required_horizon = formula_horizon(prepared.formula)

        legend,series,times,ap_values = process_simulation_results(
            raw_results, 
            ap_count=len(prepared.atomic_propositions), 
            required_horizon=required_horizon,
            )

        monitor = BLTLMonitor(times,ap_values,verbose=verbose_monitor)

        holds = monitor.holds(prepared.formula)

        return TraceResult(
            holds=holds,
            times=times,
            ap_values=ap_values,
            atomic_propositions=(
                prepared.atomic_propositions
            ),
            legend=legend,
            series=series,
        )


    def run_smc(self, prepared: PreparedProperty, config: SMCConfig) -> SMCResult:
        """
        Estimate the satisfaction probability from independent
        stochastic simulations.
        """
        if config.runs <= 0:
            raise ValueError("The number of SMC runs must be positive.")

        if config.progress_every <= 0:
            raise ValueError("progress_every must be positive.")

        successes = 0
        start_time = time.perf_counter()

        for run_index in range(1,config.runs + 1):
            trace_result = self.run_trace(prepared,verbose_monitor=False)

            successes += int(trace_result.holds)

            report_progress = (
                run_index % config.progress_every == 0
                or run_index == config.runs
            )

            if report_progress:
                estimate = successes / run_index

                print(
                    f"[{run_index}/{config.runs}] "
                    f"p_hat = {estimate:.4f}"
                )

        elapsed_seconds = time.perf_counter() - start_time

        return SMCResult(
            successes=successes,
            runs=config.runs,
            elapsed_seconds=elapsed_seconds,
            confidence_level=config.confidence_level,
        )