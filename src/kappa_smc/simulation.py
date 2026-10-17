from dataclasses import dataclass

import kappy as kp
import numpy as np


@dataclass(frozen=True)
class SimulationConfig:
    max_time: float = 100.0
    plot_period: float = 1.0

    @property
    def pause_condition(self) -> str:
        return f"[T] > {self.max_time}"


class KappaSimulator:
    """
    Thin wrapper around kappy for executing one Kappa simulation.
    """

    def __init__(self, kappa_bin_path: str, config: SimulationConfig) -> None:
        self.kappa_bin_path = kappa_bin_path
        self.config = config

    def simulate(self, model: str) -> dict:
        """
        Execute one stochastic simulation and return the KaSim plot.
        """
        client = kp.KappaStd(kappa_bin_path=self.kappa_bin_path)

        try:
            client.add_model_string(model)
            client.project_parse()

            parameters = kp.SimulationParameter(
                pause_condition=self.config.pause_condition,
                plot_period=self.config.plot_period,
            )

            client.simulation_start(parameters)
            client.wait_for_simulation_stop()

            results = client.simulation_plot()

        finally:
            client.shutdown()

        return results


def process_simulation_results(results: dict, ap_count: int, required_horizon: float | None = None) -> tuple[tuple[str, ...],np.ndarray,np.ndarray,np.ndarray]:
    """
    Sort a KaSim trajectory by time and extract the values needed for monitoring.

    If KaSim stops before the required formula horizon because the model has reached a terminal 
    state, extend the trajectory by repeating its final state at the required final time.

    Returns:
        legend: names of all simulation columns;

        series: complete simulation table, ordered by time;

        times: values of the [T] column;

        ap_values: boolean matrix containing p_0, ..., p_n.
    """

    legend = tuple(results["legend"])
    series = np.asarray(results["series"],dtype=float)

    if series.ndim != 2 or len(series) == 0:
        raise ValueError(
            "The simulation returned an empty or invalid trajectory."
        )

    time_index = legend.index("[T]")
    
    ap_indices = []

    for index in range(ap_count):
        observable_name = f"p_{index}"

        try:
            ap_indices.append(legend.index(observable_name))
        except ValueError as error:
            raise ValueError(
                f"Observable {observable_name!r} is missing.\n"
                f"Available columns: {list(legend)}"
            ) from error

    # Sort the complete table using the actual location of [T].
    order = np.argsort(series[:, time_index],kind="stable")
    series = series[order]

    # KaSim may stop early when no rule is enabled anymore.
    # In that case, the final state persists indefinitely.
    if required_horizon is not None:
        initial_time = series[0, time_index]
        required_final_time = initial_time + required_horizon
        actual_final_time = series[-1, time_index]

        if actual_final_time < required_final_time:
            final_row = series[-1].copy()
            final_row[time_index] = required_final_time

            series = np.vstack([series, final_row])

    times = series[:, time_index]
    ap_values = series[:, ap_indices] != 0

    return legend, series, times, ap_values
