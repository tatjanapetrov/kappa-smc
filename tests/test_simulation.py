"""Unit tests for KaSim plot post-processing (no KaSim process needed)."""

import numpy as np
import pytest

from kappa_smc.simulation import process_simulation_results


def test_extracts_aps_by_name_not_position_and_sorts_time():
    results = {
        "legend": ["p_1", "[T]", "count", "p_0"],
        "series": [
            [1, 2.0, 7, 0],
            [0, 0.0, 3, 1],
            [1, 1.0, 5, 1],
        ],
    }
    legend, series, times, ap_values = process_simulation_results(results, ap_count=2)

    assert legend == ("p_1", "[T]", "count", "p_0")
    np.testing.assert_array_equal(times, [0, 1, 2])
    np.testing.assert_array_equal(series[:, 2], [3, 5, 7])
    np.testing.assert_array_equal(
        ap_values, [[True, False], [True, True], [False, True]]
    )
    assert ap_values.dtype == np.bool_


def test_extends_terminal_state_to_required_horizon():
    results = {
        "legend": ["[T]", "p_0", "count"],
        "series": [[0.0, 0, 2], [1.0, 1, 4]],
    }
    _, series, times, ap_values = process_simulation_results(
        results, ap_count=1, required_horizon=3.0
    )

    np.testing.assert_array_equal(times, [0, 1, 3])
    np.testing.assert_array_equal(series[-1, 1:], series[-2, 1:])
    np.testing.assert_array_equal(ap_values[:, 0], [False, True, True])


def test_does_not_extend_a_sufficiently_long_trace():
    results = {"legend": ["[T]", "p_0"], "series": [[0, 0], [5, 1]]}
    _, series, times, _ = process_simulation_results(
        results, ap_count=1, required_horizon=3.0
    )
    assert len(series) == 2
    np.testing.assert_array_equal(times, [0, 5])


@pytest.mark.parametrize("series", [[], [1, 2, 3]])
def test_rejects_empty_or_non_matrix_trace(series):
    with pytest.raises(ValueError, match="empty or invalid"):
        process_simulation_results({"legend": ["[T]", "p_0"], "series": series}, 1)


def test_missing_ap_column_is_reported():
    results = {"legend": ["[T]", "p_0"], "series": [[0, 1], [1, 0]]}
    with pytest.raises(ValueError, match="p_1.*missing"):
        process_simulation_results(results, ap_count=2)
