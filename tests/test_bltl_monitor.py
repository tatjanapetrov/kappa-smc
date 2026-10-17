"""Unit tests for the Kappa-SMC BLTL monitor, independent of KaSim.

Place at tests/test_bltl_monitor.py and run: python -m pytest -q tests/test_bltl_monitor.py
Trace values are sampled at the timestamps listed; each AP is one matrix column.
"""

import numpy as np
import pytest

from kappa_smc.bltl import (
    AP, And, BLTLMonitor, F_leq, G_leq, Not, Or, Until_leq,
    bind_ap_indices, extract_atomic_propositions, formula_horizon, parse_formula,
)


P = AP("P>0", idx=0)
Q = AP("Q>0", idx=1)


def monitor(times, *columns):
    """Build one monitor; each column lists an AP's truth values over time."""
    ts = np.asarray(times, dtype=float)
    values = np.column_stack(columns).astype(bool)
    assert values.shape[0] == len(ts)
    return BLTLMonitor(ts, values)


@pytest.mark.parametrize(
    "values, expected",
    [([False, True, False], [False, True, False]),
     ([True, True, False], [True, True, False])],
)
def test_atomic_proposition(values, expected):
    m = monitor([0, 1, 2], values)
    np.testing.assert_array_equal(m.evaluate(P), expected)
    assert m.holds(P) == expected[0]


def test_boolean_operators():
    m = monitor([0, 1, 2], [True, False, True], [False, False, True])
    np.testing.assert_array_equal(m.evaluate(Not(P)), [False, True, False])
    np.testing.assert_array_equal(m.evaluate(And(P, Q)), [False, False, True])
    np.testing.assert_array_equal(m.evaluate(Or(P, Q)), [True, False, True])


@pytest.mark.parametrize(
    "values, expected_f, expected_g",
    [
        ([False, False, True], True, False),  # witness exactly at inclusive deadline
        ([True, True, True], True, True),
        ([False, False, False], False, False),
        ([True, False, True], True, False),
    ],
)
def test_finally_globally(values, expected_f, expected_g):
    m = monitor([0, 1, 2], values)
    assert m.holds(F_leq(2, P)) is expected_f
    assert m.holds(G_leq(2, P)) is expected_g


def test_deadline_is_inclusive_but_later_events_are_excluded():
    m = monitor([0, 1, 1.99, 2, 2.01, 3], [False, False, False, True, False, False])
    assert m.holds(F_leq(2, P)) is True
    m2 = monitor([0, 1, 2, 2.01, 3], [False, False, False, True, False])
    assert m2.holds(F_leq(2, P)) is False


def test_irregular_physical_times_not_sample_indices():
    m = monitor([0.0, 0.1, 0.2, 3.0], [False, False, True, False])
    assert m.holds(F_leq(1, P)) is True
    assert m.holds(G_leq(1, Not(P))) is False


def test_zero_bound_checks_current_state_only():
    m = monitor([0, 1], [True, False])
    assert m.holds(F_leq(0, P)) is True
    assert m.holds(G_leq(0, P)) is True
    assert m.holds(F_leq(0, Not(P))) is False


def test_until_succeeds_when_left_holds_before_right():
    m = monitor([0, 1, 2, 3],
                [True, True, False, False],
                [False, False, True, False])
    assert m.holds(Until_leq(2, P, Q)) is True


def test_until_fails_if_left_becomes_false_before_witness():
    m = monitor([0, 1, 2, 3],
                [True, False, True, False],
                [False, False, True, False])
    assert m.holds(Until_leq(2, P, Q)) is False


def test_until_right_true_immediately_does_not_require_left():
    m = monitor([0, 1], [False, False], [True, False])
    assert m.holds(Until_leq(0, P, Q)) is True


def test_until_deadline_exclusion():
    m = monitor([0, 1, 2, 3],
                [True, True, True, True],
                [False, False, False, True])
    assert m.holds(Until_leq(2, P, Q)) is False
    assert m.holds(Until_leq(3, P, Q)) is True


def test_nested_eventually_globally_success():
    # F<=2 G<=2 P: find a start <=2 with P continuously true for 2 units.
    m = monitor([0, 1, 2, 3, 4], [False, False, True, True, True])
    phi = F_leq(2, G_leq(2, P))
    assert formula_horizon(phi) == 4
    assert m.holds(phi) is True


def test_nested_eventually_globally_failure():
    m = monitor([0, 1, 2, 3, 4], [False, True, True, False, True])
    assert m.holds(F_leq(2, G_leq(2, P))) is False


def test_duality_of_eventually_and_globally():
    # On a complete window, not(F P) = G(not P).
    m = monitor([0, 1, 2, 3], [False, False, True, False])
    left = Not(F_leq(2, P))
    right = G_leq(2, Not(P))
    assert m.holds(left) == m.holds(right)


def test_incomplete_trace_rejected():
    m = monitor([0, 1], [True, True])
    with pytest.raises(ValueError, match="too short"):
        m.holds(F_leq(2, P))
    with pytest.raises(ValueError, match="too short"):
        m.holds(F_leq(1, G_leq(1, P)))


def test_unbound_or_invalid_atomic_proposition():
    m = monitor([0, 1], [True, False])
    with pytest.raises(ValueError, match="not been associated"):
        m.holds(AP("P>0"))
    with pytest.raises(IndexError, match="outside"):
        m.holds(AP("R>0", idx=1))


def test_parser_extraction_and_binding():
    phi = parse_formula('And(F_leq(2, AP("P>0")), G_leq(1, AP("Q>0")))')
    names = extract_atomic_propositions(phi)
    assert names == ["P>0", "Q>0"]
    bound = bind_ap_indices(phi, names)
    m = monitor([0, 1, 2], [False, True, False], [True, True, True])
    assert m.holds(bound) is True
