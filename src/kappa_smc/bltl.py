from dataclasses import dataclass
from typing import Optional
import numpy as np
import ast
import re

# BLTL AST nodes
class Formula:
    """Base class for BLTL formulas."""

@dataclass(frozen=True)
class AP(Formula):
    # Before binding: idx=None and name="I<2"
    # After binding:  idx is 0..m-1 matching ap_vals columns
    name: str
    idx: Optional[int] = None

@dataclass(frozen=True)
class Not(Formula):
    phi: Formula

@dataclass(frozen=True)
class And(Formula):
    left: Formula
    right: Formula

@dataclass(frozen=True)
class Or(Formula):
    left: Formula
    right: Formula

@dataclass(frozen=True)
class G_leq(Formula):
    bound: int
    phi: Formula

@dataclass(frozen=True)
class F_leq(Formula):
    bound: int
    phi: Formula

@dataclass(frozen=True)
class Until_leq(Formula):
    bound: int
    left: Formula
    right: Formula


# BLTL Monitor
class BLTLMonitor:
    """
    Evaluate a BLTL formula over one sampled stochastic trace.

    Parameters
    ----------
    times: ordered simulation times.

    ap_values: boolean matrix with shape:
        number of time points x number of atomic propositions

    verbose: print a summary of each evaluated subformula.
    """

    TIME_TOLERANCE = 1e-10

    def __init__(self,times: np.ndarray,ap_values: np.ndarray,*,verbose: bool = False) -> None:
        self.times = np.asarray(times, dtype=float)
        self.ap_values = np.asarray(ap_values, dtype=bool)
        self.verbose = verbose
        self._cache: dict[Formula, np.ndarray] = {}

    def holds(self, formula: Formula) -> bool:
        """
        Return whether the formula holds at the initial trace time.
        """
        required_horizon = formula_horizon(formula)
        available_horizon = self.times[-1] - self.times[0]

        if (available_horizon + self.TIME_TOLERANCE < required_horizon):
            raise ValueError(
                "The simulation trace is too short for this formula.\n"
                f"Required horizon: {required_horizon:g}\n"
                f"Available horizon: {available_horizon:g}"
            )

        satisfaction = self.evaluate(formula)
        result = bool(satisfaction[0])

        if self.verbose:
            print(f"\nFormula holds at t={self.times[0]:g}: {result}")

        return result

    def evaluate(self, formula: Formula) -> np.ndarray:
        """
        Return a Boolean vector whose entry at index i says whether
        the formula holds starting at times[i].
        """
        cached = self._cache.get(formula)

        if cached is not None:
            return cached

        if isinstance(formula, AP):
            result = self._evaluate_ap(formula)

        elif isinstance(formula, Not):
            result = np.logical_not(self.evaluate(formula.phi))

        elif isinstance(formula, And):
            result = np.logical_and(
                self.evaluate(formula.left),
                self.evaluate(formula.right),
            )

        elif isinstance(formula, Or):
            result = np.logical_or(
                self.evaluate(formula.left),
                self.evaluate(formula.right),
            )

        elif isinstance(formula, G_leq):
            result = self._evaluate_globally(formula)

        elif isinstance(formula, F_leq):
            result = self._evaluate_finally(formula)

        elif isinstance(formula, Until_leq):
            result = self._evaluate_until(formula)

        else:
            raise TypeError(
                f"Unsupported formula node: {type(formula).__name__}"
            )
        
        self._cache[formula] = result

        if self.verbose:
            true_count = int(np.count_nonzero(result))

            print(
                f"{type(formula).__name__}: "
                f"{true_count}/{len(result)} time points true"
            )

        return result

    def _evaluate_ap(self, formula: AP) -> np.ndarray:
        if formula.idx is None:
            raise ValueError(
                f"Atomic proposition {formula.name!r} "
                "has not been associated with a trace column."
            )

        if not 0 <= formula.idx < self.ap_values.shape[1]:
            raise IndexError(
                f"Atomic proposition index {formula.idx} "
                "is outside the monitoring matrix."
            )

        return self.ap_values[:, formula.idx].copy()

    def _window_stop(
        self,
        start_index: int,
        bound: float,
    ) -> int | None:
        """
        Return the exclusive index delimiting the interval

            [times[start_index], times[start_index] + bound].

        Return None if the trace ends before the full interval.
        """
        if bound < 0:
            raise ValueError("Temporal bounds must be non-negative.")

        deadline = self.times[start_index] + bound

        if (self.times[-1] + self.TIME_TOLERANCE < deadline):
            return None

        return int(
            np.searchsorted(
                self.times,
                deadline + self.TIME_TOLERANCE,
                side="right",
            )
        )

    def _evaluate_globally(self,formula: G_leq) -> np.ndarray:
        inner = self.evaluate(formula.phi)
        result = np.zeros(len(self.times), dtype=bool)

        for start_index in range(len(self.times)):
            stop_index = self._window_stop(
                start_index,
                formula.bound,
            )

            if stop_index is None:
                result[start_index] = False
                continue

            result[start_index] = bool(
                np.all(inner[start_index:stop_index])
            )

        return result

    def _evaluate_finally(self,formula: F_leq) -> np.ndarray:
        inner = self.evaluate(formula.phi)
        result = np.zeros(len(self.times), dtype=bool)

        for start_index in range(len(self.times)):
            stop_index = self._window_stop(
                start_index,
                formula.bound,
            )

            if stop_index is None:
                result[start_index] = False
                continue

            result[start_index] = bool(
                np.any(inner[start_index:stop_index])
            )

        return result

    def _evaluate_until(self,formula: Until_leq) -> np.ndarray:
        left = self.evaluate(formula.left)
        right = self.evaluate(formula.right)

        result = np.zeros(len(self.times), dtype=bool)

        for start_index in range(len(self.times)):
            stop_index = self._window_stop(
                start_index,
                formula.bound,
            )

            if stop_index is None:
                # A witness may still make the formula true before
                # the end of the incomplete interval.
                stop_index = len(self.times)

            for witness_index in range(start_index,stop_index):
                right_holds = right[witness_index]

                left_holds_before_witness = (
                    witness_index == start_index
                    or np.all(
                        left[start_index:witness_index]
                    )
                )

                if (right_holds and left_holds_before_witness):
                    result[start_index] = True
                    break

        return result


# Formula parsing + AP extraction

# supports e.g. "I<2", "S>50", "R<10"
_AP_RE = re.compile(r"^\s*([A-Za-z_]\w*)\s*(<|>)\s*([0-9]+)\s*$")

def parse_ap_string(ap_str: str) -> tuple[str, str, str]:
    m = _AP_RE.match(ap_str)
    if not m:
        raise ValueError(f"Bad AP string '{ap_str}'. Use like 'I<2' or 'S>50' (only < or >, int threshold).")
    var, op, c = m.group(1), m.group(2), m.group(3)
    return var, op, c

def parse_formula(formula_str: str):
    """
    Parse a formula string like:
      F_leq(50, AP("I<2"))
      G_leq(50, F_leq(20, AP("R>10")))
      And(F_leq(10, AP("I>10")), F_leq(50, AP("I<2")))

    Allowed constructors: AP, Not, And, Or, G_leq, F_leq, Until_leq
    """
    allowed = {"AP", "Not", "And", "Or", "G_leq", "F_leq", "Until_leq"}

    def conv(node):
        if isinstance(node, ast.Constant):
            return node.value

        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            fname = node.func.id
            if fname not in allowed:
                raise ValueError(f"Constructor '{fname}' not allowed.")
            args = [conv(a) for a in node.args]

            if fname == "AP":
                if len(args) != 1 or not isinstance(args[0], str):
                    raise ValueError('AP must be AP("I<2") (string).')
                # validate early
                parse_ap_string(args[0])
                return AP(idx=None, name=args[0])

            if fname == "Not":
                return Not(args[0])

            if fname == "And":
                return And(args[0], args[1])

            if fname == "Or":
                return Or(args[0], args[1])

            if fname == "G_leq":
                return G_leq(int(args[0]), args[1])

            if fname == "F_leq":
                return F_leq(int(args[0]), args[1])

            if fname == "Until_leq":
                return Until_leq(int(args[0]), args[1], args[2])

        raise ValueError(f"Unsupported syntax: {ast.dump(node)}")

    tree = ast.parse(formula_str, mode="eval")
    return conv(tree.body)


# Formula traversal utilities

def extract_atomic_propositions(formula: Formula) -> list[str]:
    """
    Extract the atomic propositions occurring in a formula.

    Duplicate propositions are removed while preserving the order
    of their first occurrence.
    """
    found: list[str] = []

    def visit(node: Formula) -> None:
        if isinstance(node, AP):
            found.append(node.name)
            return

        if isinstance(node, Not):
            visit(node.phi)
            return

        if isinstance(node, (And, Or)):
            visit(node.left)
            visit(node.right)
            return

        if isinstance(node, (G_leq, F_leq)):
            visit(node.phi)
            return

        if isinstance(node, Until_leq):
            visit(node.left)
            visit(node.right)
            return

        raise TypeError(
            f"Unsupported formula node: {type(node).__name__}"
        )

    visit(formula)

    unique: list[str] = []
    seen: set[str] = set()

    for proposition in found:
        if proposition not in seen:
            seen.add(proposition)
            unique.append(proposition)

    return unique


def bind_ap_indices(
    formula: Formula,
    atomic_propositions: list[str],
) -> Formula:
    """
    Return a copy of the formula where each AP is associated with
    the index of the corresponding column in the monitoring matrix.
    """
    index_by_name = {
        proposition: idx
        for idx, proposition in enumerate(atomic_propositions)
    }

    def bind(node: Formula) -> Formula:
        if isinstance(node, AP):
            try:
                idx = index_by_name[node.name]
            except KeyError as error:
                raise ValueError(
                    f"Atomic proposition {node.name!r} was not registered."
                ) from error

            return AP(name=node.name,idx=idx)

        if isinstance(node, Not):
            return Not(phi=bind(node.phi))

        if isinstance(node, And):
            return And(left=bind(node.left),right=bind(node.right))

        if isinstance(node, Or):
            return Or(left=bind(node.left),right=bind(node.right))

        if isinstance(node, G_leq):
            return G_leq(bound=node.bound,phi=bind(node.phi))

        if isinstance(node, F_leq):
            return F_leq(bound=node.bound,phi=bind(node.phi))

        if isinstance(node, Until_leq):
            return Until_leq(bound=node.bound,left=bind(node.left),right=bind(node.right))

        raise TypeError(f"Unsupported formula node: {type(node).__name__}")

    return bind(formula)


def formula_horizon(formula: Formula) -> float:
    """
    Compute the minimum trace horizon needed to evaluate the formula
    at the initial time.

    Examples:
        F<=50 AP       -> horizon 50
        F<=35 G<=40 AP -> horizon 75
    """
    if isinstance(formula, AP):
        return 0.0

    if isinstance(formula, Not):
        return formula_horizon(formula.phi)

    if isinstance(formula, (And, Or)):
        return max(
            formula_horizon(formula.left),
            formula_horizon(formula.right),
        )

    if isinstance(formula, (G_leq, F_leq)):
        return (
            formula.bound
            + formula_horizon(formula.phi)
        )

    if isinstance(formula, Until_leq):
        return (
            formula.bound
            + max(
                formula_horizon(formula.left),
                formula_horizon(formula.right),
            )
        )

    raise TypeError(
        f"Unsupported formula node: {type(formula).__name__}"
    )