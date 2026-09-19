"""Gauss-Legendre integration rules on the interval ``[-1, 1]``."""

import numpy as np


class GaussQuadrature:
    """Provide the one-, two-, and three-point Gauss-Legendre rules."""
    
    @staticmethod
    def get_rule(n_points: int) -> tuple[list[float], list[float]]:
        """Return integration points and weights for ``n_points``."""
        if n_points == 1:
            return [0.0], [2.0]

        if n_points == 2:
            value = 1.0 / np.sqrt(3.0)
            return [-value, value], [1.0, 1.0]

        if n_points == 3:
            value = np.sqrt(0.6)
            return [-value, 0.0, value], [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]

        raise NotImplementedError(
            f"Gauss rule of order {n_points} is not implemented; "
            "only one to three integration points are available."
        )
