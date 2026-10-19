"""Tests for the Gauss-Legendre integration rules used by the solver."""

import numpy as np
import pytest

from src.solver.quadrature import GaussQuadrature


class TestGaussQuadrature:
    def test_order_1_exact_linear(self):
        """One point integrates a linear polynomial exactly over [-1, 1]."""
        n_gauss_integration_points = 1
        points, weights = GaussQuadrature.get_rule(n_gauss_integration_points)

        assert points == [0.]
        assert weights == [2.]

        integral = 0.0
        for xi, w in zip(points, weights):
            f_xi = xi + 2.0
            integral += f_xi * w

        assert integral == pytest.approx(4.0)

    def test_order_2_exact_cubic(self):
        """Two points integrate the odd cubic x**3 exactly over [-1, 1]."""
        points, weights = GaussQuadrature.get_rule(2)
        assert points == pytest.approx([-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)])
        assert weights == pytest.approx([1.0, 1.0])

        integral = 0.0
        for xi, w in zip(points, weights):
            f_xi = xi**3
            integral += f_xi * w

        assert integral == pytest.approx(0.0)

    def test_order_3_exact_quintic(self):
        """Three points integrate 1 + x**5 exactly over [-1, 1]."""
        points, weights = GaussQuadrature.get_rule(3)

        expected_points = [-np.sqrt(0.6), 0.0, np.sqrt(0.6)]
        assert sorted(points) == pytest.approx(expected_points)
        expected_weights = [5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0]
        assert sorted(weights) == pytest.approx(expected_weights)

        integral = 0.0
        for xi, w in zip(points, weights):
            f_xi = 1 + xi**5
            integral += f_xi * w

        assert integral == pytest.approx(2.0)

    def test_not_implemented_order(self):
        with pytest.raises(NotImplementedError):
            GaussQuadrature.get_rule(10)

