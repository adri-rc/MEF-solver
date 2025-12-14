import pytest
import numpy as np
from src.solver.quadrature import GaussQuadrature

"""
@file test_quadrature.py
@brief Unit tests for the GaussQuadrature numerical integration module.
@details
This test suite validates the accuracy of the Gauss-Legendre rules provided by the solver.
Validation strategy:
    1. Check that points and weights match theoretical tabular values.
    2. Integrate known polynomials analytically vs numerically.
    3. Verify error handling for unsupported integration orders.

Covered cases:
    - Order 1: Linear functions.
    - Order 2: Cubic functions.
    - Order 3: Quintic functions.
"""

class TestGaussQuadrature:
    
    def test_order_1_exact_linear(self):
        """
        @test 1 Gauss Integration Point
        @brief Exact integration up to polynomial degree 1. (linear)
        @details 
        Test function: f(x) = x + 2 over [-1, 1].
        Analytical result: [x^2/2 + 2x]_{-1}^{1} = (0.5 + 2) - (0.5 - 2) = 4.0
        """
        n_gauss_integration_points=1
        points, weights = GaussQuadrature.get_rule(n_gauss_integration_points)
        
        assert points == [0.]
        assert weights == [2.]
        
        integral = 0.0 # Initialize
        for xi, w in zip(points, weights):
            f_xi = xi + 2.0
            integral += f_xi * w
            
        assert integral == pytest.approx(4.0)

    def test_order_2_exact_cubic(self):
        """
        @test 2 Gauss Integration Points
        @brief Exact integration up to polynomial degree 3.
        @details
        Exact integration of: f(x) = x^3 en [-1, 1].
        Analytical result: 0.0 (symmetric odd function).
        """
        points, weights = GaussQuadrature.get_rule(2)
        # Check xi points and weigths
        assert points == pytest.approx([-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)])
        assert weights == pytest.approx([1.0, 1.0])
        
        integral = 0.0 # Initialize
        for xi, w in zip(points, weights):
            f_xi = xi**3
            integral += f_xi * w
            
        assert integral == pytest.approx(0.0)
    
    def test_order_3_exact_quintic(self):
        """
        @test 3 Gauss Integration Points
        @brief Exact integration up to polynomial degree 5.
        @details
        Test function: f(x) = 1 + x^5 (Degree 4)
        Analytical result: 2.
        """
        points, weights = GaussQuadrature.get_rule(3)
        
        # Check xi points and weigths
        expected_points = [-np.sqrt(0.6), 0.0, np.sqrt(0.6)] 
        assert sorted(points) == pytest.approx(expected_points)
        expected_weights = [5.0/9.0, 5.0/9.0, 8.0/9.0]
        assert sorted(weights) == pytest.approx(expected_weights)

        integral = 0.0 # Initialize
        for xi, w in zip(points, weights):
            f_xi = 1+xi**5
            integral += f_xi * w
            
        assert integral == pytest.approx(2.)

    def test_not_implemented_order(self):
        """
        @test Error GaussLegendre integration order not suported.
        """
        with pytest.raises(NotImplementedError):
            GaussQuadrature.get_rule(10)









