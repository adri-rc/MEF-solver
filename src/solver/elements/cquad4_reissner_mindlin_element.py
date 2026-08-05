"""Classical isoparametric Reissner-Mindlin CQUAD4 formulation.

Membrane, bending and transverse-shear energies are integrated with a 2x2
Gauss rule.  This is the consistent displacement formulation of the theory;
as expected for a low-order Reissner-Mindlin element, very thin plates can
exhibit transverse-shear locking.
"""

import numpy as np

from ..quadrature import GaussQuadrature
from .cquad4_element import CQUAD4ElementBase


class CQUAD4ReissnerMindlinElementSolver(CQUAD4ElementBase):
    """Four-node shell using the classical Reissner-Mindlin formulation."""

    formulation = "reissner-mindlin"

    def _build_local_stiffness(self) -> np.ndarray:
        """Integrate the classical Reissner-Mindlin energy with a 2x2 rule."""
        stiffness = np.zeros((24, 24))
        membrane = np.zeros((8, 8))
        plate = np.zeros((12, 12))

        points, weights = GaussQuadrature.get_rule(2)
        membrane_constitutive = self._get_D_membrane()
        bending_constitutive = self._get_D_bending()
        shear_constitutive = self._get_D_shear()

        thickness = self.prop.t
        bending_inertia = thickness**3 / 12.0
        shear_area = (5.0 / 6.0) * thickness

        for i, xi in enumerate(points):
            for j, eta in enumerate(points):
                B_membrane, B_bending, B_shear, det_jacobian = (
                    self._get_B_matrices(xi, eta)
                )
                differential_area = det_jacobian * weights[i] * weights[j]
                membrane += (
                    B_membrane.T @ membrane_constitutive @ B_membrane
                ) * (differential_area * thickness)
                if thickness > 0.0:
                    plate += (
                        B_bending.T @ bending_constitutive @ B_bending
                    ) * (differential_area * bending_inertia)
                    plate += (
                        B_shear.T @ shear_constitutive @ B_shear
                    ) * (differential_area * shear_area)

        membrane_dofs = [0, 1, 6, 7, 12, 13, 18, 19]
        plate_dofs = [2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22]
        self._scatter_block(stiffness, membrane, membrane_dofs)
        if thickness > 0.0:
            self._scatter_block(stiffness, plate, plate_dofs)
        self._add_drilling_stiffness(stiffness, points, weights)
        return 0.5 * (stiffness + stiffness.T)
