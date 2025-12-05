import numpy as np

class GaussQuadrature:
    """
    @class GaussQuadrature
    @brief Gauss Quadrature provider of numerical integration schemes (Points and Weights).
    """
    
    @staticmethod
    def get_rule(n_points: int) -> tuple[list[float], list[float]]:
        """
        @brief Returns the points (xi) and weights (w) for the Gauss-Legendre quadrature
        on the natural interval [-1, 1].
        @details
        The number of integration points required (n_gp) to approximate a polynomial of degree 'p' is:
        p<=2*n_gp-1.

        SEE: TABLE 4.1 POSITION OF GAUSS POINTS AND CORRESPONDING WEIGHTS.
        (Ref: J. Fish & T. Belytschko)
        
        @param n_points user wanted number of integration points.
        @return Tuple ([chi_coord_int_points], [wigths])
        """
        if n_points == 1:
            
            return ([0.0], [2.0]) # Exact solution for degree 1 polynomials (linear)
        
        elif n_points == 2: # Exact solution for degree 3 polynomials (cubic)
            val = 1.0 / np.sqrt(3.0) # +/- 0.57735...
            return ([-val, val], [1.0, 1.0])
        
        elif n_points == 3: # Exact solution for degree 5 polynomials
            val = np.sqrt(0.6) # +/- 0.77459...
            return ([-val, 0.0, val], [5.0/9.0, 8.0/9.0, 5.0/9.0])
            
        else:
            raise NotImplementedError(f"Gauss's rule of order {n_points} is not implemented."
                                      "\nCurrently, only up to 3 integration points have been defined.")
