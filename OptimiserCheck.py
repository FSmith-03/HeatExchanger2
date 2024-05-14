import numpy as np
from Hydraulic_Functions import *
from scipy.optimize import differential_evolution


def rastrigin_2d(x1, x2):
    A = 10
    y = A * 2 + (x1**2 - A * np.cos(2 * np.pi * x1)) + (x2**2 - A * np.cos(2 * np.pi * x2))
    return -y



def test_optimization_2d(bounds, num_runs):
    """
    Test optimization algorithm with 2D Rastrigin function.

    Parameters:
        optimizer (callable): Optimization algorithm function.
        bounds (list): List of tuples defining the bounds for x and y coordinates.
        num_runs (int): Number of optimization runs to perform (default=10).

    Returns:
        list: List of tuples containing the best solution found (x, y) and its corresponding value in each run.
    """
    for _ in range(num_runs):
        # Generate random initial solution within bounds
        initial_solution = [np.random.uniform(bounds[0][0], bounds[0][1]), np.random.uniform(bounds[1][0], bounds[1][1])]
        print("HI")
        # Run optimization algorithm
        result = differential_evolution(rastrigin_2d, bounds)
        optimal_params = result.x
        optimal_LMTD = -result.fun 
    
    return result

bounds = [(0,5), (0,5)]

print(test_optimization_2d(bounds, 100))