import numpy as np
from scipy.optimize import differential_evolution, basinhopping, brute

def rastrigin_2d(x):
    A = 10
    x1, x2 = x
    y = A * 2 + (x1**2 - A * np.cos(2 * np.pi * x1)) + (x2**2 - A * np.cos(2 * np.pi * x2))
    return y

def sphere_function(x):
    x1, x2 = x
    y = x1**2 +x2**2
    return y

def test_optimization_2d(bounds, num_runs):
    """
    Test optimization algorithm with 2D Rastrigin function.

    Parameters:
        bounds (list): List of tuples defining the bounds for x and y coordinates.
        num_runs (int): Number of optimization runs to perform (default=10).

    Returns:
        list: List of tuples containing the best solution found (x, y) and its corresponding value in each run.
    """
    results = []
    for _ in range(num_runs):
        # Generate random initial solution within bounds
        initial_solution = [np.random.uniform(bounds[0][0], bounds[0][1]), np.random.uniform(bounds[1][0], bounds[1][1])]
        
        # Run optimization algorithm
        result = differential_evolution(sphere_function, bounds)
        results.append(result)
    
    return results
bounds = [(-5,5), (-5,5)]
print(test_optimization_2d(bounds, 1))
#result_bh = basinhopping(rastrigin_2d, x0=np.random.uniform(0, 5, size=2), minimizer_kwargs={"bounds": bounds})
#print("Basin-Hopping Result:")
#print(result_bh)
#result_brute = brute(rastrigin_2d, ranges=bounds, Ns=1000, finish=None)
#min_value = rastrigin_2d(result_brute) 
#print(min_value, "At location:", result_brute)