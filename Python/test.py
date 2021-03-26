import numpy as np
from general import Custom_Model

# this is a 'fake' force to simulate something that we get from the AFM
def force_basic(h, R):
    return R * 5 * h ** 3

# initially define the experimental observables
times = [np.linspace(0, 0.1, 1000), np.linspace(0, 1, 200), np.linspace(0, 4, 500)]
indentation = [time ** (3 / 2) for time in times]
R = [1, 0.5, 2]
forces = [force_basic(h, r) for h, r in zip(indentation, R)]

# define the custom model and put in the experimental observables
a = Custom_Model(forces, times, indentation, R)

# define the function for the desired observable
def force_func(params):
    return a.radii * params[0] * a.indentation ** params[1]

# set the target observable as one of the observables within the custom model
a.target_observable = a.force
# set the observable function as the previously defined function
a.observable_function = force_func

# fit the function
a.fit()