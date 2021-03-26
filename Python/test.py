import numpy as np
import matplotlib.pyplot as plt
from general import Custom_Model, LR_Maxwell_Force, LR_Voigt_Indentation, LR_Maxwell, LR_Voigt, LR_PowerLaw, row2mat

relaxance = np.array([1e4, 1e5, 1e-2])
retardance = np.array([1e-10, 1e-9, 1e-5])

times = [np.linspace(0, 1, 500) for i in range(2)]
indentations = [time / max(time) * 10 * (i + 1) * 10e-9 for i, time in enumerate(times)]
radii = [20 * (i + 1) * 1e-9 for i in range(len(times))]
forces = [LR_Maxwell_Force(relaxance, t, h, r) for t, h, r in zip(times, indentations, radii)]

a = LR_Voigt(forces, times, indentations, radii)

print(a.fit_slow())

quit()
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