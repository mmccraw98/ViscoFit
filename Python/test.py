import numpy as np
import matplotlib.pyplot as plt
from general import Custom_Model, LR_Maxwell_Force, LR_Voigt_Indentation, LR_Maxwell, LR_Voigt, LR_PowerLaw, row2mat, LR_PowerLaw_Force

def afm_force_signal(h, R):  # this is a fake force to simulate something that might have been obtained from an AFM
    return R * 5 * h ** 3

# initially define the experimental observables
times = [np.linspace(0, 0.1, 1000), np.linspace(0, 1, 200), np.linspace(0, 4, 500)]
indentation = [time ** (3 / 2) for time in times]
R = [1, 0.5, 2]
forces = [afm_force_signal(h, r) for h, r in zip(indentation, R)]

# define the custom model and put in the experimental observables
model = Custom_Model(forces, times, indentation, R, target_observable=forces)  # we are attempting to emulate the force

# define the function for the desired observable, note: the observables in the function MUST come from the model
def force_func(params):
    return model.radii * params[0] * model.indentation ** params[1]

# define the function for the model
model.observable_function = force_func

# defining the bounds (n, 2) (two parameters in the force_func here so n=2)
# [[lower 1, upper 1],
#  [lower 2, upper 2]]
param_bounds = np.array([[1, 10], [1, 6]])

print(model.fit(guess=np.array([1, 10])))

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
