import numpy as np
import matplotlib.pyplot as plt
from viscofit import Custom_Model, LR_Maxwell_Force, LR_Voigt_Indentation, LR_Maxwell, LR_Voigt, LR_PowerLaw, row2mat, LR_PowerLaw_Force

import igor


def afm_force_signal(h, R):  # this is a fake force to simulate something that might have been obtained from an AFM
    return R * 5 * h ** 3

# initially define the experimental observables
times = np.linspace(0, 0.1, 1000)
indentation = times ** (3 / 2)
R = 1
forces = afm_force_signal(indentation, R)

# define the custom model and put in the experimental observables
model = Custom_Model(forces, times, indentation, R)

# define the function for the desired observable, note: the observables in the function MUST come from the model
def force_func(params):
    return model.radii * params[0] * model.indentation ** params[1]

# defining the bounds (n, 2) (two parameters in the force_func here so n=2)
# [[lower 1, upper 1],
#  [lower 2, upper 2]]
param_bounds = [1, 10, 1, 5]

# pass the training data (data to be replicated by the function being fit), the function, and the function parameter bounds
# to the fitting function
print(model.fit(training_data=forces, function=force_func, bounds=param_bounds))
