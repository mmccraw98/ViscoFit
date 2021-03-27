from numpy import array, sum, sqrt, convolve, exp, ones, zeros, tile, insert, concatenate, argmin, diff, var
from numpy.random import uniform
from scipy.optimize import minimize, dual_annealing
from time import time


def row2mat(row, n):
    '''
    stacks a row vector (numpy (m, )) n times to create a matrix (numpy (m, n)) NOTE: CAN SLOW DOWN COMPUTATION IF DONE MANY TIMES
    :param row: numpy array row vector
    :param n: int number of replications to perform
    :return: numpy matrix (m, n) replicated row vector
    '''
    # do once at the beginning of any calculation to improve performance
    return tile(row, (n, 1)).T


def tic():
    '''
    start the clock
    :return: a global var for the initial time
    '''
    global current_time
    current_time = time()


def toc(return_numeric=False):
    '''
    stop the clock and return the elapsed time
    :param return_numeric: bool whether or not to return a numeric value or a report string
    :return: numeric value of elapsed time or a report string for elapsed time
    '''
    if return_numeric:
        return time() - current_time
    else:
        print('process completed in {:.2f}s'.format(time() - current_time))


def LR_Maxwell_Force(model_params, time, indentation, radius):
    '''
    calculates the response force for a generalized maxwell model according to the lee and radok contact mechanics formulation
    :param model_params: numpy array contains the model stiffnesses and time constants in either of the two following forms:
    to incorporate steady state fluidity: array([Elastic Stiffness, Fluidity Time Constant, Arm 1 Stiffness, Arm 1 Time Constant, ... Arm N Stiffness, Arm N Time Constant])
    normal model: array([Elastic Stiffness, Arm 1 Stiffness, Arm 1 Time Constant, ... Arm N Stiffness, Arm N Time Constant])
    :param time: (n,) numpy array time signal of the loading cycle
    :param indentation: (n,) numpy array indentation signal of the loading cycle
    :param radius: float radius of the indenter tip
    :return: (n,) numpy array force signal of the loading cycle
    '''
    if model_params.size % 2 == 0:  # fluidity case:
        time_matrix = row2mat(time, model_params[0::2].size)
        relaxance = - sum(model_params[0::2] / model_params[1::2] * exp(- time_matrix / model_params[1::2]), axis=1)
        relaxance += sum(model_params[0::2]) / (time[1] - time[0])  # add the delta function to the relaxances
        # must divide by dt since the integral of dirac delta MUST be 1 by definiton
    else:  # no fluidity case
        time_matrix = row2mat(time, model_params[1::2].size)
        relaxance = - sum(model_params[1::2] / model_params[2::2] * exp(- time_matrix / model_params[2::2]), axis=1)
        relaxance += (model_params[0] + sum(model_params[1::2])) / (time[1] - time[0])  # add the delta function of the relaxances
        # must divide by dt since the integral of dirac delta MUST be 1 by definiton
    return 16 * sqrt(radius) / 3 * convolve(relaxance, indentation ** (3 / 2), mode='full')[:time.size] * (time[1] - time[0])


def LR_Voigt_Indentation(model_params, time, force, radius):
    '''
    calculates the response indentation for a kelvin-voigt model according to the lee and radok contact mechanics formulation
    :param model_params: numpy array contains the model compliances and time constants in either of the two following forms:
    to incorporate steady state fluidity: array([Elastic Compliance, Fluidity Time Constant, Arm 1 Compliance, Arm 1 Time Constant, ... Arm N Compliance, Arm N Time Constant])
    normal model: array([Elastic Compliance, Arm 1 Compliance, Arm 1 Time Constant, ... Arm N Compliance, Arm N Time Constant])
    :param time: (n,) numpy array time signal of the loading cycle
    :param force: (n,) numpy array force signal of the loading cycle
    :param radius: float radius of the indenter tip
    :return: (n,) numpy array indentation signal of the loading cycle
    '''
    dirac = zeros(time.shape)
    dirac[0] = 1
    if model_params.size % 2 == 0:  # fluidity case
        time_matrix = row2mat(time, model_params[0::2].size)
        retardance = sum(model_params[2::2] / model_params[3::2] * exp(- time_matrix / model_params[3::2]), axis=1) + model_params[1]
        retardance += (model_params[0] * dirac) / (time[1] - time[0])  # add the delta function to the relaxances
        # must divide by dt since the integral of dirac delta MUST be 1 by definiton
    else:  # no fluidity case
        time_matrix = row2mat(time, model_params[1::2].size)
        retardance = sum(model_params[1::2] / model_params[2::2] * exp(- time_matrix / model_params[2::2]), axis=1)
        retardance += (model_params[0] * dirac) / (time[1] - time[0])  # add the delta function to the relaxances
        # must divide by dt since the integral of dirac delta MUST be 1 by definiton
    return (3 / (8 * sqrt(radius)) * convolve(force, retardance, mode='full')[:time.size] * (time[1] - time[0])) ** (2 / 3)


def LR_PowerLaw_Force(model_params, time, indentation, radius):
    '''
    calculates the force respones for a generalized maxwell model according to the lee and radok contact mechanics formulation
    :param model_params: numpy array contains the model's instantaneous relaxation (E0) and power law exponent (a) in the following form
    array([E0, a])
    :param time: (n,) numpy array time signal of the loading cycle
    :param indentation: (n,) numpy array indentation signal of the loading cycle
    :param radius: float radius of the indenter tip
    :return: (n,) numpy array force signal of the loading cycle
    '''
    scaled_indentations_deriv = concatenate(([0], diff(indentation ** (3 / 2)))) / (time[1] - time[0])
    relaxation = model_params[0] * (1 + time / (time[1] - time[0])) ** (- model_params[1])
    return 16 * sqrt(radius) / 3 * convolve(relaxation, scaled_indentations_deriv, mode='full')[:time.size] * (time[1] - time[0])


class LR_Maxwell():
    def __init__(self, forces, times, indentations, radii, E_logbounds=(1, 9), T_logbounds=(-5, 0)):  #@TODO add conical and flat punch indenter options
        '''
        initializes an instance of the LR_Maxwell class
        used for generating fits, of experimentally obtained force-distance data all belonging to the same sample,
        to a maxwell model which corresponds to the sample's viscoelastic behavior
        :param forces: either list of numpy arrays or single numpy array corresponding to the force signals from an AFM
        :param times: either list of numpy arrays or single numpy array corresponding to the time signals from an AFM
        :param indentations: either list of numpy arrays or single numpy array corresponding to the indentation signals from an AFM
        :param radii: either list of floats or single float corresponding to the tip radii of an AFM
        :param E_logbounds: tuple (float, float) high and low log bound for the elastic elements in the model
        :param T_logbounds: tuple (float, float) high and low log bound for the time constants in the model
        '''
        # if there are multiple inputs
        if type(forces) is list:
            # check for any size mismatches
            if any([len(arr) != len(forces) for arr in (times, indentations, radii)]):
                exit('Error: Size Mismatch in Experimental Observables!  All experimental observables must be the same size!')
            # concatenate the lists of experimental observables to put them into a single row vector form
            self.time = times
            self.indentation = indentations
            # create a 'mask' of dt to properly integrate each experiment
            self.dts = [dt * ones(arr.shape) for dt, arr in zip([t[1] - t[0] for t in times], times)]
            # create a 'mask' of radii to scale each experiment
            self.radii = [radius * ones(arr.shape) for radius, arr in zip(radii, forces)]
            # need to concatenate the desired observable so that it can be easily compared
            self.force = concatenate(forces)
        # if there are single inputs
        else:
            # dt is a single value rather than a 'mask' array as seen above
            self.dts = [times[1] - times[0]]
            self.time = [times]
            self.indentation = [indentations]
            # radius is a single value rather than a 'mask' array as seen above
            self.radii = [radii]
            self.force = forces
        # define the boundaries
        self.E_logbounds = E_logbounds
        self.T_logbounds = T_logbounds

    def LR_force(self, model_params):
        '''
        calculates the response force for a generalized maxwell model according to the lee and radok contact mechanics formulation
        :param model_params: numpy array contains the model stiffnesses and time constants in either of the two following forms:
        to incorporate steady state fluidity: array([Elastic Stiffness, Fluidity Time Constant, Arm 1 Stiffness, Arm 1 Time Constant, ... Arm N Stiffness, Arm N Time Constant])
        normal model: array([Elastic Stiffness, Arm 1 Stiffness, Arm 1 Time Constant, ... Arm N Stiffness, Arm N Time Constant])
        :return: numpy array 'predicted' force signals for all real (experimentally obtained) indentations
        '''
        def make_relaxance(t, dt):
            if model_params.size % 2 == 0:  # fluidity case:
                time_matrix = row2mat(t, model_params[0::2].size)
                relaxance = - sum(model_params[0::2] / model_params[1::2] * exp(- time_matrix / model_params[1::2]), axis=1)
                relaxance += sum(model_params[0::2]) / dt  # add the delta function to the relaxances
                # must divide by dt since the integral of dirac delta MUST be 1 by definiton
            else:  # no fluidity case
                time_matrix = row2mat(t, model_params[1::2].size)
                relaxance = - sum(model_params[1::2] / model_params[2::2] * exp(- time_matrix / model_params[2::2]), axis=1)
                relaxance += (model_params[0] + sum(model_params[1::2])) / dt  # add the delta function of the relaxances
                # must divide by dt since the integral of dirac delta MUST be 1 by definiton
            return relaxance
        return array([16 * sqrt(r) / 3 * convolve(make_relaxance(t, dt), h ** (3 / 2), mode='full')[: t.size] * dt
                      for r, t, dt, h in zip(self.radii, self.time, self.dts, self.indentation)]).ravel()

    def get_bounds(self, model_size, fluidity=False):
        '''
        gets the boundaries for a maxwell model of a given size
        :param model_size: int number of arms in the model
        :param fluidity: bool whether or not to include the fluidity term
        :return: numpy arrays of boundaries for the maxwell model (lower bound, upper bound)
        '''
        lower = 10 ** concatenate(([self.E_logbounds[0]],
                                   concatenate([[self.E_logbounds[0], self.T_logbounds[0]]
                                                for i in range(model_size)]))).astype(float)
        upper = 10 ** concatenate(([self.E_logbounds[1]],
                                   concatenate([[self.E_logbounds[1], self.T_logbounds[1]]
                                                for i in range(model_size)]))).astype(float)
        if fluidity:
            lower = insert(lower, 1, 10 ** self.T_logbounds[0])
            upper = insert(upper, 1, 10 ** self.T_logbounds[1])
        return lower, upper

    def get_initial_guess(self, model_size, fluidity=False):
        '''
        gets random log-uniform initial guess for a maxwell model of a given size
        :param model_size: int number of arms in the model
        :param fluidity: bool whether or not to include the fluidity term
        :return: numpy array of initial guess for a maxwell model
        '''
        guess = 10 ** concatenate(([uniform(low=self.E_logbounds[0], high=self.E_logbounds[1])],
                                   concatenate([[uniform(low=self.E_logbounds[0], high=self.E_logbounds[1]),
                                                 uniform(low=self.T_logbounds[0], high=self.T_logbounds[1])]
                                                for i in range(model_size)])))
        if fluidity:
            guess = insert(guess, 1, 10 ** uniform(low=self.T_logbounds[0], high=self.T_logbounds[1]))
        return guess

    def SSE(self, model_params, lower_bounds, upper_bounds):
        '''
        gives the sum of squared errors between the 'predicted' force and real (experimentally obtained) force signals
        :param model_params: numpy array of relaxance parameters (refer to LR_force)
        :param lower_bounds: numpy array result of a single get_bounds[0] function call (lower bounds)
        :param upper_bounds: numpy array result of a single get_bounds[1] function call (upper bounds)
        :return: float sum of squared errors between the 'predicted' and real force signals
        '''
        sse = sum((self.LR_force(model_params=model_params) - self.force) ** 2, axis=0)
        if any(lower_bounds > model_params) or any(upper_bounds < model_params):
            return 1e20 * sse
        return sse

    def fit(self, maxiter=1000, max_model_size=4, fit_sequential=True, num_attempts=5):
        '''
        fit experimental force distance curve(s) to maxwell model of arbitrary size using a nelder-mead simplex which
        typically gives good fits rather quickly
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :param max_model_size: int largest number of arms per maxwell model to test (going larger tends to give poor and unphysical fits)
        :param fit_sequential: bool whether or not to fit sequentially (cascade fit from previous model as the initial guess of the next) (RECOMMENDED)
        :param num_attempts: int number of fitting attempts to make per fit, larger number will give more statistically significant results, but
        will take longer
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        data = []  # store the global data for the fits
        for model_size in range(1, max_model_size + 1):  # fit without fluidity
            current_data = []  # store the data for the current fitting attempts
            tic()
            lower_bounds, upper_bounds = self.get_bounds(model_size, fluidity=False)
            for fit_attempt in range(num_attempts):
                guess = self.get_initial_guess(model_size, fluidity=False)
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                    guess[: 2 * (model_size - 1) + 1] = data[-1][0][: 2 * (model_size - 1) + 1]
                results = minimize(self.SSE, x0=guess, args=(lower_bounds, upper_bounds),
                                   method='Nelder-Mead', options={'maxiter': maxiter,
                                                                  'maxfev': maxiter,
                                                                  'xatol': 1e-60,
                                                                  'fatol': 1e-60})
                current_data.append([results.x, results.fun])
            current_data = array(current_data, dtype='object')
            best_fit = current_data[:, 0][argmin(current_data[:, 1])]
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True), var(current_data[:, -1])])

        for model_size in range(1, max_model_size + 1):  # fit with fluidity
            current_data = []  # store the data for the current fitting attempts
            tic()
            lower_bounds, upper_bounds = self.get_bounds(model_size, fluidity=True)
            for fit_attempt in range(num_attempts):
                guess = self.get_initial_guess(model_size, fluidity=True)
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                    guess[: 2 * model_size] = data[-1][0][: 2 * model_size]
                results = minimize(self.SSE, x0=guess, args=(lower_bounds, upper_bounds),
                                   method='Nelder-Mead', options={'maxiter': maxiter,
                                                                  'maxfev': maxiter,
                                                                  'xatol': 1e-60,
                                                                  'fatol': 1e-60})
                current_data.append([results.x, results.fun])
            current_data = array(current_data, dtype='object')
            best_fit = current_data[:, 0][argmin(current_data[:, 1])]
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True), var(current_data[:, -1])])

        data = array(data, dtype='object')
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2], 'trial_variance': best_fit[3]}

    def fit_slow(self, maxiter=1000, max_model_size=4, fit_sequential=True, num_attempts=5):
        '''
        fit experimental force distance curve(s) to maxwell model of arbitrary size using simulated annealing with
        a nelder-mead simplex local search, this is very computationally costly and will take a very long time
        though typically results in much better fits
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :param max_model_size: int largest number of arms per maxwell model to test (going larger tends to give poor and unphysical fits)
        :param fit_sequential: bool whether or not to fit sequentially (cascade fit from previous model as the initial guess of the next) (RECOMMENDED)
        :param num_attempts: int number of fitting attempts to make per fit, larger number will give more statistically significant results, but
        will take longer
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        data = []  # store the global data for the fits
        for model_size in range(1, max_model_size + 1):  # fit without fluidity
            current_data = []  # store the data for the current fitting attempts
            tic()
            lower_bounds, upper_bounds = self.get_bounds(model_size, fluidity=False)
            bound = array((lower_bounds, upper_bounds)).T
            for fit_attempt in range(num_attempts):
                guess = self.get_initial_guess(model_size, fluidity=False)
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                    guess[: 2 * (model_size - 1) + 1] = data[-1][0][: 2 * (model_size - 1) + 1]
                results = dual_annealing(self.SSE, bound, args=(lower_bounds, upper_bounds), maxiter=maxiter,
                                         local_search_options={'method': 'nelder-mead'}, x0=guess)
                current_data.append([results.x, results.fun])
            current_data = array(current_data, dtype='object')
            best_fit = current_data[:, 0][argmin(current_data[:, 1])]
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True), var(current_data[:, -1])])

        for model_size in range(1, max_model_size + 1):  # fit with fluidity
            current_data = []  # store the data for the current fitting attempts
            tic()
            lower_bounds, upper_bounds = self.get_bounds(model_size, fluidity=True)
            bound = array((lower_bounds, upper_bounds)).T
            for fit_attempt in range(num_attempts):
                guess = self.get_initial_guess(model_size, fluidity=True)
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                    guess[: 2 * model_size] = data[-1][0][: 2 * model_size]
                results = dual_annealing(self.SSE, bound, args=(lower_bounds, upper_bounds), maxiter=maxiter,
                                         local_search_options={'method': 'nelder-mead'}, x0=guess)
                current_data.append([results.x, results.fun])
            current_data = array(current_data, dtype='object')
            best_fit = current_data[:, 0][argmin(current_data[:, 1])]
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True), var(current_data[:, -1])])

        data = array(data, dtype='object')
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2], 'trial_variance': best_fit[3]}


class LR_Voigt():
    def __init__(self, forces, times, indentations, radii, J_logbounds=(-9, -1), T_logbounds=(-5, 0)):  #@TODO add conical and flat punch indenter options
        '''
        initializes an instance of the LR_Voigt class
        used for generating fits, of experimentally obtained force-distance data all belonging to the same sample,
        to a kelvin-voigt model which corresponds to the sample's viscoelastic behavior
        :param forces: either list of numpy arrays or single numpy array corresponding to the force signals from an AFM
        :param times: either list of numpy arrays or single numpy array corresponding to the time signals from an AFM
        :param indentations: either list of numpy arrays or single numpy array corresponding to the indentation signals from an AFM
        :param radii: either list of floats or single float corresponding to the tip radii of an AFM
        :param J_logbounds: tuple (float, float) high and low log bound for the compliance elements in the model
        :param T_logbounds: tuple (float, float) high and low log bound for the time constants in the model
        '''
        # if there are multiple inputs
        if type(forces) is list:
            # check for any size mismatches
            if any([len(arr) != len(forces) for arr in (times, indentations, radii)]):
                exit('Error: Size Mismatch in Experimental Observables!  All experimental observables must be the same size!')
            # concatenate the lists of experimental observables to put them into a single row vector form
            self.time = times
            self.force = forces
            # create a 'mask' of dt to properly integrate each experiment
            self.dts = [dt * ones(arr.shape) for dt, arr in zip([t[1] - t[0] for t in times], times)]
            # create a 'mask' of radii to scale each experiment
            self.radii = [radius * ones(arr.shape) for radius, arr in zip(radii, forces)]
            # need to concatenate the desired observable so that it can be easily compared
            self.scaled_indentation = concatenate([indentation ** (3 / 2) for indentation in indentations])
        # if there are single inputs
        else:
            # dt is a single value rather than a 'mask' array as seen above
            self.dts = [times[1] - times[0]]
            self.time = [times]
            self.force = [forces]
            self.scaled_indentation = indentations ** (3 / 2)
            # radius is a single value rather than a 'mask' array as seen above
            self.radii = [radii]
        # define the boundaries
        self.J_logbounds = J_logbounds
        self.T_logbounds = T_logbounds

    def LR_scaled_indentation(self, model_params):
        '''
        calculates the scaled response indentation (h^3/2) for a generalized maxwell model according to the lee and radok contact mechanics formulation
        :param model_params: numpy array contains the model compliances and time constants in either of the two following forms:
        to incorporate steady state fluidity: array([Elastic Compliance, Fluidity Time Constant, Arm 1 Compliance, Arm 1 Time Constant, ... Arm N Compliance, Arm N Time Constant])
        normal model: array([Elastic Compliance, Arm 1 Compliance, Arm 1 Time Constant, ... Arm N Compliance, Arm N Time Constant])
        :return: numpy array scaled 'predicted' indentation signals (h^3/2) for all real (experimentally obtained) forces
        '''
        def make_retardance(t, dt):
            if model_params.size % 2 == 0:  # fluidity case:
                time_matrix = row2mat(t, model_params[3::2].size)
                retardance = sum(model_params[2::2] / model_params[3::2] * exp(- time_matrix / model_params[3::2]), axis=1) + model_params[1]
                retardance += model_params[0] / dt  # add the delta function to the relaxances
                # must divide by dt since the integral of dirac delta MUST be 1 by definiton
            else:  # no fluidity case
                time_matrix = row2mat(t, model_params[1::2].size)
                retardance = sum(model_params[1::2] / model_params[2::2] * exp(- time_matrix / model_params[2::2]), axis=1)
                retardance += model_params[0] / dt  # add the delta function of the relaxances
                # must divide by dt since the integral of dirac delta MUST be 1 by definiton
            return retardance
        return array([3 / (sqrt(r) * 8) * convolve(make_retardance(t, dt), f, mode='full')[: t.size] * dt
                      for r, t, dt, f in zip(self.radii, self.time, self.dts, self.force)]).ravel()

    def get_bounds(self, model_size, fluidity=False):
        '''
        gets the boundaries for a maxwell model of a given size
        :param model_size: int number of arms in the model
        :param fluidity: bool whether or not to include the fluidity term
        :return: numpy arrays of boundaries for the kelvin-voigt model (lower bound, upper bound)
        '''
        lower = 10 ** concatenate(([self.J_logbounds[0]],
                                   concatenate([[self.J_logbounds[0], self.T_logbounds[0]]
                                                for i in range(model_size)]))).astype(float)
        upper = 10 ** concatenate(([self.J_logbounds[1]],
                                   concatenate([[self.J_logbounds[1], self.T_logbounds[1]]
                                                for i in range(model_size)]))).astype(float)
        if fluidity:
            lower = insert(lower, 1, 10 ** self.T_logbounds[0])
            upper = insert(upper, 1, 10 ** self.T_logbounds[1])
        return lower, upper

    def get_initial_guess(self, model_size, fluidity=False):
        '''
        gets random log-uniform initial guess for a kelvin-voigt model of a given size
        :param model_size: int number of arms in the model
        :param fluidity: bool whether or not to include the fluidity term
        :return: numpy array of initial guess for a kelvin-voigt model
        '''
        guess = 10 ** concatenate(([uniform(low=self.J_logbounds[0], high=self.J_logbounds[1])],
                                   concatenate([[uniform(low=self.J_logbounds[0], high=self.J_logbounds[1]),
                                                 uniform(low=self.T_logbounds[0], high=self.T_logbounds[1])]
                                                for i in range(model_size)])))
        if fluidity:
            guess = insert(guess, 1, 10 ** uniform(low=self.T_logbounds[0], high=self.T_logbounds[1]))
        return guess

    def SSE(self, model_params, lower_bounds, upper_bounds):
        '''
        gives the sum of squared errors between the scaled 'predicted' indentation and real scaled (experimentally obtained) indentation signals (h^3/2)
        :param model_params: numpy array of retardance parameters (refer to LR_force)
        :return: float sum of squared errors between the scaled 'predicted' and real indentation signals (h^3/2)
        '''
        sse = sum((self.LR_scaled_indentation(model_params=model_params) - self.scaled_indentation) ** 2, axis=0)
        if any(lower_bounds > model_params) or any(upper_bounds < model_params):
            return 1e20 * sse
        return sse

    def fit(self, maxiter=1000, max_model_size=4, fit_sequential=True, num_attempts=5):
        '''
        fit experimental force distance curve(s) to kelvin-voigt model of arbitrary size using a nelder-mead simplex which
        typically gives good fits rather quickly
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :param max_model_size: int largest number of arms per maxwell model to test (going larger tends to give poor and unphysical fits)
        :param fit_sequential: bool whether or not to fit sequentially (cascade fit from previous model as the initial guess of the next) (RECOMMENDED)
        :param num_attempts: int number of fitting attempts to make per fit, larger number will give more statistically significant results, but
        will take longer
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        data = []  # store the global data for the fits
        for model_size in range(1, max_model_size + 1):  # fit without fluidity
            current_data = []  # store the data for the current fitting attempts
            tic()
            lower_bounds, upper_bounds = self.get_bounds(model_size, fluidity=False)
            for fit_attempt in range(num_attempts):
                guess = self.get_initial_guess(model_size, fluidity=False)
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                    guess[: 2 * (model_size - 1) + 1] = data[-1][0][: 2 * (model_size - 1) + 1]
                results = minimize(self.SSE, x0=guess, args=(lower_bounds, upper_bounds),
                                   method='Nelder-Mead', options={'maxiter': maxiter,
                                                                  'maxfev': maxiter,
                                                                  'xatol': 1e-60,
                                                                  'fatol': 1e-60})
                current_data.append([results.x, results.fun])
            current_data = array(current_data, dtype='object')
            best_fit = current_data[:, 0][argmin(current_data[:, 1])]
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True), var(current_data[:, -1])])

        for model_size in range(1, max_model_size + 1):  # fit with fluidity
            current_data = []  # store the data for the current fitting attempts
            tic()
            lower_bounds, upper_bounds = self.get_bounds(model_size, fluidity=True)
            for fit_attempt in range(num_attempts):
                guess = self.get_initial_guess(model_size, fluidity=True)
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                    guess[: 2 * model_size] = data[-1][0][: 2 * model_size]
                results = minimize(self.SSE, x0=guess, args=(lower_bounds, upper_bounds),
                                   method='Nelder-Mead', options={'maxiter': maxiter,
                                                                  'maxfev': maxiter,
                                                                  'xatol': 1e-60,
                                                                  'fatol': 1e-60})
                current_data.append([results.x, results.fun])
            current_data = array(current_data, dtype='object')
            best_fit = current_data[:, 0][argmin(current_data[:, 1])]
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True), var(current_data[:, -1])])

        data = array(data, dtype='object')
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2],
                'trial_variance': best_fit[3]}

    def fit_slow(self, maxiter=1000, max_model_size=4, fit_sequential=True, num_attempts=5):
        '''
        fit experimental force distance curve(s) to kelvin-voigt model of arbitrary size using simulated annealing with
        a nelder-mead simplex local search, this is very computationally costly and will take a very long time
        though typically results in much better fits
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :param max_model_size: int largest number of arms per maxwell model to test (going larger tends to give poor and unphysical fits)
        :param fit_sequential: bool whether or not to fit sequentially (cascade fit from previous model as the initial guess of the next) (RECOMMENDED)
        :param num_attempts: int number of fitting attempts to make per fit, larger number will give more statistically significant results, but
        will take longer
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        data = []  # store the global data for the fits
        for model_size in range(1, max_model_size + 1):  # fit without fluidity
            current_data = []  # store the data for the current fitting attempts
            tic()
            lower_bounds, upper_bounds = self.get_bounds(model_size, fluidity=False)
            bound = array((lower_bounds, upper_bounds)).T
            for fit_attempt in range(num_attempts):
                guess = self.get_initial_guess(model_size, fluidity=False)
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                    guess[: 2 * (model_size - 1) + 1] = data[-1][0][: 2 * (model_size - 1) + 1]
                results = dual_annealing(self.SSE, bound, args=(lower_bounds, upper_bounds), maxiter=maxiter,
                                         local_search_options={'method': 'nelder-mead'}, x0=guess)
                current_data.append([results.x, results.fun])
            current_data = array(current_data, dtype='object')
            best_fit = current_data[:, 0][argmin(current_data[:, 1])]
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True), var(current_data[:, -1])])

        for model_size in range(1, max_model_size + 1):  # fit with fluidity
            current_data = []  # store the data for the current fitting attempts
            tic()
            lower_bounds, upper_bounds = self.get_bounds(model_size, fluidity=True)
            bound = array((lower_bounds, upper_bounds)).T
            for fit_attempt in range(num_attempts):
                guess = self.get_initial_guess(model_size, fluidity=True)
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                    guess[: 2 * model_size] = data[-1][0][: 2 * model_size]
                results = dual_annealing(self.SSE, bound, args=(lower_bounds, upper_bounds), maxiter=maxiter,
                                         local_search_options={'method': 'nelder-mead'}, x0=guess)
                current_data.append([results.x, results.fun])
            current_data = array(current_data, dtype='object')
            best_fit = current_data[:, 0][argmin(current_data[:, 1])]
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True), var(current_data[:, -1])])

        data = array(data, dtype='object')
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2],
                'trial_variance': best_fit[3]}


class LR_PowerLaw():
    def __init__(self, forces, times, indentations, radii, E0_logbounds=(1, 9), a_logbounds=(-5, 0)):  #@TODO add conical and flat punch indenter options
        '''
        initializes an instance of the Custom_Model class
        used for generating fits, of experimentally obtained force-distance data all belonging to the same sample,
        to a power law rheology model which corresponds to the sample's viscoelastic behavior
        :param forces: either list of numpy arrays or single numpy array corresponding to the force signals from an AFM
        :param times: either list of numpy arrays or single numpy array corresponding to the time signals from an AFM
        :param indentations: either list of numpy arrays or single numpy array corresponding to the indentation signals from an AFM
        :param radii: either list of floats or single float corresponding to the tip radii of an AFM
        :param E0_logbounds: tuple (float, float) high and low log bound for the compliance elements in the model
        :param a_logbounds: tuple (float, float) high and low log bound for the time constants in the model
        '''
        # if there are multiple inputs
        if type(forces) is list:
            # check for any size mismatches
            if any([len(arr) != len(forces) for arr in (times, indentations, radii)]):
                exit('Error: Size Mismatch in Experimental Observables!  All experimental observables must be the same size!')
            # concatenate the lists of experimental observables to put them into a single row vector form
            self.time = times
            self.force = concatenate(forces)
            # create a 'mask' of dt to properly integrate each experiment
            self.dts = [dt * ones(arr.shape) for dt, arr in zip([t[1] - t[0] for t in times], times)]
            # create a 'mask' of radii to scale each experiment
            self.radii = [radius * ones(arr.shape) for radius, arr in zip(radii, forces)]
            # calculate the numerical derivatives of the scaled indentations (h^3/2)
            self.scaled_indentations_deriv = [concatenate(([0], diff(indentation ** (3 / 2)))) / dt for indentation, dt in zip(indentations, self.dts)]
        # if there are single inputs
        else:
            self.time = [times]
            # dt is a single value rather than a 'mask' array as seen above
            self.dts = [times[1] - times[0]]
            # radius is a single value rather than a 'mask' array as seen above
            self.radii = [radii]
            self.force = forces
            # calculate the numerical derivatives of the scaled indentations (h^3/2)
            self.scaled_indentations_deriv = [concatenate(([0], diff(indentations ** (3 / 2)))) / self.dts]
        # define the boundaries
        self.E0_logbounds = E0_logbounds
        self.a_logbounds = a_logbounds

    def LR_force(self, model_params):
        '''
        calculates the force respones for a generalized maxwell model according to the lee and radok contact mechanics formulation
        :param model_params: numpy array contains the model's instantaneous relaxation (E0) and power law exponent (a) in the following form
        array([E0, a])
        :return: numpy array scaled 'predicted' force signals for all real (experimentally obtained) forces
        '''
        def get_relaxation(t, dt):
            return model_params[0] * (1 + t / dt) ** (- model_params[1])
        return concatenate([16 * sqrt(r) / 3 * convolve(get_relaxation(t, dt), scaled_dh, mode='full')[: t.size] * dt
                            for r, t, dt, scaled_dh in zip(self.radii, self.time, self.dts, self.scaled_indentations_deriv)])

    def get_bounds(self):
        '''
        gets the boundaries for a power law rheology model of a given size
        :return: numpy arrays of boundaries for the power law rheology model (lower bound, upper bound)
        '''
        lower = 10 ** array([self.E0_logbounds[0], self.a_logbounds[0]]).astype(float)
        upper = 10 ** array([self.E0_logbounds[1], self.a_logbounds[1]]).astype(float)
        return lower, upper

    def get_initial_guess(self):
        '''
        gets random log-uniform initial guess for a power law rheology model of a given size
        :return: numpy array of initial guess for a power law rheology model
        '''
        return 10 ** array([uniform(low=self.E0_logbounds[0], high=self.E0_logbounds[1]),
                            uniform(low=self.a_logbounds[0], high=self.a_logbounds[1])])

    def SSE(self, model_params, lower_bounds, upper_bounds):
        '''
        gives the sum of squared errors between the scaled 'predicted' indentation and real scaled (experimentally obtained) force signals
        :param model_params: numpy array of model parameters (refer to LR_force)
        :param lower_bounds: numpy array result of a single get_bounds[0] function call (lower bounds)
        :param upper_bounds: numpy array result of a single get_bounds[1] function call (upper bounds)
        :return: float sum of squared errors between the scaled 'predicted' and real force signals
        '''
        sse = sum((self.LR_force(model_params=model_params) - self.force) ** 2, axis=0)
        if any(lower_bounds > model_params) or any(upper_bounds < model_params):
            return 1e20 * sse
        return sse

    def fit(self, maxiter=1000, num_attempts=5):
        '''
        fit experimental force distance curve(s) to power law rheology model using a nelder-mead simplex which typically gives good fits rather quickly
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :param num_attempts: int number of fitting attempts to make per fit, larger number will give more statistically significant results, but
        will take longer
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        data = []  # store the global data for the fits
        tic()
        for fit_attempt in range(num_attempts):
            lower_bounds, upper_bounds = self.get_bounds()
            guess = self.get_initial_guess()
            results = minimize(self.SSE, x0=guess, args=(lower_bounds, upper_bounds),
                               method='Nelder-Mead', options={'maxiter': maxiter,
                                                                  'maxfev': maxiter,
                                                                  'xatol': 1e-60,
                                                                  'fatol': 1e-60})
            data.append([results.x, results.fun])
        data = array(data, dtype='object')
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': toc(True), 'trial_variance': var(data[:, -1])}

    def fit_slow(self, maxiter=1000, num_attempts=5):
        '''
        fit experimental force distance curve(s) to power law rheology model using simulated annealing with
        a nelder-mead simplex local search, this is very computationally costly and will take a very long time
        though typically results in much better fits
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :param num_attempts: int number of fitting attempts to make per fit, larger number will give more statistically significant results, but
        will take longer
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        data = []  # store the global data for the fits
        tic()
        for fit_attempt in range(num_attempts):
            lower_bounds, upper_bounds = self.get_bounds()
            bound = array((lower_bounds, upper_bounds)).T
            guess = self.get_initial_guess()
            results = dual_annealing(self.SSE, bound, args=(lower_bounds, upper_bounds), maxiter=maxiter,
                                     local_search_options={'method': 'nelder-mead'}, x0=guess)
            data.append([results.x, results.fun])
        data = array(data, dtype='object')
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': toc(True), 'trial_variance': var(data[:, -1])}


class Custom_Model():
    def __init__(self, forces, times, indentations, radii, target_observable):
        # if there are multiple inputs
        if type(forces) is list:
            # check for any size mismatches
            if any([len(arr) != len(forces) for arr in (times, indentations, radii, target_observable)]):
                exit(
                    'Error: Size Mismatch in Experimental Observables!  All experimental observables must be the same size!')
            # concatenate the lists of experimental observables to put them into a single row vector form
            self.force = forces
            self.time = times
            self.indentation = indentations
            # create a 'mask' of dt to properly integrate each experiment
            self.dts = [dt * ones(arr.shape) for dt, arr in zip([t[1] - t[0] for t in times], times)]
            # create a 'mask' of radii to scale each experiment
            self.radii = [radius * ones(arr.shape) for radius, arr in zip(radii, forces)]
            self.target_observable = concatenate(target_observable)
        # if there are single inputs
        else:
            self.force = [forces]
            self.time = [times]
            self.indentation = [indentations]
            # dt is a single value rather than a 'mask' array as seen above
            self.dts = [times[1] - times[0]]
            # radius is a single value rather than a 'mask' array as seen above
            self.radii = radii
            self.target_observable = target_observable
        self.observable_function = None

    def SSE(self, model_params, upper_bounds, lower_bounds):
        sse = sum((self.observable_function(model_params) - self.target_observable) ** 2, axis=0)
        if any(lower_bounds > model_params) or any(upper_bounds < model_params):
            return 1e20 * sse
        return sse

    def SSE(self, model_params):
        return sum((self.observable_function(model_params) - self.target_observable) ** 2, axis=0)

    def fit(self, bounds, maxiter=1000, num_attempts=5):
        '''
        fit experimental observable of your choice to a custom model for the observable using a nelder-mead simplex which typically gives good fits rather quickly
        :param function: function for the desired observable to be predicted
        :param bounds: (n, 2) numpy array of upper and lower bounds: [[lower1, upper1], ... [lowerN, upperN]]
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :param num_attempts: int number of fitting attempts to make per fit, larger number will give more statistically significant results, but
        will take longer
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        if self.observable_function is None:
            exit('Error: the model\'s observable_function argument is undefined!  define by setting model.observable_function = func')
        data = []  # store the global data for the fits
        tic()
        lower_bounds, upper_bounds = bounds[:, 0], bounds[:, 1]
        for fit_attempt in range(num_attempts):
            guess = [uniform(low=low, high=high) for low, high in zip(lower_bounds, upper_bounds)]
            results = minimize(self.observable_function, x0=guess, args=(lower_bounds, upper_bounds),
                               method='Nelder-Mead', options={'maxiter': maxiter,
                                                              'maxfev': maxiter,
                                                              'xatol': 1e-60,
                                                              'fatol': 1e-60})
            data.append([results.x, results.fun])
        data = array(data, dtype='object')
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': toc(True), 'trial_variance': var(data[:, -1])}

    def fit(self, guess, maxiter=1000, num_attempts=5):
        '''
        fit experimental observable of your choice to a custom model for the observable using a nelder-mead simplex which typically gives good fits rather quickly
        :param function: function for the desired observable to be predicted
        :param bounds: (n, 2) numpy array of upper and lower bounds: [[lower1, upper1], ... [lowerN, upperN]]
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :param num_attempts: int number of fitting attempts to make per fit, larger number will give more statistically significant results, but
        will take longer
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        if self.observable_function is None:
            exit('Error: the model\'s observable_function argument is undefined!  define by setting model.observable_function = func')
        data = []  # store the global data for the fits
        tic()
        for fit_attempt in range(num_attempts):
            results = minimize(self.observable_function, x0=guess,
                               method='Nelder-Mead', options={'maxiter': maxiter,
                                                              'maxfev': maxiter,
                                                              'xatol': 1e-60,
                                                              'fatol': 1e-60})
            data.append([results.x, results.fun])
        data = array(data, dtype='object')
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': toc(True), 'trial_variance': var(data[:, -1])}

    def fit_slow(self, bounds, maxiter=1000, num_attempts=5):
        '''
        fit experimental observable of your choice to a custom model for the observable using simulated annealing with
        a nelder-mead simplex local search, this is very computationally costly and will take a very long time
        though typically results in much better fits
        :param function: function for the desired observable to be predicted
        :param bounds: (n, 2) numpy array of upper and lower bounds: [[lower1, upper1], ... [lowerN, upperN]]
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :param num_attempts: int number of fitting attempts to make per fit, larger number will give more statistically significant results, but
        will take longer
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        if self.observable_function is None:
            exit('Error: the model\'s observable_function argument is undefined!  define by setting model.observable_function = func')
        data = []  # store the global data for the fits
        tic()
        lower_bounds, upper_bounds = bounds[:, 0], bounds[:, 1]
        for fit_attempt in range(num_attempts):
            guess = [uniform(low=low, high=high) for low, high in zip(lower_bounds, upper_bounds)]
            results = dual_annealing(self.observable_function, bounds, args=(lower_bounds, upper_bounds), maxiter=maxiter,
                                     local_search_options={'method': 'nelder-mead'}, x0=guess)
            data.append([results.x, results.fun])
        data = array(data, dtype='object')
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': toc(True), 'trial_variance': var(data[:, -1])}

#@TODO suppress warnings
#@TODO add conical and flat punch indenter options
#@TODO add the ibw reader
