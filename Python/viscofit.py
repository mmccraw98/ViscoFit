from numpy import array, sum, sqrt, convolve, exp, ones, zeros, insert, concatenate, argmin, diff, var, seterr
from numpy.random import uniform
from scipy.optimize import minimize, dual_annealing
from general import tic, toc, row2mat


def forceMaxwell_LeeRadok(model_params, time, indentation, radius):
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
    # ABSOLUTELY MUST START TIME FROM 0
    time = time - time[0]
    if model_params.size % 2 == 0:  # fluidity case:
        time_matrix = row2mat(time, model_params[0::2].size)
        relaxance = - sum(model_params[0::2] / model_params[1::2] * exp(- time_matrix / model_params[1::2]), axis=1)
        relaxance[0] += sum(model_params[0::2]) / (time[1] - time[0])  # add the delta function to the relaxances
        # must divide by dt since the integral of dirac delta MUST be 1 by definiton
    else:  # no fluidity case
        time_matrix = row2mat(time, model_params[1::2].size)
        relaxance = - sum(model_params[1::2] / model_params[2::2] * exp(- time_matrix / model_params[2::2]), axis=1)
        relaxance[0] += (model_params[0] + sum(model_params[1::2])) / (time[1] - time[0])  # add the delta function of the relaxances
        # must divide by dt since the integral of dirac delta MUST be 1 by definiton
    return 16 * sqrt(radius) / 3 * convolve(relaxance, indentation ** (3 / 2), mode='full')[:time.size] * (time[1] - time[0])


def indentationKelvinVoigt_LeeRadok(model_params, time, force, radius):
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
    # ABSOLUTELY MUST START TIME FROM 0
    time = time - time[0]
    if model_params.size % 2 == 0:  # fluidity case
        time_matrix = row2mat(time, model_params[0::2].size)
        retardance = sum(model_params[2::2] / model_params[3::2] * exp(- time_matrix / model_params[3::2]), axis=1) + model_params[1]
        retardance[0] += (model_params[0]) / (time[1] - time[0])  # add the delta function to the relaxances
        # must divide by dt since the integral of dirac delta MUST be 1 by definiton
    else:  # no fluidity case
        time_matrix = row2mat(time, model_params[1::2].size)
        retardance = sum(model_params[1::2] / model_params[2::2] * exp(- time_matrix / model_params[2::2]), axis=1)
        retardance[0] += (model_params[0]) / (time[1] - time[0])  # add the delta function to the relaxances
        # must divide by dt since the integral of dirac delta MUST be 1 by definiton
    return (3 / (8 * sqrt(radius)) * convolve(force, retardance, mode='full')[:time.size] * (time[1] - time[0])) ** (2 / 3)


def forcePowerLaw_LeeRadok(model_params, time, indentation, radius):
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


class maxwellModel():
    def __init__(self, forces, times, indentations, radii, E_logbounds=(1, 9), T_logbounds=(-5, 0)):  #@TODO add conical and flat punch indenter options
        '''
        initializes an instance of the maxwellModel class
        used for generating fits, of experimentally obtained force-distance data all belonging to the same sample,
        to a maxwell model which corresponds to the sample's viscoelastic behavior
        :param forces: either list of numpy arrays or single numpy array corresponding to the force signals from an AFM
        :param times: either list of numpy arrays or single numpy array corresponding to the time signals from an AFM
        :param indentations: either list of numpy arrays or single numpy array corresponding to the indentation signals from an AFM
        :param radii: either list of floats or single float corresponding to the tip radii of an AFM
        :param E_logbounds: tuple (float, float) high and low log bound for the elastic elements in the model
        :param T_logbounds: tuple (float, float) high and low log bound for the time constants in the model
        '''
        seterr(divide='ignore', under='ignore', over='ignore')  # ignore div0 and under/overflow warnings in numpy
        # if there are multiple inputs, they need to be treated differently than single inputs
        if type(forces) is list:
            # check for any size mismatches
            if any([len(arr) != len(forces) for arr in (times, indentations, radii)]):
                exit('Error: Size Mismatch in Experimental Observables!  All experimental observables must be the same size!')
            # time and indentation stay in their original forms
            # ABSOLUTELY MUST START TIME FROM 0
            self.time = [t - t[0] for t in times]
            self.indentation = indentations
            # create a dt valued vector for each experiment's dt value
            self.dts = [t[1] - t[0] for t in times]
            # create a radius valued vector for each experiment's radius value
            self.radii = [radius * ones(arr.shape) for radius, arr in zip(radii, forces)]
            self.force = forces
        # if there are single inputs, they must be formatted to a list for compatibility with the multiple input code
        else:
            # dt is a list with a single value rather than a vector as seen above
            self.dts = [times[1] - times[0]]
            # ABSOLUTELY MUST START TIME FROM 0
            self.time = [times - times[0]]
            self.indentation = [indentations]
            self.force = [forces]
            # radius is a list with a single value rather than a vector as seen above
            self.radii = [radii]
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
        # makes the relaxance time signal for a given experiment
        def make_relaxance(t, dt):
            if model_params.size % 2 == 0:  # fluidity case:
                time_matrix = row2mat(t, model_params[0::2].size)  # needed for multiplication with model parameters
                # Q = Eg * delta( t ) - sum_i_N+1( Ei / Ti * exp( -t / Ti ) ) <- sum over all arms
                relaxance = - sum(model_params[0::2] / model_params[1::2] * exp(- time_matrix / model_params[1::2]), axis=1)
                relaxance[0] += sum(model_params[0::2]) / dt  # add the delta function of magnitude Eg to the relaxances
                # must divide by dt since the integral of dirac delta MUST be 1 by definiton
            else:  # no fluidity case
                time_matrix = row2mat(t, model_params[1::2].size)  # needed for multiplication with model parameters
                # Q = Eg * delta( t ) - sum_i_N( Ei / Ti * exp( -t / Ti ) ) <- sum over all but the elastic ('restoring') arm
                relaxance = - sum(model_params[1::2] / model_params[2::2] * exp(- time_matrix / model_params[2::2]), axis=1)
                relaxance[0] += (model_params[0] + sum(model_params[1::2])) / dt  # add the delta function of magnitude Eg to the relaxances
                # must divide by dt since the integral of dirac delta MUST be 1 by definiton
            return relaxance
        # lee and radok viscoelastic contact force
        # F = 16 * sqrt(R) / 3 * integral_0_t( Q( t - u ) * d( u )^3/2 ) du
        # apply for each set of experimental data and turn to a single vector for comparison with the vectorized force
        return [16 * sqrt(r) / 3 * convolve(make_relaxance(t, dt), h ** (3 / 2), mode='full')[: t.size] * dt
                for r, t, dt, h in zip(self.radii, self.time, self.dts, self.indentation)]


    def get_bounds(self, model_size, fluidity=False):
        '''
        gets the boundaries for a maxwell model of a given size
        :param model_size: int number of arms in the model
        :param fluidity: bool whether or not to include the fluidity term
        :return: numpy arrays of boundaries for the maxwell model (lower bound, upper bound)
        '''
        # form the lower bounds by exponentiating the log lower bounds
        lower = 10 ** concatenate(([self.E_logbounds[0]],
                                   concatenate([[self.E_logbounds[0], self.T_logbounds[0]]
                                                for i in range(model_size)]))).astype(float)
        # form the upper bounds by exponentiating the log upper bounds
        upper = 10 ** concatenate(([self.E_logbounds[1]],
                                   concatenate([[self.E_logbounds[1], self.T_logbounds[1]]
                                                for i in range(model_size)]))).astype(float)
        # if there is fluidity in the model, insert a time constant bound after the first element
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
        # get a random guess within the log bounds for each parameter and exponentiate the result
        guess = 10 ** concatenate(([uniform(low=self.E_logbounds[0], high=self.E_logbounds[1])],
                                   concatenate([[uniform(low=self.E_logbounds[0], high=self.E_logbounds[1]),
                                                 uniform(low=self.T_logbounds[0], high=self.T_logbounds[1])]
                                                for i in range(model_size)])))
        # if there is fluidity, insert another exponentiated time constant guess after the first element
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
        # calculate the sum of squared errors between the predicted force vector and the real force vector
        # sse = sum_ti_tf( ( F_pred( ti ) - F_real( ti ) )^2 )
        sse = sum([sum((pred - real) ** 2) for pred, real in zip(self.LR_force(model_params=model_params), self.force)])
        # if any parameters are outside of their associated boundaries, penalize the sse with an arbitrarily large error term
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
        # attempt a trial of fits for each hypothetical model in the desired range of model sizes
        for model_size in range(1, max_model_size + 1):  # fit without fluidity
            print('{}%'.format(100 * model_size / max_model_size))#, end='\r')
            current_data = []  # store the data for the current fitting attempts
            tic()  # start the timer for the trial
            lower_bounds, upper_bounds = self.get_bounds(model_size, fluidity=False)  # get lower and upper bounds
            for fit_attempt in range(num_attempts):  # attempt a number of attempts for each trial
                guess = self.get_initial_guess(model_size, fluidity=False)  # get an initial guess
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess if sequential fitting is desired
                    guess[: 2 * (model_size - 1) + 1] = data[-1][0][: 2 * (model_size - 1) + 1]
                # minimize the SSE within the given bounds using nelder-mead running for a specified number of iterations
                # with the target cost and guess-sensitivity of 1e-60, using the random initial guess
                results = minimize(self.SSE, x0=guess, args=(lower_bounds, upper_bounds),
                                   method='Nelder-Mead', options={'maxiter': maxiter,
                                                                  'maxfev': maxiter,
                                                                  'xatol': 1e-60,
                                                                  'fatol': 1e-60})
                current_data.append([results.x, results.fun])  # save the final parameters and cost to the trial data
            current_data = array(current_data, dtype='object')  # convert the trial data to a numpy array for easier analysis
            best_fit = current_data[:, 0][argmin(current_data[:, 1])]  # get the trial data with the lowest cost -> best fit
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True), var(current_data[:, -1])])  # add the best fit trial data to the global fitting attempts data

        # attempt a trial of fits for model sizes within the specified range using steady state fluidity
        for model_size in range(1, max_model_size + 1):  # fit with fluidity
            current_data = []  # store the data for the current fitting attempts
            tic()  # start the timer for the fit
            lower_bounds, upper_bounds = self.get_bounds(model_size, fluidity=True)  # get the lower and upper bounds for the model size
            for fit_attempt in range(num_attempts):  # attempt a specified number of fits
                guess = self.get_initial_guess(model_size, fluidity=True)  # get a random initial guess within the specified bounds
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess if sequentially fitting
                    guess[: 2 * model_size] = data[-1][0][: 2 * model_size]
                # minimize the SSE within the specified bounds using nelder-mead and the random initial guess
                results = minimize(self.SSE, x0=guess, args=(lower_bounds, upper_bounds),
                                   method='Nelder-Mead', options={'maxiter': maxiter,
                                                                  'maxfev': maxiter,
                                                                  'xatol': 1e-60,
                                                                  'fatol': 1e-60})
                current_data.append([results.x, results.fun])  # append the parameters and their cost
            current_data = array(current_data, dtype='object')  # turn the trial data into a numpy array for easier analysis
            best_fit = current_data[:, 0][argmin(current_data[:, 1])]  # get the trial with the lowest cost
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True), var(current_data[:, -1])])  # save the trial with the lowest cost's data

        data = array(data, dtype='object')  # convert the global fit data into a numpy array for easier analysis
        best_fit = data[argmin(data[:, 1])]  # get the fit with the lowest cost across all model variations
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2], 'trial_variance': best_fit[3]}  # return the best fit parameters

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
            print('{}%'.format(100 * model_size / max_model_size))#, end='\r')
            current_data = []  # store the data for the current fitting attempts
            tic()  # start the timer
            lower_bounds, upper_bounds = self.get_bounds(model_size, fluidity=False)  # get the lower and upper bounds for the given model size
            bound = array((lower_bounds, upper_bounds)).T  # convert the bounds to the desired format for the scipy model
            for fit_attempt in range(num_attempts):  # attempt the desired number of fit attempts
                guess = self.get_initial_guess(model_size, fluidity=False)  # get an initial guess within the specified bounds
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess if using sequential fitting
                    guess[: 2 * (model_size - 1) + 1] = data[-1][0][: 2 * (model_size - 1) + 1]
                # minimize the SSE within the bounds using dual annealing with nelder-mead for local search
                results = dual_annealing(self.SSE, bound, args=(lower_bounds, upper_bounds), maxiter=maxiter,
                                         local_search_options={'method': 'nelder-mead'}, x0=guess)
                current_data.append([results.x, results.fun])  # save the parameters and the cost
            current_data = array(current_data, dtype='object')  # format the trial data to a numpy array for easier analysis
            best_fit = current_data[:, 0][argmin(current_data[:, 1])]  # get the best fit data
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True), var(current_data[:, -1])])  # save the best fit from the trials

        for model_size in range(1, max_model_size + 1):  # fit with fluidity
            current_data = []  # store the data for the current fitting attempts
            tic()  # start the timer
            lower_bounds, upper_bounds = self.get_bounds(model_size, fluidity=True)  # get the lower and upper bounds for the given model parameters
            bound = array((lower_bounds, upper_bounds)).T  # format the bounds for the fitting function
            for fit_attempt in range(num_attempts):  # perform a specified number of fitting attempts
                guess = self.get_initial_guess(model_size, fluidity=True)  # give an initial guess within the specified parameter bounds
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess if sequentially fitting
                    guess[: 2 * model_size] = data[-1][0][: 2 * model_size]
                # minimize the SSE within the bounds using dual annealing with a nelder mead local search
                results = dual_annealing(self.SSE, bound, args=(lower_bounds, upper_bounds), maxiter=maxiter,
                                         local_search_options={'method': 'nelder-mead'}, x0=guess)
                current_data.append([results.x, results.fun])  # save the data for each fit trial
            current_data = array(current_data, dtype='object')  # convert the trial data to a numpy array for easier analysis
            best_fit = current_data[:, 0][argmin(current_data[:, 1])]  # get the trial data with the lowest cost
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True), var(current_data[:, -1])])  # store the data from the best fit trial

        data = array(data, dtype='object')  # convert the overall data to a numpy array for easier analysis
        best_fit = data[argmin(data[:, 1])]  # get the best overall fit
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2], 'trial_variance': best_fit[3]}  # return the best overall fit data


class kelvinVoigtModel():
    def __init__(self, forces, times, indentations, radii, J_logbounds=(-9, -1), T_logbounds=(-5, 0)):  #@TODO add conical and flat punch indenter options
        '''
        initializes an instance of the kelvinVoigtModel class
        used for generating fits, of experimentally obtained force-distance data all belonging to the same sample,
        to a kelvin-voigt model which corresponds to the sample's viscoelastic behavior
        :param forces: either list of numpy arrays or single numpy array corresponding to the force signals from an AFM
        :param times: either list of numpy arrays or single numpy array corresponding to the time signals from an AFM
        :param indentations: either list of numpy arrays or single numpy array corresponding to the indentation signals from an AFM
        :param radii: either list of floats or single float corresponding to the tip radii of an AFM
        :param J_logbounds: tuple (float, float) high and low log bound for the compliance elements in the model
        :param T_logbounds: tuple (float, float) high and low log bound for the time constants in the model
        '''
        seterr(divide='ignore', under='ignore', over='ignore')  # ignore div0 and under/overflow warnings in numpy
        # if there are multiple inputs
        if type(forces) is list:
            # check for any size mismatches
            if any([len(arr) != len(forces) for arr in (times, indentations, radii)]):
                exit('Error: Size Mismatch in Experimental Observables!  All experimental observables must be the same size!')
            # store the time and force data as arrays (as they were given)
            # ABSOLUTELY MUST START TIME FROM 0
            self.time = [t - t[0] for t in times]
            self.force = forces
            # create a dt valued vector for each experiment's dt
            self.dts = [t[1] - t[0] for t in times]
            # create a radius valued vector for each experiment's radius
            self.radii = [radius * ones(arr.shape) for radius, arr in zip(radii, forces)]
            self.scaled_indentation = [indentation ** (3 / 2) for indentation in indentations]
        # if there are single inputs
        else:
            # dt is a single value rather than a 'mask' array as seen above
            self.dts = [times[1] - times[0]]
            # ABSOLUTELY MUST START TIME FROM 0
            self.time = [times - times[0]]
            self.force = [forces]
            self.scaled_indentation = [indentations ** (3 / 2)]
            # radius is a single value rather than a 'mask' array as seen above
            self.radii = [radii]
        # define the boundaries
        self.J_logbounds = J_logbounds
        self.T_logbounds = T_logbounds

    def LR_scaled_indentation(self, model_params):
        '''
        calculates the scaled response indentation (d^3/2) for a generalized maxwell model according to the lee and radok contact mechanics formulation
        :param model_params: numpy array contains the model compliances and time constants in either of the two following forms:
        to incorporate steady state fluidity: array([Elastic Compliance, Fluidity Time Constant, Arm 1 Compliance, Arm 1 Time Constant, ... Arm N Compliance, Arm N Time Constant])
        normal model: array([Elastic Compliance, Arm 1 Compliance, Arm 1 Time Constant, ... Arm N Compliance, Arm N Time Constant])
        :return: numpy array scaled 'predicted' indentation signals (d^3/2) for all real (experimentally obtained) forces
        '''
        def make_retardance(t, dt):  # function to make the retardance signal for an experiment given its parameters
            if model_params.size % 2 == 0:  # fluidity case:
                time_matrix = row2mat(t, model_params[3::2].size)  # needed for vectorized math with model parameter vectors
                # calculate the retardance
                # U = Jg + sum_i_N+1( Ji / Ti * exp( -t / Ti ) ) <- sum over all arms
                retardance = sum(model_params[2::2] / model_params[3::2] * exp(- time_matrix / model_params[3::2]), axis=1) + model_params[1]
                retardance[0] += model_params[0] / dt  # add the delta function to the relaxances
                # must divide by dt since the integral of dirac delta MUST be 1 by definiton
            else:  # no fluidity case
                time_matrix = row2mat(t, model_params[1::2].size)  # needed for vectorized math with model parameter vectors
                # calculate the retardance
                # U = Jg + sum_i_N( Ji / Ti * exp( -t / Ti ) ) <- sum over all but the elastic arm
                retardance = sum(model_params[1::2] / model_params[2::2] * exp(- time_matrix / model_params[2::2]), axis=1)
                retardance[0] += model_params[0] / dt  # add the delta function of the relaxances
                # must divide by dt since the integral of dirac delta MUST be 1 by definiton
            return retardance
        # calculate the prediction for d^3/2 according to the lee and radok viscoelastic contact mechanics for each experiment
        # and convert to a single row vector for easier comparison with the specified target d^3/2
        # d( t )^3/2 = 3 / ( 8 * sqrt( R ) ) * integral_0_t( U( t - u ) * F( u ) ) du
        return [3 / (sqrt(r) * 8) * convolve(make_retardance(t, dt), f, mode='full')[: t.size] * dt
                for r, t, dt, f in zip(self.radii, self.time, self.dts, self.force)]

    def get_bounds(self, model_size, fluidity=False):
        '''
        gets the boundaries for a maxwell model of a given size
        :param model_size: int number of arms in the model
        :param fluidity: bool whether or not to include the fluidity term
        :return: numpy arrays of boundaries for the kelvin-voigt model (lower bound, upper bound)
        '''
        # form the lower bounds by exponentiating the log lower bounds
        lower = 10 ** concatenate(([self.J_logbounds[0]],
                                   concatenate([[self.J_logbounds[0], self.T_logbounds[0]]
                                                for i in range(model_size)]))).astype(float)
        # form the upper bounds by exponentiating the log upper bounds
        upper = 10 ** concatenate(([self.J_logbounds[1]],
                                   concatenate([[self.J_logbounds[1], self.T_logbounds[1]]
                                                for i in range(model_size)]))).astype(float)
        # if there is fluidity in the model, insert a time constant bound after the first element
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
        # get a random guess within the log bounds for each parameter and exponentiate the result
        guess = 10 ** concatenate(([uniform(low=self.J_logbounds[0], high=self.J_logbounds[1])],
                                   concatenate([[uniform(low=self.J_logbounds[0], high=self.J_logbounds[1]),
                                                 uniform(low=self.T_logbounds[0], high=self.T_logbounds[1])]
                                                for i in range(model_size)])))
        # if there is fluidity, insert another exponentiated time constant guess after the first element
        if fluidity:
            guess = insert(guess, 1, 10 ** uniform(low=self.T_logbounds[0], high=self.T_logbounds[1]))
        return guess

    def SSE(self, model_params, lower_bounds, upper_bounds):
        '''
        gives the sum of squared errors between the scaled 'predicted' indentation and real scaled (experimentally obtained) indentation signals (d^3/2)
        :param model_params: numpy array of retardance parameters (refer to LR_force)
        :param lower_bounds: numpy array result of a single get_bounds[0] function call (lower bounds)
        :param upper_bounds: numpy array result of a single get_bounds[1] function call (upper bounds)
        :return: float sum of squared errors between the scaled 'predicted' and real indentation signals (d^3/2)
        '''
        # calculate the sum of squared errors between the predicted force vector and the real force vector
        # sse = sum_ti_tf( ( d^3/2_pred( ti ) - d^3/2_real( ti ) )^2 )
        sse = sum([sum((pred - real) ** 2) for pred, real in zip(self.LR_scaled_indentation(model_params=model_params), self.scaled_indentation)])
        # if any parameters are outside of their associated boundaries, penalize the sse with an arbitrarily large error term
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
        # attempt a trial of fits for each hypothetical model in the desired range of model sizes
        for model_size in range(1, max_model_size + 1):  # fit without fluidity
            print('{}%'.format(100 * model_size / max_model_size))#, end='\r')
            current_data = []  # store the data for the current fitting attempts
            tic()  # start the timer for the trial
            lower_bounds, upper_bounds = self.get_bounds(model_size, fluidity=False)  # get lower and upper bounds
            for fit_attempt in range(num_attempts):  # attempt a number of attempts for each trial
                guess = self.get_initial_guess(model_size, fluidity=False)  # get an initial guess
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess if sequential fitting is desired
                    guess[: 2 * (model_size - 1) + 1] = data[-1][0][: 2 * (model_size - 1) + 1]
                # minimize the SSE within the given bounds using nelder-mead running for a specified number of iterations
                # with the target cost and guess-sensitivity of 1e-60, using the random initial guess
                results = minimize(self.SSE, x0=guess, args=(lower_bounds, upper_bounds),
                                   method='Nelder-Mead', options={'maxiter': maxiter,
                                                                  'maxfev': maxiter,
                                                                  'xatol': 1e-60,
                                                                  'fatol': 1e-60})
                current_data.append([results.x, results.fun])  # save the final parameters and cost to the trial data
            current_data = array(current_data,
                                 dtype='object')  # convert the trial data to a numpy array for easier analysis
            best_fit = current_data[:, 0][
                argmin(current_data[:, 1])]  # get the trial data with the lowest cost -> best fit
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True),
                         var(current_data[:, -1])])  # add the best fit trial data to the global fitting attempts data

        # attempt a trial of fits for model sizes within the specified range using steady state fluidity
        for model_size in range(1, max_model_size + 1):  # fit with fluidity
            current_data = []  # store the data for the current fitting attempts
            tic()  # start the timer for the fit
            lower_bounds, upper_bounds = self.get_bounds(model_size,
                                                         fluidity=True)  # get the lower and upper bounds for the model size
            for fit_attempt in range(num_attempts):  # attempt a specified number of fits
                guess = self.get_initial_guess(model_size,
                                               fluidity=True)  # get a random initial guess within the specified bounds
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess if sequentially fitting
                    guess[: 2 * model_size] = data[-1][0][: 2 * model_size]
                # minimize the SSE within the specified bounds using nelder-mead and the random initial guess
                results = minimize(self.SSE, x0=guess, args=(lower_bounds, upper_bounds),
                                   method='Nelder-Mead', options={'maxiter': maxiter,
                                                                  'maxfev': maxiter,
                                                                  'xatol': 1e-60,
                                                                  'fatol': 1e-60})
                current_data.append([results.x, results.fun])  # append the parameters and their cost
            current_data = array(current_data,
                                 dtype='object')  # turn the trial data into a numpy array for easier analysis
            best_fit = current_data[:, 0][argmin(current_data[:, 1])]  # get the trial with the lowest cost
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True),
                         var(current_data[:, -1])])  # save the trial with the lowest cost's data

        data = array(data, dtype='object')  # convert the global fit data into a numpy array for easier analysis
        best_fit = data[argmin(data[:, 1])]  # get the fit with the lowest cost across all model variations
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2],
                'trial_variance': best_fit[3]}  # return the best fit parameters

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
            print('{}%'.format(100 * model_size / max_model_size))#, end='\r')
            current_data = []  # store the data for the current fitting attempts
            tic()  # start the timer
            lower_bounds, upper_bounds = self.get_bounds(model_size,
                                                         fluidity=False)  # get the lower and upper bounds for the given model size
            bound = array(
                (lower_bounds, upper_bounds)).T  # convert the bounds to the desired format for the scipy model
            for fit_attempt in range(num_attempts):  # attempt the desired number of fit attempts
                guess = self.get_initial_guess(model_size,
                                               fluidity=False)  # get an initial guess within the specified bounds
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess if using sequential fitting
                    guess[: 2 * (model_size - 1) + 1] = data[-1][0][: 2 * (model_size - 1) + 1]
                # minimize the SSE within the bounds using dual annealing with nelder-mead for local search
                results = dual_annealing(self.SSE, bound, args=(lower_bounds, upper_bounds), maxiter=maxiter,
                                         local_search_options={'method': 'nelder-mead'}, x0=guess)
                current_data.append([results.x, results.fun])  # save the parameters and the cost
            current_data = array(current_data,
                                 dtype='object')  # format the trial data to a numpy array for easier analysis
            best_fit = current_data[:, 0][argmin(current_data[:, 1])]  # get the best fit data
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True),
                         var(current_data[:, -1])])  # save the best fit from the trials

        for model_size in range(1, max_model_size + 1):  # fit with fluidity
            current_data = []  # store the data for the current fitting attempts
            tic()  # start the timer
            lower_bounds, upper_bounds = self.get_bounds(model_size,
                                                         fluidity=True)  # get the lower and upper bounds for the given model parameters
            bound = array((lower_bounds, upper_bounds)).T  # format the bounds for the fitting function
            for fit_attempt in range(num_attempts):  # perform a specified number of fitting attempts
                guess = self.get_initial_guess(model_size,
                                               fluidity=True)  # give an initial guess within the specified parameter bounds
                if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess if sequentially fitting
                    guess[: 2 * model_size] = data[-1][0][: 2 * model_size]
                # minimize the SSE within the bounds using dual annealing with a nelder mead local search
                results = dual_annealing(self.SSE, bound, args=(lower_bounds, upper_bounds), maxiter=maxiter,
                                         local_search_options={'method': 'nelder-mead'}, x0=guess)
                current_data.append([results.x, results.fun])  # save the data for each fit trial
            current_data = array(current_data,
                                 dtype='object')  # convert the trial data to a numpy array for easier analysis
            best_fit = current_data[:, 0][argmin(current_data[:, 1])]  # get the trial data with the lowest cost
            data.append([best_fit, self.SSE(best_fit, lower_bounds, upper_bounds), toc(True),
                         var(current_data[:, -1])])  # store the data from the best fit trial

        data = array(data, dtype='object')  # convert the overall data to a numpy array for easier analysis
        best_fit = data[argmin(data[:, 1])]  # get the best overall fit
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2],
                'trial_variance': best_fit[3]}  # return the best overall fit data


class powerLawModel():
    def __init__(self, forces, times, indentations, radii, E0_logbounds=(1, 9), a_logbounds=(-5, 0)):  #@TODO add conical and flat punch indenter options
        '''
        initializes an instance of the customModel class
        used for generating fits, of experimentally obtained force-distance data all belonging to the same sample,
        to a power law rheology model which corresponds to the sample's viscoelastic behavior
        :param forces: either list of numpy arrays or single numpy array corresponding to the force signals from an AFM
        :param times: either list of numpy arrays or single numpy array corresponding to the time signals from an AFM
        :param indentations: either list of numpy arrays or single numpy array corresponding to the indentation signals from an AFM
        :param radii: either list of floats or single float corresponding to the tip radii of an AFM
        :param E0_logbounds: tuple (float, float) high and low log bound for the compliance elements in the model
        :param a_logbounds: tuple (float, float) high and low log bound for the time constants in the model
        '''
        seterr(divide='ignore', under='ignore', over='ignore')  # ignore div0 and under/overflow warnings in numpy
        # if there are multiple inputs
        if type(forces) is list:
            # check for any size mismatches
            if any([len(arr) != len(forces) for arr in (times, indentations, radii)]):
                exit('Error: Size Mismatch in Experimental Observables!  All experimental observables must be the same size!')
            # concatenate the lists of experimental observables to put them into a single row vector form
            # ABSOLUTELY MUST START TIME FROM 0
            self.time = [t - t[0] for t in times]
            self.force = forces
            # create a list of dt and radius valued vectors for each experiment
            self.dts = [(t[1] - t[0]) * ones(t.shape) for t in times]
            self.radii = [radius * ones(arr.shape) for radius, arr in zip(radii, forces)]
            # calculate the numerical derivatives of the scaled indentations d/dt (d(t)^3/2)
            self.scaled_indentations_deriv = [concatenate(([0], diff(indentation ** (3 / 2)))) / dt for indentation, dt in zip(indentations, self.dts)]
        # if there are single inputs put them each into arrays except for the scaled indentation derivatives
        else:
            # ABSOLUTELY MUST START TIME FROM 0
            self.time = [times - times[0]]
            # dt is a single value rather than a 'mask' array as seen above
            self.dts = [times[1] - times[0]]
            # radius is a single value rather than a 'mask' array as seen above
            self.radii = [radii]
            self.force = [forces]
            # calculate the numerical derivatives of the scaled indentations (d^3/2)
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
        def get_relaxation(t, dt):  # create the relaxation signal for a given model
            # E( t ) = E0 * ( 1 + t / dt )^-a
            return model_params[0] * (1 + t / dt) ** (- model_params[1])
        # calculate the 'predicted' force for each experiment and combine the total signal to single row vector
        # for easier comparison with the experimental force
        # F( t ) = 16 * sqrt( R ) / 3 * integral_0_t( E( t - u ) * d/dt d( u )^3/2 ) du
        return [16 * sqrt(r) / 3 * convolve(get_relaxation(t, dt), scaled_dh, mode='full')[: t.size] * dt
                for r, t, dt, scaled_dh in zip(self.radii, self.time, self.dts, self.scaled_indentations_deriv)]

    def get_bounds(self):
        '''
        gets the boundaries for a power law rheology model of a given size
        :return: numpy arrays of boundaries for the power law rheology model (lower bound, upper bound)
        '''
        # calculate the lower bounds by exponentiating a lower end log bound for E0 and a
        lower = 10 ** array([self.E0_logbounds[0], self.a_logbounds[0]]).astype(float)
        # calculate the upper bounds by exponentiating a upper end log bound for E0 and a
        upper = 10 ** array([self.E0_logbounds[1], self.a_logbounds[1]]).astype(float)
        return lower, upper

    def get_initial_guess(self):
        '''
        gets random log-uniform initial guess for a power law rheology model of a given size
        :return: numpy array of initial guess for a power law rheology model
        '''
        # calculate an initial guess by exponentiating a log bound guess within the specified ranges
        # for E0 and a
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
        # calculate the sum of squared errors between the 'predicted' force and the real force
        # sse = sum_ti_tf( ( F_predicted( ti ) - F_real( ti ) )^2 )
        sse = sum([sum((pred - real) ** 2) for pred, real in zip(self.LR_force(model_params=model_params), self.force)])
        # if the current model parameters are outside of the specified boundaries, return an arbitrarily largely scaled sse
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
        tic()  # start the timer
        for fit_attempt in range(num_attempts):  # perform a specified number of fit attempts
            print('{}%'.format(100 * fit_attempt / num_attempts))#, end='\r')
            lower_bounds, upper_bounds = self.get_bounds()  # get the lower and upper bounds for the fit
            guess = self.get_initial_guess()  # get an initial guess within the specified bounds
            # minimize the SSE within the specified bounds using nelder mead
            results = minimize(self.SSE, x0=guess, args=(lower_bounds, upper_bounds),
                               method='Nelder-Mead', options={'maxiter': maxiter,
                                                                  'maxfev': maxiter,
                                                                  'xatol': 1e-60,
                                                                  'fatol': 1e-60})
            data.append([results.x, results.fun])  # save the final parameters and final cost
        data = array(data, dtype='object')  # convert the global data to a numpy array for easier analysis
        best_fit = data[argmin(data[:, 1])]  # get the trial with the lowest cost
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': toc(True), 'trial_variance': var(data[:, -1])}  # return the data from the best fit

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
        tic()  # start the timer
        for fit_attempt in range(num_attempts):  # perform a specified number of fits
            print('{}%'.format(100 * fit_attempt / num_attempts))#, end='\r')
            lower_bounds, upper_bounds = self.get_bounds()  # get the lower and upper bounds
            bound = array((lower_bounds, upper_bounds)).T  # format the bounds for the optimization code
            guess = self.get_initial_guess()  # get an initial guess within the current bounds
            # minimize the SSE within the specified bounds using dual annealing with nelder mead local search
            results = dual_annealing(self.SSE, bound, args=(lower_bounds, upper_bounds), maxiter=maxiter,
                                     local_search_options={'method': 'nelder-mead'}, x0=guess)
            data.append([results.x, results.fun])  # save the final parameters and the cost
        data = array(data, dtype='object')  # format the trial data as a numpy array for convenient analysis
        best_fit = data[argmin(data[:, 1])]  # get the trial with the lowest cost
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': toc(True), 'trial_variance': var(data[:, -1])}  # return the best fit data


class customModel():
    def __init__(self, forces, times, indentations, radii):
        '''
        initializes an instance of the customModel class used for generating fits, of experimentally obtained force-distance data all belonging to the same sample, to a custom defined rheology model which corresponds to the sample's viscoelastic behavior
        :param forces: either list of numpy arrays or single numpy array corresponding to the force signals from an AFM
        :param times: either list of numpy arrays or single numpy array corresponding to the time signals from an AFM
        :param indentations: either list of numpy arrays or single numpy array corresponding to the indentation signals from an AFM
        :param radii: either list of floats or single float corresponding to the tip radii of an AFM
        '''
        seterr(divide='ignore', under='ignore', over='ignore')  # ignore div0 and under/overflow warnings in numpy
        # if there are multiple inputs, then they must all be single row vectors
        if type(forces) is list:
            if any([len(arr) != len(forces) for arr in (times, indentations, radii)]):  # check to make sure that the experimental data is all the same size
                exit('Error: Size Mismatch in Experimental Observables!  All experimental observables must be the same size!')
            self.force = concatenate(forces)  # put the forces into a single row vector
            # ABSOLUTELY MUST START TIME FROM 0
            times = [t - t[0] for t in times]
            self.time = concatenate(times)  # do the same for time
            self.indentation = concatenate(indentations)  # and indentation
            # make row vectors containing the radius and dt values for each experiment
            self.radii = concatenate([radius * ones(arr.shape) for radius, arr in zip(radii, forces)])
            self.dts = concatenate([(t[1] - t[0]) * ones(t.shape) for t in times])
        else:  # if only a single experimental trial is given, store all the values as vectors
            self.force = forces
            # ABSOLUTELY MUST START TIME FROM 0
            self.time = times - times[0]
            self.indentation = indentations
            self.radii = radii * ones(forces.shape)
            self.dts = (times[1] - times[0]) * ones(forces.shape)
        self.observable_function = None  # the custom function being fit - initially is undefined until the fit takes place
        self.target_observable = None  # the 'training data' for the function being fit - will also be defined once the fit takes place

    def SSE(self, model_params, lower_bounds, upper_bounds):
        '''
        gives the sum of squared errors between the scaled 'predicted' indentation and real scaled (experimentally obtained) indentation signals (d^3/2)
        :param model_params: numpy array of retardance parameters (refer to LR_force)
        :param lower_bounds: numpy array result of a single get_bounds[0] function call (lower bounds)
        :param upper_bounds: numpy array result of a single get_bounds[1] function call (upper bounds)
        :return: float sum of squared errors between the scaled 'predicted' and real indentation signals (d^3/2)
        '''
        # calculate the sum of squared errors between the current prediction of the target observable and the real target observable
        # sse = sum_ti_tf( ( Custom_Function_Prediction( ti ) - Real( ti ) )^2 )
        sse = sum((self.observable_function(model_params) - self.target_observable) ** 2, axis=0)
        # if any of the parameters are outside of the specified bounds, return an sse that is scaled to an arbitrarily high degree
        if any(lower_bounds > model_params) or any(upper_bounds < model_params):
            return 1e20 * sse
        return sse

    def fit(self, function, training_data, bounds, maxiter=1000, num_attempts=5):
        '''
        fit experimental observable of your choice to a custom model for the observable using a nelder-mead simplex which typically gives good fits rather quickly
        :param function: function for the desired observable to be predicted
        :param training_data: either numpy array or list of numpy arrays the experimental data to be replicated by the function being trained
        :param bounds: (n, 2) numpy array of upper and lower bounds: [[lower1, upper1], ... [lowerN, upperN]]
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :param num_attempts: int number of fitting attempts to make per fit, larger number will give more statistically significant results, but
        will take longer
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        self.observable_function = function  # define the observable function for optimization
        if type(training_data) is list:  # if there are multiple entries for the training data, put them into row vector form
            training_data = concatenate(training_data)
        self.target_observable = training_data  # define the observable as the target for the objective function to replicate
        if type(bounds) is list:  # convert the function bounds to a numpy array if they aren't already
            bounds = array(bounds)
        try:  # put the bounds into the proper shape
            bounds = bounds.reshape(-1, 2)
        except:  # if this is impossible, the bounds must be incorrectly defined - warn the user
            exit('Error: bounds is not correctly dimensioned! size given: {} size needed: either (2*n,) or (n, 2)'.format(bounds.shape))
        data = []  # store the global data for the fits
        tic()  # start the timer
        lower_bounds, upper_bounds = bounds[:, 0], bounds[:, 1]  # define the lower and upper bounds as the first and second columns of the bounds
        for fit_attempt in range(num_attempts):  # perform a specified number of fit attempts
            print('{}%'.format(100 * fit_attempt / num_attempts))#, end='\r')
            guess = [uniform(low=low, high=high) for low, high in zip(lower_bounds, upper_bounds)]  # create an initial guess within the given parameter bounds
            # minimize the SSE of the custom function using nelder-mead with the given initial guess and parameter bounds
            results = minimize(self.SSE, x0=guess, args=(lower_bounds, upper_bounds),
                               method='Nelder-Mead', options={'maxiter': maxiter,
                                                              'maxfev': maxiter,
                                                              'xatol': 1e-60,
                                                              'fatol': 1e-60})
            data.append([results.x, results.fun])  # save the final parameters and their cost
        data = array(data, dtype='object')  # convert the fit trial data to a numpy array for easier analysis
        best_fit = data[argmin(data[:, 1])]  # get the trial with the lowest cost
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': toc(True), 'trial_variance': var(data[:, -1])}  # return the data for the best fitting trial

    def fit_slow(self, function, training_data, bounds, maxiter=1000, num_attempts=5):
        '''
        fit experimental observable of your choice to a custom model for the observable using simulated annealing with
        a nelder-mead simplex local search, this is very computationally costly and will take a very long time
        though typically results in much better fits
        :param function: function for the desired observable to be predicted
        :param training_data: either numpy array or list of numpy arrays the experimental data to be replicated by the function being trained
        :param bounds: (n, 2) numpy array of upper and lower bounds: [[lower1, upper1], ... [lowerN, upperN]]
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :param num_attempts: int number of fitting attempts to make per fit, larger number will give more statistically significant results, but
        will take longer
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        self.observable_function = function  # define the observable function for optimization
        if type(
                training_data) is list:  # if there are multiple entries for the training data, put them into row vector form
            training_data = concatenate(training_data)
        self.target_observable = training_data  # define the observable as the target for the objective function to replicate
        if type(bounds) is list:  # convert the function bounds to a numpy array if they aren't already
            bounds = array(bounds)
        try:  # put the bounds into the proper shape
            bounds = bounds.reshape(-1, 2)
        except:  # if this is impossible, the bounds must be incorrectly defined - warn the user
            exit(
                'Error: bounds is not correctly dimensioned! size given: {} size needed: either (2*n,) or (n, 2)'.format(bounds.shape))
        data = []  # store the global data for the fits
        tic()  # start the timer
        lower_bounds, upper_bounds = bounds[:, 0], bounds[:, 1]  # define the lower and upper bounds as the first and second columns of the bounds
        for fit_attempt in range(num_attempts):  # perform a specified number of fit attempts
            print('{}%'.format(100 * fit_attempt / num_attempts))#, end='\r')
            guess = [uniform(low=low, high=high) for low, high in zip(lower_bounds, upper_bounds)]  # create an initial guess within the given parameter bounds
            # minimize the SSE of the custom function within the specified bounds using dual annealing with a nelder mead local search
            results = dual_annealing(self.SSE, bounds, args=(lower_bounds, upper_bounds), maxiter=maxiter,
                                     local_search_options={'method': 'nelder-mead'}, x0=guess)
            data.append([results.x, results.fun])  # store the final parameters and their cost
        data = array(data, dtype='object')  # convert the trial data to a numpy array for easier analysis
        best_fit = data[argmin(data[:, 1])]  # get the trial with the lowest cost
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': toc(True), 'trial_variance': var(data[:, -1])}  # return the trial data of the best fitting parameter set

#@TODO test without simultaneous fits
#@TODO test with retract
#@TODO benchmark all against nonlinear least squares gradient descent methods
#@TODO test log fitting against standard fitting i.e. guessing the order of magnitude of each parameter (10**a rather than a)
#@TODO add the ibw reader

#@TODO THESE ARE FOR LATER IMPLEMENTATIONS OF THE CODE
#@TODO add conical and flat punch indenter options
#@TODO make the map reader and add it to the how-to-guide
#@TODO cuda parallelization of map reader?
#@TODO add public package requirements and get a license
