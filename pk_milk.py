# coding=utf-8
# !/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy import stats
import warnings
from scipy import integrate
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# todo ( create function for lactation schedule and remove init)


class pk_milk():
    def __init__(self,
                 gens=1,
                 y_start=0,
                 y_end=100,
                 lifespan=80,
                 brth_age=25,
                 average_lact_time=[0.5] * 5,
                 k_lac=[1e-2] * 5,
                 k_elim=[np.log(2) / 5] * 5,
                 odes_in_each_generation=2):

        self.gens = gens
        self.y_start = y_start
        self.y_end = y_end
        self.lifespan = lifespan
        self.brth_age = brth_age
        self.average_lact_time = average_lact_time
        self.k_lac = k_lac
        self.k_elim = k_elim
        self.odes_in_each_generation = odes_in_each_generation

    def dt_from_timesteps_(self, timestep_variable, method='timesteps_per_month'):

        if method == 'set_delta_t':
            self.delta_t = timestep_variable
            self.timesteps_per_month = float(1 / self.delta_t)
            print("new delta_t:", self.delta_t)
            print("new timesteps/month:", self.timesteps_per_month)

        if method == 'timesteps_per_month':
            self.timesteps_per_month = timestep_variable
            print("new timesteps/month:", self.timesteps_per_month)
            self.delta_t = np.float(1 / self.timesteps_per_month)
            print("new delta_t:", self.delta_t)

        return self

    def n_steps_(self):
        self.n_steps = np.floor((self.y_end - self.y_start) / self.delta_t) + 1
        print("calculated n_steps:", self.n_steps)
        return self

    def intake_intensity_curve_(self, intake_intensity_data=False, method='points2spline', peak_intensity=False):

        if method == 'points2spline':
            col_year = np.array(intake_intensity_data[:, 0])
            col_congener_intake = np.array(intake_intensity_data[:, 1])
            year_one = col_year[0]

            # Shift the years to zero
            adj_year = col_year - year_one
            if (col_year[0] != self.y_start) or (col_year[-1] != self.y_end):
                warnings.warn('simulation time misaligned with intake curve.')

            intake_intensity_spline = interpolate.InterpolatedUnivariateSpline(adj_year, col_congener_intake)

        if method == 'asymettric_exp_up_and_down':
            num_steps_up = (peak_intensity - self.y_start) / self.delta_t

            x_up = np.array([self.y_start, peak_intensity - self.y_start])
            y_up = np.log([1e-10, peak_intensity])
            (slope_up, intercept_up, r_value, p_value,
             std_err) = stats.linregress(x_up, y_up)

            # create the regression line
            x_up_interp = np.linspace(self.y_start, peak_intensity, num_steps_up)
            y_reg_line_up = np.polyval([slope_up, intercept_up], x_up_interp)
            # take the exp of the reg line to return it to an exp fit
            upswing_reg_lin = np.exp(y_reg_line_up)

            # the remaining intensity is set to zero
            num_steps_d = (self.y_end - (peak_intensity + self.delta_t)) / self.delta_t
            x_d_interp = np.linspace(peak_intensity + self.delta_t, self.y_end, num_steps_d)
            x_d = np.array([peak_intensity + self.delta_t, self.y_end])
            y_d = np.log([peak_intensity, 1e-10])
            (slope_d, intercept_d, r_value, p_value, std_err) = stats.linregress(x_d, y_d)

            y_reg_line_down = np.polyval([slope_d, intercept_d], x_d_interp)
            dswing_reg_lin = np.exp(y_reg_line_down)

            # concatenate the up and down swing sides
            y_up_down = np.concatenate((upswing_reg_lin, dswing_reg_lin[1:]))
            x_up_down = np.concatenate((x_up_interp, x_d_interp[1:]))

            intake_intensity_spline = interpolate.InterpolatedUnivariateSpline(x_up_down, y_up_down)

        self.intake_intensity_curve = intake_intensity_spline

        return self

    def biomonitoring_eval_(self, biomonitoring_data,
                            exponential_fit=lambda x, a, c, d: a * np.exp(-c * x) + d, p0=(1e-6, 1e-6, 1),
                            assumed_kinetic_order=1, method='lin2exp'):
        """
        Evaluating biomonitoring data for an individual congener.

        :param biomonitoring_data: year column and congener blood concentration column e.g.
                1996,13.92
                1997,17.43
                1998,11.33
        :param exponential_fit: fit function
        :param p0: initial conditions for non-linear fit
        :param assumed_kinetic_order: first or second order
        :param methodology: lin2exp or exp
        :return: year array, fitted y value array, slope
        """

        col_year = np.array(biomonitoring_data[:, 0]).flatten('C')
        col_congener = np.array(biomonitoring_data[:, 1]).flatten('C')

        year_one = col_year[0]

        if assumed_kinetic_order == 1:
            # Shift the years to zero
            adj_year = col_year - year_one

            if method == 'lin2exp':
                log_congener = np.log(col_congener)
                slope, intercept, r_value, p_value, std_err = stats.linregress(adj_year, log_congener)
                return adj_year, np.polyval([slope, intercept], adj_year), -slope, -r_value

            elif method == 'exp':
                (popt, pcov) = curve_fit(exponential_fit, col_year, col_congener, p0=p0)
                return adj_year, exponential_fit(adj_year, *popt), np.exp(-popt[1])

        elif assumed_kinetic_order == 2:
            log_year = np.log(col_year)
            log_congener = np.log(col_congener)
            (slope, intercept, r_value, p_value, std_err) = stats.linregress(log_year, log_congener)
            return log_year, np.polyval([slope, intercept], log_year), -slope, -r_value

    def age_timeline_(self, brth_sched=False):
        '''
        Generate generational critical age definitions

        :param  brth_sched: user provided birth schedule
                  e.g., brth_sched=[0, 15, 25, 35]
        :return: self
        '''
        print('default interval age at which mother gives birth: ', self.brth_age, 'years')
        self.brth_sched = np.linspace(0, self.brth_age * self.gens, self.gens + 1)
        print('calculated birth time-table based on birth age:', self.brth_sched, 'years')

        self.cbtg_child = list(map(lambda x: x + self.brth_age, self.brth_sched))
        print('generational birth schedule:', self.cbtg_child, 'years')

        self.aigd_mother = list(map(lambda x: x + self.lifespan, self.brth_sched))
        print('generational death schedule:', self.aigd_mother, 'years')

        if brth_sched:
            self.brth_sched = brth_sched
            print('provided birth time-table:', self.brth_sched, 'years')

            if self.gens < len(self.brth_sched):
                print('more birth scheduling info was provided than gens calculated.')
                print('increase number of gens else birth scheduling info will be ignored.')

                self.brth_sched = self.brth_sched[:(self.gens + 1)]
                warnings.warn('only the following section of the birth schedule will be used:', self.brth_sched)
            if self.gens > len(self.brth_sched):
                warnings.warn('insufficient birth scheduling information. increase gens to match.')

        return self

    # calculate mass of lipids for each month
    def cal_m_lip(row):
        mass_lipids = row.avg_monthly_bodyweight_kg * row.lipid_fraction
        return mass_lipids

    def lipid_mass_from_bw_and_lipid_fraction(self, bodyweight_and_lipid_fraction_data):
        """
        :param  bodyweight_and_lipid_fraction_data:
        matrix format:
        bodyweight_kg,  lipid_fraction
        3.79233266687,  0.213947530915
        4.67272833393,  0.225984883013
        5.46991188708,  0.237058988856

        :return: array
        """
        bodyweight_in_kg = np.array(bodyweight_and_lipid_fraction_data[:, 0]).flatten('C')
        lipid_fraction = np.array(bodyweight_and_lipid_fraction_data[:, 1]).flatten('C')

        # multiply every bodyweight by its corresponding lipid fraction
        lipid_mass = np.multiply(bodyweight_in_kg, lipid_fraction)
        return lipid_mass

    def age_splines_(self, lipid_mass_array):
        """
        create a list of age splines for each generation.
        :param bodyweight_and_lipid_fraction_matrix:
        :return: two arrays of spline objects describing each generation's lipid mass or
        lipid mass change (derivative)
        """

        age_spline_derivative = pd.DataFrame()
        age_spline = pd.DataFrame()

        num_steps_before = []
        num_steps_after = []
        for gen in range(0, self.gens):
            num_steps_before.append(np.int((self.brth_sched[gen] - self.y_start) / self.delta_t))
            num_steps_after.append(np.int((self.y_end - self.aigd_mother[gen]) / self.delta_t))

        for gen in range(0, self.gens):
            y_before = []
            x_before = []
            y_after = []
            x_after = []
            if num_steps_before[gen] == 0:
                y_gen = lipid_mass_array
                x_gen = np.linspace(self.brth_sched[gen], self.aigd_mother[gen], len(y_gen))

                if num_steps_after[gen] > 0:
                    y_after = np.zeros(np.int(num_steps_after[gen]))
                    x_after = np.linspace(self.aigd_mother[gen], self.y_end, num_steps_after[gen])

                # concatenate everything, but remove overlaps.
                y_all = np.concatenate([y_before[:-1], y_gen, y_after[1:-1]])
                x_all = np.concatenate([x_before[:-1], x_gen, x_after[1:-1]])

                age_spline_n = interpolate.InterpolatedUnivariateSpline(x_all, y_all, k=1)
                age_spline_d = age_spline_n.derivative()


            elif num_steps_before[gen] > 0:
                y_before = np.zeros(np.int(num_steps_before[gen]))
                x_before = np.linspace(self.brth_sched[gen - 1], self.brth_sched[gen], num_steps_before[gen])

                y_gen = lipid_mass_array
                x_gen = np.linspace(self.brth_sched[gen], self.aigd_mother[gen], len(y_gen))

                if num_steps_after[gen] > 0:
                    y_after = np.zeros(np.int(num_steps_after[gen]))
                    x_after = np.linspace(self.aigd_mother[gen], self.y_end, num_steps_after[gen])

                # concatenate everything, but remove overlaps.
                y_all = np.concatenate([y_before[:-1], y_gen, y_after[1:-1]])
                x_all = np.concatenate([x_before[:-1], x_gen, x_after[1:-1]])

                age_spline_n = interpolate.InterpolatedUnivariateSpline(x_all, y_all, k=1)
                age_spline_d = age_spline_n.derivative()

            age_spline = age_spline.append([age_spline_n], 1)
            age_spline_derivative = age_spline_derivative.append([age_spline_d], 1)

        self.age_spline = age_spline[0].ravel()
        self.age_spline_derivative = age_spline_derivative[0].ravel()

        return self

    def body_mass(self, t, y):
        '''
        The first generation's mass balance should be specified above. Every other generations balance,
        assuming it's the same, can be specified here.

        Note that 'cntr' variable is a loop tracking tool. To specify that the previous box's mass balance
        should be employed, use itr_mtrx[0][cntr]-X), where X is the total number of mass balances - 1.
        This is because the array goes from 0 to 3, with length 4.

        Use the np.int() to surround the itr_mtrx calls because arrays should be tracked with
        integers, not floats.

        You can add more mass balances, but do not change dydt_matrix[0][gen] label.
        This is because these are a placeholder variables that are reorganized from an array to
        a matrix.

        For notes:
        aig_mother[gen]  # age the mother gives birth
        cbtg_child[gen]  # year child is born from previous gen
        aigd_mother[gen]  # age of death of the mother

       '''
        cntr = 0

        # auto setup for multi-generational ode
        odes_per_gen = range(0, self.odes_in_each_generation)
        dydt_matrix = np.zeros(shape=(len(odes_per_gen), self.gens), dtype=object)

        order_array_counter = np.array(range(0, self.gens * len(odes_per_gen)))
        itr_mtrx = order_array_counter.reshape((len(odes_per_gen), self.gens), order='F')

        def k_lac_m2c_(t, gen):

            if np.all((t >= self.cbtg_child[gen]) & (t <= self.cbtg_child[gen] + self.average_lact_time[gen])):
                k_lac_m2c = self.k_lac[gen]

            else:
                k_lac_m2c = 0
            return k_lac_m2c

        for gen in range(0, self.gens):

            k_lac_mother_from_daughter = k_lac_m2c_(t, gen)

            if gen == 0:
                # print(self.age_spline_derivative[1](t), t)
                dydt_matrix[0][cntr] = self.age_spline_derivative[0](t)
                dydt_matrix[1][cntr] = self.intake_intensity_curve(t) * y[0] \
                                       - self.k_elim[gen] * y[1] \
                                       - k_lac_mother_from_daughter * y[1]
                cntr = np.int(cntr + 1)


            elif gen > 0:
                k_lac_mother_from_daughter = k_lac_m2c_(t, gen)
                k_lac_child_from_mother = k_lac_m2c_(t, gen - 1)
                dydt_matrix[0][cntr] = self.age_spline_derivative[gen](t)
                dydt_matrix[1][cntr] = self.intake_intensity_curve(t) * y[np.int(itr_mtrx[0][cntr])] \
                                       + k_lac_mother_from_daughter * y[np.int(itr_mtrx[1][cntr - 1])] \
                                       - k_lac_child_from_mother * y[np.int(itr_mtrx[1][cntr])] \
                                       - self.k_elim[gen] * y[np.int(itr_mtrx[1][cntr])]

                cntr = np.int(cntr + 1)

        dydt = np.ravel(dydt_matrix, order='F')

        return dydt

    def generation_mass_balance(self, quickplot=False):

        t = np.zeros((np.int(self.n_steps), 1))

        # use ``vode`` with "backward differentiation formula" or 'bdf'
        r = integrate.ode(self.body_mass).set_integrator('vode',
                                                         order=4,
                                                         nsteps=self.n_steps,
                                                         min_step=1e-15,
                                                         method='bdf')

        y0 = np.zeros((np.int(self.gens * self.odes_in_each_generation), 1))
        r.set_initial_value(y0, self.y_start)

        # create vectors to store trajectories
        ode_init = np.zeros((np.int(self.n_steps * self.gens * self.odes_in_each_generation), 1))
        ode_init_matrix = ode_init.reshape((self.odes_in_each_generation * self.gens,
                                            np.int(self.n_steps)), order='F')

        # initialize k for while loop
        iter_odes = range(0, self.odes_in_each_generation * self.gens, 1)
        k = 1
        while r.successful() and k < self.n_steps:
            r.integrate(r.t + self.delta_t)
            t[k] = r.t
            for ode in iter_odes:
                ode_init_matrix[ode][k] = r.y[ode]
            k += 1

        if quickplot:
            for ode in iter_odes:
                ax1 = plt.subplot(len(iter_odes), 1, iter_odes[ode] + 1)
                plt.plot(t, ode_init_matrix[ode][:])
                ax1.plot(t, ode_init_matrix[ode][:])
                ax1.set_xlim(self.y_start, self.y_end)
                ax1.grid('on')

            plt.xlim(self.y_start, self.y_end)
            plt.show()

        self.dydt_solution = ode_init_matrix
        return self
