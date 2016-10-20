import os
import numpy as np
import pandas as pd
import math
from scipy import stats
import warnings

import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

"""
             bm_data=BIOMONITORING_DATA_ELIMINATION_,
             cong_spag=_cong_spag,
             cte=cte_,
             kinetics_order=KINETICS_ORDER_,
             gens=1,
             plot_kinetics=False,
             study_start_year=STUDY_START_YEAR_,
             study_end_year=STUDY_END_YEAR_,
             lifespan=AGE_MAX_IN_YEARS_,
             prg_intrvl=AGE_MOM_IN_YEARS_,
             bw_frac_lip=BODYWEIGHT_AND_LIPID_FRACTION_,
             timesteps_per_month=TIMESTEPS_PER_MONTH,
             average_lact_time=AVERAGE_LACT_TIME_,
             k_lac=K_LAC_GEN,
             k_elim=K_ELIM,
             peak_dose_intensity=PEAK_DOSING_INTENSITY_)
"""
SOURCE_PATH_MAT = os.path.join('golden', 'P_ritter.mat')
DATA_PATH_MAT = os.path.join('golden', 'data_all_mean.mat')
BODYWEIGHT_PATH_MAT = os.path.join('golden', 'bw_cz.mat')
LIPID_FRACTION_PATH_MAT = os.path.join('golden', 'flip_cz.mat')

BODYWEIGHT_AND_LIPID_FRACTION_ = os.path.join(
    'monthly_bodyweight_and_lipid_fraction.csv')
BIOMONITORING_DATA_ELIMINATION_ = os.path.join('mikes_et_al_2012.csv')
BIOMONITORING_DATA_ELIMINATION_UNITS = [
    'ng/kg_bw/d', 'years']
_cong_spag = os.path.join(
    'cogener_start_peak_age_group.csv')
cte_ = ['hcb', 'ppddt',
        'ppdde', 'betahch', 'gammahch', 'pcb138']

_CONGENERS_HALFLIFE_DICT = {'hcb': [10, np.log(2) / 5, np.log(2) / 5],
                            'ppddt': [5, np.log(2) / 5, np.log(2) / 5]}

cte_ = list(_CONGENERS_HALFLIFE_DICT.keys())

TIMESTEPS_PER_MONTH_ = 10.0

POP_START_YEAR_ = 0
STUDY_START_YEAR_ = 0
STUDY_END_YEAR_ = 2100 - 1921
AGE_MAX_IN_YEARS_ = 80  # [years]
ABSORP_FACTOR_ = 0.9  # absorbption factor [-]
AVERAGE_LACT_TIME_ = [3] * 5  # months breastfeeding (generic source?)
K_LAC_GEN = [1e-2] * 5
AGE_MOM_IN_YEARS_ = 25

PEAK_DOSING_INTENSITY_ = 80  # [ng/kg/d]
REF_FEMALE_BODYWEIGHT_KG_ = 69
AVERAGE_MOTHER_LIPID_FRACTION_ = 0.35  # I MADE THIS UP - DEFAULT CHANGE
KINETICS_ORDER_ = 1
AVERAGE_BIRTH_LIPID_FRACTION_ = 0.215
body_w_ave = 3.1744
K_ELIM = np.log(2) / 5
PEAK_YEAR = [35]


class pk_milk():
    # print (TIMESTEPS_PER_MONTH)
    def __init__(self,
                 bm_data=BIOMONITORING_DATA_ELIMINATION_,
                 cong_spag=_cong_spag,
                 cte=cte_,
                 kinetics_order=KINETICS_ORDER_,
                 gens=1,
                 plot_kinetics=False,
                 y_start=0,
                 y_end=100,
                 lifespan=80,
                 brth_age=25,
                 brth_sched=[0, 25, 50, 75],
                 bw_frac_lip=BODYWEIGHT_AND_LIPID_FRACTION_,
                 timesteps_per_month=1,
                 average_lact_time=AVERAGE_LACT_TIME_,
                 k_lac=K_LAC_GEN,
                 k_elim=np.log(2) / 5,
                 peak_dose_intensity=80,
                 delta_t=1,
                 n_steps=1):

        self.bm_data = bm_data
        self.cong_spag = cong_spag
        self.cte = cte
        self.kinetics_order = kinetics_order
        self.gens = gens
        self.plot_kinetics = plot_kinetics
        self.y_start = y_start
        self.y_end = y_end
        self.lifespan = lifespan
        self.brth_age = brth_age
        self.bw_frac_lip = bw_frac_lip
        self.timesteps_per_month = timesteps_per_month
        self.average_lact_time = average_lact_time
        self.k_lac = k_lac
        self.k_elim = k_elim
        self.peak_dose_intensity = peak_dose_intensity
        self.delta_t = delta_t
        self.n_steps = n_steps
        self.brth_sched = brth_sched

    def dt_from_timesteps_(self, timestep_variable, method='timesteps_per_month'):
        print("default delta_t:", self.delta_t)
        print("default timesteps/month:", self.timesteps_per_month)

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
        print("default n_steps:", self.n_steps)
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
            y_up = np.log([1e-10, self.peak_dose_intensity])
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
            y_d = np.log([self.peak_dose_intensity, 1e-10])
            (slope_d, intercept_d, r_value, p_value, std_err) = stats.linregress(x_d, y_d)

            y_reg_line_down = np.polyval([slope_d, intercept_d], x_d_interp)
            dswing_reg_lin = np.exp(y_reg_line_down)

            # concatenate the up and down swing sides
            y_up_down = np.concatenate((upswing_reg_lin, dswing_reg_lin[1:]))
            x_up_down = np.concatenate((x_up_interp, x_d_interp[1:]))

            intake_intensity_spline = interpolate.InterpolatedUnivariateSpline(x_up_down, y_up_down)

        return intake_intensity_spline

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

        col_year = np.array(biomonitoring_data[:, 0])
        col_congener = np.array(biomonitoring_data[:, 1])
        year_one = col_year[0]

        if assumed_kinetic_order == 1:
            # Shift the years to zero
            adj_year = col_year - year_one

            if method == 'lin2exp':
                log_congener = np.log(col_congener)
                (slope, intercept, r_value, p_value, std_err) = stats.linregress(adj_year, log_congener)
                return adj_year, np.polyval([slope, intercept], adj_year), slope

            elif method == 'exp':
                (popt, pcov) = curve_fit(exponential_fit, col_year, col_congener, p0=p0)
                return adj_year, exponential_fit(adj_year, *popt), np.exp(-popt[1])

        elif assumed_kinetic_order == 2:
            log_year = np.log(col_year)
            log_congener = np.log(col_congener)
            (slope, intercept, r_value, p_value, std_err) = stats.linregress(log_year, log_congener)
            return log_year, np.polyval([slope, intercept], log_year), slope

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
        self.aigd_mother = list(map(lambda x: x + self.lifespan, self.brth_sched))

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

    def k_lac_m2c_(self, t, gen):
        self.k_lac_m2c = 0
        if (t >= self.cbtg_child[gen]) & (t <= self.cbtg_child[gen] + self.average_lact_time[gen]):
            self.k_lac_m2c = self.k_lac[gen]

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
        bodyweight_in_kg = np.array(bodyweight_and_lipid_fraction_data[:, 0])
        lipid_fraction = np.array(bodyweight_and_lipid_fraction_data[:, 1])

        # multiply every bodyweight by its corresponding lipid fraction
        lipid_mass = np.multiply.outer(bodyweight_in_kg, lipid_fraction).ravel()
        return lipid_mass




    def create_age_splines_(self, lipid_mass_array):
        """
        create a list of age splines for each generation.
        :param bodyweight_and_lipid_fraction_matrix:
        :return: array of spline objects for each generation's life.
        """

        num_steps_before = []
        num_steps_after = []
        for gen in range(0, self.gens):
            num_steps_before.append(np.int((self.brth_sched[gen] - self.y_start) / self.delta_t))
            num_steps_after.append(np.int((self.y_end - self.aigd_mother[gen]) / self.delta_t))

        age_spline_df = pd.DataFrame()
        for gen in range(0, self.gens):
            y_before = []
            x_before = []
            y_after = []
            x_after = []
            y_gen = lipid_mass_array
            if num_steps_before[gen] == 0:

                x_gen = np.linspace(self.aigd_mother[gen], self.aigd_mother[gen], len(y_gen))

                if num_steps_after[gen] > 0:
                    y_after = np.zeros(np.int(num_steps_after[gen]))
                    x_after = np.linspace(self.aigd_mother[gen], self.y_end, num_steps_after[gen])

                # concatenate everything, but remove overlaps.
                y_all = np.concatenate([y_before[:-1], y_gen, y_after[1:-1]])
                x_all = np.concatenate([x_before[:-1], x_gen, x_after[1:-1]])

                age_spline = interpolate.InterpolatedUnivariateSpline(x_all, y_all).derivative()

            elif num_steps_before[gen] > 0:
                y_before = np.zeros(np.int(num_steps_before[gen]))
                x_before = np.linspace(self.brth_sched[gen - 1], self.brth_sched[gen], num_steps_before[gen])

                x_gen = np.linspace(self.brth_sched[gen], self.aigd_mother[gen], len(y_gen))

                if num_steps_after[gen] > 0:
                    y_after = np.zeros(np.int(num_steps_after[gen]))
                    x_after = np.linspace(self.aigd_mother[gen], self.y_end, num_steps_after[gen])

                # concatenate everything, but remove overlaps.
                y_all = np.concatenate([y_before[:-1], y_gen, y_after[1:-1]])
                x_all = np.concatenate([x_before[:-1], x_gen, x_after[1:-1]])

                age_spline = interpolate.InterpolatedUnivariateSpline(x_all, y_all).derivative()

            age_spline_df = age_spline_df.append([age_spline], 1)

        return age_spline_df[0].ravel()

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

        for gen in range(0, self.gens):

            k_lac_m2c_r = k_lac_m2c_(self, t, gen)
            k_lac_2cfm_r = k_lac_m2c_(self, t, gen)

            if gen == 0:
                dydt_matrix[0][cntr] = bw_spl_der[0](t)
                dydt_matrix[1][cntr] = intakeCall(t) * y[0] \
                                       - k_elim * y[1] \
                                       - k_lac_m2c_r * y[1]
                cntr = np.int(cntr + 1)

            elif gen > 0:
                dydt_matrix[0][cntr] = bw_spl_der[gen](t)
                dydt_matrix[1][cntr] = intakeCall(t) * y[np.int(itr_mtrx[0][cntr])] \
                                       + k_lac_2cfm_r * y[np.int(itr_mtrx[1][cntr - 1])] \
                                       - k_lac_m2c_r * y[np.int(itr_mtrx[1][cntr])] \
                                       - k_elim * y[np.int(itr_mtrx[1][cntr])]

                cntr = np.int(cntr + 1)

        dydt = np.ravel(dydt_matrix, order='F')

        return dydt

    def generation_mass_balance(self, y,
                                congener,
                                gens,
                                simulation_start,
                                simulation_end,
                                delta_t,
                                intakeCall,
                                prg_intrvl,
                                lifespan,
                                num_steps,
                                bw_frac_lip,
                                k_lac,
                                average_lact_time,
                                k_elim):

        # dynamically set the number of odes as a function of gens
        n = len(y) * gens
        num_odes_in_gen = 2
        n = num_odes_in_gen * gens
        dydt = np.zeros((n, 1))

        # simulation bounds
        t_start = np.float(simulation_start)
        t_final = np.float(simulation_end)

        start_array_dydt = linspace(0, num_odes_in_gen * (gens - 1), gens)

        # automatically generate generational critical age definitions
        aig_mother = linspace(0,
                              prg_intrvl * (gens),
                              gens + 1)

        cbtg_child = linspace(prg_intrvl,
                              prg_intrvl * (gens + 1),
                              gens + 1)

        aigd_mother = linspace(lifespan,
                               (lifespan + (25 * gens)),
                               gens + 1)

        print("aig_mother", aig_mother)  # age the mother gives birth
        print("cbtg_child", cbtg_child)  # year child is born from previous gen
        print("aigd_mother", aigd_mother)  # age of death of the mother

        for gen in range(0, gens + 1):
            if np.all(gens >= 1):
                start_dydt_gen = []
                start_dydt_gen.append(start_array_dydt)
                for i in range(1, gens - 1):
                    start_dydt_gen.append([x + i for x in start_array_dydt])
                start_dydt_gen = np.array(start_dydt_gen)

        odes_per_gen = range(0, num_odes_in_gen)
        dydt_matrix = np.zeros(shape=(len(odes_per_gen),
                                      gens),
                               dtype=object)

        order_array_counter = np.array(range(0, gens * len(odes_per_gen)))
        itr_mtrx = order_array_counter.reshape((len(odes_per_gen), gen),
                                               order='F')

        bw_spl_der = age_splines(gens,
                                 aig_mother,
                                 aigd_mother,
                                 t_start,
                                 t_final,
                                 delta_t,
                                 bw_frac_lip)

        t = np.zeros((np.int(num_steps), 1))

        # use ``vode`` with "backward differentiation formula" or 'bdf'
        r = integrate.ode(body_mass).set_integrator('vode',
                                                    order=4,
                                                    nsteps=num_steps,
                                                    min_step=1e-10,
                                                    method='bdf')

        y0 = np.zeros((np.int(gens * num_odes_in_gen), 1))
        r.set_initial_value(y0, t_start)

        # create vectors to store trajectories
        ode_init = np.zeros((np.int(num_steps) * gens * num_odes_in_gen))
        ode_init_matrix = ode_init.reshape((num_odes_in_gen * gens,
                                            np.int(num_steps)),
                                           order='F')

        iter_odes = range(0, num_odes_in_gen * gens, 1)

        # initialize k for while loop
        k = 1
        while r.successful() and k < num_steps:
            r.integrate(r.t + delta_t)
            t[k] = r.t
            for ode in iter_odes:
                ode_init_matrix[ode][k] = r.y[ode]
            k += 1

        for ode in iter_odes:
            ax1 = plt.subplot(len(iter_odes), 1, iter_odes[ode] + 1)
            plt.plot(t, ode_init_matrix[ode][:])
            ax1.plot(t, ode_init_matrix[ode][:])
            ax1.set_xlim(t_start, t_final)
            ax1.grid('on')

        plt.xlim(t_start, t_final)
        plt.legend(iter_odes)
        plt.show()

        return (y, t)


A = pk_milk(gens=4)
A.dt_from_timesteps_(timestep_variable=5,method='timesteps_per_month')
A.n_steps_()
# A.bm_eval_(use_scipy_exp=True, plt_kin=False, out_kin=True, kin_order=1)
print(A.intake_intensity_curve_(method='asymettric_exp_up_and_down', peak_intensity=35)(25))
# A.age_timeline_()  #



# determine lipid mass from bodyweight and lipid fraction.
# import matrix 'bodyweight_and_lipid_fraction_data'
# example
loc = './monthly_bodyweight_and_lipid_fraction.csv'
bodyweight_and_lipid_fraction_dataframe = pd.read_csv(loc)
bodyweight_and_lipid_fraction_matrix =bodyweight_and_lipid_fraction_dataframe.as_matrix()

lipid_mass_array =  (A.lipid_mass_from_bw_and_lipid_fraction(
    bodyweight_and_lipid_fraction_data=bodyweight_and_lipid_fraction_matrix))
print(lipid_mass_array)

# create age spline array that holds the continuous functions describing the mass of fat in the body.
# rerun age_timeline_ function.
A.age_timeline_()
print(A.create_age_splines_(lipid_mass_array=lipid_mass_array))
