import os
import numpy as np
import pandas as pd
import math
from scipy import stats
import warnings

import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from itertools import accumulate

con_list = ['1', '2', '3']
#
#
# class pk_milk(type):
#
#     def con_list(self, con_list):
#         self.con_list = con_list

# An example of a class
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
BODYWEIGHT_PATH_MAT = os.path.join('golden', 'bw_cz.mat')  # bodyweight data
LIPID_FRACTION_PATH_MAT = os.path.join('golden', 'flip_cz.mat')

BODYWEIGHT_AND_LIPID_FRACTION_ = os.path.join(
    'monthly_bodyweight_and_lipid_fraction.csv')
BIOMONITORING_DATA_ELIMINATION_ = os.path.join('mikes_et_al_2012.csv')
BIOMONITORING_DATA_ELIMINATION_UNITS = [
    'ng/kg_bw/d', 'years']  # just to keep track
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
        self.con_list = con_list
        self.delta_t = delta_t
        self.n_steps = n_steps
        self.brth_sched = brth_sched

    def dt_from_timesteps_(self, dt=False, tpm=False):
        print("default dt:", self.delta_t)
        print("default timesteps/month:", self.timesteps_per_month)

        if not all([dt, tpm]):
            if dt:
                self.delta_t = dt
                self.timesteps_per_month = float(1 / dt)
                print("new dt:", self.delta_t)
                print("new timesteps/month:", self.timesteps_per_month)

            if tpm:
                self.timesteps_per_month = tpm
                print("new timesteps/month:", self.timesteps_per_month)
                self.delta_t = np.float(1 / self.timesteps_per_month)
                print("new dt:", self.delta_t)
        if all([dt, tpm]):
            warnings.warn('specify dt or timesteps, not both')
        return self

    def n_steps_(self):
        print("default n_steps:", self.n_steps)
        self.n_steps = np.floor((self.y_end - self.y_start) / self.delta_t) + 1
        print("calculated n_steps:", self.n_steps)
        return self

    def intake_call_(self, y_peak):
        num_steps_up = (y_peak - self.y_start) / self.delta_t
        x_up_t_step = np.linspace(self.y_start, y_peak, num=num_steps_up)

        x_up = np.array([self.y_start, y_peak - self.y_start])
        y_up = np.log([1e-10, self.peak_dose_intensity])
        (slope_up, intercept_up, r_value, p_value,
         std_err) = stats.linregress(x_up, y_up)

        # interpolate using the provided slope and intercept to create the
        # regression line
        x_up_interp = np.linspace(self.y_start, y_peak, num_steps_up)
        y_reg_line_up = np.polyval([slope_up, intercept_up], x_up_interp)
        # take the exp of the reg line to return it to an exp fit
        upswing_reg_lin = np.exp(y_reg_line_up)

        # the remaining intensity is set to zero
        num_steps_d = (self.y_end - (y_peak + self.delta_t)) / self.delta_t
        x_d_interp = np.linspace(y_peak + self.delta_t, self.y_end, num_steps_d)
        x_d = np.array([y_peak + self.delta_t, self.y_end])
        y_d = np.log([self.peak_dose_intensity, 1e-10])
        (slope_d, intercept_d, r_value, p_value, std_err) = stats.linregress(x_d, y_d)

        y_reg_line_down = np.polyval([slope_d, intercept_d], x_d_interp)
        dswing_reg_lin = np.exp(y_reg_line_down)

        # concatenate the up and down swing sides
        y_up_down = np.concatenate((upswing_reg_lin, dswing_reg_lin[1:]))
        x_up_down = np.concatenate((x_up_interp, x_d_interp[1:]))

        peak_intensity_spline = interpolate.InterpolatedUnivariateSpline(x_up_down, y_up_down)

        return peak_intensity_spline

    def bm_eval_(self, kin_order=1, plt_kin=True, out_kin=False, use_scipy_exp=False, p0=(1e-6, 1e-6, 1),
                 use_scipy_lin2exp=False):

        kinetics_from_biomonitoring_data = []
        bm_data = pd.read_csv(self.bm_data)
        bm_data_t = bm_data.ix[:, 0]
        colors = plt.cm.rainbow(np.linspace(0, 1, len(self.cte)))

        def exp_fit(x, a, c, d):
            return a * np.exp(-c * x) + d

        x = np.array(bm_data_t)
        yr1 = x[0]

        for congener, c in zip(self.cte, colors):

            if kin_order == 1:

                x = x - x[0]
                y = np.array(bm_data[congener])

                if not all([use_scipy_exp, use_scipy_lin2exp]):

                    if use_scipy_exp:

                        popt, pcov = curve_fit(exp_fit, x, y, p0=p0)

                        if out_kin:
                            print('congoner:', congener)
                            print('exp reg stats: ')
                            print("y = (", popt[0], ") * (e^(", -popt[1], " x )) + ", popt[2])
                            print('note: x can be a range of years (e.g., 2000-2010),')
                            print('or sequential (e.g.,''0-9). exp reg is performed by shifting x')
                            print('such that x[0] is zero. When plotted, it is shifted back.')

                        if plt_kin:
                            # shift the x values back to their respective years
                            x_shift = list(map(lambda x: x + yr1, x))
                            plt.plot(x_shift, y, color=c, marker='o', linestyle='',
                                     label=str('biomonitoring data; ' + congener))
                            plt.plot(x_shift, exp_fit(x, *popt),
                                     color=c, marker='', linestyle='-')

                    if use_scipy_lin2exp:
                        y = np.log(np.array(bm_data[congener]))
                        (slope, intercept, r_value, p_value,
                         std_err) = stats.linregress(x, y)
                        y_regression_line = np.polyval([slope, intercept], x)

                        if out_kin:
                            print('congoner:', congener)
                            print('log-lin reg stats: ')
                            print("ln(y) = (", intercept, ") + (", slope, ") x")

                        if plt_kin:
                            # shift the x values back to their respective years
                            x_shift = list(map(lambda x: x + yr1, x))
                            plt.plot(x_shift, np.exp(y), color=c, marker='o', linestyle='',
                                     label=str('biomonitoring data; ' + congener))
                            plt.plot(x_shift, np.exp(y_regression_line),
                                     color=c, marker='', linestyle='-')


            elif kin_order == 2:
                x = np.log(np.array(bm_data_t))
                y = np.log(np.array(bm_data[congener]))
                (slope, intercept, r_value,
                 p_value, std_err) = stats.linregress(x, y)

                y_regression_line = np.polyval([slope, intercept], x)
                kinetics_from_biomonitoring_data.append(slope)
                if plt_kin:
                    plt.plot(np.exp(x), np.exp(y), color=c,
                             marker='o', linestyle='', label=str(
                            'biomonitoring data; ' + congener))
                    plt.plot(np.exp(x), np.exp(y_regression_line),
                             color=c, marker='', linestyle='-')

        if plt_kin:
            plt.show()

    # automatically generate generational critical age definitions
    def age_timeline(self, brth_sched=False):

        '''
        aig_mother [  0.  25.  50.  75.]
        cbtg_child [  25.   50.   75.  100.]
        aigd_mother [  80.  105.  130.  155.]
        :param brth_sched:
        :return:
        '''
        print('default interval age at which mother gives birth: ', self.brth_age, 'years')
        if brth_sched:
            self.brth_sched = brth_sched
            print('provided birth time-table:', self.brth_sched, 'years')

            if self.gens < len(self.brth_sched):
                print('more birth scheduling info was provided than gens calculated.')
                print('increase number of gens else birth scheduling info will be ignored.')

                self.brth_sched = self.brth_sched[:(self.gens+1)]
                warnings.warn('only the following section of the birth schedule will be used:',self.brth_sched)
            if self.gens > len(self.brth_sched):
                warnings.warn('insufficient birth scheduling information. increase gens to match.')

        if not brth_sched:
            self.brth_sched = np.linspace(0, self.brth_age * self.gens, self.gens + 1)
            print('calculated birth time-table based on birth age:', self.brth_sched, 'years')

        return self

        # aig_mother = linspace(0,
        #                       prg_intrvl * (gens),
        #                       gens + 1)
        #
        # cbtg_child = linspace(prg_intrvl,
        #                       prg_intrvl * (gens + 1),
        #                       gens + 1)
        #
        # aigd_mother = linspace(lifespan,
        #                        (lifespan + (25 * gens)),
        #                        gens + 1)


A = pk_milk(gens=5)
A.dt_from_timesteps_(tpm=5)
A.n_steps_()
# A.bm_eval_(use_scipy_exp=True, plt_kin=False, out_kin=True, kin_order=1)
A.intake_call_(y_peak=35)
A.age_timeline(brth_sched=[0, 15, 25, 35]) #


# pk_milk()


# print (A.gen_delta_t())
# A = Shape.con_list = con_list
# print(A.con_list)
