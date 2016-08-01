# coding=utf-8
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd
from scipy import polyval, stats
import numpy as np
import os
import matplotlib.cm as cm
import scipy.interpolate as interpolate
from scipy.integrate import odeint
from scipy.ndimage.interpolation import shift
from scipy import integrate
from matplotlib.pylab import *

"""
Notes: It doesn't make sense that the first order up/down coefficients are statically defined. The upgradient curve
would have eq: y = e^(kdec_up*x). If you fit a first order coefficient, such that the curve started at 0 and went
through the peak intensity, it would result in only one curve. Regarldess of the kdec_up value.

"""
def mat_to_df(path):
    mat = scipy.io.loadmat(path,
                           matlab_compatible=True,
                           variable_names=None)
    mat_heads_init = list(mat)[0]
    mat_matrix = mat[mat_heads_init]
    df = pd.DataFrame(mat_matrix)
    return df


# These will be removed after all have been converted to .csv files.
# TODO (Figure out what some of these are...).
SOURCE_PATH_MAT = os.path.join('golden', 'P_ritter.mat')
DATA_PATH_MAT = os.path.join('golden', 'data_all_mean.mat')
BODYWEIGHT_PATH_MAT = os.path.join('golden', 'bw_cz.mat')  # bodyweight data
LIPID_FRACTION_PATH_MAT = os.path.join('golden', 'flip_cz.mat')

BODYWEIGHT_AND_LIPID_FRACTION_ = os.path.join(
    'monthly_bodyweight_and_lipid_fraction.csv')
BIOMONITORING_DATA_ELIMINATION_ = os.path.join('mikes_et_al_2012.csv')
BIOMONITORING_DATA_ELIMINATION_UNITS = ['ng/kg_bw/d', 'years']  # just to keep track
_CONGENER_START_PEAK_AGE_GROUP = os.path.join('cogener_start_peak_age_group.csv')
_CONGENERS_TO_EVALUATE = ['hcb', 'ppddt', 'ppdde', 'betahch', 'gammahch', 'pcb138']

_CONGENERS_HALFLIFE_DICT = {'hcb': [10, np.log(2) / 5, np.log(2) / 5],
                            'ppddt': [5, np.log(2) / 5, np.log(2) / 5]}

_CONGENERS_TO_EVALUATE = list(_CONGENERS_HALFLIFE_DICT.keys())

TIMESTEPS_PER_MONTH = 1
N_MALES_ = 100
N_FEMALES_ = 100
STUDY_START_YEAR_ = 1994
STUDY_END_YEAR_ = 2009
AGE_MAX_IN_YEARS_ = 80  # [years]
ABSORP_FACTOR_ = 0.9  # absorbption factor [-]
AVERAGE_LACT_TIME_MONTHS_ = 6.  # months breastfeeding (generic source?)
AGE_MOM_IN_YEARS_ = 25
AGE_GROUPINGS_ = range(15, 45, 5)  # start, stop, step
PEAK_DOSING_INTENSITY_ = 80  # [ng/kg/d]
N_OPTIMIZATION_RUNS_ = 1
REF_FEMALE_BODYWEIGHT_KG_ = 69
AVERAGE_MOTHER_LIPID_FRACTION_ = 0.35  # I MADE THIS UP - DEFAULT CHANGE
KINETICS_ORDER_ = 1
AVERAGE_BIRTH_LIPID_FRACTION_ = 0.215
AVERAGE_BIRTH_BODYWEIGHT_KG_ = 3.1744


def biomonitoring_data_error_handling(a, b):
    eval_bio_diff = list(set(a) - set(b))

    if not eval_bio_diff:
        print '''The following chemicals with biomonitoring data
         will not be evaluated'''.replace(
            '/t', ''), eval_bio_diff
    else:
        print ''' Biomonitoring data may be available, but not
        evaluating chemicals:'''.replace(
            '/t', ''), eval_bio_diff


def main(n_optimization_runs=N_OPTIMIZATION_RUNS_,
         biomonitoring_data=BIOMONITORING_DATA_ELIMINATION_,
         congener_start_peak_age_group=_CONGENER_START_PEAK_AGE_GROUP,
         congeners_to_evaluate=_CONGENERS_TO_EVALUATE,  # list
         kinetics_order=KINETICS_ORDER_,
         plot_kinetics=False,
         n_male=N_MALES_,
         n_female=N_FEMALES_,
         study_start_year=STUDY_START_YEAR_,
         study_end_year=STUDY_END_YEAR_,
         average_lact_time_months_=AVERAGE_LACT_TIME_MONTHS_,
         age_max_in_years=AGE_MAX_IN_YEARS_,
         age_mom_in_years=AGE_MOM_IN_YEARS_,
         average_mother_bodyweight=REF_FEMALE_BODYWEIGHT_KG_,
         average_mother_lipid_fraction=AVERAGE_MOTHER_LIPID_FRACTION_,
         absorbtion_factor=ABSORP_FACTOR_,
         average_birth_lipid_fraction_=AVERAGE_BIRTH_LIPID_FRACTION_,
         average_birth_bodyweight_kg_=AVERAGE_BIRTH_BODYWEIGHT_KG_,
         bodyweight_and_lipid_fraction=BODYWEIGHT_AND_LIPID_FRACTION_,
         timesteps_per_month=TIMESTEPS_PER_MONTH):

    # kinetics - import biomonitoring data if needed
    # TODO(build in argparse here for csv extension)
    kinetics_from_biomonitoring_data = []
    if biomonitoring_data:
        biomonitoring_data = pd.read_csv(biomonitoring_data)
        biomonitoring_data_time = biomonitoring_data.ix[:, 0]
        congeners_to_evaluate = congeners_to_evaluate
        biomonitoring_data_congoners = list(biomonitoring_data[-1:])

        biomonitoring_data_error_handling(
            congeners_to_evaluate, biomonitoring_data_congoners)

        colors = cm.rainbow(np.linspace(0, 1, len(congeners_to_evaluate)))
        for congener, c in zip(congeners_to_evaluate, colors):
            if kinetics_order == 1:
                x = np.array(biomonitoring_data_time)
                y = np.log(np.array(biomonitoring_data[congener]))
                (slope, intercept, r_value,
                 p_value, std_err) = stats.linregress(x, y)

                # TODO(output the kinetic rates and fit values to a table)
                y_regression_line = polyval([slope, intercept], x)
                kinetics_from_biomonitoring_data.append(slope)
                if plot_kinetics:
                    plt.plot(x, np.exp(y), color=c,
                             marker='o', linestyle='', label=str(
                            'biomonitoring data; ' + congener))
                    plt.plot(x, np.exp(y_regression_line),
                             color=c, marker='', linestyle='-')

            elif kinetics_order == 2:
                x = np.log(np.array(biomonitoring_data_time))
                y = np.log(np.array(biomonitoring_data[congener]))
                (slope, intercept, r_value,
                 p_value, std_err) = stats.linregress(x, y)

                # TODO(output the kinetic rates and fit values to a table)
                y_regression_line = polyval([slope, intercept], x)
                kinetics_from_biomonitoring_data.append(slope)
                if plot_kinetics:
                    plt.plot(np.exp(x), np.exp(y), color=c,
                             marker='o', linestyle='', label=str(
                            'biomonitoring data; ' + congener))
                    plt.plot(np.exp(x), np.exp(yy_regression_line),
                             color=c, marker='', linestyle='-')

        if plot_kinetics:
            plt.show()

    # assigning chemical specific parameters
    congener_start_peak_age_group_df = pd.read_csv(
        congener_start_peak_age_group)
    biomonitoring_data_error_handling(congeners_to_evaluate,
                                      congener_start_peak_age_group_df)


    # calculate mass of lipids for each month
    def calculate_mass_lipids(row):
        mass_lipids = row.avg_monthly_bodyweight_kg * row.lipid_fraction
        return mass_lipids

    if bodyweight_and_lipid_fraction:
        bodyweight_and_lipid_fraction = pd.read_csv(
            bodyweight_and_lipid_fraction)
        bodyweight_and_lipid_fraction[
            'mass_lipids_kg'] = bodyweight_and_lipid_fraction.apply(
            calculate_mass_lipids, 1)

        # calculate growth kinetics at timestep.
        # step 1: paramterize growth with spline interpretation
        y_bm = np.array(bodyweight_and_lipid_fraction['mass_lipids_kg'])
        x_bm = np.linspace(1, len(y_bm), len(y_bm))
        spl_int_bw_lip_mass = interpolate.InterpolatedUnivariateSpline(
            x_bm, y_bm)

        # step 2: get derivative of spline (for kinetics)
        spl_int_bw_lip_mass_deriv = spl_int_bw_lip_mass.derivative()

        # interpolate in between spline points based on timestamp
        timesteps = timesteps_per_month * len(bodyweight_and_lipid_fraction)
        timesteps_array = np.linspace(1, len(y_bm), timesteps)
        interpolated_bw_lip_mass = spl_int_bw_lip_mass(timesteps_array)
        delta_bw_lip_mass = shift(interpolated_bw_lip_mass, 1, cval=np.NaN)

        # replace the first value (which was shifted).
        # calculate the first delta mass then calculate the kinetics
        delta_bw_lip_mass[0] = bodyweight_and_lipid_fraction['mass_lipids_kg'][
                                   0] - average_birth_bodyweight_kg_ * average_birth_lipid_fraction_

        interpolated_bw_lip_mass_deriv = spl_int_bw_lip_mass_deriv(timesteps_array)

        growth_kinetics_df = pd.DataFrame()
        growth_kinetics_df['mass_lip_kg'] = interpolated_bw_lip_mass
        growth_kinetics_df['growth_kinetics_dKg/dt'] = interpolated_bw_lip_mass_deriv

        # TODO(include function to print bodyweight and derivative (kinetics))
        # plt.plot(timesteps_array, interpolated_bw_lip_mass, 'b')
        # plt.plot(timesteps_array, interpolated_bw_lip_mass_deriv, 'r')
        # plt.show()

    else:
        print('''please include .csv file containing the'average_body_weight_in_kg'
            in the first column and the 'lipid_fraction'
            in the second column''').replace(
            '/t', '')

    # calculating breastfeeding kinetics

    # [months] age of mother at birth in months
    age_mom_in_months = age_mom_in_years * 12

    # [years] age in years with interval of 1 year
    age_in_years_array = range(1, age_max_in_years)  # start, stop, step=1

    # [years] age in years with interval of 1 month
    age_years_in_months = np.linspace(start=1.0 / 12,
                                      stop=age_max_in_years,
                                      num=age_max_in_years * 12,
                                      endpoint=True)
    lact_kinetics = pd.DataFrame()
    lact_kinetics['age_index_years'] = age_years_in_months
    lact_kinetics['age_index_months'] = lact_kinetics['age_index_years'] * 12

    # # print lact_kinetics.head()
    # p_source = mat_to_df(SOURCE_PATH_MAT).T
    # all_mean_data = mat_to_df(DATA_PATH_MAT).to_csv('test.csv')


    def set_peak_intesity_distribution(peak_intensity, peak_year, year_begin, year_end):
        """
        In this function, we create a callable function from the peak intensity distribution.
        The intensity distribution function is assumed to be symmetric around the peak year.
        Further, the upswing portion of intensity distribution function is assumed to first order increase.
        As such, it is parameterized as an expontial fit between a small input, here assumed to be ~ 1ng/kg/day.
        the known peak intensity @ the peak year
        I want to input the peak distributions and output the intensity distribution.
        If the intensity distribution is symettric, but the tail is clipped, I want to
        export the array of y values corresponding to the end.
        y_peak, y_begin, y_end are all in years.
        """
        # set the upswing and downswing times as a function of the given years.
        if (peak_year - year_begin) * 2 < (year_end - year_begin):
            print ("unclipped symmetric intensity distribution")
            upswing_time = (peak_year - year_begin) * 12

            # linearize the known datapoints and then create exponential function:
            x_up = np.array([0, upswing_time])
            y_up = np.log([1, peak_intensity])
            (slope_up, intercept_up, r_value, p_value, std_err) = stats.linregress(x_up, y_up)

            # interpolate using the provided slope and intercept to create the regression line
            x_up_interp = linspace(0, upswing_time, upswing_time)
            y_reg_line_up = polyval([slope_up, intercept_up], x_up_interp)
            upswing_reg_lin = np.exp(y_reg_line_up)  # take the exp of the reg line to return it to an exp fit

            # downswing side of symmetric intensity distribution
            x_down = linspace(upswing_time, upswing_time * 2, upswing_time)
            y_down = upswing_reg_lin[::-1]  # reversed array

            # concatenate the up and down swing sides
            y_up_down = np.append(upswing_reg_lin, y_down)
            x_up_down = linspace(0, upswing_time * 2, upswing_time * 2)

            # create a callable spline function for the up and down functions
            upswing_to_peak_intensity_spline = interpolate.InterpolatedUnivariateSpline(
                x_up_interp, upswing_reg_lin)

            downswing_to_peak_intensity_spline = interpolate.InterpolatedUnivariateSpline(
                x_down, y_down)

            return upswing_to_peak_intensity_spline, downswing_to_peak_intensity_spline


        else:
            print ("clipped symettric intensity distribution")
            # todo(add symettric intensity distribution when clipped)


    (upswing_callable_intake, downswing_callable_intake) = set_peak_intesity_distribution(80/1e9, 1977, 1930, 2100)

    #todo(sort out why for values around 6 does it oscillate between one and 2)

    intake_curve_track_mother = []
    def intake_amount(t, year_begin,year_peak,year_end,mass_mother,days_per_month, timesteps_per_month):

        return_value = []
        if t < (year_peak - year_begin)*12:
            return_value = upswing_callable_intake(t)*mass_mother*days_per_month/timesteps_per_month

        elif t >= year_peak*12 and t <2*year_peak*12:
            return_value = downswing_callable_intake(t)*mass_mother*days_per_month/timesteps_per_month

        else:
            return_value = 0

        intake_curve_track_mother.append(return_value)
        return return_value

    intake_curve_track_child = []

    def intake_amount_child(t, year_begin, year_peak, year_end, mass_child, days_per_month, timesteps_per_month):

        return_value = []
        if t < average_lact_time_months_:
            return_value = 0

        elif t > average_lact_time_months_ and t < (year_peak - year_begin) * 12:
            return_value = upswing_callable_intake(t) * mass_child * days_per_month / timesteps_per_month

        elif t >= (year_peak - year_begin) * 12 and t < 2 * (year_peak - year_begin) * 12:
            return_value = downswing_callable_intake(t) * mass_child * days_per_month / timesteps_per_month
        else:
            return_value = 0

        intake_curve_track_child.append(return_value)

        return return_value


    def mother_breastfeeding_window(k_lactation, t):
        if t < average_lact_time_months_:
            x = k_lactation
        else:
            x = 0

        return x

    # body mass balances
    # Set the time range
    t_start = 0.
    t_final = (80-25)*12 #age_max_in_years * 12
    print t_final

    # set the timestep & number of steps
    delta_t = 1. / timesteps_per_month
    num_steps = np.floor((t_final - t_start) / delta_t) + 1

    # set initial conditions
    # initial mass of fat in mother
    d_M_m_i = average_mother_bodyweight * average_mother_lipid_fraction

    # initial mass of chemical in mother
    d_M_chem_m_i = 0

    # initial mass of chemical in child
    d_C_chem_m_i = 0

    # initial mass of the child
    d_C_m_m_i = average_birth_lipid_fraction_ * average_birth_bodyweight_kg_

    def body_balance(t, y):
        #todo(put into function the k_elimination as a generic)

        k_elimination = np.log(2)/(5*12)
        k_lactation = 0.005

        n = len(y)  # 2: implies we have n ODEs
        dydt = np.zeros((n, 1))

        # mother mass, mass balance
        dydt[0] = 0  # mass of the mother is constant

        # mother chemical mass, mass balance
        dydt[1] = intake_amount(t, 1930,1977,2100,y[0],365/12,timesteps_per_month) - k_elimination * y[1] - mother_breastfeeding_window(k_lactation, t) * y[0]
        print t

        # child chemical mass, mass balance
        dydt[2] = intake_amount_child(t, 1930,1977,2100,y[3],365/12,timesteps_per_month) + k_elimination * y[3] + mother_breastfeeding_window(k_lactation, t) * y[0] - k_elimination * y[2]

        # child mass, mass balance (initial weight = 3.17kg)
        dydt[3] = spl_int_bw_lip_mass_deriv(t)

        return dydt

    print intake_curve_track_child

    # use ``vode`` with "backward differentiation formula"
    r = integrate.ode(body_balance).set_integrator('vode', method='bdf')
    r.set_initial_value([d_M_m_i, d_M_chem_m_i, d_C_chem_m_i, d_C_m_m_i], t_start)

    # Additional Python step: create vectors to store trajectories
    t = np.zeros((np.int(num_steps), 1))
    mother_mass = np.zeros((np.int(num_steps), 1))
    mass_chemical_in_mother = np.zeros((np.int(num_steps), 1))
    mass_chemical_in_child = np.zeros((np.int(num_steps), 1))
    child_mass = np.zeros((np.int(num_steps), 1))

    '''
    Mass balance equations:
    Mass of mother (fat mass only)
    # interpolated_bw_lip_mass_deriv
    d(M_mother)/dt =  spl_int_bw_lip_mass_deriv(timesteps_array)

    Mass of Chemical in Mother
    d(M_chemical_in_mother)/dt = Intake(t) - K_elimination * (M_mother) - K_lactation * (M_mother)

    Mass of child (fat mass only)
    d(M_child)/dt = interpreted mass.

    Mass of chemical in child during breastfeeding
    d(M_chemical_in_child)/dt = K_lactation * (M_mother) # for 6 months

    Mass of chemical in child after breastfeeding
    d(M_chemical_in_child/dt) = Intake(t) - K_elimination * (M_child)
    '''

    # Integrate the ODE(s) across each delta_t timestep
    #todo(figure out why the time, t is weird near 6 months)
    k = 1
    while r.successful() and k < num_steps:
        r.integrate(r.t + delta_t)

        t[k] = r.t
        mother_mass[k] = r.y[0]
        mass_chemical_in_mother[k] = r.y[1]
        mass_chemical_in_child[k] = r.y[2]
        child_mass[k] = r.y[3]
        k += 1

    # All done!  Plot the trajectories in two separate plots:
    ax1 = subplot(411)
    ax1.plot(t, mother_mass)
    ax1.set_xlim(t_start, t_final)
    ax1.set_xlabel('Time [minutes]')
    ax1.set_ylabel('Mother mass [kg]')
    ax1.grid('on')

    ax2 = plt.subplot(412)
    ax2.plot(t, mass_chemical_in_mother, 'r')
    ax2.set_xlim(t_start, t_final)
    ax2.set_xlabel('Time [minutes]')
    ax2.set_ylabel('Mass of chemical in mother [Kg]')
    ax2.grid('on')

    ax3 = plt.subplot(413)
    ax3.plot(t, mass_chemical_in_child, 'r')
    ax3.set_xlim(t_start, t_final)
    ax3.set_xlabel('Time [minutes]')
    ax3.set_ylabel('Mass of chemical in child [Kg]')
    ax3.grid('on')

    ax4 = plt.subplot(414)
    ax4.plot(t, child_mass, 'r')
    ax4.set_xlim(t_start, t_final)
    ax4.set_xlabel('Time [minutes]')
    ax4.set_ylabel('Mass of child [Kg]')
    ax4.grid('on')



    plt.show()
    plt.close('all')

    plt.plot(intake_curve_track_child)
    plt.show()
main()
