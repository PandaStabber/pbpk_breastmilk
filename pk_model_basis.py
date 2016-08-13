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

I've also made is such that you can dynamically assign the number of generations and year of childbirth.
This will automatically make the ODE's for the 4 compartments and calculate them. There's plenty yet to do, but
the basic form of the model is complete.

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
BIOMONITORING_DATA_ELIMINATION_UNITS = [
    'ng/kg_bw/d', 'years']  # just to keep track
_CONGENER_START_PEAK_AGE_GROUP = os.path.join(
    'cogener_start_peak_age_group.csv')
_CONGENERS_TO_EVALUATE = ['hcb', 'ppddt',
                          'ppdde', 'betahch', 'gammahch', 'pcb138']

_CONGENERS_HALFLIFE_DICT = {'hcb': [10, np.log(2) / 5, np.log(2) / 5],
                            'ppddt': [5, np.log(2) / 5, np.log(2) / 5]}

_CONGENERS_TO_EVALUATE = list(_CONGENERS_HALFLIFE_DICT.keys())

TIMESTEPS_PER_MONTH = 1
N_MALES_ = 100
N_FEMALES_ = 100
POP_START_YEAR_ = 1921
STUDY_START_YEAR_ = 1921
STUDY_END_YEAR_ = 2100
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
         congeners_half_life=_CONGENERS_HALFLIFE_DICT,
         kinetics_order=KINETICS_ORDER_,
         generations=1,
         plot_kinetics=False,
         study_start_year=STUDY_START_YEAR_,
         study_end_year=STUDY_END_YEAR_,
         average_lact_time_months_=AVERAGE_LACT_TIME_MONTHS_,
         age_max_in_years=AGE_MAX_IN_YEARS_,
         age_mom_in_years=AGE_MOM_IN_YEARS_,
         average_mother_bodyweight=REF_FEMALE_BODYWEIGHT_KG_,
         average_mother_lipid_fraction=AVERAGE_MOTHER_LIPID_FRACTION_,
         average_birth_lipid_fraction_=AVERAGE_BIRTH_LIPID_FRACTION_,
         average_birth_bodyweight_kg_=AVERAGE_BIRTH_BODYWEIGHT_KG_,
         bodyweight_and_lipid_fraction=BODYWEIGHT_AND_LIPID_FRACTION_,
         timesteps_per_month=TIMESTEPS_PER_MONTH):

    # TODO(build in argparse here for csv extension)

    kinetics_from_biomonitoring_data = []

    # set the timeline and step parameters for the study for the study
    t_start = study_start_year
    t_final = study_end_year
    delta_t = 1.  # / timesteps_per_month
    num_steps = np.floor((t_final - t_start) / delta_t) + 1

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

        def mass_spline(study_start_year, study_end_year, bodyweight_and_lipid_fraction):
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
            timesteps = timesteps_per_month * \
                len(bodyweight_and_lipid_fraction)
            timesteps_array = np.linspace(1, len(y_bm), timesteps)
            interpolated_bw_lip_mass = spl_int_bw_lip_mass(timesteps_array)
            delta_bw_lip_mass = shift(interpolated_bw_lip_mass, 1, cval=np.NaN)

            # replace the first value (which was shifted).
            # calculate the first delta mass then calculate the kinetics
            delta_bw_lip_mass[0] = bodyweight_and_lipid_fraction['mass_lipids_kg'][
                0] - average_birth_bodyweight_kg_ * average_birth_lipid_fraction_

            interpolated_bw_lip_mass_deriv = spl_int_bw_lip_mass_deriv(
                timesteps_array)

            growth_kinetics_df = pd.DataFrame()
            growth_kinetics_df['mass_lip_kg'] = interpolated_bw_lip_mass
            growth_kinetics_df[
                'growth_kinetics_dKg/dt'] = interpolated_bw_lip_mass_deriv

            return spl_int_bw_lip_mass_deriv
            # TODO(include function to print bodyweight and derivative (kinetics))
            # plt.plot(timesteps_array, interpolated_bw_lip_mass, 'b')
            # plt.plot(timesteps_array, interpolated_bw_lip_mass_deriv, 'r')
            # plt.show()

    else:
        print('''please include .csv file containing the'average_body_weight_in_kg'
            in the first column and the 'lipid_fraction'
            in the second column''').replace(
            '/t', '')
    spl_int_bw_lip_mass_deriv = mass_spline(study_start_year, study_end_year,
                                            bodyweight_and_lipid_fraction)

    # [years] age in years with interval of 1 month
    age_years_in_months = np.linspace(start=1.0 / 12,
                                      stop=age_max_in_years,
                                      num=age_max_in_years * 12,
                                      endpoint=True)
    lact_kinetics = pd.DataFrame()
    lact_kinetics['age_index_years'] = age_years_in_months
    lact_kinetics['age_index_months'] = lact_kinetics['age_index_years'] * 12

    def set_peak_intesity_distribution(peak_intensity, peak_year, year_begin, year_end, delta_t):
        # todo(remove peak_year,year_begin,year_end, and replace it with
        # default names in function)
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
            print("unclipped symmetric intensity distribution")
            total_time_years = (year_end - year_begin)
            total_time_months = total_time_years * 12
            num_steps = total_time_months / delta_t

            # set total step space for simulation
            x_up_down = linspace(
                start=0, stop=total_time_months, num=num_steps)

            upswing_time_years = peak_year - year_begin
            upswing_time_months = upswing_time_years * 12
            upswing_steps = upswing_time_months / delta_t
            print upswing_steps

            downswing_steps = upswing_steps  # if symmetric

            # linearize the known datapoints and then create exponential
            # function:
            x_up = np.array([0, upswing_time_months])
            y_up = np.log([1, peak_intensity])
            (slope_up, intercept_up, r_value, p_value,
             std_err) = stats.linregress(x_up, y_up)

            # interpolate using the provided slope and intercept to create the
            # regression line
            x_up_interp = linspace(0, upswing_time_months, upswing_steps)
            y_reg_line_up = polyval([slope_up, intercept_up], x_up_interp)
            # take the exp of the reg line to return it to an exp fit
            upswing_reg_lin = np.exp(y_reg_line_up)

            # downswing side of symmetric intensity distribution
            x_down = linspace(upswing_time_months,
                              upswing_time_months * 2, downswing_steps)
            y_down = upswing_reg_lin[::-1]  # reversed array

            # the remaining intensity is set to zero
            remaining_steps = num_steps - upswing_steps - downswing_steps
            y_remaining = np.zeros(np.int(remaining_steps), dtype=np.int)

            # concatenate the up and down swing sides
            y_up_down = np.concatenate((upswing_reg_lin, y_down, y_remaining))

            # create a callable spline function for the up and down functions
            peak_intensity_spline = interpolate.InterpolatedUnivariateSpline(
                x_up_down, y_up_down)

            return peak_intensity_spline

        else:
            print("clipped symettric intensity distribution")
            # todo(add symettric intensity distribution when clipped)

    # todo(sort out why for values oscillate above and below 6)

    intake_curve_track_mother = []
    intake_curve_track_child = []

    (upswing_callable_intake) = set_peak_intesity_distribution(
        80, 1977, 1921, 2100, delta_t)

    def generation_mass_balance(y, congener, generations, simulation_start, simulation_end):

        # the number of ODEs is equal to 4 * generations
        n = len(y) * generations  # 2: implies we have n ODEs
        num_odes_in_gen = 4
        n = num_odes_in_gen * generations
        dydt = np.zeros((n, 1))
        # print dydt

        # simulation bounds
        t_start = np.float(simulation_start)
        t_final = np.float(simulation_end)
        simulation_time_total = t_final - t_start
        print t_start, t_final, simulation_time_total

        intake = 1

        # create simulation timestamps for important events
        # when people are born and when they die.
        # this defines the start array dydt numbering system.
        start_array_dydt = []
        array_mothers_age_in_each_gen = []
        array_mothers_death_in_each_gen = []
        ode_birth_lactation_death_limits = []

        # for i in range(0,generations+1):
        start_array_dydt = linspace(0, 4 * (generations - 1), generations)
        array_mothers_age_in_each_gen = linspace(
            0, age_mom_in_years * (generations), generations + 1)
        array_children_birth_time_in_each_gen = linspace(age_mom_in_years, age_mom_in_years * (generations + 1),
                                                         generations + 1)
        array_mothers_death_in_each_gen = linspace(age_max_in_years,
                                                   age_max_in_years + (age_mom_in_years * generations), generations + 1)

        print "start_array_dydt", start_array_dydt
        print "array_mothers_age_in_each_gen", array_mothers_age_in_each_gen
        print "array_children_birth_time_in_each_gen", array_children_birth_time_in_each_gen
        print "array_mothers_death_in_each_gen", array_mothers_death_in_each_gen

        if generations >= 1:
            start_array_dydt_gen = []
            start_array_dydt_gen.append(start_array_dydt)
            for i in range(1, generations + 1):
                start_array_dydt_gen.append([x + i for x in start_array_dydt])
            start_array_dydt_gen = np.array(start_array_dydt_gen)

        odes_per_gen = range(0, num_odes_in_gen)
        dydt_matrix = np.zeros(
            shape=(len(odes_per_gen), generations), dtype=object)

        order_array_counter = np.array(
            range(0, generations * len(odes_per_gen)))
        iteration_matrix = order_array_counter.reshape(
            (len(odes_per_gen), generations), order='F')

        #
        # print iteration_matrix
        # print 'generations', range(0, generations)

        # def mother_timeline(t, counter):
        #
        #     if t < current_gen_death_time:
        #         print 'option a'
        #         return spl_int_bw_lip_mass_deriv(t - next_gen_birth_time)
        #     else:
        #         print 'shutoff'
        #         y[iteration_matrix[0][counter]] = 0

        def body_mass(t, y):
            # print t
            counter = 0
            for generation in range(0, generations):

                k_elimination = 1e-3
                k_lactation = 1e-4
                ode_numbers = start_array_dydt_gen[:, generation]
                current_gen_birth_time = array_mothers_age_in_each_gen[
                    generation]
                next_gen_birth_time = array_children_birth_time_in_each_gen[
                    generation]
                current_gen_death_time = array_mothers_death_in_each_gen[
                    generation]

                # print 'generation:', generation, 'counter:',counter
                if generation == 0:
                    for ode in odes_per_gen:
                        # print next_gen_birth_time
                        print spl_int_bw_lip_mass_deriv(t)
                        dydt_matrix[0][counter] = spl_int_bw_lip_mass_deriv(t)
                        dydt_matrix[1][counter] = - \
                            k_elimination * y[0] - k_lactation * y[0]
                        dydt_matrix[2][counter] = k_lactation * \
                            y[0] - k_elimination * y[2]
                        dydt_matrix[3][counter] = next_gen_birth_time
                        # next_gen_birth_time#spl_int_bw_lip_mass_deriv(t -
                        # next_gen_birth_time)

                elif generation > 0:
                    counter = counter + 1
                    dydt_matrix[0][
                        counter] = current_gen_birth_time  # iteration_matrix[0][counter]#spl_int_bw_lip_mass_deriv(t)
                    dydt_matrix[1][counter] = intake - k_elimination * y[
                        iteration_matrix[1][counter] - 1] - k_lactation * y[iteration_matrix[1][counter] - 1]
                    dydt_matrix[2][counter] = intake + k_lactation * y[
                        iteration_matrix[2][counter] - 2] - k_elimination * y[iteration_matrix[2][counter]]
                    # iteration_matrix[3][counter]
                    dydt_matrix[3][counter] = next_gen_birth_time

            dydt = np.ravel(dydt_matrix, order='F')

            # print dydt
            return dydt

        t = np.zeros((np.int(num_steps), 1))
        # use ``vode`` with "backward differentiation formula" or 'bdf'

        r = integrate.ode(body_mass).set_integrator('vode', order=5, nsteps=num_steps, max_step=delta_t, min_step=0,
                                                    method='bdf')
        r.set_initial_value([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], t_start)

        # create vectors to store trajectories
        t = np.zeros((np.int(num_steps), 1))

        mm_0 = np.zeros((np.int(num_steps), 1))
        mcm_0 = np.zeros((np.int(num_steps), 1))
        mcc_0 = np.zeros((np.int(num_steps), 1))
        mc_0 = np.zeros((np.int(num_steps), 1))

        mm_1 = np.zeros((np.int(num_steps), 1))
        mcm_1 = np.zeros((np.int(num_steps), 1))
        mcc_1 = np.zeros((np.int(num_steps), 1))
        mc_1 = np.zeros((np.int(num_steps), 1))

        mm_2 = np.zeros((np.int(num_steps), 1))
        mcm_2 = np.zeros((np.int(num_steps), 1))
        mcc_2 = np.zeros((np.int(num_steps), 1))
        mc_2 = np.zeros((np.int(num_steps), 1))

        k = 1
        while r.successful() and k < num_steps:
            r.integrate(r.t + delta_t)

            t[k] = r.t
            mm_0[k] = r.y[0]
            mcm_0[k] = r.y[1]
            mcc_0[k] = r.y[2]
            mc_0[k] = r.y[3]

            mm_1[k] = r.y[4]
            mcm_1[k] = r.y[5]
            mcc_1[k] = r.y[6]
            mc_1[k] = r.y[7]

            mm_2[k] = r.y[8]
            mcm_2[k] = r.y[9]
            mcc_2[k] = r.y[10]
            mc_2[k] = r.y[11]

            # todo(create these automatically based on the number of
            # generations)

            k += 1

        ax1 = subplot(511)
        ax1.plot(t, mm_0)
        ax1.set_xlim(t_start, t_final)
        ax1.set_xlabel('Age [years]')
        ax1.set_ylabel('Mother mass [kg]')
        ax1.grid('on')

        ax2 = subplot(512)
        ax2.plot(t, mcm_0)
        ax2.set_xlim(t_start, t_final)
        ax2.set_xlabel('Age [years]')
        ax2.set_ylabel('Mother mass [kg]')
        ax2.grid('on')

        plt.show()

    generation_mass_balance(y, 'hcb', 3, t_start, t_final)

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

    # All done!  Plot the trajectories in two separate plots:
    ax1 = subplot(511)
    ax1.plot(t / 12, mother_mass)
    ax1.set_xlim(t_start, t_final / 12)
    ax1.set_xlabel('Age [years]')
    ax1.set_ylabel('Mother mass [kg]')
    ax1.grid('on')

    ax2 = plt.subplot(512)
    ax2.plot(t / 12, mass_chemical_in_mother, 'r')
    ax2.set_xlim(t_start, t_final / 12)
    ax2.set_xlabel('Age [years]')
    ax2.set_ylabel('Mass of chemical in mother [Kg]')
    ax2.grid('on')

    ax3 = plt.subplot(513)
    ax3.plot(t / 12, mass_chemical_in_child, 'r')
    ax3.set_xlim(t_start, t_final / 12)
    ax3.set_xlabel('Age [years]')
    ax3.set_ylabel('Mass of chemical in child [Kg]')
    ax3.grid('on')

    ax4 = plt.subplot(514)
    ax4.plot(t / 12, child_mass, 'r')
    ax4.set_xlim(t_start, t_final / 12)
    ax4.set_xlabel('Age [years]]')
    ax4.set_ylabel('Mass of child [Kg]')
    ax4.grid('on')

    ax5 = plt.subplot(515)
    ax5.plot(intake_curve_track_child, 'r')
    # ax5.set_xlim(t_start, t_final/12)
    ax5.set_xlabel('Age [years]]')
    ax5.set_ylabel('Intake Intensity [ng/kg/m]')
    ax5.grid('on')

    # plt.show()
    # plt.close('all')

    # plt.plot(intake_curve_track_child)
    # plt.show()


main()
