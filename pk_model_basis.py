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
from functions import *

"""
Notes: It doesn't make sense that the first order up/down coefficients are statically defined. The upgradient curve
would have eq: y = e^(kdec_up*x). If you fit a first order coefficient, such that the curve started at 0 and went
through the peak intensity, it would result in only one curve. Regarldess of the kdec_up value.

I've also made is such that you can dynamically assign the number of gens and year of childbirth.
This will automatically make the ODE's for the 4 compartments and calculate them. There's plenty yet to do, but
the basic form of the model is complete.

"""

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
_cong_spag = os.path.join(
    'cogener_start_peak_age_group.csv')
_cte = ['hcb', 'ppddt',
        'ppdde', 'betahch', 'gammahch', 'pcb138']

_CONGENERS_HALFLIFE_DICT = {'hcb': [10, np.log(2) / 5, np.log(2) / 5],
                            'ppddt': [5, np.log(2) / 5, np.log(2) / 5]}

_cte = list(_CONGENERS_HALFLIFE_DICT.keys())

TIMESTEPS_PER_MONTH = 10
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
body_w_ave = 3.1744


def main(n_optimization_runs=N_OPTIMIZATION_RUNS_,
         bm_data=BIOMONITORING_DATA_ELIMINATION_,
         cong_spag=_cong_spag,
         cte=_cte,  # congoners to evaluate
         congeners_half_life=_CONGENERS_HALFLIFE_DICT,
         kinetics_order=KINETICS_ORDER_,
         gens=1,
         plot_kinetics=False,
         study_start_year=STUDY_START_YEAR_,
         study_end_year=STUDY_END_YEAR_,
         average_lact_time_months_=AVERAGE_LACT_TIME_MONTHS_,
         lifespan=AGE_MAX_IN_YEARS_,
         prg_intrvl=AGE_MOM_IN_YEARS_,
         average_mother_bodyweight=REF_FEMALE_BODYWEIGHT_KG_,
         average_mother_lipid_fraction=AVERAGE_MOTHER_LIPID_FRACTION_,
         average_birth_lipid_fraction_=AVERAGE_BIRTH_LIPID_FRACTION_,
         body_w_ave=body_w_ave,
         bw_frac_lip=BODYWEIGHT_AND_LIPID_FRACTION_,
         timesteps_per_month=TIMESTEPS_PER_MONTH):


    t_start = study_start_year
    t_final = study_end_year
    delta_t = 1  # / timesteps_per_month (e.g., 1 = 1./1. timesteps per month)
    num_steps = np.floor((t_final - t_start) / delta_t) + 1

    if bm_data:
        bm_data_eval(bm_data,
                     cte,
                     kinetics_order,
                     plot_kinetics)

    cong_spag_df = pd.read_csv(cong_spag)

    bm_data_err(cte, cong_spag_df)

    if bw_frac_lip:
        def mass_spline(study_start_year,
                        study_end_year,
                        bw_frac_lip,
                        t_step,
                        body_w_ave,
                        average_birth_lipid_fraction_):
            bw_frac_lip = pd.read_csv(bw_frac_lip)
            bw_frac_lip['mass_lipids_kg'] = bw_frac_lip.apply(cal_m_lip, 1)

            # step 1: paramterize growth with spline interpretation
            y_bm = np.array(bw_frac_lip['mass_lipids_kg'])

            x_bm = np.linspace(study_start_year, study_end_year, num=len(y_bm))
            spl_bw = interpolate.InterpolatedUnivariateSpline(
                x_bm, y_bm)

            # step 2: get derivative of spline (for dydt kinetics)
            spl_bw_deriv = spl_bw.derivative()

            return spl_bw_deriv, spl_bw

        (spl_bw_deriv, spl_bw) = mass_spline(study_start_year,
                                             study_end_year,
                                             bw_frac_lip,
                                             timesteps_per_month,
                                             body_w_ave,
                                             average_birth_lipid_fraction_)

    # TODO(Automate entry of peak year and peak max)
    intakeCall = intake_int_dist(80, 1977, t_start, t_final+average_lact_time_months_, delta_t)

    def generation_mass_balance(y,
                                congener,
                                gens,
                                simulation_start,
                                simulation_end,
                                delta_t,
                                intakeCall,
                                bw_spl_der,
                                prg_intrvl,
                                lifespan,
                                num_steps):

        # dynamically set the number of odes as a function of gens
        n = len(y) * gens
        num_odes_in_gen = 4
        n = num_odes_in_gen * gens
        dydt = np.zeros((n, 1))  # initialize dydt array

        # simulation bounds
        t_start = np.float(simulation_start)
        t_final = np.float(simulation_end)
        simulation_time_total = t_final - t_start

        # create simulation timestamps for important events
        # when people are born and when they die.
        # this defines the start array dydt numbering system.
        start_array_dydt = []
        aig_mother = []
        aigd_mother = []
        ode_birth_lactation_death_limits = []

        start_array_dydt = linspace(0, 4 * (gens - 1), gens)
        # print start_array_dydt

        # age that a mother gives birth in each gen
        aig_mother = linspace(0,
                              prg_intrvl * (gens),
                              gens + 1)

        # array of child births for each gen
        cbtg_child = linspace(prg_intrvl, prg_intrvl * (gens + 1),
                              gens + 1)

        # array of ages of death of the mother for each gen
        aigd_mother = linspace(
            lifespan, (lifespan + (25 * gens)), gens + 1)

        print "start_array_dydt", start_array_dydt
        print "aig_mother", aig_mother
        print "cbtg_child", cbtg_child
        print "aigd_mother", aigd_mother

        """
        Note1: Ok. I've figured out what to do. I need to make a loop that iteratively
        offsets the start of each new generations mass. Effectively, this means 
        making a new mass spline for each generation that is 0 until the next 
        generation is born. To be honest, I'm not sure how to do that. Perhaps, 
        I can initialize a matrix with each row as each generations subsequent spline. 
        Then I make a loop where 0s are added to the front and back of each spline.

        Note2: I've figured out Note1. However, this creates time-scale problems in the plot.
        As such, the plots will need to be removed from their original timescales and re-plotted.
        """

        for gen in range(0, gens + 1):
            intakeCall = intake_int_dist(80, 1977, t_start, t_final+ aig_mother[gen], delta_t)
            print intakeCall(1)
            (spl_bw_deriv_child, spl_bw_child) = mass_spline(study_start_year + aig_mother[gen],
                                                             study_end_year + aig_mother[gen],
                                                             bw_frac_lip,
                                                             timesteps_per_month,
                                                             body_w_ave,
                                                             average_birth_lipid_fraction_)

            if np.all(gens >= 1):
                start_dydt_gen = []
                start_dydt_gen.append(start_array_dydt)
                for i in range(1, gens - 1):
                    start_dydt_gen.append([x + i for x in start_array_dydt])
                start_dydt_gen = np.array(start_dydt_gen)

        odes_per_gen = range(0, num_odes_in_gen)
        dydt_matrix = np.zeros(shape=(len(odes_per_gen), gens), dtype=object)

        order_array_counter = np.array(
            range(0, gens * len(odes_per_gen)))

        itr_mtrx = start_dydt_gen

        def body_mass(t, y):
            cntr = 0
            for gen in range(0, gens):

                k_elim = np.log(2) / 5
                k_lac = 1e-1

                # Note; unused, but provided for tracking purposes
                ode_numbers = start_dydt_gen[:, gen]
                cgbt = aig_mother[gen]  # age the mother gives birth
                ngbt = cbtg_child[gen]
                cgdt = aigd_mother[gen]

                if gen == 0:
                    for ode in odes_per_gen:
                        dydt_matrix[0][gen] = bw_spl_der(t)
                        dydt_matrix[1][gen] = intakeCall(t) * y[0] - k_elim * y[0] - k_lac * y[0]
                        dydt_matrix[2][gen] = k_lac * y[0] - k_elim * y[2]

                        # todo(create a function that starts the mass at t+birthtime*gen and is 0 before)
                        dydt_matrix[3][gen] = bw_spl_der(t)

                elif gen >= 0:
                    t_r = t + gen * 25

                    cntr = np.int(cntr + 1)
                    # print np.int(itr_mtrx[0][cntr] - 1), 'cntr'
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

                    '''
                    dydt_matrix[0][gen] = bw_spl_der(t_r)
                    dydt_matrix[1][gen] = intakeCall(t_r) - k_elim * y[np.int(itr_mtrx[0]
                                                                              [cntr] - 1)] - k_lac * y[
                        np.int(itr_mtrx[0][cntr] - 1)]
                    dydt_matrix[2][gen] = k_lac * y[np.int(itr_mtrx[0]
                                                           [cntr] - 2)] - k_elim * y[np.int(itr_mtrx[0][cntr] - 1)]
                    dydt_matrix[3][gen] = bw_spl_der(t_r)

            dydt = np.ravel(dydt_matrix, order='F')

            return dydt

        t = np.zeros((np.int(num_steps), 1))

        # use ``vode`` with "backward differentiation formula" or 'bdf'
        r = integrate.ode(body_mass).set_integrator('vode', order=4, nsteps=num_steps, min_step=1e-14,
                                                    method='bdf')
        y0 = np.zeros((np.int(gens * 4), 1))
        r.set_initial_value(y0, t_start)

        # create vectors to store trajectories
        ode_init = np.zeros((np.int(num_steps) * gens * 4))
        ode_init_matrix = ode_init.reshape((4 * gens, np.int(num_steps)), order='F')
        iter_odes = range(0, 4 * gens, 1)

        # initialize k for while loop
        k = 1
        while r.successful() and k < num_steps:
            r.integrate(r.t + delta_t)
            t[k] = r.t
            for ode in iter_odes:
                ode_init_matrix[ode][k] = r.y[ode]
            k += 1



        for ode in iter_odes:
            ax1 = plt.subplot(len(iter_odes),1, iter_odes[ode]+1)
            plt.plot(t, ode_init_matrix[ode][:])
            ax1.plot(t, ode_init_matrix[ode][:])
            # ax1.legend(iter_odes[ode])
            ax1.set_xlim(t_start, t_final)
            ax1.grid('on')
            # plt.show()

        plt.xlim(t_start, t_final)
        plt.legend(iter_odes)
        plt.show()

        return (y, t)

    y = []

    (y, t) = generation_mass_balance(y=y,
                                     congener='hcb',
                                     gens=1,
                                     simulation_start=1921,
                                     simulation_end=2100,
                                     delta_t=delta_t,
                                     intakeCall=intakeCall,
                                     bw_spl_der=spl_bw_deriv,
                                     prg_intrvl=prg_intrvl,
                                     lifespan=lifespan,
                                     num_steps=num_steps)

    '''
    Mass balance equations:
    Mass of mother (fat mass only)
    # interpolated_bw_lip_mass_deriv
    d(M_mother)/dt =  spl_bw_deriv(timesteps_array)

    Mass of Chemical in Mother
    d(M_chemical_in_mother)/dt = Intake(t) - K_elimination * (M_mother) - K_lactation * (M_mother)

    Mass of child (fat mass only)
    d(M_child)/dt = interpreted mass.

    Mass of chemical in child during breastfeeding
    d(M_chemical_in_child)/dt = K_lactation * (M_mother) # for 6 months

    Mass of chemical in child after breastfeeding
    d(M_chemical_in_child/dt) = Intake(t) - K_elimination * (M_child)
    '''

    # # All done!  Plot the trajectories in two separate plots:
    # ax1 = subplot(511)
    # ax1.plot(t / 12, y[0])
    # ax1.set_xlim(t_start, t_final / 12)
    # ax1.set_xlabel('Age [years]')
    # ax1.set_ylabel('Mother mass [kg]')
    # ax1.grid('on')


main()
