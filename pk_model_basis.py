# coding=utf-8

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
cte_ = ['hcb', 'ppddt',
        'ppdde', 'betahch', 'gammahch', 'pcb138']

_CONGENERS_HALFLIFE_DICT = {'hcb': [10, np.log(2) / 5, np.log(2) / 5],
                            'ppddt': [5, np.log(2) / 5, np.log(2) / 5]}

cte_ = list(_CONGENERS_HALFLIFE_DICT.keys())

TIMESTEPS_PER_MONTH = 1
N_MALES_ = 100
N_FEMALES_ = 100
POP_START_YEAR_ = 0
STUDY_START_YEAR_ = 0
STUDY_END_YEAR_ = 2100 - 1921
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


def main(bm_data=BIOMONITORING_DATA_ELIMINATION_,
         cong_spag=_cong_spag,
         cte=cte_,
         kinetics_order=KINETICS_ORDER_,
         gens=1,
         plot_kinetics=False,
         study_start_year=STUDY_START_YEAR_,
         study_end_year=STUDY_END_YEAR_,
         lifespan=AGE_MAX_IN_YEARS_,
         prg_intrvl=AGE_MOM_IN_YEARS_,
         bw_frac_lip=BODYWEIGHT_AND_LIPID_FRACTION_):
    t_start = study_start_year
    t_final = study_end_year
    delta_t = .1  # / timesteps_per_month (e.g., 1 = 1./1. timesteps per month)
    num_steps = np.floor((t_final - t_start) / delta_t) + 1

    if bm_data:
        bm_data_eval(bm_data,
                     cte,
                     kinetics_order,
                     plot_kinetics)

    cong_spag_df = pd.read_csv(cong_spag)

    bm_data_err(cte, cong_spag_df)

    intakeCall = asym_int_dist(80, 50, t_start, t_final, delta_t)

    def generation_mass_balance(y,
                                congener,
                                gens,
                                simulation_start,
                                simulation_end,
                                delta_t,
                                intakeCall,
                                # bw_spl_der,
                                prg_intrvl,
                                lifespan,
                                num_steps):

        # dynamically set the number of odes as a function of gens
        n = len(y) * gens
        num_odes_in_gen = 2
        n = num_odes_in_gen * gens
        dydt = np.zeros((n, 1))  # initialize dydt array

        # simulation bounds
        t_start = np.float(simulation_start)
        t_final = np.float(simulation_end)

        aig_mother = []
        aigd_mother = []

        start_array_dydt = linspace(0, num_odes_in_gen * (gens - 1), gens)
        aig_mother = linspace(0,
                              prg_intrvl * (gens),
                              gens + 1)

        # array of child births for each gen
        cbtg_child = linspace(prg_intrvl,
                              prg_intrvl * (gens + 1),
                              gens + 1)

        # array of ages of death of the mother for each gen
        # array of ages of death of the mother for each gen
        aigd_mother = linspace(
            lifespan, (lifespan + (25 * gens)), gens + 1)

        print("start_array_dydt", start_array_dydt)
        print("aig_mother", aig_mother)
        print("cbtg_child", cbtg_child)
        print("aigd_mother", aigd_mother)

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

        order_array_counter = np.array(range(0,
                                             gens * len(odes_per_gen)))

        itr_mtrx = order_array_counter.reshape((len(odes_per_gen),
                                                gen),
                                               order='F')

        print(bw_frac_lip)

        bw_spl_der = age_splines(gens,
                                 aig_mother,
                                 aigd_mother,
                                 t_start,
                                 t_final,
                                 delta_t,
                                 bw_frac_lip)
        print(bw_spl_der)

        def body_mass(t, y):
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
            for gen in range(0, gens):

                k_elim = np.log(2) / 5
                k_lac = 0
                # if (t >= cbtg_child[gen]) & (t <= cbtg_child[gen] + 0.5):
                #     k_lac = 1e-1

                if gen == 0:

                    dydt_matrix[0][cntr] = bw_spl_der[0](t)
                    dydt_matrix[1][cntr] = intakeCall(t) * y[0] \
                                           - k_elim * y[1] \
                                           - k_lac * y[1]
                    cntr = np.int(cntr + 1)

                elif gen > 0:
                    dydt_matrix[0][cntr] = bw_spl_der[gen](t)
                    dydt_matrix[1][cntr] = intakeCall(t) * y[np.int(itr_mtrx[0][cntr])] + \
                                           k_lac * y[np.int(itr_mtrx[1][cntr - 1])] \
                                           - k_elim * y[
                        np.int(itr_mtrx[1][cntr])]

                    cntr = np.int(cntr + 1)

            dydt = np.ravel(dydt_matrix, order='F')

            return dydt

        t = np.zeros((np.int(num_steps), 1))

        # use ``vode`` with "backward differentiation formula" or 'bdf'
        r = integrate.ode(body_mass).set_integrator('vode',
                                                    order=4,
                                                    nsteps=num_steps,
                                                    min_step=1e-8,
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

    y = []

    (y, t) = generation_mass_balance(y=y,
                                     congener='hcb',
                                     gens=4,
                                     simulation_start=t_start,
                                     simulation_end=t_final,
                                     delta_t=delta_t,
                                     intakeCall=intakeCall,
                                     prg_intrvl=prg_intrvl,
                                     lifespan=lifespan,
                                     num_steps=num_steps)


main()
